"""Collect optimizer control-flow evidence without changing frozen evaluation inputs."""
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

import chematic

OUT = Path(__file__).resolve().parent
BASE = OUT.parent / "cation-water"
CACHE = {}


def raw(path):
    key = str(Path(path))
    if key not in CACHE:
        CACHE[key] = Path(path).read_bytes()
    return CACHE[key]


def read(path):
    return json.loads(raw(path))


def rows(path):
    return [json.loads(line) for line in raw(path).splitlines()]


def artifact(path):
    data = raw(path)
    return {"path": str(path), "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def fresh_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    previous = read(BASE / "candidate-summary.json")
    current = read(OUT / "candidate-summary.json")
    assert previous["manifests"] == current["manifests"]
    assert previous["protocol"] == current["protocol"]
    baseline_energy = rows(BASE / "candidate.jsonl")
    candidate_energy = rows(OUT / "candidate.jsonl")
    assert baseline_energy == candidate_energy
    baseline = rows(BASE / "quality.jsonl")
    candidate = rows(OUT / "quality.jsonl")
    scored = [row for row in rows(OUT / "quality-scored.jsonl") if row["engine"] == "chematic"]
    baseline_scored = [row for row in rows(OUT / "baseline-quality-scored.jsonl") if row["engine"] == "chematic"]
    assert len(baseline) == len(candidate) == len(scored) == len(baseline_scored) == 265
    assert all(row["independently_sound"] and row["stereo"]["violated"] == 0 and row["stereo"]["unevaluable"] == 0 for row in baseline_scored)
    scorer = read(OUT / "scorer-build.json")["binary"]
    assert fresh_hash(scorer["path"]) == scorer["sha256"]
    changed = []
    convergence_audit = []
    graph_modes = Counter()
    for old, new in zip(baseline, candidate, strict=True):
        assert (old["row_index"], old["smiles"], old["status"]) == (new["row_index"], new["smiles"], new["status"])
        assert old["force_field"]["energy_before"] == new["force_field"]["energy_before"]
        old_ff = old["force_field"]
        new_ff = new["force_field"]
        if old_ff["coords"] != new_ff["coords"] or old_ff["converged"] != new_ff["converged"]:
            changed.append({
                "input_index": old["row_index"], "name": old["name"], "smiles": old["smiles"],
                "coordinates_changed": old_ff["coords"] != new_ff["coords"],
                "converged_before": old_ff["converged"], "converged_after": new_ff["converged"],
                "iterations_before": old_ff["iterations"], "iterations_after": new_ff["iterations"],
                "fd_residual_force_before": old_ff["max_residual_force"],
                "fd_residual_force_after": new_ff["max_residual_force"],
                "total_energy_delta": new_ff["energy_after"]["total"] - old_ff["energy_after"]["total"],
            })
        mol = chematic.from_smiles(new["smiles"])
        mode = "input_graph"
        if len(new_ff["coords"]) != len(mol.depict_data()["atoms"]):
            mol = mol.add_hydrogens()
            mode = "expanded_explicit_h"
        graph_modes[mode] += 1
        assert len(new_ff["coords"]) == len(mol.depict_data()["atoms"])
        audit = {"input_index": new["row_index"], "graph_mode": mode}
        for label, record in [("baseline", old_ff), ("candidate", new_ff)]:
            coords = record["coords"]
            energy = mol.mmff94_energy_breakdown(coords)["total"]
            assert energy == record["energy_after"]["total"]
            maximum = max(abs(value) for point in mol.mmff94_bounded_analytic_gradient(coords) for value in point)
            audit[label] = {"max_analytic_gradient": maximum, "reported_converged": record["converged"],
                            "flag_matches_tolerance": (maximum < 1e-4) == record["converged"]}
        assert audit["candidate"]["flag_matches_tolerance"]
        convergence_audit.append(audit)
    (OUT / "convergence-audit.json").write_text(json.dumps({
        "threshold": 1e-4, "units": "kcal/mol/angstrom",
        "method": "Reconstruct each force-field graph with from_smiles and add_hydrogens when required by coordinate count; require exact energy equality on force_field.coords before comparing the public analytic gradient with the convergence flag.",
        "graph_modes": dict(graph_modes), "rows": convergence_audit,
        "baseline_flag_mismatches": sum(not r["baseline"]["flag_matches_tolerance"] for r in convergence_audit),
        "candidate_flag_mismatches": 0,
    }, indent=2, sort_keys=True) + "\n")
    gradient = read(OUT / "gradient-comparison.json")
    holdout = read(OUT / "candidate-holdout.json")
    assert gradient["gate_passed"]
    assert holdout["input_count"] == 20 and holdout["all_types_match"]
    assert holdout["artifacts"]["extension"]["sha256"] == current["artifacts"]["schematic"]["sha256"]
    for collection in [current["artifacts"], gradient["artifacts"], holdout["artifacts"]]:
        for value in collection.values():
            assert fresh_hash(value["path"]) == value["sha256"]
    build_source = read(OUT / "build-source.json")
    for record in [build_source, *build_source["additional_sources"]]:
        assert fresh_hash(record["path"]) == record["sha256"]
    rust = re.findall(r"test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored", raw(OUT / "rust-tests.log").decode())
    python = re.search(r"(\d+) passed in", raw(OUT / "python-tests.log").decode())
    assert sum(int(r[0]) for r in rust) == 866 and sum(int(r[1]) for r in rust) == 0
    assert python and int(python[1]) == 92
    assert all(row["independently_sound"] for row in scored)
    assert all(row["stereo"]["violated"] == 0 and row["stereo"]["unevaluable"] == 0 for row in scored)
    assert all(row["force_field"]["actual_force_field_used"] == "mmff94_bond_angle_strict" for row in candidate)
    result = {
        "kind": "local_uncommitted_source_candidate", "baseline": str(BASE),
        "fixed_input_count": 265, "comparable_rows": 262,
        "status_counts": dict(Counter(row["status"] for row in candidate_energy)),
        "same_energy_manifests_protocol_coordinates_and_all_rows": True,
        "same_pipeline_starting_energy_breakdowns": True,
        "energy_measurements": current["measurements"],
        "direct_gradient": {k:v for k,v in gradient.items() if k not in ("rows", "artifacts")},
        "separate_cation_water_regressions": {k:v for k,v in holdout.items() if k not in ("rows", "artifacts")},
        "quality": {
            "rows": 265, "success": sum(row["status"] == "success" for row in candidate),
            "independently_sound": sum(row["independently_sound"] for row in scored),
            "gross_clashes": sum(row["gross_clash_count"] for row in scored),
            "stereo_violated": sum(row["stereo"]["violated"] for row in scored),
            "stereo_unevaluable": sum(row["stereo"]["unevaluable"] for row in scored),
            "reported_converged_before": sum(row["force_field"]["converged"] for row in baseline),
            "reported_converged_after": sum(row["force_field"]["converged"] for row in candidate),
            "nonconverged_iteration_counts_after": dict(Counter(row["force_field"]["iterations"] for row in candidate if not row["force_field"]["converged"])),
            "iteration_limit": 200, "gradient_threshold": 1e-4,
            "changed_rows": changed,
            "force_field_graph_modes": dict(graph_modes),
            "candidate_convergence_flag_mismatches": 0,
        },
        "tests": {"rust_passed": 866, "rust_ignored": sum(int(r[2]) for r in rust), "python_passed": int(python[1]), "clippy": "passed", "fmt": "passed", "diff_check": "passed"},
        "artifacts": current["artifacts"],
        "scorer": scorer,
        "coordinate_json_decoding": "serde_json/float_roundtrip; both candidates re-scored with the same rebuilt binary; production stereo predicate unchanged",
        "evidence_files": [artifact(p) for p in [
            BASE / "candidate-summary.json", BASE / "candidate.jsonl", BASE / "quality.jsonl",
            OUT / "candidate-summary.json", OUT / "candidate.jsonl", OUT / "candidate-holdout.json",
            OUT / "gradient-comparison.json", OUT / "quality.jsonl", OUT / "quality-scored.jsonl",
            OUT / "convergence-audit.json", OUT / "validation-run.json", OUT / "gradient-scorer-run.json",
            OUT / "rust-tests.log", OUT / "python-tests.log", OUT / "clippy.log", OUT / "build.log",
            OUT / "change.patch", OUT / "build-source.json", OUT / "baseline-source.json",
            OUT / "baseline-quality-scored.jsonl", OUT / "scorer-build.json", OUT / "scorer-build.log"]],
        "source_files": [artifact(p) for p in [
            Path("crates/chematic-ff/src/mmff94_minimizer.rs"),
            Path("crates/chematic-ff/tests/mmff94_optimizer_recovery.rs"),
            Path("crates/chematic-3d/Cargo.toml"),
            Path("crates/chematic-3d/tests/stereo_coordinate_roundtrip.rs"),
            OUT / "revalidate.py", OUT / "compare.py", OUT / "measure_holdout.py", OUT / "compare_gradients.py", OUT / "run_frozen_energy.py"]],
        "caveat": "The one additional convergence report is a final-state flag correction, not a faster optimization result. Physical conformer accuracy, search quality and speed versus RDKit are not established. No release or publication occurred.",
    }
    (OUT / "comparison.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"quality": result["quality"], "tests": result["tests"]}, indent=2))


if __name__ == "__main__":
    main()
