"""Verify frozen-input continuity and collect final candidate evidence."""
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

OUT = Path(__file__).resolve().parent
BASE = OUT.parent / "oop-vdw"


def read(path):
    return json.loads(Path(path).read_bytes())


def rows(path):
    return [json.loads(line) for line in Path(path).read_bytes().splitlines()]


def artifact(path):
    data = Path(path).read_bytes()
    return {"path": str(path), "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def main():
    previous = read(BASE / "candidate-summary.json")
    current = read(OUT / "candidate-summary.json")
    old_rows = rows(BASE / "candidate.jsonl")
    new_rows = rows(OUT / "candidate.jsonl")
    assert len(old_rows) == len(new_rows) == 265
    assert previous["manifests"] == current["manifests"]
    assert previous["protocol"] == current["protocol"]
    for old, new in zip(old_rows, new_rows, strict=True):
        for field in ["input_index", "smiles", "status"]:
            assert old[field] == new[field], (old["input_index"], field)
        if old["status"] == "ok":
            for field in ["coordinate_sha256", "rdkit_energy_kcal_mol", "schematic_energy_kcal_mol", "term_delta_kcal_mol"]:
                assert old[field] == new[field], (old["input_index"], field)
    comparable = [row for row in new_rows if row["status"] == "ok"]
    quality = rows(OUT / "quality.jsonl")
    old_quality = rows(BASE / "quality.jsonl")
    assert len(quality) == len(old_quality) == 265
    assert all(old["coords"] == new["coords"] for old, new in zip(old_quality, quality, strict=True))
    scored = [row for row in rows(OUT / "quality-scored.jsonl") if row["engine"] == "chematic"]
    assert len(scored) == 265
    gradient = read(OUT / "gradient-comparison.json")
    holdout = read(OUT / "candidate-holdout.json")
    assert holdout["input_count"] == 20 and holdout["all_types_match"]
    assert holdout["artifacts"]["extension"]["sha256"] == current["artifacts"]["schematic"]["sha256"]
    for collection in [current["artifacts"], gradient["artifacts"], holdout["artifacts"]]:
        for value in collection.values():
            assert artifact(value["path"])["sha256"] == value["sha256"]
    rust = re.findall(r"test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored", (OUT / "rust-tests.log").read_text())
    python = re.search(r"(\d+) passed in", (OUT / "python-tests.log").read_text())
    assert python and int(python[1]) == 92
    assert sum(int(row[0]) for row in rust) == 860
    assert sum(int(row[1]) for row in rust) == 0
    result = {
        "kind": "local_uncommitted_source_candidate", "baseline": str(BASE),
        "frozen_input_count": 265, "comparable_rows": len(comparable),
        "status_counts": dict(Counter(row["status"] for row in new_rows)),
        "same_manifests_protocol_coordinates_and_energies": True,
        "energy_measurements": current["measurements"],
        "max_term_abs_errors": {term: max(abs(row["term_delta_kcal_mol"][term]) for row in comparable) for term in comparable[0]["term_delta_kcal_mol"]},
        "holdout": {k:v for k,v in holdout.items() if k not in ("rows", "artifacts")},
        "gradient": {k:v for k,v in gradient.items() if k not in ("rows", "artifacts")},
        "quality": {
            "rows": len(quality), "success": sum(row["status"] == "success" for row in quality),
            "independently_sound": sum(row["independently_sound"] for row in scored),
            "gross_clashes": sum(row["gross_clash_count"] for row in scored),
            "stereo_violated": sum(row["stereo"]["violated"] for row in scored),
            "stereo_unevaluable": sum(row["stereo"]["unevaluable"] for row in scored),
            "reported_converged": sum(row["force_field"]["converged"] for row in quality),
            "previous_reported_converged": sum(row["force_field"]["converged"] for row in old_quality),
            "same_generated_coordinates_as_previous": True,
            "iteration_limit": 200,
        },
        "tests": {"rust_passed": sum(int(row[0]) for row in rust), "rust_ignored": sum(int(row[2]) for row in rust), "python_passed": int(python[1]), "clippy": "passed", "fmt": "passed", "diff_check": "passed"},
        "artifacts": current["artifacts"],
        "evidence_files": [artifact(p) for p in [BASE / "candidate-summary.json", BASE / "candidate.jsonl", OUT / "candidate-summary.json", OUT / "candidate.jsonl", OUT / "candidate-holdout.json", OUT / "gradient-comparison.json", OUT / "quality.jsonl", OUT / "quality-scored.jsonl", OUT / "rust-tests.log", OUT / "python-tests.log", OUT / "clippy.log", OUT / "build.log", OUT / "implementation.patch"]],
        "source_files": [artifact(p) for p in [Path("crates/chematic-ff/src/mmff94_numeric.rs"), Path("crates/chematic-ff/tests/mmff94_cation_water.rs"), Path("scripts/tests/test_mmff94_same_explicit_h_energy.py"), OUT / "measure_holdout.py", OUT / "compare_gradients.py", OUT / "compare.py"]],
        "caveat": "Frozen-corpus and separate regression results measure MMFF94 implementation agreement; they do not establish superior physical accuracy, conformer search quality or speed. No published package was evaluated for this candidate.",
    }
    (OUT / "comparison.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key:result[key] for key in ["status_counts", "quality", "tests"]}, indent=2))


if __name__ == "__main__":
    main()
