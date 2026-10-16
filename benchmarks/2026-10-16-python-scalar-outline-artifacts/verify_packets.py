#!/usr/bin/env python3
"""Audit the four frozen packets without producing new timing samples."""

from pathlib import Path
import gzip
import hashlib
import json
import re
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
STEMS = [
    "2026-10-16-python-aromatic-edge-proof",
    "2026-10-16-python-aromatic-normalization",
    "2026-10-16-python-scalar-outline",
    "2026-10-16-python-scalar-property-cache",
]


def digest(data):
    return hashlib.sha256(data).hexdigest()


def read_json(path):
    return json.loads(path.read_text())


def old_file(base, name):
    p = subprocess.run(
        ["git", "show", base + ":" + name], cwd=ROOT, capture_output=True
    )
    return p.stdout if p.returncode == 0 else None


def main():
    packets = []
    for number, stem in enumerate(STEMS, 1):
        folder = ROOT / "benchmarks" / (stem + "-artifacts")
        manifest = read_json(folder / "source-manifest.json")
        decision = read_json(folder / "decision.json")
        base = manifest["base_revision"]
        patch = gzip.decompress((folder / "candidate.patch.gz").read_bytes())
        assert digest(patch) == manifest["patch_sha256"]
        for name, expected in manifest["baseline_source_sha256"].items():
            assert digest(old_file(base, name)) == expected, name
        with tempfile.TemporaryDirectory(prefix="chematic-packet-audit-") as tmp:
            tree = Path(tmp)
            for name in manifest["candidate_source_sha256"]:
                data = old_file(base, name)
                if data is not None:
                    p = tree / name
                    p.parent.mkdir(parents=True, exist_ok=True)
                    p.write_bytes(data)
            subprocess.run(
                ["git", "apply", "--unsafe-paths", "-"],
                cwd=tree, input=patch, check=True,
            )
            for name, expected in manifest["candidate_source_sha256"].items():
                assert digest((tree / name).read_bytes()) == expected, name
        if number == 3:
            for name, expected in manifest["candidate_source_sha256"].items():
                assert digest((ROOT / name).read_bytes()) == expected, name
        parity = read_json(folder / "parity.json")
        assert parity["baseline"] == parity["candidate"]
        assert sum(parity["counts"].values()) == 197008
        assert all(v["rows"] == 14072 and v["errors"] == 0
                   for v in parity["candidate"].values())
        assert parity["candidate"]["descriptor"]["nonfinite_float_outcomes"] == 6
        cold_pass = True
        prepared_pass = True
        floor_pass = True
        negative_count = 0
        for offset in ("0", "5000"):
            data = json.loads(gzip.decompress(
                (folder / f"frontier-{offset}.json.gz").read_bytes()
            ))
            assert data["complete"] and len(data["blocks"]) == 20
            assert len(data["summary"]) == 21
            below_floor = {}
            negatives = {}
            for key, lanes in data["summary"].items():
                ratio = lanes["baseline"]["speedup"]
                assert ratio["blocks"] == 20
                if ratio["median"] < .99:
                    below_floor[key] = ratio
                if ratio["ci95"][1] < 1:
                    negatives[key] = ratio
                for engine in ("rdkit", "cosmolkit"):
                    lane = lanes[engine]
                    if lane["demonstrated_win"]:
                        assert lane["speedup"]["ci95"][0] > 1
                        assert lane["speedup"]["winning_blocks"] == 20
            assert below_floor == decision["source_lanes_below_floor"][offset]
            assert negatives == decision["all_source_negative_intervals"][offset]
            cold_pass &= all(data["summary"][k]["baseline"]["speedup"]["ci95"][0] > 1
                             for k in ("tpsa/first_use", "labute_asa/first_use"))
            prepared_pass &= data["summary"]["tpsa/repeated_use"]["baseline"]["speedup"]["ci95"][0] > 1
            floor_pass &= not below_floor
            negative_count += len(negatives)
        gate = cold_pass and floor_pass and (prepared_pass if number >= 3 else True)
        assert gate == decision["local_adoption_gate_passed"]
        assert not gate
        report = ROOT / "benchmarks" / (stem + ".md")
        for target in re.findall(r"\[[^\]]+\]\(([^)]+)\)", report.read_text()):
            if target.startswith(("http://", "https://", "mailto:", "#")):
                continue
            assert (report.parent / target.split("#")[0]).is_file(), target
        packets.append({"candidate": number, "source_reconstruction": True,
                        "exact_output_outcomes": 197008, "preserved_nonfinite": 6,
                        "complete_profiles": 2, "blocks_per_profile": 20,
                        "source_lanes": 42, "cold_gate": bool(cold_pass),
                        "floor_gate": bool(floor_pass), "prepared_tpsa_gate": bool(prepared_pass),
                        "negative_intervals": negative_count, "adoption_gate": bool(gate)})
    selection = read_json(ROOT / "benchmarks" / (STEMS[2] + "-artifacts") /
                          "final-source-selection.json")
    for name, expected in selection["python_binding_restored_from_base_sha256"].items():
        assert digest((ROOT / name).read_bytes()) == expected
    assert not (ROOT / "crates/chematic-py/tests/test_scalar_property_cache.py").exists()
    print(json.dumps({"passed": True, "retained_candidate": 3, "packets": packets}, indent=2))


if __name__ == "__main__":
    main()
