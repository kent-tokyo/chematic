#!/usr/bin/env python3
"""Same-input Open Babel #3000 check; original HIN fixture, no timing claims."""
import argparse
import hashlib
import importlib
import json
from pathlib import Path
import subprocess
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bond_types(block):
    rows = block.splitlines()
    atoms, bonds = int(rows[3][:3]), int(rows[3][3:6])
    return [int(row[6:9]) for row in rows[4 + atoms:4 + atoms + bonds]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--obabel", default="/opt/homebrew/bin/obabel")
    parser.add_argument("--output", type=Path, default=ROOT / "validation/results/openbabel-hin-aromatic-2026-10-15.json")
    args = parser.parse_args()
    import chematic
    fixture = ROOT / "validation/openbabel_issues/hin/aromatic-ring.hin"
    source = fixture.read_text(encoding="utf-8")
    record_json = chematic.parse_hin_record_json(source)
    candidate = chematic.from_hin(source).to_mol_block()
    command = [args.obabel, "-ihin", "-", "-osdf"]
    oracle = subprocess.run(command, input=source, text=True, capture_output=True, timeout=20, check=True)
    observed = bond_types(oracle.stdout)
    expected = bond_types(candidate)
    assert expected == [4] * 6, expected
    assert len(observed) == 6, observed
    assert chematic.roundtrip_hin_record(source) == source
    assert json.loads(record_json)["atoms"][0]["partial_charge"] == -0.75
    artifacts = args.output.with_suffix(".artifacts")
    artifacts.mkdir(parents=True, exist_ok=True)
    payloads = {"openbabel.sdf": oracle.stdout, "openbabel.stderr.txt": oracle.stderr,
                "chematic.mol": candidate, "chematic-record.json": record_json}
    artifact_hashes = {}
    for name, payload in payloads.items():
        path = artifacts / name
        path.write_text(payload, encoding="utf-8")
        artifact_hashes[name] = sha(path)
    sources = ["crates/chematic-mol/src/hin.rs", "crates/chematic-mol/src/hin/reader.rs",
               "crates/chematic-mol/src/hin/json.rs", "crates/chematic-py/src/formats/hin.rs",
               "crates/chematic-wasm/src/hin.rs"]
    report = {
        "schema_version": 1, "checked_at": datetime.now(timezone.utc).isoformat(),
        "scope": "Original synthetic aromatic HIN to MDL bond-type preservation; no general chemistry or speed claim",
        "issue": {"url": "https://github.com/openbabel/openbabel/issues/3000", "state_checked": "open"},
        "candidate": {"base_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                      "source_candidate": True, "published_artifact": False,
                      "source_sha256": {name: sha(ROOT / name) for name in sources},
                      "extension_sha256": sha(Path(importlib.import_module("chematic.chematic").__file__)),
                      "mdl_bond_types": expected, "exact_source_retention": True},
        "oracle": {"command": command, "version": subprocess.check_output([args.obabel, "-V"], text=True).strip(),
                   "binary_sha256": sha(Path(args.obabel).resolve()), "exit_code": oracle.returncode,
                   "mdl_bond_types": observed, "reported_defect_reproduced": observed == [5] * 6},
        "fixture": {"path": str(fixture.relative_to(ROOT)), "sha256": sha(fixture), "license": "MIT OR Apache-2.0"},
        "artifact_sha256": artifact_hashes,
    }
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"candidate_bond_types": expected, "openbabel_bond_types": observed, "report": str(args.output)}))


if __name__ == "__main__":
    main()
