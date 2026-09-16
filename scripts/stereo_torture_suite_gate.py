#!/usr/bin/env python3
"""Check semantic preservation and canonical reparse stability for T5.1."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path

from rdkit import Chem, rdBase


ROOT = Path(__file__).resolve().parents[1]
SUITE = ROOT / "validation" / "stereo_torture_suite_development.jsonl"
CLI = ROOT / "target" / "debug" / "chematic"


def identity(smiles: str) -> str | None:
    molecule = Chem.MolFromSmiles(smiles)
    return None if molecule is None else Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)


def run(cli: Path, smiles: list[str], directory: Path) -> list[dict[str, object]]:
    input_path = directory / "input.smi"
    input_path.write_text("\n".join(smiles) + "\n", encoding="utf-8")
    completed = subprocess.run([str(cli), "batch-canonicalize", "--input", str(input_path)], cwd=ROOT, text=True, capture_output=True)
    if completed.returncode:
        raise RuntimeError(completed.stderr.strip() or "batch-canonicalize failed")
    return json.loads(completed.stdout)["records"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, default=CLI)
    parser.add_argument("--suite", type=Path, default=SUITE)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cli, suite = args.cli.resolve(), args.suite.resolve()
    rows = [json.loads(line) for line in suite.read_text(encoding="utf-8").splitlines()]
    manifest, cases = rows[0], rows[1:]
    if not manifest.get("_manifest") or len(cases) != 300:
        raise ValueError("expected a 300-case suite with a manifest")
    with tempfile.TemporaryDirectory(prefix="chematic-stereo-suite-") as directory:
        first = run(cli, [row["smiles"] for row in cases], Path(directory))
        second = run(cli, [str(row.get("canonical_smiles", "")) for row in first], Path(directory))
    failures = []
    for case, output, reparsed in zip(cases, first, second, strict=True):
        expected = identity(case["smiles"])
        actual = identity(str(output.get("canonical_smiles", ""))) if output.get("error") is None else None
        stable = output.get("error") is None and reparsed.get("error") is None and output.get("canonical_smiles") == reparsed.get("canonical_smiles")
        if expected != actual or not stable:
            failures.append({"id": case["id"], "semantic_equal": expected == actual, "canonical_reparse_stable": stable, "error": output.get("error") or reparsed.get("error")})
    result = {
        "schema_version": 1, "suite_sha256": hashlib.sha256(suite.read_bytes()).hexdigest(),
        "suite_cases": len(cases), "rdkit_version": rdBase.rdkitVersion,
        "semantic_equal": len(cases) - sum(not row["semantic_equal"] for row in failures),
        "canonical_reparse_stable": len(cases) - sum(not row["canonical_reparse_stable"] for row in failures),
        "failure_count": len(failures), "failures": failures, "gate_passed": not failures,
        "scope": "development suite only; RDKit semantic identity and chematic canonical reparse stability",
        "not_claimed": ["independent accuracy", "CIP label adjudication", "sealed challenge performance"],
    }
    output = args.output if args.output.is_absolute() else ROOT / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in ("suite_cases", "semantic_equal", "canonical_reparse_stable", "failure_count", "gate_passed")}, sort_keys=True))
    return 0 if result["gate_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
