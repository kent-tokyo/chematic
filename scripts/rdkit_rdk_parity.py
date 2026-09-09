#!/usr/bin/env python3
"""Compare v1.0.8 Rust ``rdkit_rdk`` fingerprints with RDKit's RDKFingerprint."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_chematic(corpus: Path) -> dict:
    command = [
        "cargo", "run", "-p", "chematic-cli", "--release", "--",
        "batch-fingerprints", "--input", str(corpus), "--algorithm", "rdkit_rdk",
    ]
    result = subprocess.run(command, cwd=ROOT, check=True, capture_output=True, text=True)
    return json.loads(result.stdout)


def rdkit_bits(smiles: str) -> list[int] | None:
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return sorted(Chem.RDKFingerprint(mol).GetOnBits())


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("corpus", type=Path)
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()

    batch = run_chematic(args.corpus)
    records = batch["records"]
    matches = 0
    mismatches = 0
    chematic_failures = 0
    rdkit_failures = 0
    examples: list[dict[str, object]] = []
    for record in records:
        smiles = record["input_smiles"]
        if record.get("error") is not None:
            chematic_failures += 1
            continue
        expected = rdkit_bits(smiles)
        if expected is None:
            rdkit_failures += 1
            continue
        actual = record["fingerprint"]["set_bits"]
        if actual == expected:
            matches += 1
        else:
            mismatches += 1
            if len(examples) < 25:
                examples.append({
                    "smiles": smiles,
                    "chematic_set_bits": actual,
                    "rdkit_set_bits": expected,
                    "extra": sorted(set(actual) - set(expected)),
                    "missing": sorted(set(expected) - set(actual)),
                })

    report = {
        "schema_version": 1,
        "target_version": "1.0.8",
        "corpus": str(args.corpus),
        "corpus_sha256": sha256(args.corpus),
        "operation": "rdkit_rdk",
        "comparison_boundary": "same SMILES inputs; v1.0.8 Rust CLI versus RDKit Chem.RDKFingerprint",
        "configuration": {
            "minPath": 1,
            "maxPath": 7,
            "useHs": True,
            "branchedPaths": True,
            "useBondOrder": True,
            "fpSize": 2048,
            "numBitsPerFeature": 2,
        },
        "rows": len(records),
        "chematic_valid": batch["valid_count"],
        "chematic_failures": chematic_failures,
        "rdkit_failures": rdkit_failures,
        "exact_matches": matches,
        "exact_mismatches": mismatches,
        "exact_match_pct_of_common_success": (100.0 * matches / (matches + mismatches)) if matches + mismatches else None,
        "rdkit_version": __import__("rdkit").rdBase.rdkitVersion,
        "mismatch_examples": examples,
    }
    args.json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    raise SystemExit(main())
