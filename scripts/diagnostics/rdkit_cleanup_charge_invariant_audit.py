#!/usr/bin/env python3
"""Audit formal-charge preservation in a pinned Cleanup source comparison."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def charge(smiles: str | None) -> int | None:
    if smiles is None:
        return None
    mol = Chem.MolFromSmiles(smiles)
    return Chem.GetFormalCharge(mol) if mol is not None else None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dump", required=True, type=Path)
    parser.add_argument("--comparison", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    RDLogger.DisableLog("rdApp.*")

    comparison = json.loads(args.comparison.read_text())
    if rdBase.rdkitVersion != comparison["rdkit_version"]:
        parser.error(
            f"expected RDKit {comparison['rdkit_version']}; found {rdBase.rdkitVersion}"
        )
    dump_hash = sha256(args.dump)
    if dump_hash != comparison["source_dump_sha256"]:
        parser.error("source dump hash differs from the pinned comparison")
    rows = [json.loads(line) for line in args.dump.read_text().splitlines()]
    by_index = {row["index"]: row for row in rows}
    if len(rows) != comparison["source_rows"] or len(by_index) != len(rows):
        parser.error("source dump row count or indices differ from the comparison")
    excluded = {row["index"] for row in comparison["excluded"]}
    oracle_errors = {row["index"] for row in comparison["oracle_errors"]}
    mismatches = {row["index"]: row for row in comparison["mismatches"]}
    if len(mismatches) != len(comparison["mismatches"]):
        parser.error("comparison has duplicate mismatch indices")

    counts = {
        "both_preserve": 0,
        "source_only_preserves": 0,
        "rdkit_only_preserves": 0,
        "neither_preserves": 0,
    }
    asymmetric = []
    parser_failures = []
    for index, row in sorted(by_index.items()):
        if index in excluded or index in oracle_errors:
            continue
        source_smiles = row.get("cleanup", {}).get("Ok")
        recorded = mismatches.get(index)
        if recorded is not None and (
            recorded["input"] != row["input"]
            or recorded["chematic_cleanup"] != source_smiles
        ):
            parser.error(f"source mismatch at row {index}")
        # The frozen comparator records only differing RDKit outputs. On an
        # exact row its RDKit output is identical to the source output.
        rdkit_smiles = recorded["rdkit_cleanup"] if recorded else source_smiles
        charges = {
            "input": charge(row["input"]),
            "source": charge(source_smiles),
            "rdkit": charge(rdkit_smiles),
        }
        if any(value is None for value in charges.values()):
            parser_failures.append({"index": index, "charges": charges})
            continue
        source_preserves = charges["source"] == charges["input"]
        rdkit_preserves = charges["rdkit"] == charges["input"]
        category = (
            "both_preserve" if source_preserves and rdkit_preserves
            else "source_only_preserves" if source_preserves
            else "rdkit_only_preserves" if rdkit_preserves
            else "neither_preserves"
        )
        counts[category] += 1
        if source_preserves != rdkit_preserves:
            asymmetric.append({"index": index, "input": row["input"],
                               "charges": charges, "source_output": source_smiles,
                               "rdkit_output": rdkit_smiles})

    if sum(counts.values()) + len(parser_failures) != comparison["compatible_rows"]:
        parser.error("charge audit did not account for every comparable row")
    report = {
        "profile": "cleanup_formal_charge_preservation",
        "cleanup_oracle_version": comparison["rdkit_version"],
        "formal_charge_parser_version": rdBase.rdkitVersion,
        "comparison_sha256": sha256(args.comparison),
        "source_dump_sha256": dump_hash,
        "compatible_rows": comparison["compatible_rows"],
        "counts": counts,
        "parser_failures": parser_failures,
        "asymmetric_rows": asymmetric,
        "limits": [
            "selected stress corpus",
            "formal charge preservation alone does not prove chemical correctness",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(f"Charge audit: {counts}; parser failures: {len(parser_failures)}")


if __name__ == "__main__":
    main()
