#!/usr/bin/env python3
"""Reproducible core-descriptor parity lane on a large SMILES corpus.

This deliberately measures the descriptors exposed in the shared binding
contract. The broader descriptor census remains a separate, longer lane.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("smiles_csv", type=Path)
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()

    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
    import chematic

    with args.smiles_csv.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    column = "SMILES" if rows and "SMILES" in rows[0] else (next(iter(rows[0])) if rows else "")

    fields = {
        "molecular_weight": {"tolerance": 1e-6, "published_tolerance": 0.01, "matches": 0, "mismatches": 0, "published_matches": 0},
        "tpsa": {"tolerance": 1e-6, "published_tolerance": 0.01, "matches": 0, "mismatches": 0, "published_matches": 0},
        "hbd": {"tolerance": 0, "matches": 0, "mismatches": 0, "published_matches": 0},
        "hba": {"tolerance": 0, "matches": 0, "mismatches": 0, "published_matches": 0},
        "heavy_atoms": {"tolerance": 0, "matches": 0, "mismatches": 0, "published_matches": 0},
    }
    parsed = 0
    parse_failures = 0
    examples: list[dict[str, object]] = []

    for row in rows:
        smiles = row[column].strip()
        if not smiles:
            continue
        rd_mol = Chem.MolFromSmiles(smiles)
        try:
            ch_mol = chematic.from_smiles(smiles)
        except (TypeError, ValueError, RuntimeError):
            ch_mol = None
        if rd_mol is None or ch_mol is None:
            parse_failures += 1
            continue
        parsed += 1
        expected = {
            "molecular_weight": Descriptors.MolWt(rd_mol),
            "tpsa": rdMolDescriptors.CalcTPSA(rd_mol, includeSandP=True),
            "hbd": rdMolDescriptors.CalcNumHBD(rd_mol),
            "hba": rdMolDescriptors.CalcNumHBA(rd_mol),
            "heavy_atoms": rd_mol.GetNumHeavyAtoms(),
        }
        actual = {
            "molecular_weight": ch_mol.mw,
            "tpsa": ch_mol.tpsa,
            "hbd": ch_mol.hbd,
            "hba": ch_mol.hba,
            "heavy_atoms": ch_mol.heavy_atoms,
        }
        row_mismatches = {}
        for name, stats in fields.items():
            delta = abs(float(actual[name]) - float(expected[name]))
            if delta <= stats["tolerance"]:
                stats["matches"] += 1
            else:
                stats["mismatches"] += 1
                row_mismatches[name] = {"chematic": actual[name], "rdkit": expected[name], "delta": delta}
            if delta <= stats.get("published_tolerance", stats["tolerance"]):
                stats["published_matches"] += 1
        if row_mismatches and len(examples) < 25:
            examples.append({"smiles": smiles, "fields": row_mismatches})

    result = {
        "schema_version": 1,
        "corpus": str(args.smiles_csv),
        "rows": len(rows),
        "parsed": parsed,
        "parse_failures": parse_failures,
        "fields": fields,
        "mismatch_examples": examples,
        "comparison_boundary": "same SMILES inputs; chematic Python binding versus RDKit Python API",
    }
    args.json.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
