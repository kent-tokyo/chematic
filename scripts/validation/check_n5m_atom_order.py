#!/usr/bin/env python3
"""Check frozen RDKit atom types and partial charges on permuted SMILES."""
import json
from pathlib import Path

import chematic


def main():
    fixture = json.loads((Path(__file__).resolve().parents[1] / "fixtures" / "mmff94-n5m-atom-order-2026-11-22.json").read_text(encoding="utf-8"))
    if len(fixture["molecules"]) != 30:
        raise ValueError("unexpected fixture size")
    for row in fixture["molecules"]:
        mol = chematic.from_smiles(row["smiles"]).add_hydrogens()
        types = list(mol.mmff94_numeric_atom_types())
        charges = mol.mmff94_charges_numeric()
        if types != row["types"]:
            raise AssertionError(f"{row['name']}: atom type mismatch")
        if len(charges) != len(row["charges"]):
            raise AssertionError(f"{row['name']}: charge length mismatch")
        for index, (actual, expected) in enumerate(zip(charges, row["charges"])):
            if abs(actual - expected) >= 1e-10:
                raise AssertionError(f"{row['name']} atom {index}: {actual} != {expected}")
    print(f"{len(fixture['molecules'])} exposed atom-order variants matched the RDKit oracle")


if __name__ == "__main__":
    main()
