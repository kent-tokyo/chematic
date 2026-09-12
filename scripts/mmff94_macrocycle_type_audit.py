#!/usr/bin/env python3
"""Compare RDKit and schematic MMFF numeric atom types on macrocycle fixtures."""

from __future__ import annotations

import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

import chematic


ROOT = Path(__file__).resolve().parents[1]
MANIFESTS = (
    ROOT / "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_a.json",
    ROOT / "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_b.json",
)


def main() -> int:
    rows = []
    for path in MANIFESTS:
        rows.extend(json.loads(path.read_text(encoding="utf-8"))["molecules"])
    for row in rows:
        if row.get("primary_category") != "macrocycle":
            continue
        smiles = row["smiles"]
        rdkit = Chem.MolFromSmiles(smiles)
        AllChem.MMFFSanitizeMolecule(rdkit)
        props = AllChem.MMFFGetMoleculeProperties(rdkit)
        schematic = chematic.from_smiles(smiles)
        schematic_types = schematic.mmff94_numeric_atom_types()
        mismatches = [
            {
                "atom": index,
                "rdkit": props.GetMMFFAtomType(index),
                "chematic": schematic_types[index],
            }
            for index in range(len(schematic_types))
            if props.GetMMFFAtomType(index) != schematic_types[index]
        ]
        print(json.dumps({"name": row["name"], "mismatches": mismatches}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
