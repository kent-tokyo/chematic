#!/usr/bin/env python3
"""Issue #227 Phase 0.4: RDKit oracle for MMFF-specific atom + bond
aromaticity, corpus-wide (the same 265-molecule Wave 1 corpus used
elsewhere), not just the 12 curated fixtures in
`mmff94_aromaticity_bond_parity_227.py`.

Ground truth: `Chem.SetAromaticity(mol, Chem.AROMATICITY_MMFF94)` on a
freshly-Kekulized copy -- RDKit's real `setMMFFAromaticity`.

Run: .venv/bin/python scripts/mmff94_aromaticity_corpus_parity_227.py \
  > validation/results/mmff94_aromaticity_corpus_parity_227_rdkit.jsonl
"""

import json
import sys

from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

MANIFESTS = [
    ("A", "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_a.json"),
    ("B", "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_b.json"),
]


def process(name, smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"name": name, "smiles": smiles, "status": "parse_failure"}

    amol = Chem.Mol(mol)
    try:
        Chem.Kekulize(amol, clearAromaticFlags=True)
    except Exception as e:  # noqa: BLE001 -- typed failure, recorded not swallowed
        return {"name": name, "smiles": smiles, "status": "kekulize_failure", "error": str(e)}
    Chem.SetAromaticity(amol, Chem.AROMATICITY_MMFF94)

    atom_aromatic = {a.GetIdx(): a.GetIsAromatic() for a in amol.GetAtoms()}
    bond_aromatic = {}
    for b in amol.GetBonds():
        key = f"{min(b.GetBeginAtomIdx(), b.GetEndAtomIdx())}-{max(b.GetBeginAtomIdx(), b.GetEndAtomIdx())}"
        bond_aromatic[key] = b.GetIsAromatic()

    return {
        "name": name,
        "smiles": smiles,
        "status": "ok",
        "n_heavy_atoms": mol.GetNumAtoms(),
        "atom_aromatic": atom_aromatic,
        "bond_aromatic": bond_aromatic,
    }


def main():
    for tier, path in MANIFESTS:
        with open(path) as f:
            manifest = json.load(f)
        for m in manifest["molecules"]:
            row = process(m["name"], m["smiles"])
            row["tier"] = tier
            print(json.dumps(row))
            sys.stdout.flush()


if __name__ == "__main__":
    main()
