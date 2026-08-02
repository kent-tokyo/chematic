#!/usr/bin/env python3
"""Issue #227 Phase 0.3: RDKit oracle for the real hybridization gate that
`compute_mmff94_aromatic_view` approximates as `total_degree(atom) > 3`.

RDKit's `setMMFFAromaticity` (`Code/GraphMol/Aromaticity.cpp` line 1023,
pinned commit -- see `scripts/mmff94_provenance/PROVENANCE.md`) rejects a
ring C/N atom's ring from aromaticity when
`atom->getHybridization() != Atom::SP2`. `GetHybridization()` is RDKit's
general valence-model hybridization perception, set during standard
`Chem.MolFromSmiles` sanitization -- not itself MMFF-specific, so it is
queried directly on the plain (non-Kekulized-copy) molecule.

For every ring C/N heavy atom in the 265-molecule Wave 1 corpus, records
`element`, `total_degree` (matches chematic's own definition:
GetTotalDegree(), heavy-atom bonds + implicit/explicit Hs), and whether
RDKit's real hybridization gate rejects the atom
(`hybridization != SP2`).

Run: .venv/bin/python scripts/mmff94_hybridization_gate_gap_227.py \
  > validation/results/mmff94_hybridization_gate_gap_227_rdkit.jsonl
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

    ri = mol.GetRingInfo()
    ring_atoms = {a for ring in ri.AtomRings() for a in ring}

    atoms = []
    for idx in sorted(ring_atoms):
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetSymbol() not in ("C", "N"):
            continue
        atoms.append(
            {
                "index": idx,
                "element": atom.GetSymbol(),
                "total_degree": atom.GetTotalDegree(),
                "hybridization": str(atom.GetHybridization()),
                "gate_fires_reject": atom.GetHybridization() != Chem.HybridizationType.SP2,
            }
        )

    return {"name": name, "smiles": smiles, "status": "ok", "atoms": atoms}


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
