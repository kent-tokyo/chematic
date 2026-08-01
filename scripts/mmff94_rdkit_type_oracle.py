#!/usr/bin/env python3
"""Phase 1A audit for issue #227 (MMFF94 strict coverage gap).

RDKit oracle: for every molecule in the 265-molecule Wave 1 corpus
(validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_{a,b}.json), records
RDKit's real MMFF94 atom type per heavy atom, MMFF properties availability,
force-field construction success, single-point energy, and Minimize() return
code. Atom order is preserved as RDKit parses it (Chem.MolFromSmiles does not
reorder/canonicalize atoms by default) -- index-aligned with chematic's own
dump (mmff94_numeric_type_dump.rs) per PR #226's already-verified 265/265
heavy-atom element-sequence atom mapping between the two engines.

Run: .venv/bin/python scripts/mmff94_rdkit_type_oracle.py \
  > validation/results/mmff94_rdkit_type_oracle.jsonl
"""

import json
import sys

from rdkit import Chem
from rdkit.Chem import AllChem
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

    n_heavy = mol.GetNumAtoms()
    molH = Chem.AddHs(mol)
    seed = 20260801
    embed_ret = AllChem.EmbedMolecule(molH, randomSeed=seed, useRandomCoords=True)
    if embed_ret != 0:
        # Retry with a different seed before giving up -- ETKDG can fail
        # stochastically on some seeds without meaning MMFF typing is broken.
        embed_ret = AllChem.EmbedMolecule(molH, randomSeed=seed + 1, useRandomCoords=True)

    props = AllChem.MMFFGetMoleculeProperties(molH) if embed_ret == 0 else None
    mmff_properties_available = props is not None

    atom_types = []
    if props is not None:
        for atom in molH.GetAtoms():
            if atom.GetIdx() < n_heavy:  # heavy atoms are always indexed first by AddHs
                atom_types.append(
                    {
                        "index": atom.GetIdx(),
                        "element": atom.GetSymbol(),
                        "aromatic": atom.GetIsAromatic(),
                        "rdkit_mmff_type": props.GetMMFFAtomType(atom.GetIdx()),
                    }
                )

    ff_construction_ok = False
    energy = None
    optimize_return_code = None
    if props is not None:
        try:
            ff = AllChem.MMFFGetMoleculeForceField(molH, props)
            ff_construction_ok = ff is not None
            if ff is not None:
                energy = ff.CalcEnergy()
                optimize_return_code = ff.Minimize(maxIts=200)
        except Exception as e:  # noqa: BLE001 -- typed failure, recorded not swallowed
            ff_construction_ok = False
            energy = None
            optimize_return_code = f"exception:{e}"

    return {
        "name": name,
        "smiles": smiles,
        "status": "ok",
        "n_heavy_atoms": n_heavy,
        "embed_return_code": embed_ret,
        "mmff_properties_available": mmff_properties_available,
        "atom_types": atom_types,
        "ff_construction_ok": ff_construction_ok,
        "energy_kcal_mol": energy,
        "optimize_return_code": optimize_return_code,
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
