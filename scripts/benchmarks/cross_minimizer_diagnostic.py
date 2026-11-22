#!/usr/bin/env python3
"""Compare both MMFF minimizers from each engine's ETKDG starting coordinates."""

import argparse
import json
import subprocess
from pathlib import Path

import chematic
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

from platinum_development_probe import graph_matches, rms_values
from platinum_stereo_strict import stereo_violations_strict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--sdf", type=Path, required=True)
    parser.add_argument("--rows", type=Path, required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--disable-macrocycle-torsions", action="store_true")
    parser.add_argument("--disable-macrocycle-14", action="store_true")
    args = parser.parse_args()
    row = next(r for r in json.loads(args.rows.read_text())["results"] if r["name"] == args.name)
    reference = Chem.SDMolSupplier(str(args.sdf), removeHs=False)[row["index"]]
    heavy_reference = Chem.RemoveHs(reference)
    mol = Chem.AddHs(Chem.MolFromSmiles(row["smiles"]))
    other = chematic.from_smiles(row["smiles"]).add_hydrogens()
    if not graph_matches(other, mol):
        raise ValueError("atom-indexed molecular graphs differ")
    result = subprocess.run(
        [str(args.binary.resolve())], input=json.dumps({
            "smiles": row["smiles"],
            "macrocycle_torsions": not args.disable_macrocycle_torsions,
            "macrocycle_14": not args.disable_macrocycle_14,
        }) + "\n",
        text=True, capture_output=True, check=True,
    )
    source = json.loads(result.stdout)
    expected_atoms = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    expected_bonds = sorted(
        (min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
         str(b.GetBondType()).upper()) for b in mol.GetBonds()
    )
    actual_bonds = sorted((a, b, "SINGLE" if order in ("UP", "DOWN") else order)
                          for a, b, order in source["bonds"])
    if source["atomic_numbers"] != expected_atoms or actual_bonds != expected_bonds:
        raise ValueError("diagnostic binary produced a different indexed graph")
    generated = Chem.Mol(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.numThreads = 1
    ids = AllChem.EmbedMultipleConfs(generated, numConfs=20, params=params)
    rdkit_coords = [
        [list(generated.GetConformer(cid).GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
        for cid in ids
    ]
    output = {
        "name": args.name, "rdkit": rdBase.rdkitVersion, "seed": 42,
        "macrocycle_torsions": not args.disable_macrocycle_torsions,
        "macrocycle_14": not args.disable_macrocycle_14, "sources": {},
    }
    for label, coordinates in (("chematic", source["coords"]), ("rdkit", rdkit_coords)):
        records = []
        for initial in coordinates:
            cfield = other.mmff_force_field(initial)
            cstatus, cenergy = cfield.minimize_(max_iters=200)
            ccoords = cfield.positions
            rdmol = Chem.Mol(mol)
            rdmol.RemoveAllConformers()
            conformer = Chem.Conformer(mol.GetNumAtoms())
            for index, point in enumerate(initial):
                conformer.SetAtomPosition(index, point)
            rdmol.AddConformer(conformer)
            props = AllChem.MMFFGetMoleculeProperties(rdmol, mmffVariant="MMFF94")
            rfield = AllChem.MMFFGetMoleculeForceField(rdmol, props)
            rstatus = rfield.Minimize(maxIts=200)
            renergy = rfield.CalcEnergy()
            flat = rfield.Positions()
            rcoords = [list(flat[3 * i:3 * i + 3]) for i in range(mol.GetNumAtoms())]
            rms = rms_values(mol, heavy_reference, [initial, ccoords, rcoords])
            records.append({
                "initial_rmsd": rms[0], "chematic_rmsd": rms[1], "rdkit_rmsd": rms[2],
                "chematic_status": cstatus, "rdkit_status": rstatus,
                "chematic_energy": cenergy, "rdkit_energy": renergy,
                "stereo": stereo_violations_strict(mol, [initial, ccoords, rcoords]),
            })
        output["sources"][label] = records
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
