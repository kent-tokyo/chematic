#!/usr/bin/env python3
"""Compare MMFF minimizers from the same RDKit ETKDG starting coordinates."""

import argparse
import json
from pathlib import Path

import chematic
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

from platinum_development_probe import graph_matches, rms_values, stereo_violations


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sdf", required=True, type=Path)
    parser.add_argument("--rows", required=True, type=Path)
    parser.add_argument("--name", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    measured = json.loads(args.rows.read_text())
    row = next(r for r in measured["results"] if r["name"] == args.name)
    reference = Chem.SDMolSupplier(str(args.sdf), removeHs=False)[row["index"]]
    ref_heavy = Chem.RemoveHs(reference)
    mol = Chem.AddHs(Chem.MolFromSmiles(row["smiles"]))
    other = chematic.from_smiles(row["smiles"]).add_hydrogens()
    assert graph_matches(other, mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.numThreads = 1
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=20, params=params))
    output = {"name": args.name, "rdkit": rdBase.rdkitVersion, "seed": 42, "max_iters": 200, "results": []}
    for cid in ids:
        initial = [list(mol.GetConformer(cid).GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
        field = other.mmff_force_field(initial)
        cstatus, cenergy = field.minimize_(max_iters=200)
        ccoords = field.positions
        rdmol = Chem.Mol(mol)
        props = AllChem.MMFFGetMoleculeProperties(rdmol, mmffVariant="MMFF94")
        rdfield = AllChem.MMFFGetMoleculeForceField(rdmol, props, confId=cid)
        rstatus = rdfield.Minimize(maxIts=200)
        renergy = rdfield.CalcEnergy()
        flat = rdfield.Positions()
        rcoords = [list(flat[3*i:3*i+3]) for i in range(mol.GetNumAtoms())]
        rms = rms_values(mol, ref_heavy, [initial, ccoords, rcoords])
        stereo = stereo_violations(mol, [ccoords, rcoords])
        output["results"].append({"cid": cid, "initial_rmsd": rms[0], "chematic_rmsd": rms[1], "rdkit_rmsd": rms[2], "chematic_status": cstatus, "rdkit_status": rstatus, "chematic_energy": cenergy, "rdkit_energy": renergy, "stereo": stereo})
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    print(json.dumps(output["results"], indent=2))


if __name__ == "__main__":
    main()
