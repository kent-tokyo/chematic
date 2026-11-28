#!/usr/bin/env python3
"""Reassign the selenium CIP center in a fixed square-planar development result."""

import argparse
import hashlib
import json
from pathlib import Path

import chematic
from rdkit import Chem, rdBase

SMILES = "Brc1ccc2c3c(cccc13)[Se@@](c1ccccc1)=[Se@SP3]2Br"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    records = json.loads(Path(args.input).read_text(encoding="utf-8"))["records"]
    row = next(r for r in records if r["smiles"] == SMILES and r["seed"] == 42)
    coords = row["retry"]["coordinates"]
    source = Chem.AddHs(Chem.MolFromSmiles(SMILES))
    native = chematic.from_smiles(SMILES).add_hydrogens()
    if ([a.GetAtomicNum() for a in source.GetAtoms()] != [a[1] for a in native.atom_table]
            or {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in source.GetBonds()}
            != {tuple(sorted((b[0], b[1]))) for b in native.bond_table}):
        raise ValueError("indexed graph mismatch")
    center = 11
    expected = source.GetAtomWithIdx(center).GetProp("_CIPCode")
    observed = []
    for mirror in (False, True):
        ref = Chem.Mol(source)
        conformer = Chem.Conformer(ref.GetNumAtoms())
        for i, (x, y, z) in enumerate(coords):
            conformer.SetAtomPosition(i, [-x if mirror else x, y, z])
        ref.AddConformer(conformer, assignId=True)
        Chem.AssignStereochemistryFrom3D(ref, replaceExistingTags=True)
        atom = ref.GetAtomWithIdx(center)
        observed.append(atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None)
    if observed[0] != expected or observed[1] is None or observed[1] == expected:
        raise ValueError(f"RDKit could not distinguish original and mirror: {expected} {observed}")
    result = {"smiles": SMILES, "center": center, "coords": coords,
              "expected_cip": expected, "rdkit_3d_cip": observed[0], "rdkit_mirror_cip": observed[1],
              "rdkit_version": rdBase.rdkitVersion, "input_sha256": digest(args.input),
              "script_sha256": digest(__file__)}
    Path(args.output).write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
