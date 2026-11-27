#!/usr/bin/env python3
"""Recheck one exposed Se=Se square-planar development case after stereo support."""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import chematic
from rdkit import Chem, rdBase

SMILES = "Brc1ccc2c3c(cccc13)[Se@@](c1ccccc1)=[Se@SP3]2Br"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--previous", required=True)
    parser.add_argument("--wheel", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    previous = json.loads(Path(args.previous).read_text(encoding="utf-8"))["records"]
    before = next(r["checked"] for r in previous if r["smiles"] == SMILES and r["seed"] == 42)
    native = chematic.from_smiles(SMILES).add_hydrogens()
    source = Chem.AddHs(Chem.MolFromSmiles(SMILES))
    native_bonds = {tuple(sorted((b[0], b[1]))) for b in native.bond_table}
    ref_bonds = {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in source.GetBonds()}
    if [a[1] for a in native.atom_table] != [a.GetAtomicNum() for a in source.GetAtoms()] or native_bonds != ref_bonds:
        raise ValueError("indexed graph mismatch")
    expected = source.GetAtomWithIdx(11).GetProp("_CIPCode")
    try:
        candidate = native.rdkit_embed_square_planar_checked(random_seed=42, max_unit_triple_product=0.3)
        ref = Chem.Mol(source)
        conformer = Chem.Conformer(ref.GetNumAtoms())
        for index, coords in enumerate(candidate["coordinates"]):
            conformer.SetAtomPosition(index, coords)
        ref.AddConformer(conformer, assignId=True)
        Chem.AssignStereochemistryFrom3D(ref, replaceExistingTags=True)
        center = ref.GetAtomWithIdx(11)
        observed = center.GetProp("_CIPCode") if center.HasProp("_CIPCode") else None
        after = {"status": "returned", **candidate, "rdkit_3d_cip": observed}
    except (RuntimeError, ValueError) as error:
        after = {"status": "error", "error": str(error)}
    result = {"smiles": SMILES, "seed": 42, "planarity_limit": 0.3, "center": 11,
              "expected_cip": expected, "before": before, "after": after,
              "rdkit_version": rdBase.rdkitVersion, "previous_sha256": digest(args.previous),
              "wheel_sha256": digest(args.wheel),
              "extension_sha256": digest(importlib.util.find_spec("chematic.chematic").origin),
              "script_sha256": digest(__file__)}
    Path(args.output).write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
