#!/usr/bin/env python3
"""Recheck stored square-planar angles with RDKit at identical coordinates."""

import argparse
import hashlib
import json
from pathlib import Path

import chematic
from rdkit import Chem, rdBase
from rdkit.Chem import rdMolTransforms


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    source = json.loads(Path(args.input).read_text(encoding="utf-8"))
    checks = []
    for row in source["records"]:
        if row["retry"]["status"] != "returned":
            continue
        mol = chematic.from_smiles(row["smiles"]).add_hydrogens()
        reference = Chem.AddHs(Chem.MolFromSmiles(row["smiles"]))
        atoms_match = [atom[1] for atom in mol.atom_table] == [a.GetAtomicNum() for a in reference.GetAtoms()]
        bonds = {tuple(sorted((bond[0], bond[1]))) for bond in mol.bond_table}
        ref_bonds = {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in reference.GetBonds()}
        if not atoms_match or bonds != ref_bonds:
            raise ValueError(f"indexed graph mismatch for {row['smiles']}")
        conf = Chem.Conformer(reference.GetNumAtoms())
        for index, point in enumerate(row["retry"]["coordinates"]):
            conf.SetAtomPosition(index, point)
        for center in row["retry"]["angles"]:
            for category in ("trans", "cis"):
                for a, b, native in center[category]:
                    oracle = rdMolTransforms.GetAngleDeg(conf, a, center["center"], b)
                    checks.append({"group": row["group"], "smiles": row["smiles"], "seed": row["seed"],
                                   "center": center["center"], "atoms": [a, b], "category": category,
                                   "native_degrees": native, "rdkit_degrees": oracle,
                                   "absolute_error": abs(native - oracle)})
    max_error = max((check["absolute_error"] for check in checks), default=0.0)
    result = {"input_sha256": digest(args.input), "script_sha256": digest(__file__),
              "rdkit_version": rdBase.rdkitVersion, "comparison": "same coordinates and indexed graph",
              "max_absolute_error_degrees": max_error, "checks": checks}
    Path(args.output).write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    if not checks or max_error > 1e-10:
        raise SystemExit("missing checks or angle tolerance exceeded")


if __name__ == "__main__":
    main()
