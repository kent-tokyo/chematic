#!/usr/bin/env python3
"""Compare declared selenium stereo to RDKit's 3D reassignment."""

import argparse
import hashlib
import importlib.util
import json
import platform
import sys
from pathlib import Path

import chematic
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

DEVELOPMENT = ["C[Se{}](=O)F", "F[Se{}](C)=O", "O=[Se{}](F)C",
               "C[Se{}](=N)F", "C[Se{}](=S)F", "C[Se{}](=[Se])F"]
CONTROLS = ["C[Se{}](=O)Cl", "Br[Se{}](C)=O", "N=[Se{}](Cl)C",
            "CC[Se{}](=S)Cl", "Cl[Se{}](C)=[Se]", "S=[Se{}](F)CC"]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def graph_matches(native, reference):
    atoms = [a[1] for a in native.atom_table] == [a.GetAtomicNum() for a in reference.GetAtoms()]
    bonds = {tuple(sorted((b[0], b[1]))) for b in native.bond_table}
    rdkit_bonds = {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in reference.GetBonds()}
    return atoms and bonds == rdkit_bonds


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wheel", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    rdBase.DisableLog("rdApp.*")
    rows = []
    for group, patterns, seeds in (("development", DEVELOPMENT, (7, 42)),
                                   ("controls", CONTROLS, (13, 99))):
        for pattern in patterns:
            for tag in ("@", "@@"):
                smiles = pattern.format(tag)
                native = chematic.from_smiles(smiles).add_hydrogens()
                reference = Chem.AddHs(Chem.MolFromSmiles(smiles))
                if not graph_matches(native, reference):
                    raise ValueError(f"indexed graph mismatch: {smiles}")
                centers = [a.GetIdx() for a in reference.GetAtoms()
                           if a.GetAtomicNum() == 34 and a.HasProp("_CIPCode")]
                if len(centers) != 1:
                    raise ValueError(f"expected one selenium CIP center: {smiles}")
                center = centers[0]
                expected = reference.GetAtomWithIdx(center).GetProp("_CIPCode")
                for seed in seeds:
                    molecule = Chem.Mol(reference)
                    if AllChem.EmbedMolecule(molecule, randomSeed=seed) != 0:
                        rows.append({"group": group, "smiles": smiles, "seed": seed,
                                     "status": "rdkit_embed_failed"})
                        continue
                    coordinates = [list(molecule.GetConformer().GetAtomPosition(i))
                                   for i in range(molecule.GetNumAtoms())]
                    for mirror in (False, True):
                        xyz = [[-x, y, z] if mirror else [x, y, z] for x, y, z in coordinates]
                        trial = Chem.Mol(reference)
                        conformer = Chem.Conformer(trial.GetNumAtoms())
                        for i, point in enumerate(xyz):
                            conformer.SetAtomPosition(i, point)
                        trial.AddConformer(conformer, assignId=True)
                        Chem.AssignStereochemistryFrom3D(trial, replaceExistingTags=True)
                        atom = trial.GetAtomWithIdx(center)
                        observed = atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None
                        rows.append({"group": group, "smiles": smiles, "seed": seed, "mirror": mirror,
                                     "center": center, "expected_cip": expected, "rdkit_3d_cip": observed,
                                     "rdkit_matches_declared": None if observed is None else observed == expected,
                                     "chematic_stereo_valid": native._timed_stereo_valid(xyz),
                                     "coordinates": xyz})
    extension = importlib.util.find_spec("chematic.chematic").origin
    report = {"selection": "authored development and synthetic controls, not an independent holdout",
              "rdkit_version": rdBase.rdkitVersion, "rdkit_warnings_suppressed": True,
              "platform": platform.platform(), "python": sys.version,
              "wheel_sha256": digest(args.wheel), "extension_sha256": digest(extension),
              "script_sha256": digest(__file__), "rows": rows}
    Path(args.output).write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
