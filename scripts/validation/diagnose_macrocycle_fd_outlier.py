#!/usr/bin/env python3
"""Reproduce a finite-difference outlier on an exposed macrocycle record."""
import argparse
import hashlib
import json
import platform
from pathlib import Path

import chematic
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Previously scored field-shard input")
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--axis", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    raw = args.input.read_bytes()
    matches = [r for r in json.loads(raw)["results"] if r["index"] == args.index]
    if len(matches) != 1:
        raise ValueError("index must identify exactly one record")
    row = matches[0]
    mol = Chem.AddHs(Chem.MolFromSmiles(row["smiles"]))
    params = AllChem.ETKDGv3()
    params.randomSeed, params.numThreads = 42, 1
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise ValueError("RDKit embedding failed")
    points = [list(mol.GetConformer().GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    positions = [value for point in points for value in point]
    if not 0 <= args.axis < len(positions):
        raise ValueError("axis out of range")
    other = chematic.from_smiles(row["smiles"]).add_hydrogens()
    result = {"name": row["name"], "index": args.index, "axis": args.axis,
              "input_sha256": hashlib.sha256(raw).hexdigest(),
              "coordinates_sha256": hashlib.sha256(json.dumps(points, separators=(",", ":")).encode()).hexdigest(),
              "platform": platform.platform(), "rdkit_version": rdBase.rdkitVersion,
              "seed": 42, "fields": {}}
    for variant in ("MMFF94", "UFF"):
        graph = Chem.Mol(mol)
        if variant == "UFF":
            reference = AllChem.UFFGetMoleculeForceField(graph)
            candidate = other.uff_force_field(points)
        else:
            properties = AllChem.MMFFGetMoleculeProperties(graph, mmffVariant=variant)
            reference = AllChem.MMFFGetMoleculeForceField(graph, properties)
            candidate = other.mmff_force_field(points, variant=variant)
        analytic = reference.CalcGrad(positions)[args.axis]
        candidate_gradient = [v for atom in candidate.gradient() for v in atom][args.axis]
        steps = []
        for step in (1e-4, 1e-5, 1e-6, 1e-7, 1e-8):
            plus, minus = positions.copy(), positions.copy()
            plus[args.axis] += step
            minus[args.axis] -= step
            finite = (reference.CalcEnergy(plus) - reference.CalcEnergy(minus)) / (2 * step)
            steps.append({"step": step, "finite_difference": finite,
                          "absolute_error": abs(finite - analytic)})
        result["fields"][variant] = {"rdkit_analytic": analytic,
                                     "chematic_analytic": candidate_gradient,
                                     "steps": steps}
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
