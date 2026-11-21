#!/usr/bin/env python3
"""Reproduce a near-linear MMFF minimum and independently inspect its gradient."""

import argparse
import json
import math
from pathlib import Path

import chematic
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem


SMILES = "CC(C)([NH3+])CNc1nc(-c2ccncc2)cc2cnccc12"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    other = chematic.from_smiles(SMILES).add_hydrogens()
    mol = Chem.AddHs(Chem.MolFromSmiles(SMILES))
    assert [atom[1] for atom in other.atom_table] == [
        atom.GetAtomicNum() for atom in mol.GetAtoms()
    ]
    cbonds = sorted(
        (min(i, j), max(i, j), order) for i, j, order, _ in other.bond_table
    )
    rbonds = sorted(
        (
            min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
            max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
            str(b.GetBondType()).upper(),
        )
        for b in mol.GetBonds()
    )
    assert cbonds == rbonds
    results = []
    for index, (status, energy, coords) in enumerate(
        other.etkdg_embed_multiple_confs_mmff_corrected(20, random_seed=42)
    ):
        field = other.mmff_force_field(coords)
        cgradient = [v for atom in field.gradient() for v in atom]
        reference = Chem.Mol(mol)
        conformer = Chem.Conformer(reference.GetNumAtoms())
        for i, point in enumerate(coords):
            conformer.SetAtomPosition(i, point)
        reference.AddConformer(conformer)
        props = AllChem.MMFFGetMoleculeProperties(reference)
        rfield = AllChem.MMFFGetMoleculeForceField(reference, props)
        rgradient = rfield.CalcGrad()
        assert len(cgradient) == len(rgradient)
        a = [coords[4][axis] - coords[5][axis] for axis in range(3)]
        b = [coords[6][axis] - coords[5][axis] for axis in range(3)]
        cosine = sum(x * y for x, y in zip(a, b)) / math.sqrt(
            sum(x * x for x in a) * sum(x * x for x in b)
        )
        record = {
            "index": index,
            "status": status,
            "energy": energy,
            "angle_degrees": math.degrees(math.acos(max(-1.0, min(1.0, cosine)))),
            "chematic_gradient_max": max(map(abs, cgradient)),
            "rdkit_gradient_max": max(map(abs, rgradient)),
            "max_component_difference": max(
                abs(x - y) for x, y in zip(cgradient, rgradient)
            ),
        }
        if record["chematic_gradient_max"] > 10000:
            axis = max(range(len(cgradient)), key=lambda i: abs(cgradient[i]))
            terms = ("Bond", "Angle", "StretchBend", "Oop", "Torsion", "VdW", "Ele")
            parts = {}
            for term in terms:
                properties = AllChem.MMFFGetMoleculeProperties(reference)
                for name in terms:
                    getattr(properties, "SetMMFF" + name + "Term")(name == term)
                part = AllChem.MMFFGetMoleculeForceField(reference, properties)
                positions = list(part.Positions())
                derivative = part.CalcGrad()[axis]
                numeric = {}
                for step in (1e-8, 1e-10):
                    positions[axis] += step
                    plus = part.CalcEnergy(positions)
                    positions[axis] -= 2 * step
                    minus = part.CalcEnergy(positions)
                    positions[axis] += step
                    numeric[str(step)] = (plus - minus) / (2 * step)
                parts[term] = {"analytic": derivative, "central_difference": numeric}
            record["largest_gradient_component_index"] = axis
            record["rdkit_terms"] = parts
        results.append(record)
    args.output.write_text(
        json.dumps(
            {
                "smiles": SMILES,
                "rdkit": rdBase.rdkitVersion,
                "seed": 42,
                "results": results,
            },
            indent=2,
            allow_nan=False,
        )
        + "\n"
    )
    print(
        json.dumps(
            [row for row in results if row["chematic_gradient_max"] > 0.01], indent=2
        )
    )


if __name__ == "__main__":
    main()
