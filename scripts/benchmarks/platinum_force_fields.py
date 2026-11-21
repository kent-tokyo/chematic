#!/usr/bin/env python3
"""Compare force fields on common explicit-H coordinates for development rows."""

import argparse
import hashlib
import json
import math
import platform
from pathlib import Path

import chematic
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

from platinum_development_probe import graph_matches


def compare(mol, other, coords, variant):
    result = {}
    field_mol = Chem.Mol(mol)
    try:
        if variant == "UFF":
            field = AllChem.UFFGetMoleculeForceField(field_mol)
        else:
            props = AllChem.MMFFGetMoleculeProperties(field_mol, mmffVariant=variant)
            field = (
                AllChem.MMFFGetMoleculeForceField(field_mol, props) if props else None
            )
        result["rdkit_supported"] = field is not None
    except Exception as exc:
        field = None
        result["rdkit_supported"] = False
        result["rdkit_error"] = repr(exc)
    try:
        ours = (
            other.uff_force_field(coords)
            if variant == "UFF"
            else other.mmff_force_field(coords, variant=variant)
        )
        result["chematic_supported"] = True
    except Exception as exc:
        ours = None
        result["chematic_supported"] = False
        result["chematic_error"] = repr(exc)
    if field is None or ours is None:
        return result
    energy_delta = ours.energy() - field.CalcEnergy()
    ours_gradient = [value for atom in ours.gradient() for value in atom]
    reference_gradient = list(field.CalcGrad())
    if len(ours_gradient) != len(reference_gradient):
        raise ValueError("gradient size mismatch")
    errors = [abs(a - b) for a, b in zip(ours_gradient, reference_gradient)]
    if not all(math.isfinite(x) for x in [energy_delta, *errors]):
        raise ValueError("nonfinite energy or gradient")
    positions = [point[:] for point in coords]
    finite_difference = []
    for axis in sorted({round(i * (len(errors) - 1) / 5) for i in range(6)}):
        atom, component = divmod(axis, 3)
        positions[atom][component] += 1e-5
        plus = ours.set_positions_(positions).energy()
        positions[atom][component] -= 2e-5
        minus = ours.set_positions_(positions).energy()
        positions[atom][component] += 1e-5
        finite_difference.append(abs((plus - minus) / 2e-5 - ours_gradient[axis]))
    if not all(math.isfinite(value) for value in finite_difference):
        raise ValueError("nonfinite finite-difference comparison")
    result.update(
        energy_absolute_error=abs(energy_delta),
        gradient_absolute_errors=errors,
        finite_difference_step=1e-5,
        finite_difference_absolute_errors=finite_difference,
    )
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--development", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.development.read_text())
    output = {
        "development_sha256": hashlib.sha256(args.development.read_bytes()).hexdigest(),
        "rdkit_version": rdBase.rdkitVersion,
        "chematic_version": getattr(chematic, "__version__", "unknown"),
        "platform": platform.platform(),
        "seed": 42,
        "results": [],
    }
    for row in source["results"]:
        result = {"name": row["name"], "index": row["index"], "fields": {}}
        try:
            if "smiles" not in row:
                raise ValueError("no development SMILES")
            mol = Chem.AddHs(Chem.MolFromSmiles(row["smiles"]))
            other = chematic.from_smiles(row["smiles"]).add_hydrogens()
            if not graph_matches(other, mol):
                raise ValueError("atom or bond graph mismatch")
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            params.numThreads = 1
            if AllChem.EmbedMolecule(mol, params) != 0:
                raise ValueError("RDKit embedding failed")
            coords = [
                list(mol.GetConformer().GetAtomPosition(i))
                for i in range(mol.GetNumAtoms())
            ]
            for variant in ("MMFF94", "MMFF94s", "UFF"):
                result["fields"][variant] = compare(mol, other, coords, variant)
        except Exception as exc:
            result["error"] = repr(exc)
        output["results"].append(result)
        checkpoint = args.output.with_suffix(args.output.suffix + ".tmp")
        checkpoint.write_text(json.dumps(output, allow_nan=False) + "\n")
        checkpoint.replace(args.output)
        print(row["name"], result.get("error", "ok"), flush=True)


if __name__ == "__main__":
    main()
