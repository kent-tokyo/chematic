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


def compare(mol, other, coords, variant, term_breakdown=False):
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
    if term_breakdown and variant in ("MMFF94", "MMFF94s"):
        properties = AllChem.MMFFGetMoleculeProperties(field_mol, mmffVariant=variant)
        types = other.mmff94_numeric_atom_types()
        charges = other.mmff94_charges_numeric()
        rdkit_types = [properties.GetMMFFAtomType(i) for i in range(field_mol.GetNumAtoms())]
        rdkit_charges = [properties.GetMMFFPartialCharge(i) for i in range(field_mol.GetNumAtoms())]
        if len(types) != len(rdkit_types) or len(charges) != len(rdkit_charges):
            raise ValueError("MMFF atom typing or charge count mismatch")
        result["atom_type_mismatches"] = [i for i, (a, b) in enumerate(zip(types, rdkit_types)) if a != b]
        result["charge_absolute_errors"] = [abs(a - b) for a, b in zip(charges, rdkit_charges)]
        if not all(math.isfinite(value) for value in result["charge_absolute_errors"]):
            raise ValueError("nonfinite MMFF charge comparison")
        setters = {"bond": "SetMMFFBondTerm", "angle": "SetMMFFAngleTerm",
                   "stretch_bend": "SetMMFFStretchBendTerm", "oop": "SetMMFFOopTerm",
                   "torsion": "SetMMFFTorsionTerm", "vdw": "SetMMFFVdWTerm",
                   "ele": "SetMMFFEleTerm"}
        result["terms"] = {}
        for term, setter in setters.items():
            properties = AllChem.MMFFGetMoleculeProperties(field_mol, mmffVariant=variant)
            for name in setters.values():
                getattr(properties, name)(False)
            getattr(properties, setter)(True)
            reference_field = AllChem.MMFFGetMoleculeForceField(field_mol, properties)
            energy, gradient = other._rdkit_mmff_terms(coords, terms=[term], variant=variant)
            reference_energy = reference_field.CalcEnergy()
            reference_gradient = reference_field.CalcGrad()
            if len(gradient) != len(reference_gradient):
                raise ValueError(f"gradient size mismatch for {term}")
            errors = [abs(a - b) for a, b in zip(gradient, reference_gradient)]
            energy_error = abs(energy - reference_energy)
            if not all(math.isfinite(value) for value in [energy_error, *errors]):
                raise ValueError(f"nonfinite energy or gradient for {term}")
            result["terms"][term] = {
                "energy_absolute_error": energy_error,
                "gradient_absolute_errors": errors,
            }
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--development", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--checkpoint-every", type=int, default=1)
    parser.add_argument("--store-coordinates", action="store_true",
                        help="include the explicit-hydrogen common coordinates")
    parser.add_argument("--term-breakdown", action="store_true",
                        help="compare MMFF types, charges and each force-field term")
    args = parser.parse_args()
    if args.checkpoint_every < 1:
        parser.error("checkpoint-every must be positive")
    source = json.loads(args.development.read_text())
    output = {
        "development_sha256": hashlib.sha256(args.development.read_bytes()).hexdigest(),
        "rdkit_version": rdBase.rdkitVersion,
        "chematic_version": getattr(chematic, "__version__", "unknown"),
        "platform": platform.platform(),
        "seed": 42,
        "results": [],
    }
    if args.store_coordinates:
        output["coordinates_recorded"] = True
    if args.term_breakdown:
        output["term_breakdown"] = True
    if args.store_coordinates or args.term_breakdown:
        output["script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.resume:
        existing = json.loads(args.output.read_text())
        if existing.get("coordinates_recorded", False) != args.store_coordinates:
            parser.error("resume coordinate recording mode differs")
        if existing.get("term_breakdown", False) != args.term_breakdown:
            parser.error("resume term breakdown mode differs")
        if (args.store_coordinates or args.term_breakdown) and existing.get("script_sha256") != output["script_sha256"]:
            parser.error("resume scorer version differs")
        for key in ("development_sha256", "rdkit_version", "chematic_version", "platform", "seed"):
            if existing[key] != output[key]:
                parser.error(f"resume metadata differs: {key}")
        previous = existing["results"]
        expected_indices = [r["index"] for r in source["results"][:len(previous)]]
        if [r["index"] for r in previous] != expected_indices:
            parser.error("resume checkpoint does not match the input prefix")
        output["results"] = previous
    for row in source["results"][len(output["results"]):]:
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
            if args.store_coordinates:
                result["coordinates"] = coords
            for variant in ("MMFF94", "MMFF94s", "UFF"):
                result["fields"][variant] = compare(mol, other, coords, variant, args.term_breakdown)
        except Exception as exc:
            result["error"] = repr(exc)
        output["results"].append(result)
        if len(output["results"]) % args.checkpoint_every == 0 or len(output["results"]) == len(source["results"]):
            checkpoint = args.output.with_suffix(args.output.suffix + ".tmp")
            checkpoint.write_text(json.dumps(output, allow_nan=False) + "\n")
            checkpoint.replace(args.output)
        print(row["name"], result.get("error", "ok"), flush=True)


if __name__ == "__main__":
    main()
