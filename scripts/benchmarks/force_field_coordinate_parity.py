#!/usr/bin/env python3
"""Compare published force fields on identical, explicitly hydrogenated coordinates."""

import argparse
import hashlib
import json
import math
import platform
from pathlib import Path

import chematic
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

from etkdg_quality_seed_budget_probe import MANIFEST, select_rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    results = []
    for row in select_rows(start=12, count=21):
        mol = Chem.AddHs(Chem.MolFromSmiles(row["smiles"]))
        other = chematic.from_smiles(row["smiles"]).add_hydrogens()
        atoms_match = [a[1] for a in other.atom_table] == [a.GetAtomicNum() for a in mol.GetAtoms()]
        bonds = sorted((min(a, b), max(a, b), order) for a, b, order, _ in other.bond_table)
        expected = sorted((min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), str(b.GetBondType()).upper()) for b in mol.GetBonds())
        result = {"name": row["name"], "smiles": row["smiles"], "graph_matches": atoms_match and bonds == expected, "fields": {}}
        if not result["graph_matches"]:
            results.append(result)
            continue
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        if AllChem.EmbedMolecule(mol, params) != 0:
            result["embedding_failed"] = True
            results.append(result)
            continue
        coords = [list(mol.GetConformer().GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
        for variant in ("MMFF94", "MMFF94s", "UFF"):
            try:
                # MMFF property preparation can change aromaticity on its input.
                # Keep every force field on a fresh copy of the same graph.
                field_mol = Chem.Mol(mol)
                if variant == "UFF":
                    rf = AllChem.UFFGetMoleculeForceField(field_mol)
                    cf = other.uff_force_field(coords)
                else:
                    props = AllChem.MMFFGetMoleculeProperties(field_mol, mmffVariant=variant)
                    rf = AllChem.MMFFGetMoleculeForceField(field_mol, props) if props else None
                    cf = other.mmff_force_field(coords, variant=variant)
                if rf is None:
                    result["fields"][variant] = {"rdkit_supported": False}
                    continue
                energy_error = abs(cf.energy() - rf.CalcEnergy())
                cg = [component for atom in cf.gradient() for component in atom]
                rg = list(rf.CalcGrad())
                assert len(cg) == len(rg)
                errors = [abs(a - b) for a, b in zip(cg, rg)]
                assert all(math.isfinite(x) for x in [energy_error, *errors])
                numeric_errors = []
                positions = [point[:] for point in coords]
                # Check axes spanning the molecule, including its first and last atoms.
                for axis in sorted({round(i * (len(cg) - 1) / 5) for i in range(6)}):
                    atom, component = divmod(axis, 3)
                    positions[atom][component] += 1e-5
                    plus = cf.set_positions_(positions).energy()
                    positions[atom][component] -= 2e-5
                    minus = cf.set_positions_(positions).energy()
                    positions[atom][component] += 1e-5
                    numeric_errors.append(abs((plus - minus) / 2e-5 - cg[axis]))
                result["fields"][variant] = {"energy_absolute_error": energy_error, "gradient_max_absolute_error": max(errors), "gradient_absolute_errors": errors, "finite_difference_step_angstrom": 1e-5, "sampled_gradient_fd_absolute_errors": numeric_errors}
            except Exception as exc:
                result["fields"][variant] = {"error": repr(exc)}
        results.append(result)
        print(row["name"], flush=True)
    output = {"rdkit_version": rdBase.rdkitVersion, "chematic_version": getattr(chematic, "__version__", "unknown"), "host": platform.platform(), "seed": 42, "selector": "existing development tier B eligibility and SHA-256 order, positions 12:33", "manifest_sha256": hashlib.sha256(MANIFEST.read_bytes()).hexdigest(), "results": results}
    args.output.write_text(json.dumps(output, indent=2) + "\n")


if __name__ == "__main__":
    main()
