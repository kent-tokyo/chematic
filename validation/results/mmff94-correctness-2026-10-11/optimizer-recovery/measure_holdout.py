"""Separate regression holdout: never appended to the frozen 265-input corpus.

Run from repository root with the development Python environment:
python .../measure_holdout.py candidate
"""
import hashlib
import json
import sys
from pathlib import Path

import chematic
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
from scripts.mmff94_same_explicit_h_energy import gradient_diagnostic, rdkit_term_energy_breakdown

OUT = Path(__file__).resolve().parent
SMILES = [
    "O", "NC(=[NH2+])N", "[NH2+]=C(N)N", "CNC(=[NH2+])N",
    "CN(C)C(=[N+](C)C)N(C)C", "NNC(=[NH2+])N", "N=C(N)N",
    "CC(=[NH2+])N", "CNC=[NH2+]", "NC(=[NH2+])NC(=[NH2+])N",
    "N=C(N)NC(=[NH2+])N", "CO", "COC", "[OH-]", "[OH3+]",
    "C=[OH+]", "C=[NH2+]", "C[NH3+]", "NC(=[NH2+])NCCCCC",
    "NC(=[NH2+])[N+](=O)[O-]",
]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    assert rdBase.rdkitVersion == "2026.03.6"
    extension = Path(chematic._core.__file__) if hasattr(chematic, "_core") else Path("crates/chematic-py/python/chematic/chematic.cpython-313-darwin.so")
    artifacts = {"extension": {"path": str(extension), "sha256": digest(extension)},
                 "rdBase": {"path": rdBase.__file__, "sha256": digest(rdBase.__file__)}}
    rows = []
    for smiles in SMILES:
        reference = Chem.AddHs(Chem.MolFromSmiles(smiles))
        assert AllChem.EmbedMolecule(reference, randomSeed=43) == 0
        properties = AllChem.MMFFGetMoleculeProperties(reference, mmffVariant="MMFF94")
        mol = chematic.from_smiles(smiles).add_hydrogens()
        assert [a["element"] for a in mol.depict_data()["atoms"]] == [a.GetSymbol() for a in reference.GetAtoms()]
        assert {(min(a, b), max(a, b), order) for a, b, order, *_ in mol.bond_table} == {
            (min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), str(b.GetBondType()))
            for b in reference.GetBonds()
        }
        coords = reference.GetConformer().GetPositions().tolist()
        for index, point in enumerate(coords):
            point[2] += 0.03 * ((index % 5) - 2)
            reference.GetConformer().SetAtomPosition(index, point)
        types = list(mol.mmff94_numeric_atom_types())
        charges = list(mol.mmff94_charges_numeric())
        expected_types = [properties.GetMMFFAtomType(i) for i in range(reference.GetNumAtoms())]
        expected_charges = [properties.GetMMFFPartialCharge(i) for i in range(reference.GetNumAtoms())]
        force_field = AllChem.MMFFGetMoleculeForceField(reference, properties, confId=0)
        terms = mol.mmff94_energy_breakdown(coords)
        expected_terms = rdkit_term_energy_breakdown(reference)
        gradient = [v for point in mol.mmff94_bounded_analytic_gradient(coords) for v in point]
        expected_gradient = force_field.CalcGrad()
        errors = [abs(a-b) for a,b in zip(gradient, expected_gradient, strict=True)]
        rows.append({
            "smiles": smiles, "coordinate_sha256": hashlib.sha256(json.dumps(coords, separators=(",", ":")).encode()).hexdigest(),
            "types": types, "rdkit_types": expected_types, "types_match": types == expected_types,
            "charges": charges, "rdkit_charges": expected_charges,
            "max_charge_error": max(abs(a-b) for a,b in zip(charges, expected_charges, strict=True)),
            "charge_sum": sum(charges), "formal_charge": Chem.GetFormalCharge(reference),
            "terms": terms, "rdkit_terms": expected_terms,
            "rdkit_total": force_field.CalcEnergy(), "abs_total_error": abs(terms["total"] - force_field.CalcEnergy()),
            "term_abs_errors": {term: abs(terms[term]-value) for term,value in expected_terms.items()},
            "gradient_components": len(errors), "max_gradient_abs_error": max(errors),
            "max_gradient_scaled_error": max(error/(1.0+abs(value)) for error,value in zip(errors,expected_gradient,strict=True)),
            "internal_gradient": gradient_diagnostic(mol, coords, 1e-5),
        })
    for artifact in artifacts.values():
        assert digest(artifact["path"]) == artifact["sha256"]
    result = {"kind": "separate_cation_water_regression_holdout", "outside_frozen_265": True,
              "rdkit_runtime": rdBase.rdkitVersion, "artifacts": artifacts,
              "protocol": {"seed": 43, "explicit_hydrogens": True, "coordinate_z_perturbation": "0.03 * ((atom_index % 5) - 2)", "source_tree": "dirty_candidate"},
              "input_count": len(rows), "rows": rows,
              "all_types_match": all(row["types_match"] for row in rows),
              "max_charge_error": max(row["max_charge_error"] for row in rows),
              "max_total_energy_error": max(row["abs_total_error"] for row in rows),
              "max_gradient_abs_error": max(row["max_gradient_abs_error"] for row in rows),
              "max_gradient_scaled_error": max(row["max_gradient_scaled_error"] for row in rows),
              "max_internal_gradient_scaled_error": max(row["internal_gradient"]["max_scaled_error"] for row in rows)}
    (OUT / (sys.argv[1] + "-holdout.json")).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k:v for k,v in result.items() if k not in ("rows", "artifacts")}, indent=2))
    if sys.argv[1] == "candidate":
        assert result["all_types_match"]
        assert result["max_charge_error"] < 1e-12
        assert result["max_total_energy_error"] < 1e-8
        assert result["max_gradient_scaled_error"] < 1e-7
        assert result["max_internal_gradient_scaled_error"] < 1e-5


if __name__ == "__main__":
    main()
