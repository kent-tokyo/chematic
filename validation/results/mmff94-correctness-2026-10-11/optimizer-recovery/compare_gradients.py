"""Compare analytic gradients with RDKit on the energy diagnostic's exact coordinates."""

import hashlib
import json
from pathlib import Path

import chematic
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

OUT = Path("validation/results/mmff94-correctness-2026-10-11/optimizer-recovery")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    summary = json.loads((OUT / "candidate-summary.json").read_bytes())
    assert rdBase.rdkitVersion == summary["versions"]["rdkit_runtime"]
    artifacts = summary["artifacts"]
    for artifact in artifacts.values():
        assert digest(artifact["path"]) == artifact["sha256"]
    raw = (OUT / "candidate.jsonl").read_bytes()
    rows = [json.loads(line) for line in raw.splitlines()]
    results = []
    for row in rows:
        if row["status"] != "ok":
            results.append({"input_index": row["input_index"], "status": "not_compared", "energy_status": row["status"]})
            continue
        rdkit_mol = Chem.AddHs(Chem.MolFromSmiles(row["smiles"]))
        assert AllChem.EmbedMolecule(rdkit_mol, randomSeed=summary["protocol"]["seed"] + row["input_index"]) == 0
        coords = rdkit_mol.GetConformer().GetPositions().tolist()
        coordinate_hash = hashlib.sha256(json.dumps(coords, separators=(",", ":")).encode()).hexdigest()
        assert coordinate_hash == row["coordinate_sha256"]
        properties = AllChem.MMFFGetMoleculeProperties(rdkit_mol, mmffVariant="MMFF94")
        force_field = AllChem.MMFFGetMoleculeForceField(rdkit_mol, properties, confId=0)
        assert force_field.CalcEnergy() == row["rdkit_energy_kcal_mol"]
        mol = chematic.from_smiles(row["smiles"]).add_hydrogens()
        assert [a["element"] for a in mol.depict_data()["atoms"]] == [a.GetSymbol() for a in rdkit_mol.GetAtoms()]
        assert {tuple(sorted((a, b))) for a, b, *_ in mol.bond_table} == {
            tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in rdkit_mol.GetBonds()
        }
        actual = [value for point in mol.mmff94_bounded_analytic_gradient(coords) for value in point]
        expected = force_field.CalcGrad()
        errors = [abs(a - b) for a, b in zip(actual, expected, strict=True)]
        scaled = [error / (1.0 + abs(b)) for error, b in zip(errors, expected, strict=True)]
        results.append({
            "input_index": row["input_index"], "status": "compared", "components": len(errors),
            "coordinate_sha256": coordinate_hash,
            "max_abs_error_kcal_mol_angstrom": max(errors), "max_scaled_error": max(scaled),
        })
    for artifact in artifacts.values():
        assert digest(artifact["path"]) == artifact["sha256"]
    comparable = [r for r in results if r["status"] == "compared"]
    result = {
        "kind": "direct_rdkit_vs_chematic_analytic_gradient",
        "artifacts": artifacts,
        "energy_rows_sha256": hashlib.sha256(raw).hexdigest(),
        "rows": results,
        "input_rows": len(results), "comparable_rows": len(comparable),
        "components": sum(r["components"] for r in comparable),
        "max_abs_error_kcal_mol_angstrom": max(r["max_abs_error_kcal_mol_angstrom"] for r in comparable),
        "max_scaled_error": max(r["max_scaled_error"] for r in comparable),
        "absolute_error_threshold_kcal_mol_angstrom": 1e-6,
        "scaled_error_threshold": 1e-7,
    }
    result["gate_passed"] = result["max_abs_error_kcal_mol_angstrom"] < 1e-6 and result["max_scaled_error"] < 1e-7
    (OUT / "gradient-comparison.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key not in ("rows", "artifacts")}, indent=2))
    assert result["gate_passed"]


if __name__ == "__main__":
    main()
