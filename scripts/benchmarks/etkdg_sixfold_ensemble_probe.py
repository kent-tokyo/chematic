"""Independent RDKit 2026.09.1 MMFF oracle for the sixfold ETKDG probe.

Run the Rust example to a JSONL file, then pass that file as the sole argument.
The six development SMILES were selected because GetExperimentalTorsions has
nonzero V6 for each; no frozen holdout rows are used.
"""

import json
import sys

from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom


def summarize(rows):
    return {
        "count": len(rows),
        "converged": sum(row["status"] == 0 for row in rows),
        "best_energy": min((row["energy"] for row in rows), default=None),
    }


def max_energy_error(mol, rows):
    """Re-score CheMatic coordinates with RDKit's own MMFF94 force field."""
    properties = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
    assert properties is not None
    errors = []
    for row in rows:
        assert len(row["coords"]) == mol.GetNumAtoms()
        conformer = Chem.Conformer(mol.GetNumAtoms())
        for atom_index, point in enumerate(row["coords"]):
            conformer.SetAtomPosition(atom_index, point)
        mol.RemoveAllConformers()
        mol.AddConformer(conformer)
        field = AllChem.MMFFGetMoleculeForceField(mol, properties)
        assert field is not None
        errors.append(abs(field.CalcEnergy() - row["energy"]))
    return max(errors, default=0.0)


def main():
    assert Chem.rdBase.rdkitVersion == "2026.09.1"
    for line in open(sys.argv[1], encoding="utf-8"):
        source = json.loads(line)
        mol = Chem.AddHs(Chem.MolFromSmiles(source["smiles"]))
        sixfold = [
            float(torsion["V"][5])
            for torsion in rdDistGeom.GetExperimentalTorsions(mol)
            if abs(float(torsion["V"][5])) > 0
        ]
        assert sixfold, source["smiles"]
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        params.numThreads = 0
        ids = AllChem.EmbedMultipleConfs(mol, numConfs=20, params=params)
        optimized = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0, maxIters=200)
        assert len(ids) == len(optimized)
        oracle = [{"status": status, "energy": energy} for status, energy in optimized]
        mol_40 = Chem.AddHs(Chem.MolFromSmiles(source["smiles"]))
        ids_40 = AllChem.EmbedMultipleConfs(mol_40, numConfs=40, params=params)
        optimized_40 = AllChem.MMFFOptimizeMoleculeConfs(
            mol_40, numThreads=0, maxIters=200
        )
        assert len(ids_40) == len(optimized_40)
        oracle_40 = [
            {"status": status, "energy": energy} for status, energy in optimized_40
        ]
        scoring_mol = Chem.AddHs(Chem.MolFromSmiles(source["smiles"]))
        energy_error = max(
            max_energy_error(scoring_mol, source["legacy"]),
            max_energy_error(scoring_mol, source["corrected"]),
        )
        assert energy_error < 1e-6, (source["smiles"], energy_error)
        compatibility = summarize(source["legacy"])
        corrected = summarize(source["corrected"])
        chosen = min(
            (compatibility, corrected),
            key=lambda summary: summary["best_energy"],
        )
        print(
            json.dumps(
                {
                    "smiles": source["smiles"],
                    "sixfold_coefficients": sixfold,
                    "max_same_coordinate_mmff_energy_error": energy_error,
                    "rdkit": summarize(oracle),
                    "rdkit_40": summarize(oracle_40),
                    "chematic_compatibility": compatibility,
                    "chematic_corrected": corrected,
                    "chematic_best_gradient": chosen,
                },
                sort_keys=True,
            )
        )


if __name__ == "__main__":
    main()
