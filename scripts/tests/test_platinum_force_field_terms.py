import sys
from pathlib import Path

import chematic
from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks"))
from platinum_force_fields import compare


def test_term_comparison_records_types_charges_and_components():
    reference = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    assert AllChem.EmbedMolecule(reference, randomSeed=42) == 0
    coords = [list(reference.GetConformer().GetAtomPosition(i))
              for i in range(reference.GetNumAtoms())]
    mol = chematic.from_smiles("CCO").add_hydrogens()
    result = compare(reference, mol, coords, "MMFF94", term_breakdown=True)
    assert result["atom_type_mismatches"] == []
    assert max(result["charge_absolute_errors"]) < 1e-9
    assert set(result["terms"]) == {"bond", "angle", "stretch_bend", "oop", "torsion", "vdw", "ele"}
    for value in result["terms"].values():
        assert value["energy_absolute_error"] < 1e-7
        assert max(value["gradient_absolute_errors"], default=0) < 1e-5
