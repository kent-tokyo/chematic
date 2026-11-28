import sys
from pathlib import Path

from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "validation"))
from classify_conformer_holdout import classify


def test_classification_uses_ring_and_rotatable_bond_thresholds():
    assert classify(Chem.MolFromSmiles("c1ccccc1")) == "ordinary"
    assert classify(Chem.MolFromSmiles("CCCCCCCCCCC")) == "flexible"
    assert classify(Chem.MolFromSmiles("C1CCCCCCCCCCC1")) == "macrocycle"
