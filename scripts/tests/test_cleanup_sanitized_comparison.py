import importlib.util
from pathlib import Path
import sys

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

SCRIPTS = Path(__file__).parents[1]
sys.path.insert(0, str(SCRIPTS))
SPEC = importlib.util.spec_from_file_location("cleanup_gate", SCRIPTS / "rdkit_cleanup_source_gate.py")
gate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gate)


def test_compares_sanitized_copy_without_changing_oracle():
    molecule = Chem.MolFromSmiles("CCCN(CCC)C1CCC2=C(C1)C([C](=O)[Zn])=CC=C2")
    oracle = rdMolStandardize.Cleanup(molecule)
    before = Chem.MolToSmiles(oracle)
    actual = "CCCN(CCC)C1CCc2cccc([C]=O)c2C1.[Zn]"
    assert gate.equivalent_after_sanitization(oracle, actual) is True
    assert Chem.MolToSmiles(oracle) == before
    assert gate.equivalent_after_sanitization(oracle, "CC") is False
    assert gate.equivalent_after_sanitization(oracle, None) is None
    assert gate.equivalent_after_sanitization(oracle, "not-smiles") is None
