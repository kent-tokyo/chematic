import importlib.util
import copy
import json
from pathlib import Path
import sys

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

SCRIPTS = Path(__file__).parents[1]
sys.path.insert(0, str(SCRIPTS))
SPEC = importlib.util.spec_from_file_location("cleanup_gate", SCRIPTS / "rdkit_cleanup_source_gate.py")
gate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gate)
CHECK_SPEC = importlib.util.spec_from_file_location(
    "cleanup_evidence", SCRIPTS / "checks/evidence/check_standardization_source_evidence.py"
)
checker = importlib.util.module_from_spec(CHECK_SPEC)
CHECK_SPEC.loader.exec_module(checker)


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


def test_ferrocene_charge_distribution_remains_a_chemical_difference():
    molecule = Chem.MolFromSmiles("C1=CC=CC1[Fe]C1C=CC=C1")
    oracle = rdMolStandardize.Cleanup(molecule)
    actual = "C1=C[CH-]C=C1.C1=C[CH-]C=C1.[Fe+2]"
    assert gate.equivalent_after_sanitization(oracle, actual) is False


def test_record_checker_rejects_inconsistent_semantic_classification():
    import pytest

    record = json.loads(checker.RECORD.read_text())
    checker.check_current_record(record)
    altered = copy.deepcopy(record)
    altered["mismatches"][0]["semantic_classification"] = "representation_after_sanitization"
    with pytest.raises(AssertionError):
        checker.check_current_record(altered)
