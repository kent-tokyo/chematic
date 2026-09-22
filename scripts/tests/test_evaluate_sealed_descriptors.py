import importlib.util
import math
from pathlib import Path


ROOT = Path(__file__).parents[2]
SCRIPT = ROOT / "scripts/evaluate_sealed_descriptors.py"
SPEC = importlib.util.spec_from_file_location("sealed_descriptors", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def expected_values():
    return {
        "molecular_weight": 100.0,
        "hba": 2,
        "hbd": 1,
        "tpsa": 30.0,
        "logp": 1.0,
        "molar_refractivity": 20.0,
        "fsp3": 0.5,
        "aromatic_ring_count": 1,
    }


def test_compare_fields_accepts_the_complete_strict_profile():
    compared = MODULE.compare_fields(expected_values(), expected_values())
    assert set(compared) == set(MODULE.FIELD_TOLERANCES)
    assert all(result["strict_match"] for result in compared.values())


def test_compare_fields_rejects_missing_nonfinite_and_fractional_integer_values():
    actual = expected_values()
    actual.pop("tpsa")
    actual["logp"] = math.nan
    actual["hba"] = 1.5
    compared = MODULE.compare_fields(actual, expected_values())
    assert compared["tpsa"]["status"] == "unsupported"
    assert compared["logp"]["status"] == "unsupported"
    assert compared["hba"]["status"] == "mismatch"
