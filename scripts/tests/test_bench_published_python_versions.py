"""Unit tests for the paired Python-wheel benchmark worker."""

import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "bench_published_python_versions.py"
SPEC = importlib.util.spec_from_file_location("bench_published_python_versions", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class StubMol:
    tpsa = 12.5
    labute_asa = 27.25


def test_numeric_descriptor_digests_use_exact_float_encoding():
    mol = StubMol()
    assert MODULE.evaluate(mol, "tpsa") == float(12.5).hex()
    assert MODULE.evaluate(mol, "labute_asa") == float(27.25).hex()


def test_paired_ratio_name_matches_the_calculation_direction():
    result = MODULE.ci([2.0] * 20)
    assert result["median_a_over_b"] == 2.0
    assert "median_b_over_a" not in result
