import importlib.util
from pathlib import Path


ROOT = Path(__file__).parents[2]
MODULE_PATH = ROOT / "validation/indigo_workflows/benchmark.py"
SPEC = importlib.util.spec_from_file_location("indigo_workflow_benchmark", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_ratio_summary_is_deterministic_and_requires_every_block_for_strict_win():
    result = MODULE.ratio_summary([2_000] * 21, [1_000] * 21)
    assert result == {
        "indigo_over_chematic_median": 2.0,
        "ci95": [2.0, 2.0],
        "winning_blocks": 21,
        "blocks": 21,
    }


def test_rgroup_corpus_has_parallel_expected_rows():
    import json

    cases = json.loads((MODULE.RGROUP_CORPUS).read_text())["cases"]
    assert sum(len(case["molecules"]) for case in cases) == sum(
        len(case["expected"]) for case in cases
    )


def test_comparator_versions_are_pinned():
    assert MODULE.EXPECTED_CHEMATIC == "1.1.0"
    assert MODULE.EXPECTED_INDIGO == "1.46.0"
    assert MODULE.EXPECTED_RDKIT == "2026.03.6"
