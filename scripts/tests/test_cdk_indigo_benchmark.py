import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).parents[2]
MODULE_PATH = ROOT / "validation/cdk_indigo_comparison/benchmark.py"
SPEC = importlib.util.spec_from_file_location("cdk_indigo_benchmark", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_ratio_ci_is_deterministic_and_requires_every_block_for_strict_win():
    result = MODULE.ratio_ci([2_000] * 21, [1_000] * 21, seed=17)
    assert result == {
        "median": 2.0,
        "ci95": [2.0, 2.0],
        "winning_blocks": 21,
        "blocks": 21,
    }


def test_contract_pins_stable_comparator_versions_and_claim_boundary():
    contract = json.loads(
        (ROOT / "validation/cdk_indigo_comparison/contract.json").read_text()
    )
    assert contract["comparators"]["cdk"]["version"] == "2.13"
    assert contract["comparators"]["indigo"]["version"] == "1.46.0"
    assert "not universal superiority" in contract["acceptance"]["claim_boundary"]
