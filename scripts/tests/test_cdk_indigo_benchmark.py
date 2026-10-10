import importlib.util
import json
from pathlib import Path
import sys


ROOT = Path(__file__).parents[2]
MODULE_PATH = ROOT / "validation/cdk_indigo_comparison/benchmark.py"
SPEC = importlib.util.spec_from_file_location("cdk_indigo_benchmark", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)

ADAPTER_MODULE_PATH = ROOT / "validation/cdk_indigo_comparison/python_adapter.py"
ADAPTER_SPEC = importlib.util.spec_from_file_location(
    "cdk_indigo_python_adapter", ADAPTER_MODULE_PATH
)
ADAPTER_MODULE = importlib.util.module_from_spec(ADAPTER_SPEC)
assert ADAPTER_SPEC.loader is not None
sys.modules[ADAPTER_SPEC.name] = ADAPTER_MODULE
ADAPTER_SPEC.loader.exec_module(ADAPTER_MODULE)


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


def test_accuracy_operations_do_not_inherit_mutated_engine_state():
    parsed = []

    def parse(_smiles):
        molecule = {"mutated": False}
        parsed.append(molecule)
        return molecule

    def mutate(molecule):
        molecule["mutated"] = True
        return True

    engine = ADAPTER_MODULE.Engine(
        "fake",
        "1",
        parse,
        {"mutate": mutate, "observe": lambda molecule: molecule["mutated"]},
    )
    assert ADAPTER_MODULE.fresh_operation_result(engine, "C", "mutate")["value"]
    assert not ADAPTER_MODULE.fresh_operation_result(engine, "C", "observe")["value"]
    assert len(parsed) == 2
