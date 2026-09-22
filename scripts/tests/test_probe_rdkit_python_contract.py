import json
from pathlib import Path

from scripts.probe_rdkit_python_contract import (
    backend_record,
    load_smiles,
    observe,
    timed,
)


def test_load_smiles_supports_smi_and_manifest_jsonl(tmp_path: Path):
    smi = tmp_path / "rows.smi"
    smi.write_text("# comment\nCCO ethanol\n\nCCN amine\n", encoding="utf-8")
    assert load_smiles(smi, 10) == ["CCO", "CCN"]

    jsonl = tmp_path / "rows.jsonl"
    jsonl.write_text(
        "\n".join(
            [
                json.dumps({"_manifest": True}),
                json.dumps({"smiles": "CCC"}),
            ]
        ),
        encoding="utf-8",
    )
    assert load_smiles(jsonl, 10) == ["CCC"]


def test_observe_preserves_return_and_exception_contracts():
    assert observe(lambda: True) == {"outcome": "returned", "value": True}

    def fail():
        raise TypeError("typed failure")

    result = observe(fail)
    assert result["outcome"] == "raised"
    assert result["exception_type"] == "builtins.TypeError"
    assert result["message"] == "typed failure"


def test_backend_record_does_not_infer_an_unknown_callable():
    result = backend_record(lambda: None)
    assert result["value"] is None
    assert "runtime callable type=" in result["evidence"]


def test_timed_keeps_each_repetition_and_result_count():
    result = timed(lambda: 3, 3)
    assert result["repetitions"] == 3
    assert result["result_count"] == 3
    assert len(result["samples_ns"]) == 3
