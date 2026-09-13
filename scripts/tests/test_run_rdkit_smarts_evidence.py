import json

from scripts.run_rdkit_smarts_evidence import promote_completed_dump


def row(identifier: str):
    return {"record_type": "molecule", "id": identifier, "smiles": "CC", "patterns": {"[C]": {}}}


def test_runner_adds_distinct_verified_footer_for_legacy_producer(tmp_path):
    corpus = tmp_path / "corpus.smi"
    corpus.write_text("CC\n", encoding="utf-8")
    raw = tmp_path / "raw.jsonl"
    raw.write_text(json.dumps(row("hand")) + "\n" + json.dumps(row("corpus_0")) + "\n", encoding="utf-8")
    output = tmp_path / "output.jsonl"
    promote_completed_dump(raw, output, corpus)
    records = [json.loads(line) for line in output.read_text(encoding="utf-8").splitlines()]
    assert records[-1]["footer_origin"] == "runner_verified"
    assert records[-1]["input_rows"] == 1
    assert records[-1]["emitted_molecule_rows"] == 2


def test_runner_rejects_legacy_output_with_missing_corpus_rows(tmp_path):
    corpus = tmp_path / "corpus.smi"
    corpus.write_text("CC\nCCC\n", encoding="utf-8")
    raw = tmp_path / "raw.jsonl"
    raw.write_text(json.dumps(row("corpus_0")) + "\n", encoding="utf-8")
    try:
        promote_completed_dump(raw, tmp_path / "output.jsonl", corpus)
    except ValueError as exc:
        assert "corpus row accounting" in str(exc)
    else:
        raise AssertionError("missing corpus row was accepted")
