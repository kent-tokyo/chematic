import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).parents[2]
HARNESS = ROOT / "validation/cosmolkit_comparison"
CORPUS = HARNESS / "smoke_corpus.jsonl"

sys.path.insert(0, str(HARNESS))
import validate as comparison_validate


def result_file(tmp_path, engine, *, mismatch=False):
    corpus_hash = hashlib.sha256(CORPUS.read_bytes()).hexdigest()
    rows = []
    for line in CORPUS.read_text().splitlines():
        row = json.loads(line)
        rows.append({
            "schema_version": 1,
            "engine": engine,
            "engine_version": "test",
            "source_commit": "deadbeef",
            "corpus_sha256": corpus_hash,
            "id": row["id"],
            "smiles": row["smiles"],
            "status": "ok",
            "operations": {
                "formula": {"status": "ok", "value": "same" if not mismatch or row["id"] != "ethanol" else "different"},
            },
        })
    path = tmp_path / f"{engine}.jsonl"
    path.write_text("".join(json.dumps(row) + "\n" for row in rows))
    return path


def test_scorecard_separates_statuses_and_mismatches(tmp_path):
    reference = result_file(tmp_path, "rdkit")
    candidate = result_file(tmp_path, "chematic", mismatch=True)
    output = tmp_path / "scorecard.json"
    completed = subprocess.run(
        [
            sys.executable,
            str(HARNESS / "scorecard.py"),
            "--result", f"rdkit={reference}",
            "--result", f"chematic={candidate}",
            "--reference", "rdkit",
            "--output", str(output),
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    assert json.loads(output.read_text())["operations"]["formula"]["against_reference"]["chematic"] == {
        "match": 9,
        "mismatch": 1,
    }
    assert json.loads(completed.stdout)["valid"] is True


def test_validator_rejects_manifest_drift(tmp_path):
    corpus = tmp_path / "smoke_corpus.jsonl"
    corpus.write_text(CORPUS.read_text().replace('"CCO"', '"CCN"', 1))
    manifest = tmp_path / "corpus_manifest.json"
    manifest.write_text((HARNESS / "corpus_manifest.json").read_text())
    result = result_file(tmp_path, "chematic")
    old_corpus, old_manifest = comparison_validate.CORPUS, comparison_validate.MANIFEST
    comparison_validate.CORPUS, comparison_validate.MANIFEST = corpus, manifest
    try:
        errors = comparison_validate.validate(result)
    finally:
        comparison_validate.CORPUS, comparison_validate.MANIFEST = old_corpus, old_manifest
    assert any("manifest sha256" in error for error in errors)


def test_external_runner_rejects_engine_mismatch_without_output(tmp_path):
    adapter = tmp_path / "adapter.py"
    adapter.write_text(
        "import json, sys\n"
        "for line in open(sys.argv[sys.argv.index('--corpus') + 1]):\n"
        " row = json.loads(line); print(json.dumps({'schema_version': 1, 'engine': 'wrong', 'engine_version': 'test', 'source_commit': None, 'corpus_sha256': '07e0f5f7a4ac2743b54d6c7d4fa63a5c1e0f6278c0637a80fe9a0dc85ba145ec', 'id': row['id'], 'smiles': row['smiles'], 'status': 'ok', 'operations': {}}))\n"
    )
    output = tmp_path / "result.jsonl"
    completed = subprocess.run(
        [sys.executable, str(HARNESS / "run_external.py"), "--engine", "cosmolkit",
         "--adapter", f"{sys.executable} {adapter}", "--output", str(output)],
        cwd=ROOT, capture_output=True, text=True,
    )
    assert completed.returncode != 0
    assert "engine does not match" in completed.stderr
    assert not output.exists()
