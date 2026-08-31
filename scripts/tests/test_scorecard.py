import hashlib
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).parents[2]
HARNESS = ROOT / "validation/cosmolkit_comparison"
CORPUS = HARNESS / "smoke_corpus.jsonl"


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
