import hashlib
import json
import subprocess
import sys

from pathlib import Path


ROOT = Path(__file__).parents[2]
VALIDATOR = ROOT / "scripts/check_rdkit_accuracy_manifest.py"


def make_manifest(tmp_path: Path) -> tuple[Path, Path]:
    corpus = tmp_path / "corpus.smi"
    corpus.write_text("CCO\nCCN\n", encoding="utf-8")
    artifact = tmp_path / "artifact.json"
    artifact.write_text("{}\n", encoding="utf-8")
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema_version": 2,
                "protocol": "rdkit-accuracy-v2",
                "comparator": {"engine": "RDKit", "version": "test"},
                "operations": [
                    {
                        "id": "test",
                        "profile": "rdkit_compatibility",
                        "corpus": {
                            "path": str(corpus),
                            "sha256": digest(corpus),
                            "expected_rows": 2,
                        },
                        "split": "development",
                        "tolerances": {
                            "molecular_weight": 0,
                            "hba": 0,
                            "hbd": 0,
                            "tpsa": 0,
                            "logp": 0,
                            "molar_refractivity": 0,
                            "fsp3": 0,
                            "aromatic_ring_count": 0,
                        },
                        "artifacts": [{"role": "raw", "path": str(artifact), "sha256": digest(artifact)}],
                    }
                ],
                "acceptance": {"sealed_evaluation_required": True, "candidate_freeze_required": True},
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return manifest, corpus


def run_validator(manifest: Path, *extra: str):
    return subprocess.run(
        [sys.executable, str(VALIDATOR), str(manifest), *extra],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )


def test_manifest_accepts_matching_corpus_digest_and_row_count(tmp_path):
    manifest, _ = make_manifest(tmp_path)
    result = run_validator(manifest)
    assert result.returncode == 0, result.stderr


def test_manifest_rejects_corpus_row_count_drift(tmp_path):
    manifest, corpus = make_manifest(tmp_path)
    corpus.write_text("CCO\n", encoding="utf-8")
    result = run_validator(manifest)
    assert result.returncode != 0
    assert "corpus row count" in result.stderr


def test_manifest_rejects_corpus_digest_drift(tmp_path):
    manifest, corpus = make_manifest(tmp_path)
    corpus.write_text("CCO\nCCC\n", encoding="utf-8")
    result = run_validator(manifest)
    assert result.returncode != 0
    assert "corpus digest mismatch" in result.stderr


def test_manifest_rejects_empty_artifact(tmp_path):
    manifest, _ = make_manifest(tmp_path)
    artifact = tmp_path / "artifact.json"
    artifact.write_bytes(b"")
    result = run_validator(manifest)
    assert result.returncode != 0
    assert "artifact is empty" in result.stderr


def test_release_validation_rejects_development_manifest(tmp_path):
    manifest, _ = make_manifest(tmp_path)
    result = run_validator(manifest, "--require-sealed")
    assert result.returncode != 0
    assert "status=sealed" in result.stderr


def test_release_validation_requires_candidate_freeze_and_sealed_split(tmp_path):
    manifest, _ = make_manifest(tmp_path)
    data = json.loads(manifest.read_text(encoding="utf-8"))
    data["status"] = "sealed"
    manifest.write_text(json.dumps(data), encoding="utf-8")
    result = run_validator(manifest, "--require-sealed")
    assert result.returncode != 0
    assert "candidate_freeze" in result.stderr


def test_release_validation_accepts_frozen_manifest(tmp_path):
    manifest, _ = make_manifest(tmp_path)
    data = json.loads(manifest.read_text(encoding="utf-8"))
    data["status"] = "sealed"
    data["candidate_freeze"] = {
        "candidate_commit": "0123456789abcdef0123456789abcdef01234567",
        "candidate_tag": "v-test",
    }
    data["operations"][0]["split"] = "sealed_holdout"
    manifest.write_text(json.dumps(data), encoding="utf-8")
    result = run_validator(manifest, "--require-sealed")
    assert result.returncode == 0, result.stderr
