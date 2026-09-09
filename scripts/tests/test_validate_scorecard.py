import importlib.util
from pathlib import Path


ROOT = Path(__file__).parents[2]
SPEC = importlib.util.spec_from_file_location("validate_scorecard", ROOT / "scripts/validate_scorecard.py")
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def scorecard():
    return {
        "schema_version": 1,
        "target_version": "1.0.7",
        "corpus_sha256": "0" * 64,
        "corpus_records": 1,
        "configuration": {"profile": "default"},
        "engines": {
            "chematic": {"engine_version": "1.0.7", "source_commits": ["abcdef1"]},
        },
        "operations": {
            "mw": {"status_counts": {"chematic": {"ok": 1, "unsupported": 2}}},
        },
        "claims": [{"operation": "mw", "engine": "chematic", "status": "ok"}],
    }


def test_accepts_provenance_and_explicit_non_ok_rows():
    assert MODULE.validate(scorecard(), "1.0.7") == []


def test_rejects_stale_version_and_non_positive_claim():
    document = scorecard()
    document["target_version"] = "1.0.5"
    document["claims"][0]["status"] = "unsupported"
    errors = MODULE.validate(document, "1.0.7")
    assert any("target_version" in error for error in errors)
    assert any("non-positive status" in error for error in errors)


def test_rejects_claim_that_has_no_positive_engine_row():
    document = scorecard()
    document["operations"]["mw"]["status_counts"]["chematic"] = {
        "ok": 0,
        "unsupported": 2,
    }
    errors = MODULE.validate(document, "1.0.7")
    assert any("no 'ok' row" in error for error in errors)


def test_rejects_missing_corpus_and_configuration_metadata():
    document = scorecard()
    del document["corpus_records"]
    del document["configuration"]
    errors = MODULE.validate(document, "1.0.7")
    assert any("corpus_records" in error for error in errors)
    assert any("configuration" in error for error in errors)
