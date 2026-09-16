import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).parents[2]
SCRIPT = ROOT / "scripts/prepare_sealed_accuracy_cohort.py"
SPEC = importlib.util.spec_from_file_location("sealed_cohort", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def write_identity(path: Path, rows: list[tuple[str, str, str, str]]) -> None:
    keys = ("input_smiles", "canonical_smiles", "parent_smiles", "scaffold_smiles")
    path.write_text(
        "".join(json.dumps(dict(zip(keys, row))) + "\n" for row in rows),
        encoding="utf-8",
    )


def test_preparation_excludes_all_reference_identity_forms_and_never_seals_by_default(tmp_path):
    source = tmp_path / "source.smi"
    source.write_text("A\nB\nC\nD\nE\n", encoding="utf-8")
    identity = tmp_path / "identity.jsonl"
    write_identity(identity, [("A", "a", "pa", "sa"), ("B", "b", "pb", "sb"), ("C", "c", "pc", "sc"), ("D", "d", "pd", "sd"), ("E", "e", "pe", "se")])
    reference = tmp_path / "reference.jsonl"
    write_identity(reference, [("old1", "b", "x", "x"), ("old2", "x", "pc", "x"), ("old3", "y", "y", "sd")])
    metadata = tmp_path / "source.json"
    metadata.write_text(json.dumps({"source_url": "https://example.invalid/source", "source_release": "test", "license": "CC0", "acquired_at": "2026-09-13T00:00:00Z", "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest()}), encoding="utf-8")
    output = tmp_path / "out"
    result = subprocess.run([sys.executable, str(SCRIPT), "--source", str(source), "--source-metadata", str(metadata), "--identity-audit", str(identity), "--reference-identity-audit", str(reference), "--output-dir", str(output), "--development-rows", "1", "--sealed-rows", "1"], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["status"] == "prepared_not_sealed"
    assert manifest["splits"]["development"]["rows"] == 1
    assert manifest["splits"]["sealed_holdout"]["rows"] == 1
    assert manifest["deduplication"]["excluded_rows"] == {"duplicate_canonical": 0, "reference_canonical": 1, "reference_parent": 1, "reference_scaffold": 1}


def test_sealing_requires_freeze_and_unused_attestation(tmp_path):
    source = tmp_path / "source.smi"
    source.write_text("A\nB\n", encoding="utf-8")
    identity = tmp_path / "identity.jsonl"
    write_identity(identity, [("A", "a", "pa", "sa"), ("B", "b", "pb", "sb")])
    metadata = tmp_path / "source.json"
    metadata.write_text(json.dumps({"source_url": "https://example.invalid/source", "source_release": "test", "license": "CC0", "acquired_at": "2026-09-13T00:00:00Z", "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest()}), encoding="utf-8")
    result = subprocess.run([sys.executable, str(SCRIPT), "--source", str(source), "--source-metadata", str(metadata), "--identity-audit", str(identity), "--output-dir", str(tmp_path / "out"), "--development-rows", "1", "--sealed-rows", "1", "--attest-unused"], capture_output=True, text=True)
    assert result.returncode != 0
    assert "requires a candidate commit" in result.stderr


def test_freeze_rejects_a_tag_that_does_not_resolve_to_the_candidate_commit():
    try:
        MODULE.verify_candidate_freeze("HEAD", "definitely-not-a-candidate-tag")
    except ValueError as error:
        assert "does not resolve" in str(error)
    else:
        raise AssertionError("a nonexistent tag must not be accepted as a freeze")


def test_attestation_binds_to_the_exact_source_hash(tmp_path):
    attestation = tmp_path / "attestation.json"
    attestation.write_text(
        json.dumps(
            {
                "attestor": "test maintainer",
                "attested_at": "2026-09-13T00:00:00Z",
                "statement": "This source was unused during development.",
                "source_sha256": "different",
            }
        ),
        encoding="utf-8",
    )
    try:
        MODULE.read_unused_attestation(attestation, "expected")
    except ValueError as error:
        assert "does not match" in str(error)
    else:
        raise AssertionError("an attestation for another source must not seal this cohort")


def test_attestation_must_follow_the_annotated_candidate_tag():
    freeze = {
        "candidate_commit": "0" * 40,
        "candidate_tag": "candidate",
        "candidate_tagged_at": "2026-09-13T18:09:19+09:00",
    }
    try:
        MODULE.verify_attestation_timing(freeze, "2026-09-13T16:50:00+09:00")
    except ValueError as error:
        assert "predates" in str(error)
    else:
        raise AssertionError("an attestation before candidate freeze must not be accepted")

    MODULE.verify_attestation_timing(freeze, "2026-09-13T18:09:19+09:00")
