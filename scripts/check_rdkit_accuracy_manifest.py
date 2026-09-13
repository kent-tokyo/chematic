#!/usr/bin/env python3
"""Validate the versioned RDKit accuracy contract and artifact digests."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIELDS = {
    "molecular_weight", "hba", "hbd", "tpsa", "logp",
    "molar_refractivity", "fsp3", "aromatic_ring_count",
}


def fail(message: str) -> int:
    print(f"RDKit accuracy manifest invalid: {message}", file=sys.stderr)
    return 1


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def corpus_row_count(path: Path) -> int:
    """Count logical input rows for the checked-in line-oriented corpora."""
    return sum(
        1
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not (path.suffix.lower() in {".smi", ".smiles"} and line.lstrip().startswith("#"))
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path, nargs="?", default=ROOT / "validation/manifests/rdkit_accuracy_v2.json")
    parser.add_argument("--skip-digests", action="store_true")
    parser.add_argument(
        "--require-sealed",
        action="store_true",
        help="require a sealed evaluation manifest and an explicit candidate freeze",
    )
    args = parser.parse_args()
    try:
        data = json.loads(args.manifest.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        return fail(f"cannot read manifest: {error}")
    if not isinstance(data, dict) or data.get("schema_version") != 2:
        return fail("schema_version must be 2")
    if data.get("protocol") != "rdkit-accuracy-v2":
        return fail("protocol is missing or unexpected")
    comparator = data.get("comparator")
    if not isinstance(comparator, dict) or comparator.get("engine") != "RDKit" or not comparator.get("version"):
        return fail("comparator engine/version pin is required")
    operations = data.get("operations")
    if not isinstance(operations, list) or not operations:
        return fail("operations must be non-empty")
    seen_ids: set[str] = set()
    for operation in operations:
        if not isinstance(operation, dict):
            return fail("every operation must be an object")
        operation_id = operation.get("id")
        if not isinstance(operation_id, str) or not operation_id or operation_id in seen_ids:
            return fail("operation ids must be unique non-empty strings")
        seen_ids.add(operation_id)
        if operation.get("profile") not in {
            "native",
            "rdkit_compatibility",
            "binding_consistency",
            "independent_correctness",
        }:
            return fail(
                f"{operation_id}: profile must identify native, rdkit_compatibility, "
                "binding_consistency, or independent_correctness"
            )
        corpus = operation.get("corpus")
        if not isinstance(corpus, dict) or not isinstance(corpus.get("path"), str):
            return fail(f"{operation_id}: corpus path is required")
        if not isinstance(corpus.get("expected_rows"), int) or isinstance(corpus["expected_rows"], bool) or corpus["expected_rows"] <= 0:
            return fail(f"{operation_id}: expected_rows must be positive")
        corpus_path = ROOT / corpus["path"]
        if not corpus_path.is_file():
            return fail(f"{operation_id}: corpus does not exist: {corpus['path']}")
        try:
            actual_rows = corpus_row_count(corpus_path)
        except (OSError, UnicodeError) as error:
            return fail(f"{operation_id}: cannot read corpus: {error}")
        if actual_rows != corpus["expected_rows"]:
            return fail(
                f"{operation_id}: corpus row count {actual_rows} != expected {corpus['expected_rows']}"
            )
        corpus_digest = corpus.get("sha256")
        if not isinstance(corpus_digest, str) or (
            len(corpus_digest) != 64
            and corpus_digest != "not_frozen_in_v2_development_manifest"
        ):
            return fail(f"{operation_id}: corpus sha256 is invalid")
        if (
            not args.skip_digests
            and corpus_digest != "not_frozen_in_v2_development_manifest"
            and sha256(corpus_path) != corpus_digest
        ):
            return fail(f"{operation_id}: corpus digest mismatch: {corpus['path']}")
        tolerances = operation.get("tolerances")
        if set(tolerances or {}) != FIELDS or any(not isinstance(value, (int, float)) or isinstance(value, bool) or value < 0 for value in tolerances.values()):
            return fail(f"{operation_id}: all eight non-negative tolerances are required")
        if not isinstance(operation.get("split"), str) or not operation["split"]:
            return fail(f"{operation_id}: split is required")
        artifacts = operation.get("artifacts")
        if not isinstance(artifacts, list) or not artifacts:
            return fail(f"{operation_id}: artifact digests are required")
        for artifact in artifacts:
            if not isinstance(artifact, dict) or not isinstance(artifact.get("role"), str) or not isinstance(artifact.get("path"), str):
                return fail(f"{operation_id}: malformed artifact entry")
            digest = artifact.get("sha256")
            if not isinstance(digest, str) or (len(digest) != 64 and digest != "not_frozen_in_v2_development_manifest"):
                return fail(f"{operation_id}: artifact digest is invalid")
            path = ROOT / artifact["path"]
            if not path.is_file():
                return fail(f"{operation_id}: artifact does not exist: {artifact['path']}")
            if path.stat().st_size == 0:
                return fail(f"{operation_id}: artifact is empty: {artifact['path']}")
            if not args.skip_digests and digest != "not_frozen_in_v2_development_manifest" and sha256(path) != digest:
                return fail(f"{operation_id}: artifact digest mismatch: {artifact['path']}")
    acceptance = data.get("acceptance")
    if not isinstance(acceptance, dict) or acceptance.get("sealed_evaluation_required") is not True or acceptance.get("candidate_freeze_required") is not True:
        return fail("sealed evaluation and candidate freeze requirements are mandatory")
    if args.require_sealed:
        if data.get("status") != "sealed":
            return fail("release validation requires status=sealed")
        freeze = data.get("candidate_freeze")
        if not isinstance(freeze, dict):
            return fail("release validation requires candidate_freeze metadata")
        commit = freeze.get("candidate_commit")
        if not isinstance(commit, str) or len(commit) != 40 or any(c not in "0123456789abcdef" for c in commit):
            return fail("candidate_freeze.candidate_commit must be a 40-character lowercase git SHA")
        if not isinstance(freeze.get("candidate_tag"), str) or not freeze["candidate_tag"]:
            return fail("candidate_freeze.candidate_tag is required")
        for operation in operations:
            if not str(operation.get("split", "")).startswith("sealed"):
                return fail(f"{operation['id']}: release validation requires a sealed split")
            corpus_digest = operation["corpus"]["sha256"]
            if not isinstance(corpus_digest, str) or len(corpus_digest) != 64:
                return fail(f"{operation['id']}: sealed corpus must have a frozen SHA-256")
            for artifact in operation["artifacts"]:
                digest = artifact["sha256"]
                if not isinstance(digest, str) or len(digest) != 64:
                    return fail(f"{operation['id']}: sealed artifact must have a frozen SHA-256")
    print(f"RDKit accuracy manifest OK: schema v2, {len(operations)} operations, comparator RDKit {comparator['version']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
