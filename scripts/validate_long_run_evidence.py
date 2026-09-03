#!/usr/bin/env python3
"""Validate the checked-in v1.0 long-run execution evidence manifest."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "validation" / "manifests" / "v1.0.0-long-run-evidence.json"
REVISION_RE = re.compile(r"^[0-9a-f]{8,40}$")


def fail(message: str) -> int:
    print(f"long-run evidence invalid: {message}", file=sys.stderr)
    return 1


def main() -> int:
    try:
        evidence = json.loads(MANIFEST.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return fail(str(exc))

    if evidence.get("schema") != "chematic.v1.long-run-evidence.v1":
        return fail("unexpected schema")
    if evidence.get("release") != "1.0.0-rc.1":
        return fail("unexpected release")
    revision = evidence.get("candidate_revision")
    if not isinstance(revision, str) or not REVISION_RE.fullmatch(revision):
        return fail("candidate_revision must be a hexadecimal git revision")
    if not isinstance(evidence.get("runs"), list):
        return fail("runs must be an array")

    runs = {run.get("name"): run for run in evidence["runs"] if isinstance(run, dict)}
    required = {
        "experimental_3d": "cargo test -p chematic-3d --lib -- --ignored",
        "canonical_corpus_nci_and_descriptor": (
            "cargo test -p chematic-chem --test "
            "canonical_idempotency_corpus_standardized -- --ignored"
        ),
        "canonical_corpus_chembl": (
            "cargo test -p chematic-chem --test "
            "canonical_idempotency_corpus_standardized -- --ignored "
            "chembl_accuracy_corpus_standardized_is_canonically_idempotent"
        ),
    }
    for name, command in required.items():
        run = runs.get(name)
        if run is None:
            return fail(f"missing run {name}")
        if run.get("status") != "pass":
            return fail(f"run {name} is not pass")
        if run.get("command") != command:
            return fail(f"run {name} command drifted")
        if run.get("failed") != 0:
            return fail(f"run {name} reports failures")

    if runs["experimental_3d"].get("passed") != 9:
        return fail("experimental_3d must record 9 passing tests")
    if runs["canonical_corpus_chembl"].get("passed") != 1:
        return fail("canonical_corpus_chembl must record one passing test")
    for name in required:
        timeout = runs[name].get("timeout_seconds")
        if timeout != 1800:
            return fail(f"run {name} must retain the 1800 second timeout policy")

    attempt = evidence.get("initial_parallel_corpus_attempt")
    if not isinstance(attempt, dict) or attempt.get("status") != "timeout":
        return fail("initial parallel corpus attempt must remain recorded as timeout")
    if attempt.get("exit_code") != 124:
        return fail("initial parallel corpus timeout must use exit code 124")

    print(
        "Long-run evidence verified: 3D 9 passed; "
        "NCI/descriptor passed; ChEMBL serialized rerun passed"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
