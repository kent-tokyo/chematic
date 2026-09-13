#!/usr/bin/env python3
"""Validate the combined evidence packet for the three RDKit search gates."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULTS = (
    ROOT / "validation/results/rdkit-search-oracle-current-rerun.json",
    ROOT / "validation/results/rdkit-search-threshold-gate-v1.0.13.json",
    ROOT / "validation/results/rdkit-search-cross-binding-parity-v1.0.13.json",
)
REQUIRED_PROVENANCE = {
    "source_commit",
    "tracked_diff_sha256",
    "untracked_source_sha256",
    "lockfile_sha256",
    "python_executable",
    "platform",
    "machine",
}
SHA256 = re.compile(r"^[0-9a-f]{64}$")
COMMIT = re.compile(r"^[0-9a-f]{40}$")


def fail(message: str) -> int:
    print(f"RDKit search evidence invalid: {message}", file=sys.stderr)
    return 1


def load(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"{path}: cannot read JSON: {error}") from error
    if not isinstance(value, dict):
        raise ValueError(f"{path}: top-level JSON value must be an object")
    return value


def mismatch_values(value: object):
    if isinstance(value, dict):
        if "mismatches" in value:
            yield value["mismatches"]
        for child in value.values():
            yield from mismatch_values(child)
    elif isinstance(value, list):
        for child in value:
            yield from mismatch_values(child)


def validate(paths: tuple[Path, Path, Path]) -> list[str]:
    errors: list[str] = []
    reports = []
    for path in paths:
        try:
            report = load(path)
        except ValueError as error:
            errors.append(str(error))
            continue
        reports.append((path, report))
        if report.get("gate_passed") is not True:
            errors.append(f"{path}: gate_passed is not true")
        provenance = report.get("provenance")
        if not isinstance(provenance, dict) or not REQUIRED_PROVENANCE <= provenance.keys():
            errors.append(f"{path}: provenance is incomplete")
        elif not COMMIT.fullmatch(str(provenance["source_commit"])):
            errors.append(f"{path}: source_commit is not a full git SHA")
        else:
            for key in ("tracked_diff_sha256", "untracked_source_sha256"):
                if not SHA256.fullmatch(str(provenance[key])):
                    errors.append(f"{path}: {key} is not a SHA-256 digest")
        for value in mismatch_values(report):
            if value != 0:
                errors.append(f"{path}: mismatch count is {value!r}, expected 0")
    if len(reports) == len(paths):
        corpus_ids = {(r.get("corpus", {}).get("path"), r.get("corpus", {}).get("sha256")) for _, r in reports}
        if len(corpus_ids) != 1:
            errors.append("reports do not use the same corpus identity")
        provenance_ids = {
            (
                r["provenance"].get("source_commit"),
                r["provenance"].get("tracked_diff_sha256"),
                r["provenance"].get("untracked_source_sha256"),
            )
            for _, r in reports
            if isinstance(r.get("provenance"), dict)
        }
        if len(provenance_ids) != 1:
            errors.append("reports do not share the same source/diff provenance")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("reports", nargs="*", type=Path, default=list(DEFAULTS))
    args = parser.parse_args()
    if len(args.reports) != 3:
        return fail("exactly three reports are required: oracle, threshold, cross-binding")
    errors = validate(tuple(path.resolve() for path in args.reports))
    if errors:
        print("RDKit search evidence: BLOCKED")
        print("\n".join(f"- {error}" for error in errors))
        return 1
    print("RDKit search evidence: PASS (three reports share corpus, provenance, and zero mismatches)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
