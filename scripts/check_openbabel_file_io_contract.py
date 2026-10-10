#!/usr/bin/env python3
"""Validate the bounded Open Babel file-I/O competition contract."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "validation" / "openbabel_file_io_contract_v1.json"
REQUIRED_FORMATS = {
    "sdf",
    "mol_v2000",
    "mol_v3000",
    "mol2",
    "cml",
    "cdxml",
    "pdb",
    "mmcif",
}
REQUIRED_ACCOUNTING = {
    "success",
    "typed_refusal",
    "unsupported",
    "invalid_input",
    "internal_error",
}
REQUIRED_LANES = {"cold_start", "parse", "write", "round_trip", "peak_rss"}


def validate_contract(payload: dict[str, object]) -> list[str]:
    errors: list[str] = []
    if payload.get("schema_version") != 1:
        errors.append("schema_version must be 1")
    if payload.get("claim_scope") != "selected production formats; not total format-count parity":
        errors.append("claim_scope must prohibit total format-count parity claims")
    if set(payload.get("accounting_categories", [])) != REQUIRED_ACCOUNTING:
        errors.append("accounting_categories must use the complete five-way outcome partition")

    formats = payload.get("formats")
    if not isinstance(formats, list):
        return [*errors, "formats must be an array"]
    ids = [row.get("id") for row in formats if isinstance(row, dict)]
    if len(ids) != len(set(ids)):
        errors.append("format ids must be unique")
    if set(ids) != REQUIRED_FORMATS:
        errors.append(f"format ids must be exactly {sorted(REQUIRED_FORMATS)}")
    for row in formats:
        if not isinstance(row, dict):
            errors.append("each format row must be an object")
            continue
        fmt = row.get("id", "unknown")
        if row.get("tier") != "A":
            errors.append(f"{fmt}: tier must be A")
        if row.get("read") is not True or row.get("write") is not True:
            errors.append(f"{fmt}: both read and write must be in scope")
        semantics = row.get("semantics")
        if not isinstance(semantics, list) or not semantics:
            errors.append(f"{fmt}: semantic preservation fields are missing")
        fixture = ROOT / str(row.get("fixture", ""))
        if not fixture.is_file():
            errors.append(f"{fmt}: fixture does not exist: {fixture}")
            continue
        digest = hashlib.sha256(fixture.read_bytes()).hexdigest()
        if digest != row.get("fixture_sha256"):
            errors.append(f"{fmt}: fixture SHA-256 mismatch")

    criteria = payload.get("win_criteria")
    if not isinstance(criteria, dict):
        return [*errors, "win_criteria must be an object"]
    correctness = criteria.get("correctness", {})
    reliability = criteria.get("reliability", {})
    performance = criteria.get("performance", {})
    if not isinstance(correctness, dict) or correctness.get("silent_loss_rows") != 0:
        errors.append("correctness must require zero silent-loss rows")
    if not isinstance(correctness, dict) or correctness.get("wrong_confident_rows") != 0:
        errors.append("correctness must require zero wrong-confident rows")
    if not isinstance(reliability, dict) or reliability.get("internal_errors") != 0:
        errors.append("reliability must require zero internal errors")
    if not isinstance(reliability, dict) or reliability.get("panics_or_crashes") != 0:
        errors.append("reliability must require zero panics or crashes")
    if not isinstance(performance, dict) or performance.get("alternating_blocks_minimum", 0) < 21:
        errors.append("performance must require at least 21 alternating blocks")
    if not isinstance(performance, dict) or performance.get("require_equivalent_output") is not True:
        errors.append("performance wins must require equivalent output")
    if not isinstance(performance, dict) or set(performance.get("separate_lanes", [])) != REQUIRED_LANES:
        errors.append("performance lanes must separate startup, parse, write, round-trip, and memory")
    return errors


def main() -> int:
    try:
        payload = json.loads(CONTRACT.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"Open Babel file-I/O contract check failed: {exc}", file=sys.stderr)
        return 1
    errors = validate_contract(payload)
    if errors:
        print("Open Babel file-I/O contract failures:", file=sys.stderr)
        print("\n".join(errors), file=sys.stderr)
        return 1
    print("Open Babel file-I/O contract OK: 8 Tier-A formats, pinned fixtures, strict win criteria")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
