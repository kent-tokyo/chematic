#!/usr/bin/env python3
"""Validate the immutable-lane contract used for RDKit rebaselining."""

from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "validation" / "rdkit_rebaseline_manifest.json"
OPERATION_IDS = {"smiles_parse_write", "cip", "smarts", "morgan", "browser_runtime"}
AVAILABILITY = {"measured", "partial", "unavailable"}


def main() -> int:
    try:
        document = json.loads(PATH.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"RDKit rebaseline manifest invalid: {exc}", file=sys.stderr)
        return 1

    errors: list[str] = []
    if document.get("schema_version") != 1:
        errors.append("schema_version must be 1")
    policy = document.get("policy")
    if not isinstance(policy, dict) or policy.get("baseline_lanes_are_immutable") is not True:
        errors.append("policy must make baseline lanes immutable")
    operations = document.get("operations")
    operation_ids = {entry.get("id") for entry in operations if isinstance(entry, dict)} if isinstance(operations, list) else set()
    if operation_ids != OPERATION_IDS:
        errors.append("operations must declare exactly the required operation ids")

    lanes = document.get("lanes")
    if not isinstance(lanes, list) or not lanes:
        errors.append("lanes must be a non-empty array")
        lanes = []
    ids: set[str] = set()
    has_pending = False
    for lane in lanes:
        if not isinstance(lane, dict):
            errors.append("every lane must be an object")
            continue
        required = {"id", "channel", "availability", "rdkit_version", "wrapper_backend", "corpus_policy", "operation_evidence"}
        if set(lane) != required:
            errors.append("every lane must contain exactly the required keys")
            continue
        lane_id = lane["id"]
        if not isinstance(lane_id, str) or not lane_id or lane_id in ids:
            errors.append("lane ids must be unique non-empty strings")
        ids.add(lane_id)
        availability = lane["availability"]
        if availability not in AVAILABILITY:
            errors.append(f"{lane_id}: invalid availability")
        evidence = lane["operation_evidence"]
        if not isinstance(evidence, dict) or set(evidence) != OPERATION_IDS:
            errors.append(f"{lane_id}: operation_evidence must cover every operation")
            continue
        evidence_count = 0
        for operation, paths in evidence.items():
            if not isinstance(paths, list):
                errors.append(f"{lane_id}/{operation}: evidence must be an array")
                continue
            evidence_count += len(paths)
            for relative in paths:
                if not isinstance(relative, str) or not (ROOT / relative).is_file():
                    errors.append(f"{lane_id}/{operation}: missing evidence {relative}")
        if availability == "measured" and evidence_count == 0:
            errors.append(f"{lane_id}: measured lane requires evidence")
        if availability == "unavailable":
            has_pending = True
            if lane["rdkit_version"] is not None or lane["wrapper_backend"] is not None or evidence_count:
                errors.append(f"{lane_id}: unavailable lane must not claim an artifact or evidence")
    if not has_pending:
        errors.append("a next-stable unavailable lane is required until every artifact is verified")

    if errors:
        print("RDKit rebaseline manifest invalid:", file=sys.stderr)
        print("\n".join(f"- {error}" for error in errors), file=sys.stderr)
        return 1
    print(f"RDKit rebaseline manifest OK: {len(lanes)} lanes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
