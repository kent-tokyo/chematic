#!/usr/bin/env python3
"""Validate runnable RDKit rebaseline commands and explicit gaps."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "validation" / "rdkit_rebaseline_execution.json"
PLACEHOLDER = re.compile(r"\$\{([A-Z][A-Z0-9_]*)\}")
DIFFERENCE_CLASSES = {
    "chematic_regression",
    "oracle_change",
    "contract_difference",
    "unresolved",
}


def main() -> int:
    try:
        document = json.loads(MANIFEST.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"RDKit execution manifest invalid: {exc}", file=sys.stderr)
        return 1
    errors: list[str] = []
    if document.get("schema_version") != 1:
        errors.append("schema_version must be 1")
    if document.get("sealed_accuracy_cohort_reused") is not False:
        errors.append("sealed accuracy data must not be reused")
    if set(document.get("difference_classes", [])) != DIFFERENCE_CLASSES:
        errors.append("difference_classes must contain the four adjudication outcomes")
    required = set(document.get("required_operation_ids", []))
    if not required:
        errors.append("required_operation_ids must be non-empty")
    lanes = document.get("lanes")
    if not isinstance(lanes, list) or not lanes:
        errors.append("lanes must be a non-empty array")
        lanes = []
    lane_ids: set[str] = set()
    command_ids: set[str] = set()
    for lane in lanes:
        if not isinstance(lane, dict):
            errors.append("every lane must be an object")
            continue
        lane_id = lane.get("id")
        if not isinstance(lane_id, str) or not lane_id or lane_id in lane_ids:
            errors.append("lane ids must be unique non-empty strings")
            continue
        lane_ids.add(lane_id)
        corpus = lane.get("corpus")
        if not isinstance(corpus, str) or not (ROOT / corpus).is_file():
            errors.append(f"{lane_id}: corpus must exist")
        covered: set[str] = set()
        commands = lane.get("commands")
        if not isinstance(commands, list) or not commands:
            errors.append(f"{lane_id}: commands must be a non-empty array")
            commands = []
        for command in commands:
            if not isinstance(command, dict):
                errors.append(f"{lane_id}: command must be an object")
                continue
            command_id = command.get("id")
            if not isinstance(command_id, str) or not command_id or command_id in command_ids:
                errors.append(f"{lane_id}: command ids must be globally unique")
                continue
            command_ids.add(command_id)
            covers = command.get("covers")
            argv = command.get("argv")
            if not isinstance(covers, list) or not set(covers) <= required:
                errors.append(f"{command_id}: covers contains an unknown operation")
            else:
                covered.update(covers)
            if not isinstance(argv, list) or len(argv) < 2 or not all(
                isinstance(token, str) and token for token in argv
            ):
                errors.append(f"{command_id}: argv must be a non-empty string array")
                continue
            script_tokens = [token for token in argv if token.startswith("scripts/")]
            if not script_tokens or not all((ROOT / token).is_file() for token in script_tokens):
                errors.append(f"{command_id}: referenced script does not exist")
            for token in argv:
                for placeholder in PLACEHOLDER.findall(token):
                    if placeholder in {"SEALED_CORPUS", "SEALED_8K"}:
                        errors.append(f"{command_id}: sealed-data placeholder is prohibited")
        unavailable = lane.get("unavailable_operations")
        if not isinstance(unavailable, dict) or not all(
            isinstance(key, str) and isinstance(value, str) and value
            for key, value in unavailable.items()
        ):
            errors.append(f"{lane_id}: unavailable_operations must name reasons")
            unavailable = {}
        unavailable_ids = set(unavailable)
        if not unavailable_ids <= required:
            errors.append(f"{lane_id}: unknown unavailable operation")
        if covered & unavailable_ids:
            errors.append(f"{lane_id}: operation cannot be both covered and unavailable")
        if covered | unavailable_ids != required:
            missing = sorted(required - covered - unavailable_ids)
            errors.append(f"{lane_id}: unaccounted operations {missing}")
    if errors:
        print("RDKit execution manifest invalid:", file=sys.stderr)
        print("\n".join(f"- {error}" for error in errors), file=sys.stderr)
        return 1
    print(f"RDKit execution manifest OK: {len(lanes)} lanes, {len(command_ids)} commands")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
