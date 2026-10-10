#!/usr/bin/env python3
"""Validate the checked-in bounded Open Babel semantic and speed records."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
SEMANTIC = ROOT / "benchmarks/2026-10-10-openbabel-file-io-semantics-v1.0.42.json"
SPEED = ROOT / "benchmarks/2026-10-10-openbabel-file-io-cli-roundtrip-v1.0.42.json"
FORMATS = {"v3000", "mol2", "cml", "cdxml"}


def load(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path.name} must contain an object")
    return payload


def validate() -> list[str]:
    errors: list[str] = []
    try:
        semantic = load(SEMANTIC)
        speed = load(SPEED)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        return [str(exc)]

    for name, payload in (("semantic", semantic), ("speed", speed)):
        if payload.get("status") != "local-verified":
            errors.append(f"{name}: status is not local-verified")
        source = payload.get("source")
        if not isinstance(source, dict) or source.get("dirty") is not False:
            errors.append(f"{name}: source must be a clean commit")
        tool = payload.get("tool_versions", {}).get("openbabel", "")
        if "Open Babel 3.2.1" not in tool:
            errors.append(f"{name}: Open Babel 3.2.1 is not pinned")
        formats = payload.get("formats")
        if not isinstance(formats, dict) or set(formats) != FORMATS:
            errors.append(f"{name}: expected exactly {sorted(FORMATS)}")

    for fmt, row in semantic.get("formats", {}).items():
        if row.get("chematic", {}).get("status") != "exact_semantic_match":
            errors.append(f"semantic/{fmt}: chematic is not exact")
        if row.get("result") not in {"chematic", "tie"}:
            errors.append(f"semantic/{fmt}: result is not a win or tie")
    score = semantic.get("score", {})
    if score.get("openbabel_wins") != 0 or score.get("chematic_wins", 0) < 1:
        errors.append("semantic: expected at least one chematic win and zero Open Babel wins")

    if speed.get("alternating_order") is not True or speed.get("blocks", 0) < 21:
        errors.append("speed: alternating order and at least 21 blocks are required")
    for fmt, row in speed.get("formats", {}).items():
        blocks = row.get("blocks", [])
        if len(blocks) < 21:
            errors.append(f"speed/{fmt}: fewer than 21 raw blocks")
        if any(float(block.get("speedup", 0.0)) <= 1.0 for block in blocks):
            errors.append(f"speed/{fmt}: at least one block is not faster")
        if row.get("every_block_faster") is not True:
            errors.append(f"speed/{fmt}: every_block_faster is false")
        if float(row.get("paired_speedup_ci95_lower_bound", 0.0)) <= 1.0:
            errors.append(f"speed/{fmt}: paired lower bound does not exceed 1.0")
        if row.get("status") != "faster":
            errors.append(f"speed/{fmt}: status is not faster")
    return errors


def main() -> int:
    errors = validate()
    if errors:
        print("Open Babel file-I/O evidence failures:", file=sys.stderr)
        print("\n".join(errors), file=sys.stderr)
        return 1
    print("Open Babel file-I/O evidence OK: 4 semantic lanes, 4 paired CLI wins")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
