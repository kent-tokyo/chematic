"""Shared report envelope for the retained same-process format contracts.

Each runner owns its parser, semantic comparison, and malformed-input policy.
This module only records common provenance and writes the resulting report.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from benchmark_version import workspace_version


def make_report(
    root: Path,
    path: Path,
    payload: bytes,
    repeats: int,
    gate: str,
    rows: dict,
    comparison: dict,
    errors: list[str],
    *,
    extra: dict | None = None,
) -> dict:
    from rdkit import rdBase

    report = {
        "schema_version": 1,
        "target_version": workspace_version(root),
        "status": "local-verified" if not errors else "failed",
        "gate": gate,
        "fixture": {
            "path": str(path.relative_to(root)),
            "bytes": len(payload),
            "sha256": hashlib.sha256(payload).hexdigest(),
        },
        "repeats": repeats,
    }
    if extra:
        report.update(extra)
    report.update({
        "rows": rows,
        "comparison": comparison,
        "tool_versions": {"rdkit": rdBase.rdkitVersion},
        "errors": errors,
    })
    return report


def write_report(root: Path, output: Path | None, report: dict) -> int:
    encoded = json.dumps(report, indent=2) + "\n"
    if output:
        target = output if output.is_absolute() else root / output
        target.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    return 1 if report["errors"] else 0
