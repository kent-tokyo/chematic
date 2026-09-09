#!/usr/bin/env python3
"""Generate the deterministic, repository-local compatibility dashboard.

This intentionally reads only checked-in JSON. It performs no network access,
does not run an oracle, and emits no current-time field, so a clean checkout
can reproduce the same dashboard byte-for-byte.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from benchmark_version import workspace_version


ROOT = Path(__file__).resolve().parent.parent
CONTRACT = ROOT / "validation" / "cross_binding_contract.json"
DEFAULT_OUTPUT = ROOT / "docs" / "compatibility-dashboard.md"


def current_streaming_matrix() -> Path:
    """Select the checked-in matrix for the current workspace version."""
    version = workspace_version(ROOT)
    candidates = sorted(ROOT.glob(f"benchmarks/*-streaming-cross-engine-matrix-v{version}.json"))
    if not candidates:
        raise FileNotFoundError(f"no streaming matrix for workspace version {version}")
    return candidates[-1]


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def binding_rows(contract: dict) -> list[tuple[str, int]]:
    sections = (
        "fixtures",
        "descriptor_contract",
        "standardization_contract",
        "fingerprint_contract",
        "fingerprint_detail_contract",
        "reaction_application_contract",
        "adversarial",
    )
    rows = []
    for section in sections:
        value = contract[section]
        if isinstance(value, list):
            count = len(value)
        else:
            count = sum(len(item) for item in value.values() if isinstance(item, list))
        rows.append((section, count))
    return rows


def render(contract: dict, streaming: dict, streaming_path: Path) -> str:
    contract_path = "validation/cross_binding_contract.json"
    streaming_display_path = streaming_path.relative_to(ROOT).as_posix()
    lines = [
        "# Compatibility dashboard",
        "",
        "Generated from checked-in manifests by `python3 scripts/generate_compatibility_dashboard.py`.",
        "This is a compatibility-contract dashboard, not a universal RDKit parity or speed claim.",
        "",
        f"- Target version: `{streaming['target_version']}`",
        "- Regeneration: deterministic, offline, clean-checkout compatible",
        f"- Contract manifest: `{contract_path}` (SHA-256 `{digest(CONTRACT)}`)",
        f"- Streaming matrix: `{streaming_display_path}` (SHA-256 `{digest(streaming_path)}`)",
        "",
        "## Shared binding contract",
        "",
        f"Operation inventory: `{len(contract['operation_manifest']['operations'])}` currently shared operations; validate with `python3 scripts/check_cross_binding_manifest.py`.",
        "",
        "| Area | Checked-in assertions | Status |",
        "|---|---:|---|",
    ]
    for name, count in binding_rows(contract):
        lines.append(f"| `{name}` | {count} | covered by versioned fixture contract |")

    lines += [
        "",
        "## Streaming record/failure contract",
        "",
        f"Pinned matrix: {streaming['repeats']} repetitions across {len(streaming['rows'])} formats.",
        "The engine/process boundaries remain explicit; records and failures are the only cross-engine claims.",
        "",
        "| Format | Expected records | Engines with zero failures |",
        "|---|---:|---:|",
    ]
    for row in streaming["rows"]:
        zero_failure_engines = sum(item["failures"] == 0 for item in row["engines"])
        lines.append(f"| `{row['format']}` | {row['expected_records']} | {zero_failure_engines}/{len(row['engines'])} |")

    lines += [
        "",
        "## Boundaries",
        "",
        "- A `covered` row means the checked-in assertion and its local consumer tests exist; it does not mean every edge case is supported.",
        "- Missing optional engines, external review, publication, and broader corpus quality remain separate release gates.",
        "- Regenerate after changing either source manifest and review the resulting digest changes before committing.",
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output = args.output if args.output.is_absolute() else ROOT / args.output
    streaming_path = current_streaming_matrix()
    output.write_text(render(load(CONTRACT), load(streaming_path), streaming_path), encoding="utf-8")
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
