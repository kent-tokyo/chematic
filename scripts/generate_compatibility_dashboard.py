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
ACCURACY_MANIFEST = ROOT / "validation" / "manifests" / "rdkit_accuracy_v2.json"
DEFAULT_OUTPUT = ROOT / "docs" / "compatibility-dashboard.md"


def current_streaming_matrix() -> Path:
    """Select the checked-in matrix for the current workspace version."""
    version = workspace_version(ROOT)
    # Current matrix refreshes may be stored under validation/results when the
    # record is a machine-readable gate rather than a dated benchmark note.
    # Prefer the canonical validation result, then retain compatibility with
    # older dated benchmark snapshots.
    candidates = [
        ROOT / "validation" / "results" / f"cross-engine-matrix-v{version}.json",
        *sorted(ROOT.glob(f"benchmarks/*-streaming-cross-engine-matrix-v{version}.json")),
    ]
    candidates = [path for path in candidates if path.is_file()]
    if not candidates:
        raise FileNotFoundError(f"no streaming matrix for workspace version {version}")
    return candidates[0]


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


def render(
    contract: dict,
    streaming: dict,
    streaming_path: Path,
    accuracy_manifest: dict | None = None,
) -> str:
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
    if accuracy_manifest is not None:
        comparator = accuracy_manifest.get("comparator", {})
        operations = accuracy_manifest.get("operations", [])
        lines += [
            "## RDKit accuracy profiles",
            "",
            "This section is generated from `validation/manifests/rdkit_accuracy_v2.json`. "
            "It records declared evidence lanes; it does not convert development or exposed data into a sealed evaluation.",
            "",
            f"- Manifest status: `{accuracy_manifest.get('status', 'unknown')}`",
            f"- Comparator: `{comparator.get('engine', 'unknown')} {comparator.get('version', 'unknown')}`",
            "",
            "| Operation | Scope | Split | Expected rows | Profile |",
            "|---|---|---|---:|---|",
        ]
        for operation in operations:
            corpus = operation.get("corpus", {})
            lines.append(
                f"| `{operation.get('id', 'unknown')}` | `{operation.get('scope', 'unknown')}` | "
                f"`{operation.get('split', 'unknown')}` | {corpus.get('expected_rows', '—')} | "
                f"`{operation.get('profile', 'not_recorded')}` |"
            )
        lines += [
            "",
            "A missing profile is an explicit contract gap: native, RDKit-compatible, "
            "and binding-consistency lanes must not be conflated in a scorecard.",
            "",
        ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output = args.output if args.output.is_absolute() else ROOT / args.output
    streaming_path = current_streaming_matrix()
    accuracy_manifest = load(ACCURACY_MANIFEST) if ACCURACY_MANIFEST.is_file() else None
    output.write_text(
        render(load(CONTRACT), load(streaming_path), streaming_path, accuracy_manifest),
        encoding="utf-8",
    )
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
