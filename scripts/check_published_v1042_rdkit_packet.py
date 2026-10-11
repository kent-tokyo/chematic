#!/usr/bin/env python3
"""Validate a pinned Python/RDKit packet (v1.0.42 defaults retained)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def validate(
    report: dict,
    *,
    expected_schema: str = "published-v1.0.42-rdkit-packet/v1",
    expected_chematic: str = "1.0.42",
    expected_rdkit: tuple[str, ...] = ("2026.9.1", "2026.09.1"),
    expected_revision: str | None = None,
    expected_blocks: int = 21,
    minimum_artifacts: int = 2,
) -> list[str]:
    errors: list[str] = []
    if report.get("schema") != expected_schema:
        errors.append("unexpected schema")
    binding = report.get("binding")
    if binding not in {"python", "node-wasm"}:
        errors.append("unexpected binding")
    versions = report.get("versions", {})
    if versions.get("chematic") != expected_chematic:
        errors.append(f"not a chematic {expected_chematic} artifact")
    if versions.get("rdkit") not in expected_rdkit:
        errors.append(f"not an expected RDKit artifact: {versions.get('rdkit')}")
    if expected_revision is not None:
        revision = report.get("source", {}).get("chematic_revision")
        if revision != expected_revision:
            errors.append(f"unexpected chematic revision: {revision}")
    corpus = report.get("corpus", {})
    if corpus.get("rows") != 10_000 or corpus.get("timing_rows") != 1_000:
        errors.append("unexpected corpus denominator")
    if len(corpus.get("sha256", "")) != 64:
        errors.append("missing corpus SHA-256")
    artifacts = report.get("artifacts", [])
    if len(artifacts) < minimum_artifacts:
        errors.append(
            f"expected at least {minimum_artifacts} artifact records, got {len(artifacts)}"
        )
    for artifact in artifacts:
        if len(artifact.get("sha256", "")) != 64 or artifact.get("bytes", 0) <= 0:
            errors.append(f"invalid artifact record: {artifact.get('file')}")
    method = report.get("method", {})
    blocks = method.get("blocks")
    if blocks != expected_blocks:
        errors.append(f"expected {expected_blocks} alternating blocks")
    accuracy = report.get("accuracy", {})
    speed = report.get("speedup_rdkit_over_chematic", {})
    required_ops = {"labute_asa", "morgan2_chiral"}
    if binding == "python":
        required_ops |= {"tpsa", "num_rings"}
    for op in sorted(required_ops):
        gate = accuracy.get(op)
        if not gate:
            errors.append(f"missing accuracy gate for {op}")
            continue
        if gate.get("compared") != 10_000:
            errors.append(f"wrong accuracy denominator for {op}")
        lane = speed.get(op, {})
        if lane.get("eligible") and not gate.get("equivalent"):
            errors.append(f"speed eligibility bypasses output gate for {op}")
        if not lane.get("eligible") and gate.get("equivalent") and not lane.get("ineligible_reason"):
            errors.append(f"equivalent {op} lacks a speed-ineligibility reason")
        if lane.get("eligible"):
            for name in ("prepared", "pipeline"):
                interval = lane.get(name)
                if not interval or interval.get("blocks") != expected_blocks:
                    errors.append(f"missing {op} {name} interval")
        elif lane.get("prepared") is not None or lane.get("pipeline") is not None:
            errors.append(f"non-equivalent {op} reports a speed interval")
    samples = report.get("raw_samples_seconds", {})
    for engine in ("chematic", "rdkit"):
        if len(samples.get(engine, [])) != expected_blocks:
            errors.append(f"missing raw {engine} samples")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("record", type=Path)
    parser.add_argument("--expected-schema", default="published-v1.0.42-rdkit-packet/v1")
    parser.add_argument("--expected-chematic", default="1.0.42")
    parser.add_argument("--expected-rdkit", action="append", default=None)
    parser.add_argument("--expected-revision", default=None)
    parser.add_argument("--expected-blocks", type=int, default=21)
    parser.add_argument("--minimum-artifacts", type=int, default=2)
    args = parser.parse_args()
    report = json.loads(args.record.read_text())
    errors = validate(
        report,
        expected_schema=args.expected_schema,
        expected_chematic=args.expected_chematic,
        expected_rdkit=tuple(args.expected_rdkit or ("2026.9.1", "2026.09.1")),
        expected_revision=args.expected_revision,
        expected_blocks=args.expected_blocks,
        minimum_artifacts=args.minimum_artifacts,
    )
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        return 1
    print(f"OK: {args.record} ({report['binding']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
