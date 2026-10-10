#!/usr/bin/env python3
"""Validate a raw published-v1.0.42/RDKit benchmark record."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def validate(report: dict) -> list[str]:
    errors: list[str] = []
    if report.get("schema") != "published-v1.0.42-rdkit-packet/v1":
        errors.append("unexpected schema")
    binding = report.get("binding")
    if binding not in {"python", "node-wasm"}:
        errors.append("unexpected binding")
    versions = report.get("versions", {})
    if versions.get("chematic") != "1.0.42":
        errors.append("not a chematic 1.0.42 artifact")
    if versions.get("rdkit") not in {"2026.9.1", "2026.09.1"}:
        errors.append("not an RDKit 2026.09.1 artifact")
    corpus = report.get("corpus", {})
    if corpus.get("rows") != 10_000 or corpus.get("timing_rows") != 1_000:
        errors.append("unexpected corpus denominator")
    if len(corpus.get("sha256", "")) != 64:
        errors.append("missing corpus SHA-256")
    artifacts = report.get("artifacts", [])
    if len(artifacts) < 2:
        errors.append("published artifact records are missing")
    for artifact in artifacts:
        if len(artifact.get("sha256", "")) != 64 or artifact.get("bytes", 0) <= 0:
            errors.append(f"invalid artifact record: {artifact.get('file')}")
    method = report.get("method", {})
    blocks = method.get("blocks")
    if blocks != 21:
        errors.append("expected 21 alternating blocks")
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
                if not interval or interval.get("blocks") != 21:
                    errors.append(f"missing {op} {name} interval")
        elif lane.get("prepared") is not None or lane.get("pipeline") is not None:
            errors.append(f"non-equivalent {op} reports a speed interval")
    samples = report.get("raw_samples_seconds", {})
    for engine in ("chematic", "rdkit"):
        if len(samples.get(engine, [])) != 21:
            errors.append(f"missing raw {engine} samples")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("record", type=Path)
    args = parser.parse_args()
    report = json.loads(args.record.read_text())
    errors = validate(report)
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        return 1
    print(f"OK: {args.record} ({report['binding']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
