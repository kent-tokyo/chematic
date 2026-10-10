#!/usr/bin/env python3
"""Run the current source record-accounting baseline against Open Babel.

This is the first layer of the production file-I/O gate. It intentionally does
not turn CLI timings into a parser or writer speed claim; semantic round-trip
and equivalent-work performance are separate follow-up lanes.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

from benchmark_version import workspace_version


ROOT = Path(__file__).resolve().parents[1]


def run_report(command: list[str]) -> dict[str, object]:
    completed = subprocess.run(command, cwd=ROOT, check=False, text=True, capture_output=True)
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout).strip()
        raise RuntimeError(f"runner failed ({' '.join(command)}): {detail}")
    report = json.loads(completed.stdout)
    if not isinstance(report, dict):
        raise RuntimeError(f"runner returned a non-object: {' '.join(command)}")
    return report


def normalize_report_paths(report: dict[str, object]) -> None:
    fixture = report.get("fixture")
    rows = report.get("rows")
    if not isinstance(fixture, dict) or not isinstance(rows, dict):
        return
    fixture_path = fixture.get("path")
    if not isinstance(fixture_path, str):
        return
    for row in rows.values():
        if isinstance(row, dict) and "path" in row:
            row["path"] = fixture_path


def source_provenance() -> dict[str, object]:
    revision = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        check=True,
        text=True,
        capture_output=True,
    ).stdout.strip()
    status = subprocess.run(
        ["git", "status", "--porcelain", "--untracked-files=normal"],
        cwd=ROOT,
        check=True,
        text=True,
        capture_output=True,
    ).stdout
    return {"commit": revision, "dirty": bool(status.strip())}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--openbabel", default="obabel")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")

    reports: dict[str, dict[str, object]] = {}
    for fmt in ("sdf", "mol", "v3000", "mol2"):
        reports[fmt] = run_report(
            [
                "python3",
                "scripts/check_streaming_cross_engine.py",
                "--format",
                fmt,
                "--repeats",
                str(args.repeats),
                "--openbabel",
                args.openbabel,
            ]
        )
    for fmt in ("cml", "cdxml", "mmcif", "pdb"):
        reports[fmt] = run_report(
            [
                "python3",
                "scripts/check_streaming_cml_openbabel.py",
                "--format",
                fmt,
                "--repeats",
                str(args.repeats),
                "--openbabel",
                args.openbabel,
            ]
        )

    version = workspace_version(ROOT)
    errors: list[str] = []
    for fmt, report in reports.items():
        normalize_report_paths(report)
        if report.get("target_version") != version:
            errors.append(f"{fmt}: target_version does not match workspace {version}")
        rows = report.get("rows")
        if not isinstance(rows, dict):
            errors.append(f"{fmt}: rows are missing")
            continue
        for engine in ("chematic", "openbabel"):
            row = rows.get(engine)
            if not isinstance(row, dict):
                errors.append(f"{fmt}: {engine} row is missing")
                continue
            if row.get("failures") != 0:
                errors.append(f"{fmt}: {engine} reported failures")
            if row.get("records") != report.get("expected_records"):
                errors.append(f"{fmt}: {engine} record count mismatch")

    summary = {
        "schema_version": 1,
        "target_version": version,
        "source": source_provenance(),
        "status": "local-verified" if not errors else "failed",
        "gate": "openbabel_file_io_record_accounting_v1",
        "contract": "validation/openbabel_file_io_contract_v1.json",
        "repeats": args.repeats,
        "formats": reports,
        "claim_boundary": "record/failure accounting and fixture provenance only; no semantic or speed superiority claim",
        "errors": errors,
    }
    encoded = json.dumps(summary, indent=2) + "\n"
    if args.output:
        target = args.output if args.output.is_absolute() else ROOT / args.output
        target.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
