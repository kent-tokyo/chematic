#!/usr/bin/env python3
"""Validate a JSON Parent report without third-party dependencies."""

import argparse
import json
from pathlib import Path

STAGES = ["fragment", "charge", "isotope", "stereo", "tautomer"]


def validate(path: Path) -> list[str]:
    try:
        report = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        return [f"invalid JSON: {exc}"]
    errors = []
    if not isinstance(report, dict):
        return ["report must be an object"]
    if not isinstance(report.get("smiles"), str) or not report["smiles"]:
        errors.append("smiles must be a non-empty string")
    if not isinstance(report.get("status"), str) or not report["status"]:
        errors.append("status must be a non-empty string")
    stages = report.get("stages")
    if not isinstance(stages, list) or [s.get("name") for s in stages if isinstance(s, dict)] != STAGES:
        errors.append(f"stages must follow the fixed order {STAGES}")
    else:
        for index, stage in enumerate(stages):
            if not isinstance(stage.get("smiles"), str) or not stage["smiles"]:
                errors.append(f"stages[{index}].smiles must be a non-empty string")
    unknown = set(report) - {"smiles", "status", "stages"}
    if unknown:
        errors.append(f"unknown report fields: {sorted(unknown)}")
    return errors


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("report", type=Path)
    args = parser.parse_args()
    errors = validate(args.report)
    if errors:
        print(json.dumps({"valid": False, "errors": errors}, indent=2, sort_keys=True))
        raise SystemExit(1)
    print(json.dumps({"valid": True}, sort_keys=True))


if __name__ == "__main__":
    main()
