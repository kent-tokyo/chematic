#!/usr/bin/env python3
"""Check native descriptor output against the frozen cross-binding fixture.

This gate is intentionally independent from the RDKit-compatible descriptor
profile. It protects the existing native defaults from being changed while
compatibility work is developed behind named opt-in APIs.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


FIELDS = ("mw", "tpsa", "hbd", "hba", "heavy_atoms")
FLOAT_FIELDS = {"mw", "tpsa"}


def load_jsonl(path: Path) -> list[dict]:
    rows: list[dict] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        value = json.loads(line)
        if not isinstance(value, dict):
            raise ValueError(f"{path}:{line_number}: record must be an object")
        rows.append(value)
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("report", type=Path, help="JSONL from descriptor_binding_dump")
    parser.add_argument("contract", type=Path, help="cross_binding_contract.json")
    parser.add_argument("--json", type=Path, help="write a machine-readable gate result")
    args = parser.parse_args()

    report = load_jsonl(args.report)
    document = json.loads(args.contract.read_text(encoding="utf-8"))
    fixtures = document["descriptor_contract"]["fixtures"]
    expected = {row["smiles"]: row for row in fixtures}
    errors: list[str] = []
    seen: set[str] = set()
    for row in report:
        smiles = row.get("smiles")
        if not isinstance(smiles, str) or smiles not in expected:
            errors.append(f"unexpected or missing fixture identity: {smiles!r}")
            continue
        if smiles in seen:
            errors.append(f"duplicate fixture: {smiles}")
        seen.add(smiles)
        if row.get("status") != "ok":
            errors.append(f"{smiles}: status={row.get('status')!r}")
            continue
        descriptors = row.get("descriptors")
        if not isinstance(descriptors, dict):
            errors.append(f"{smiles}: descriptors must be an object")
            continue
        for field in FIELDS:
            actual = descriptors.get(field)
            if field == "mw":
                expected_value = expected[smiles]["molecular_weight"]
            elif field == "tpsa":
                expected_value = expected[smiles]["tpsa"]
            else:
                expected_value = expected[smiles][field]
            if not isinstance(actual, (int, float)) or isinstance(actual, bool):
                errors.append(f"{smiles} {field}: non-numeric output")
            elif not math.isfinite(float(actual)):
                errors.append(f"{smiles} {field}: non-finite output")
            elif field in FLOAT_FIELDS and abs(float(actual) - expected_value) > 1e-6:
                errors.append(f"{smiles} {field}: {actual} != {expected_value}")
            elif field not in FLOAT_FIELDS and actual != expected_value:
                errors.append(f"{smiles} {field}: {actual} != {expected_value}")

    missing = set(expected) - seen
    errors.extend(f"missing fixture: {smiles}" for smiles in sorted(missing))
    result = {
        "schema_version": 1,
        "profile": "native",
        "fixtures": len(expected),
        "matched": len(expected) - len(missing),
        "failed": len(errors),
        "gate_passed": not errors,
        "report": str(args.report),
        "contract": str(args.contract),
    }
    if args.json:
        args.json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    if errors:
        print("native descriptor regression: FAIL")
        print("\n".join(f"- {error}" for error in errors))
        return 1
    print(f"native descriptor regression: PASS ({len(expected)} fixed fixtures)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
