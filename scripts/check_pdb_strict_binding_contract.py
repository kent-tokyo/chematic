#!/usr/bin/env python3
"""Validate the checked-in strict-PDB binding contract and its anchors.

The strict PDB API is intentionally a source-only contract until every
published artifact is regenerated.  This gate still verifies that the shared
fixture, Rust/Python/WASM implementations, generated Node artifact, and
binding tests agree on the fixed-column error boundary.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "validation" / "cross_binding_contract.json"


def fail(message: str) -> int:
    print(f"strict PDB binding contract invalid: {message}", file=sys.stderr)
    return 1


def main() -> int:
    try:
        document = json.loads(FIXTURE.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        return fail(f"cannot read shared fixture: {error}")

    contract = document.get("pdb_strict_contract")
    if not isinstance(contract, dict) or contract.get("schema_version") != 1:
        return fail("pdb_strict_contract schema is missing")
    expected = contract.get("expected")
    if expected != {"atom_count": 1, "coords": [[10.0, 11.0, 12.0]]}:
        return fail("valid fixture expectation changed unexpectedly")
    if contract.get("malformed_error_field") != "z":
        return fail("baseline malformed field must remain z")
    cases = contract.get("fixed_column_cases")
    if not isinstance(cases, list) or [case.get("field") for case in cases] != [
        "serial", "residue sequence", "x", "y", "z"
    ]:
        return fail("fixed-column case list is incomplete or reordered")
    if any(
        not isinstance(case.get("start"), int)
        or not isinstance(case.get("end"), int)
        or not isinstance(case.get("replacement"), str)
        or case["end"] - case["start"] != len(case["replacement"])
        for case in cases
    ):
        return fail("fixed-column case offsets and replacement widths are invalid")

    checks = {
        ROOT / "crates/chematic-3d/src/pdb.rs": r"pub fn parse_pdb_atoms_strict",
        ROOT / "crates/chematic-py/src/formats.rs": r"fn from_pdb_strict",
        ROOT / "crates/chematic-py/tests/test_format_convert.py":
        r"test_pdb_strict_parser_reports_each_fixed_column",
        ROOT / "crates/chematic-wasm/src/mol_io.rs": r"pub fn mol_from_pdb_strict",
        ROOT / "crates/chematic-wasm/tests/source_only_contract.test.mjs":
        r"fixed_column_cases",
        ROOT / "crates/chematic-wasm/pkg-node/chematic_wasm.js":
        r"exports\.mol_from_pdb_strict",
        ROOT / "crates/chematic-wasm/pkg-node/chematic_wasm.d.ts":
        r"export function mol_from_pdb_strict",
    }
    for path, pattern in checks.items():
        try:
            content = path.read_text(encoding="utf-8")
        except OSError as error:
            return fail(f"cannot read {path}: {error}")
        if re.search(pattern, content) is None:
            return fail(f"missing anchor {pattern!r} in {path}")

    print("strict PDB binding contract OK: shared fixture and Rust/Python/Node/WASM anchors present")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
