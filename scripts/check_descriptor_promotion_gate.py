#!/usr/bin/env python3
"""Evaluate the opt-in RDKit descriptor profile promotion contract."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


REQUIRED_FIELDS = {
    "molecular_weight", "hba", "hbd", "tpsa", "logp",
    "molar_refractivity", "fsp3", "aromatic_ring_count",
}


def nonnegative_int(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value >= 0


def strict_field_checks(fields: dict, rows: int) -> bool:
    if set(fields) != REQUIRED_FIELDS:
        return False
    return all(
        isinstance(record, dict)
        and nonnegative_int(record.get("parsed"))
        and nonnegative_int(record.get("matches"))
        and nonnegative_int(record.get("strict_matches"))
        and nonnegative_int(record.get("mismatches"))
        and record["parsed"] == rows
        and record["matches"] == rows
        and record["strict_matches"] == rows
        and record["mismatches"] == 0
        for record in fields.values()
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("diagnostics", type=Path)
    parser.add_argument("holdout", type=Path)
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()

    diagnostics = json.loads(args.diagnostics.read_text())
    holdout = json.loads(args.holdout.read_text())
    fields = diagnostics.get("fields", {})
    rows = diagnostics.get("rows")
    parser_accounting = (
        nonnegative_int(rows)
        and nonnegative_int(diagnostics.get("parsed"))
        and nonnegative_int(diagnostics.get("parse_failures"))
        and diagnostics["parsed"] + diagnostics["parse_failures"] == rows
    )
    checks = {
        "scorecard_schema": isinstance(fields, dict) and set(fields) == REQUIRED_FIELDS,
        "rdkit_version_recorded": bool(diagnostics.get("rdkit_version")),
        "parser_accounting": parser_accounting,
        "mw_compatibility": fields.get("molecular_weight", {}).get("matches", 0)
            >= fields.get("molecular_weight", {}).get("parsed", 1),
        "all_fields_strict": strict_field_checks(fields, rows) if nonnegative_int(rows) else False,
        "finite_error_metrics": all(
            isinstance(record.get(key), (int, float))
            and not isinstance(record.get(key), bool)
            and math.isfinite(record[key])
            for record in fields.values()
            for key in ("mae", "median_abs_error", "p95_abs_error", "max_abs_error")
        ),
        "holdout_passes": holdout.get("failed") == 0 and holdout.get("rows", 0) > 0,
        "native_default_preserved": diagnostics.get("native_default_preserved") is True,
        "unclassified_mismatches_absent": diagnostics.get("mismatch_causes", {}).get("unclassified", 0) == 0,
    }
    result = {
        "schema_version": 1,
        "profile": diagnostics.get("profile"),
        "decision": "adopted_opt_in" if all(checks.values()) else "blocked",
        "default_profile_changed": False,
        "checks": checks,
        "known_residuals": diagnostics.get("mismatch_causes", {}),
        "evidence": {
            "diagnostics": str(args.diagnostics),
            "holdout": str(args.holdout),
        },
    }
    args.json.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    if result["decision"] == "blocked":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
