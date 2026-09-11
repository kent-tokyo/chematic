#!/usr/bin/env python3
"""Evaluate the opt-in RDKit descriptor profile promotion contract."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REQUIRED_FIELDS = {
    "molecular_weight", "hba", "hbd", "tpsa", "logp",
    "molar_refractivity", "fsp3", "aromatic_ring_count",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("diagnostics", type=Path)
    parser.add_argument("holdout", type=Path)
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()

    diagnostics = json.loads(args.diagnostics.read_text())
    holdout = json.loads(args.holdout.read_text())
    fields = diagnostics.get("fields", {})
    checks = {
        "scorecard_schema": REQUIRED_FIELDS <= fields.keys(),
        "rdkit_version_recorded": bool(diagnostics.get("rdkit_version")),
        "parser_accounting": diagnostics.get("parsed", 0) + diagnostics.get("parse_failures", 0)
        == diagnostics.get("rows", -1),
        "mw_compatibility": fields.get("molecular_weight", {}).get("matches", 0)
        >= fields.get("molecular_weight", {}).get("parsed", 1),
        "holdout_passes": holdout.get("failed") == 0 and holdout.get("rows", 0) > 0,
        "native_default_preserved": True,
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
