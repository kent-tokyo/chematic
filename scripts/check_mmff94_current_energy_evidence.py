#!/usr/bin/env python3
"""Validate the current-source 265-row MMFF94 same-coordinate packet."""

from __future__ import annotations

import hashlib
import json
import math
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "validation" / "results"
SUMMARY = RESULTS / "mmff94-same-explicit-h-energy-current-main-v1.0.19-2026-09-23.json"
ROWS = RESULTS / "mmff94-same-explicit-h-energy-current-main-v1.0.19-2026-09-23.jsonl"
GRADIENT_SUMMARY = (
    RESULTS / "mmff94-same-explicit-h-gradient-current-main-v1.0.19-2026-09-23.json"
)
GRADIENT_ROWS = (
    RESULTS / "mmff94-same-explicit-h-gradient-current-main-v1.0.19-2026-09-23.jsonl"
)
TERMS = {"bond", "angle", "stretch_bend", "torsion", "oop", "vdw", "electrostatic"}


def require(condition: bool, message: str, errors: list[str]) -> None:
    if not condition:
        errors.append(message)


def main() -> int:
    errors: list[str] = []
    try:
        summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"MMFF94 current-energy evidence invalid: {exc}", file=sys.stderr)
        return 1

    require(summary.get("schema_version") == 2, "schema version changed", errors)
    require(
        summary.get("profile") == "mmff94_same_explicit_h_current_source_v2",
        "profile changed",
        errors,
    )
    require(summary.get("gate_passed") is True, "evidence gate failed", errors)
    require(
        summary.get("missing_dimensions") == [],
        "provenance has missing dimensions",
        errors,
    )
    require(
        summary.get("source", {}).get("git_revision")
        == "8ca05c8b2636e5e1e52ef2c37c3c6c2ce5cff41a",
        "measured source revision changed",
        errors,
    )
    require(
        summary.get("source", {}).get("tracked_tree_dirty") is False,
        "measured source tree was dirty",
        errors,
    )
    require(
        summary.get("versions", {})
        == {
            "chematic": "1.0.19",
            "python": "3.13.6",
            "rdkit_distribution": "2026.3.6",
            "rdkit_runtime": "2026.03.6",
        },
        "version packet changed",
        errors,
    )
    require(
        summary.get("artifacts", {}).get("rdkit", {}).get("sha256")
        == "e16c467cb254a223e59a0cf81358c6b39da15a99d2909170d693e95778fddb41",
        "RDKit wheel hash changed",
        errors,
    )
    require(
        summary.get("artifacts", {}).get("schematic", {}).get("sha256")
        == "e340782b114b23bcfa7365f85182d504fbd9a1e152edcd574de74914614cd6dd",
        "CheMatic extension hash changed",
        errors,
    )

    expected = {
        "comparable_rows": 262,
        "median_abs_delta_kcal_mol": 0.21754491059398795,
        "p90_abs_delta_kcal_mol": 1.14467915306426,
        "max_abs_delta_kcal_mol": 9.874152506372653,
        "within_1_kcal_mol": 232,
        "within_5_kcal_mol": 260,
    }
    measurements = summary.get("measurements", {})
    for key, value in expected.items():
        require(measurements.get(key) == value, f"{key} changed", errors)
    require(
        measurements.get("row_accounting")
        == {
            "input_count": 265,
            "terminal_count": 265,
            "status_counts": {
                "declared_unsupported": 1,
                "embed_failure": 2,
                "ok": 262,
            },
        },
        "row accounting changed",
        errors,
    )

    rows_hash = hashlib.sha256(ROWS.read_bytes()).hexdigest()
    require(
        rows_hash == summary.get("rows_artifact", {}).get("sha256"),
        "rows artifact hash changed",
        errors,
    )
    statuses: Counter[str] = Counter()
    over_five: list[int] = []
    row_count = 0
    with ROWS.open(encoding="utf-8") as handle:
        for expected_index, line in enumerate(handle):
            row = json.loads(line)
            row_count += 1
            require(
                row.get("input_index") == expected_index,
                f"row {expected_index}: input index changed",
                errors,
            )
            statuses[row["status"]] += 1
            if row["status"] != "ok":
                continue
            require(
                isinstance(row.get("coordinate_sha256"), str)
                and len(row["coordinate_sha256"]) == 64,
                f"row {expected_index}: coordinate hash missing",
                errors,
            )
            breakdown = row.get("schematic_energy_breakdown_kcal_mol", {})
            require(
                set(breakdown) == TERMS | {"total"},
                f"row {expected_index}: energy term set changed",
                errors,
            )
            require(
                math.isclose(
                    sum(breakdown[term] for term in TERMS),
                    breakdown["total"],
                    rel_tol=1e-12,
                    abs_tol=1e-10,
                ),
                f"row {expected_index}: energy terms do not sum to total",
                errors,
            )
            if row["abs_delta_kcal_mol"] > 5.0:
                over_five.append(expected_index)
    require(row_count == 265, f"expected 265 rows, found {row_count}", errors)
    require(
        statuses == {"ok": 262, "embed_failure": 2, "declared_unsupported": 1},
        f"raw-row statuses changed: {dict(statuses)}",
        errors,
    )
    require(
        over_five == [166, 231], f">5 kcal/mol residuals changed: {over_five}", errors
    )

    try:
        gradient_summary = json.loads(GRADIENT_SUMMARY.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        errors.append(f"gradient packet unreadable: {exc}")
        gradient_summary = {}
    require(
        gradient_summary.get("source", {}).get("git_revision")
        == "a6540d1a09f41191c97274957453c4c8d05c1c8e",
        "gradient source revision changed",
        errors,
    )
    require(
        gradient_summary.get("artifacts", {}).get("schematic", {}).get("sha256")
        == "b51ce7c7d9790efd145b2194556b655d57ee8cfd5e70ae01a1afe8eef657162a",
        "gradient CheMatic extension hash changed",
        errors,
    )
    gradient_measurement = gradient_summary.get("measurements", {}).get(
        "gradient_diagnostic", {}
    )
    require(
        gradient_measurement.get("input_indices") == [166, 231]
        and gradient_measurement.get("rows") == 2,
        "gradient residual cohort changed",
        errors,
    )
    require(
        0.0 <= gradient_measurement.get("max_scaled_error", math.inf) < 1e-6,
        "gradient scaled error exceeded 1e-6",
        errors,
    )
    require(
        hashlib.sha256(GRADIENT_ROWS.read_bytes()).hexdigest()
        == gradient_summary.get("rows_artifact", {}).get("sha256"),
        "gradient rows artifact hash changed",
        errors,
    )
    gradient_rows = [
        json.loads(line) for line in GRADIENT_ROWS.read_text().splitlines()
    ]
    require(len(gradient_rows) == 265, "gradient packet row count changed", errors)
    require(
        [row["input_index"] for row in gradient_rows if "gradient_diagnostic" in row]
        == [166, 231],
        "raw gradient rows changed",
        errors,
    )

    if errors:
        print("MMFF94 current-energy evidence invalid:", file=sys.stderr)
        print("\n".join(f"- {error}" for error in errors), file=sys.stderr)
        return 1
    print(
        "MMFF94 current-energy evidence OK: 265 terminal rows; 262 comparable; "
        "p90 1.144679 kcal/mol; two >5 kcal/mol residuals; residual-gradient "
        "max scaled error <1e-6"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
