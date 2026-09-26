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
ISSUE637 = "v1.0.25-issue637-2026-09-25"
ISSUE637_BEFORE = "v1.0.25-before-issue637-2026-09-25"


def load_packet(stem: str) -> tuple[dict, list[dict]]:
    summary_path = RESULTS / f"{stem}.json"
    rows_path = RESULTS / f"{stem}.jsonl"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    rows = [json.loads(line) for line in rows_path.read_text(encoding="utf-8").splitlines()]
    if hashlib.sha256(rows_path.read_bytes()).hexdigest() != summary.get("rows_artifact", {}).get("sha256"):
        raise ValueError(f"{rows_path.name}: rows artifact hash changed")
    return summary, rows


def recount(rows: list[dict]) -> dict:
    ok = [row for row in rows if row["status"] == "ok"]
    deltas = sorted(row["abs_delta_kcal_mol"] for row in ok)
    per_term = {
        term: max(abs(row["term_delta_kcal_mol"][term]) for row in ok) for term in TERMS
    }
    return {
        "statuses": dict(Counter(row["status"] for row in rows)),
        "comparable_rows": len(ok),
        "max_abs_delta_kcal_mol": deltas[-1],
        "p90_abs_delta_kcal_mol": deltas[math.ceil(0.9 * len(deltas)) - 1],
        "within_1_kcal_mol": sum(value <= 1.0 for value in deltas),
        "over_5": [row["input_index"] for row in ok if row["abs_delta_kcal_mol"] > 5.0],
        "per_term_max": per_term,
        "term_sum_residual": max(
            abs(sum(row["rdkit_term_energies_kcal_mol"].values()) - row["rdkit_energy_kcal_mol"])
            for row in ok
        ),
        "coordinates": [row.get("coordinate_sha256") for row in rows],
    }


def validate_issue637(errors: list[str]) -> None:
    """#637: per-term RDKit oracle before/after the MMFF94 typing fixes."""
    try:
        before_summary, before_rows = load_packet(
            f"mmff94-same-explicit-h-energy-per-term-{ISSUE637_BEFORE}"
        )
        after_summary, after_rows = load_packet(
            f"mmff94-same-explicit-h-energy-per-term-{ISSUE637}"
        )
        gradient_summary, gradient_rows = load_packet(
            f"mmff94-same-explicit-h-gradient-delta1e-6-{ISSUE637}"
        )
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        errors.append(f"#637 packet unreadable: {exc}")
        return
    for name, summary, revision in (
        ("before", before_summary, "049c12e3dfa24c6142ce145462f6a39bd69be37c"),
        ("after", after_summary, "13d70a2e514f8db8abab3fa2c1a158405f7fc1f8"),
        ("gradient", gradient_summary, "13d70a2e514f8db8abab3fa2c1a158405f7fc1f8"),
    ):
        require(summary.get("schema_version") == 3, f"#637 {name}: schema changed", errors)
        require(summary.get("gate_passed") is True, f"#637 {name}: gate failed", errors)
        require(
            summary.get("source", {}).get("git_revision") == revision
            and summary.get("source", {}).get("tracked_tree_dirty") is False,
            f"#637 {name}: source revision changed or dirty",
            errors,
        )
    before = recount(before_rows)
    after = recount(after_rows)
    statuses = {"ok": 262, "embed_failure": 2, "declared_unsupported": 1}
    require(before["statuses"] == statuses and after["statuses"] == statuses, "#637 statuses changed", errors)
    require(
        before["coordinates"] == after["coordinates"],
        "#637 before/after rows were not measured on identical coordinates",
        errors,
    )
    require(before["over_5"] == [166, 231], "#637 baseline >5 kcal/mol rows changed", errors)
    require(after["over_5"] == [], "#637 candidate still has >5 kcal/mol rows", errors)
    require(after["within_1_kcal_mol"] == 262, "#637 candidate rows outside 1 kcal/mol", errors)
    require(after["max_abs_delta_kcal_mol"] < 0.5, "#637 candidate max delta >= 0.5 kcal/mol", errors)
    for term in ("bond", "electrostatic", "stretch_bend"):
        require(after["per_term_max"][term] < 0.01, f"#637 {term} term not at parity", errors)
    require(
        max(before["term_sum_residual"], after["term_sum_residual"]) < 1e-6,
        "#637 RDKit per-term energies do not sum to the RDKit total",
        errors,
    )
    for key in ("max_abs_delta_kcal_mol", "p90_abs_delta_kcal_mol", "within_1_kcal_mol", "comparable_rows"):
        require(
            after_summary.get("measurements", {}).get(key) == after[key],
            f"#637 summary {key} does not match its rows",
            errors,
        )
    gradients = {
        row["input_index"]: row["gradient_diagnostic"]
        for row in after_rows
        if "gradient_diagnostic" in row
    }
    require(sorted(gradients) == [166, 178, 231], "#637 gradient cohort changed", errors)
    require(
        all(gradients[i]["max_scaled_error"] < 1e-6 for i in (166, 231)),
        "#637 residual-row gradient scaled error >= 1e-6",
        errors,
    )
    fine = [row["gradient_diagnostic"] for row in gradient_rows if "gradient_diagnostic" in row]
    require(
        len(fine) == 1
        and fine[0]["central_difference_delta_angstrom"] == 1e-6
        and fine[0]["max_scaled_error"] < 1e-6,
        "#637 row 178 fine-step gradient check changed",
        errors,
    )
    census_before = json.loads(
        (RESULTS / f"mmff94-atom-type-census-{ISSUE637_BEFORE}.json").read_text(encoding="utf-8")
    )
    census_after = json.loads(
        (RESULTS / f"mmff94-atom-type-census-{ISSUE637}.json").read_text(encoding="utf-8")
    )
    require(
        census_before["hydrogen"]["differing_atoms"] == 3101
        and census_after["hydrogen"]["differing_atoms"] == 56
        and census_after["hydrogen"]["differing_atoms_with_agreeing_parent_type"] == 0
        and census_before["heavy"]["differing_atoms"] == census_after["heavy"]["differing_atoms"],
        "#637 atom-type census changed",
        errors,
    )


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

    validate_issue637(errors)

    if errors:
        print("MMFF94 current-energy evidence invalid:", file=sys.stderr)
        print("\n".join(f"- {error}" for error in errors), file=sys.stderr)
        return 1
    print(
        "MMFF94 current-energy evidence OK: 265 terminal rows; 262 comparable; "
        "p90 1.144679 kcal/mol; two >5 kcal/mol residuals; residual-gradient "
        "max scaled error <1e-6; #637 per-term oracle: 262/262 within 1 kcal/mol, "
        "bond/electrostatic/stretch-bend at parity, H typing 3,101 -> 56"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
