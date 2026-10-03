#!/usr/bin/env python3
"""Keep A6 quality dimensions separate; historical evidence cannot pass a current gate."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

import check_mmff94_current_energy_evidence as energy_checker
import check_mmff94_published_quality_evidence as published_checker


ROOT = Path(__file__).resolve().parents[1]
PUBLISHED = ROOT / "validation/results/mmff94-public-v1026-20260926"
TYPING = ROOT / "validation/results/mmff94-atom-type-census-v1.0.25-issue637-2026-09-25.json"
ENERGY = (
    ROOT
    / "validation/results/mmff94-same-explicit-h-energy-current-main-v1.0.19-2026-09-23.json"
)


def load_rows(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def published_dimensions(rows: list[dict]) -> dict:
    if len(rows) != 265 or any(row.get("status") != "success" for row in rows):
        raise ValueError("expected 265 successful published MMFF94 rows")
    convergence = Counter(row["force_field"].get("converged") for row in rows)
    if set(convergence) - {True, False}:
        raise ValueError("missing or non-Boolean convergence result")
    missing_by_term = {
        term: sum(len(row["force_field"]["coverage"][f"{term}_missing"]) for row in rows)
        for term in ("bonds", "angles", "torsions", "oop", "stretch_bend")
    }
    return {
        "converged": convergence[True],
        "not_converged": convergence[False],
        "not_converged_indices": [
            row["row_index"] for row in rows if not row["force_field"]["converged"]
        ],
        "missing_terms": missing_by_term,
        "fallbacks": sum(row["force_field"].get("fallback_reason") is not None for row in rows),
        "internal_sound": sum(row["final_validation"].get("sound") is True for row in rows),
        "internal_stereo_ok": sum(row["final_validation"].get("stereo_ok") is True for row in rows),
        "internal_gross_clashes": sum(
            row["final_validation"].get("gross_clash_count", 0) for row in rows
        ),
    }


def build_matrix(target_version: str) -> dict:
    # These check original package/manifest hashes, rows, independent scores,
    # matched seeds, and the two reversed acquisition orders before profiling.
    published_checker.main()
    if energy_checker.main() != 0:
        raise ValueError("historical same-coordinate energy packet failed its checker")
    runs = []
    inputs = {}
    for suffix in ("", "-r2"):
        path = PUBLISHED / f"chematic-3d{suffix}.jsonl"
        runs.append(published_dimensions(load_rows(path)))
        inputs[path.name] = digest(path)
    if runs[0] != runs[1]:
        raise ValueError("A6 quality dimensions changed between reversed run orders")
    typing = json.loads(TYPING.read_text(encoding="utf-8"))
    energy = json.loads(ENERGY.read_text(encoding="utf-8"))
    if (
        typing["chematic_version"] != "1.0.25"
        or typing["heavy"]["differing_atoms"] != 2908
        or energy["versions"]["chematic"] != "1.0.19"
    ):
        raise ValueError("historical typing or energy provenance changed")
    inputs[TYPING.name] = digest(TYPING)
    inputs[ENERGY.name] = digest(ENERGY)
    published_version = "1.0.26"
    run = runs[0]
    dimensions = {
        "typing_and_parameters": {
            "status": "open",
            "published_version": published_version,
            "published_missing_terms": run["missing_terms"],
            "published_fallbacks": run["fallbacks"],
            "separate_census_version": typing["chematic_version"],
            "separate_census_heavy_type_mismatches": typing["heavy"]["differing_atoms"],
            "note": (
                "Strict bond-angle lane permits missing torsions; the 10k "
                "typing census is a different artifact and corpus."
            ),
        },
        "convergence": {
            "status": "open" if run["not_converged"] else "pass_on_measured_artifact",
            "published_version": published_version,
            "converged": run["converged"],
            "not_converged": run["not_converged"],
            "not_converged_indices": run["not_converged_indices"],
        },
        "stereo": {
            "status": "pass_on_measured_artifact",
            "published_version": published_version,
            "independent_scorer_clean": 265,
            "internal_clean": run["internal_stereo_ok"],
        },
        "clash_and_geometry": {
            "status": "pass_on_measured_artifact",
            "published_version": published_version,
            "independent_scorer_sound_and_clash_free": 265,
            "internal_sound": run["internal_sound"],
            "internal_gross_clashes": run["internal_gross_clashes"],
        },
        "same_coordinate_energy": {
            "status": "open",
            "source_version": energy["versions"]["chematic"],
            "comparable_rows": energy["measurements"]["comparable_rows"],
            "within_5_kcal_mol": energy["measurements"]["within_5_kcal_mol"],
            "over_5_kcal_mol_indices": [166, 231],
            "note": "Historical source-built explicit-H coordinates; not the v1.0.26 wheel.",
        },
        "independent_conformer_quality": {
            "status": "not_measured",
            "note": "Between-engine RMSD is not an independent non-inferiority criterion.",
        },
    }
    return {
        "schema": "a6-separated-quality-matrix/v1",
        "target_version": target_version,
        "ready_for_target_version": all(
            item["status"] == "pass_on_measured_artifact"
            and item.get("published_version") == target_version
            for item in dimensions.values()
        ),
        "published_run_rows": 265,
        "published_run_repetitions": 2,
        "input_sha256": inputs,
        "dimensions": dimensions,
        "speed_excluded": True,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-version", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--require-ready", action="store_true")
    args = parser.parse_args()
    result = build_matrix(args.target_version)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "ready_for_target_version": result["ready_for_target_version"],
                "dimensions": {
                    key: value["status"] for key, value in result["dimensions"].items()
                },
            },
            sort_keys=True,
        )
    )
    return 2 if args.require_ready and not result["ready_for_target_version"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
