#!/usr/bin/env python3
"""Summarize common-coordinate force-field agreement and sampled differences."""

import argparse
import json
import math
from pathlib import Path


def percentile(values, fraction):
    if not values:
        return None
    ordered = sorted(values)
    return ordered[math.ceil(fraction * len(ordered)) - 1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = json.loads(args.input.read_text())
    rows = report["results"]
    summary = {
        "rows": len(rows),
        "row_errors": sum(bool(row.get("error")) for row in rows),
        "fields": {},
    }
    for variant in ("MMFF94", "MMFF94s", "UFF"):
        fields = [row.get("fields", {}).get(variant, {}) for row in rows]
        comparable = [
            field
            for field in fields
            if field.get("rdkit_supported")
            and field.get("chematic_supported")
            and "energy_absolute_error" in field
        ]
        energies = [field["energy_absolute_error"] for field in comparable]
        gradients = [
            value for field in comparable for value in field["gradient_absolute_errors"]
        ]
        differences = [
            value
            for field in comparable
            for value in field["finite_difference_absolute_errors"]
        ]
        summary["fields"][variant] = {
            "evaluated": sum(bool(field) for field in fields),
            "rdkit_supported": sum(
                bool(field.get("rdkit_supported")) for field in fields
            ),
            "chematic_supported": sum(
                bool(field.get("chematic_supported")) for field in fields
            ),
            "comparable": len(comparable),
            "energy_absolute_error_p95": percentile(energies, 0.95),
            "energy_absolute_error_max": max(energies, default=None),
            "molecules_energy_error_over_1e-4": sum(value > 1e-4 for value in energies),
            "gradient_absolute_error_p95": percentile(gradients, 0.95),
            "gradient_absolute_error_max": max(gradients, default=None),
            "molecules_gradient_error_over_1e-3": sum(
                max(field["gradient_absolute_errors"], default=0.0) > 1e-3
                for field in comparable
            ),
            "finite_difference_absolute_error_max": max(differences, default=None),
        }
    args.output.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
