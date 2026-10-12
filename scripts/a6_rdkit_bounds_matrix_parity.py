#!/usr/bin/env python3
"""Compare CheMatic and RDKit ETKDGv3 bounds on the frozen A6 corpus."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import chematic
import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import rdDistGeom

from public_package_3d_chematic import load_rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tiers", choices=("A", "B", "AB"), default="AB")
    parser.add_argument("--count", type=int)
    parser.add_argument("--no-triangle-smoothing", action="store_true")
    args = parser.parse_args()

    source_rows = load_rows(args.tiers)
    if args.count is not None:
        source_rows = source_rows[: args.count]

    exact_molecules = 0
    exact_cells = 0
    tolerance = 1.0e-12
    matching_molecules = 0
    matching_cells = 0
    total_cells = 0
    maximum_absolute_delta = 0.0
    failures: list[dict[str, object]] = []
    residuals: list[dict[str, object]] = []

    for index, row in enumerate(source_rows):
        smiles = str(row["smiles"])
        try:
            ours = np.asarray(
                chematic.from_smiles(smiles)
                .add_hydrogens()
                .rdkit_bounds_matrix(
                    use_macrocycle_14_config=True,
                    do_triangle_smoothing=not args.no_triangle_smoothing,
                ),
                dtype=np.float64,
            )
            reference_molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
            reference = np.asarray(
                rdDistGeom.GetMoleculeBoundsMatrix(
                    reference_molecule,
                    set15bounds=True,
                    scaleVDW=False,
                    doTriangleSmoothing=not args.no_triangle_smoothing,
                    useMacrocycle14config=True,
                ),
                dtype=np.float64,
            )
        except Exception as error:  # noqa: BLE001 - record both engines' boundary failures
            failures.append({"row_index": index, "id": row.get("id"), "error": str(error)})
            continue

        if ours.shape != reference.shape:
            failures.append(
                {
                    "row_index": index,
                    "id": row.get("id"),
                    "error": f"shape mismatch: {ours.shape} != {reference.shape}",
                }
            )
            continue

        equal = ours == reference
        close = np.abs(ours - reference) <= tolerance
        molecule_exact_cells = int(equal.sum())
        molecule_matching_cells = int(close.sum())
        molecule_total_cells = int(equal.size)
        exact_cells += molecule_exact_cells
        matching_cells += molecule_matching_cells
        total_cells += molecule_total_cells
        delta = float(np.max(np.abs(ours - reference), initial=0.0))
        maximum_absolute_delta = max(maximum_absolute_delta, delta)
        if molecule_exact_cells == molecule_total_cells:
            exact_molecules += 1
        if molecule_matching_cells == molecule_total_cells:
            matching_molecules += 1
        else:
            residuals.append(
                {
                    "row_index": index,
                    "id": row.get("id"),
                    "exact_cells": molecule_exact_cells,
                    "matching_cells": molecule_matching_cells,
                    "total_cells": molecule_total_cells,
                    "max_absolute_delta": delta,
                }
            )

    payload = {
        "schema": "chematic.a6-rdkit-bounds-parity.v1",
        "rdkit_version": rdBase.rdkitVersion,
        "tiers": args.tiers,
        "triangle_smoothing": not args.no_triangle_smoothing,
        "molecules": len(source_rows),
        "exact_molecules": exact_molecules,
        "exact_cells": exact_cells,
        "tolerance": tolerance,
        "matching_molecules": matching_molecules,
        "matching_cells": matching_cells,
        "total_cells": total_cells,
        "maximum_absolute_delta": maximum_absolute_delta,
        "failures": failures,
        "residuals": residuals,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in payload.items() if key not in {"residuals"}}, indent=2))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
