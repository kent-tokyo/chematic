#!/usr/bin/env python3
"""Reproducible supported-profile comparison against RDKit 2026.09.1.

Run with candidate chematic and pinned RDKit on PYTHONPATH. Synthetic inputs
are original project data; no upstream source or issue attachment is copied.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import random
from pathlib import Path

import numpy as np
import rdkit
from rdkit.Chem import rdMIF as reference
from chematic.rdkit_compat import rdMIF as candidate


def validate() -> dict:
    if rdkit.__version__ != "2026.09.1":
        raise RuntimeError("this comparison requires RDKit 2026.09.1")
    rng = random.Random(9601)
    cases = []
    comparisons = 0
    max_absolute = 0.0
    max_scaled = 0.0
    boundary_comparisons = 0
    boundary_differences = 0
    boundary_examples = []
    for size in (0, 1, 2, 7, 32, 129):
        charges = [rng.uniform(-2, 2) for _ in range(size)]
        positions = [[rng.uniform(-5, 5) for _ in range(3)] for _ in range(size)]
        queries = [[rng.uniform(-8, 8) for _ in range(3)] for _ in range(18)]
        queries += positions[:3]
        thresholds = [(x, False) for x in (0.0, 1e-9, 1.0, 4.0, 25.0, 400.0, 1e200)]
        # Include computed boundary distances and the adjacent floats.
        for point in queries[:2]:
            if positions:
                dx, dy, dz = (a - b for a, b in zip(point, positions[0]))
                d2 = dx * dx + dy * dy + dz * dz
                thresholds += [
                    (x, True)
                    for x in (math.nextafter(d2, 0), d2, math.nextafter(d2, math.inf))
                ]
        for order in (list(range(size)), list(reversed(range(size)))):
            q = [charges[i] for i in order]
            xyz = [positions[i] for i in order]
            # Copies, tuples, and non-contiguous arrays share the same values.
            q_storage = np.zeros(size * 2, dtype=np.float64)
            q_storage[::2] = q
            xyz_storage = np.zeros((size, 6), dtype=np.float64)
            xyz_storage[:, ::2] = np.asarray(xyz).reshape(size, 3)
            for form, cq, cp in (
                ("list", q, xyz),
                ("tuple", tuple(q), tuple(tuple(row) for row in xyz)),
                ("strided_numpy", q_storage[::2], xyz_storage[:, ::2]),
            ):
                for alpha, cutoff in (
                    (0.0, 0.001),
                    (0.0, 0.2),
                    (0.0, 1.0),
                    (1e-6, 1.0),
                    (0.3, 1.0),
                    (4.0, 2.0),
                ):
                    for probe, absolute in (
                        (1.0, False),
                        (-0.75, False),
                        (2.0, True),
                        (0.0, False),
                    ):
                        ours = candidate.Coulomb(cq, cp, probe, absolute, alpha, cutoff)
                        # Reference constructor receives ordinary lists; array
                        # acceptance is a separate checked candidate contract.
                        theirs = reference.Coulomb(
                            q, xyz, probe, absolute, alpha, cutoff
                        )
                        cases.append(
                            {
                                "q": q,
                                "xyz": xyz,
                                "form": form,
                                "alpha": alpha,
                                "cutoff": cutoff,
                                "probe": probe,
                                "abs": absolute,
                            }
                        )
                        for point in queries:
                            for threshold, boundary in thresholds:
                                a = ours(*point, threshold)
                                b = theirs(*point, threshold)
                                delta = abs(a - b)
                                scale = delta / max(1.0, abs(b))
                                close = math.isclose(a, b, rel_tol=2e-12, abs_tol=2e-9)
                                if boundary:
                                    boundary_comparisons += 1
                                    if not close:
                                        near_site = False
                                        for row in xyz:
                                            dx, dy, dz = (
                                                v - w for v, w in zip(point, row)
                                            )
                                            d2 = dx * dx + dy * dy + dz * dz
                                            near_site |= abs(
                                                d2 - threshold
                                            ) <= 8 * math.ulp(d2)
                                        if not near_site:
                                            raise AssertionError(
                                                (size, point, threshold, a, b)
                                            )
                                        boundary_differences += 1
                                        if len(boundary_examples) < 3:
                                            boundary_examples.append(
                                                {
                                                    "sites": size,
                                                    "point": point,
                                                    "threshold": threshold,
                                                    "candidate": a,
                                                    "reference": b,
                                                }
                                            )
                                else:
                                    if not close:
                                        raise AssertionError(
                                            (size, point, threshold, a, b)
                                        )
                                    max_absolute = max(max_absolute, delta)
                                    max_scaled = max(max_scaled, scale)
                                    comparisons += 1
    probe = reference.Coulomb([1.0], [[0.0, 0.0, 0.0]])
    # Diagnostics for deliberately refused profiles are not parity assertions.
    tiny_cutoff = reference.Coulomb([1.0], [[0.0, 0.0, 0.0]], cutoff=0.0)
    payload = json.dumps(cases, sort_keys=True, separators=(",", ":")).encode()
    return {
        "schema_version": 1,
        "scope": "explicit-charge vacuum scalar Coulomb; synthetic supported profiles",
        "seed": 9601,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "rdkit": rdkit.__version__,
        "constructor_cases": len(cases),
        "scalar_comparisons": comparisons,
        "boundary_diagnostics": {
            "comparisons": boundary_comparisons,
            "differences_above_tolerance": boundary_differences,
            "examples": boundary_examples,
            "interpretation": "inclusion at rounded squared-distance boundaries may differ across Python/C++ arithmetic; every difference was within 8 ULP of a site boundary",
        },
        "max_absolute_error_kj_mol": max_absolute,
        "max_error_scaled_by_max_1_reference": max_scaled,
        "tolerance": {"relative": 2e-12, "absolute_kj_mol": 2e-9},
        "case_parameters_sha256": hashlib.sha256(payload).hexdigest(),
        "threshold_probe": {
            "distance": 2.0,
            "at_4": probe(2.0, 0.0, 0.0, 4.0),
            "above_4": probe(2.0, 0.0, 0.0, math.nextafter(4.0, math.inf)),
        },
        "unsupported_reference_diagnostic": {
            "cutoff_0_at_origin": tiny_cutoff(0.0, 0.0, 0.0, 1.0),
            "candidate": "cutoff below 0.001 is refused, not reproduced",
        },
        "exclusions": [
            "molecule/conformer overload",
            "automatic charges",
            "dielectric fields",
            "grids",
            "performance claim",
        ],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = json.dumps(validate(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(result)
    print(result, end="")
