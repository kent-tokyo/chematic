#!/usr/bin/env python3
"""Seeded explicit-charge dielectric comparison against RDKit 2026.09.1."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import random

import numpy as np
import rdkit
from rdkit.Chem import rdMIF as reference
from chematic.rdkit_compat import rdMIF as candidate


def validate() -> dict:
    if rdkit.__version__ != "2026.09.1":
        raise RuntimeError("requires RDKit 2026.09.1")
    rng = random.Random(960102)
    cases = []
    queries_digest = hashlib.sha256()
    comparisons = 0
    max_absolute = 0.0
    max_scaled = 0.0
    boundary_count = 0
    boundary_differences = 0
    examples = []
    for size in (0, 1, 7, 8, 9, 10, 11, 12, 13, 16, 64):
        charges = [rng.uniform(-2, 2) if i % 3 else 0.0 for i in range(size)]
        extent = 0.3 if 7 <= size <= 13 else 6.0
        positions = [
            [rng.uniform(-extent, extent) for _ in range(3)] for _ in range(size)
        ]
        points = [[rng.uniform(-1, 1) for _ in range(3)] for _ in range(8)]
        points += positions[:2]
        if size:
            dx, dy, dz = (a - b for a, b in zip(points[0], positions[0]))
            edge = dx * dx + dy * dy + dz * dz
            edges = [math.nextafter(edge, 0), edge, math.nextafter(edge, math.inf)]
        else:
            edges = []
        for reverse in (False, True):
            q = list(reversed(charges)) if reverse else charges
            xyz = list(reversed(positions)) if reverse else positions
            q_store = np.zeros(2 * size)
            q_store[::2] = q
            xyz_store = np.zeros((size, 6))
            xyz_store[:, ::2] = np.asarray(xyz).reshape(size, 3)
            for form, cq, cp in [
                ("list", q, xyz),
                ("strided_numpy", q_store[::2], xyz_store[:, ::2]),
            ]:
                for epsilon, xi in (
                    (80.0, 4.0),
                    (4.0, 80.0),
                    (1.0, 1.0),
                    (1e6, 1.0),
                    (1.0, 1e6),
                    (20.0, 20.0),
                ):
                    for alpha, cutoff in (
                        (0.0, 0.001),
                        (0.0, 1.0),
                        (1e-6, 1.0),
                        (0.3, 1.0),
                    ):
                        for probe, absolute in (
                            (1.0, False),
                            (-0.75, False),
                            (2.0, True),
                        ):
                            a = candidate.CoulombDielectric(
                                cq, cp, probe, absolute, alpha, cutoff, epsilon, xi
                            )
                            b = reference.CoulombDielectric(
                                q, xyz, probe, absolute, alpha, cutoff, epsilon, xi
                            )
                            cases.append(
                                {
                                    "q": q,
                                    "xyz": xyz,
                                    "form": form,
                                    "epsilon": epsilon,
                                    "xi": xi,
                                    "alpha": alpha,
                                    "cutoff": cutoff,
                                    "probe": probe,
                                    "abs": absolute,
                                }
                            )
                            for point in points:
                                for threshold in (0.0, 1e-6, 0.25, 4.0, 1000.0):
                                    actual = a(*point, threshold)
                                    expected = b(*point, threshold)
                                    delta = abs(actual - expected)
                                    if not math.isclose(
                                        actual, expected, rel_tol=2e-12, abs_tol=2e-9
                                    ):
                                        raise AssertionError(
                                            (
                                                size,
                                                point,
                                                threshold,
                                                cases[-1],
                                                actual,
                                                expected,
                                            )
                                        )
                                    comparisons += 1
                                    max_absolute = max(max_absolute, delta)
                                    max_scaled = max(
                                        max_scaled, delta / max(1.0, abs(expected))
                                    )
                                for threshold in edges:
                                    actual = a(*point, threshold)
                                    expected = b(*point, threshold)
                                    boundary_count += 1
                                    if not math.isclose(
                                        actual, expected, rel_tol=2e-12, abs_tol=2e-9
                                    ):
                                        near = False
                                        for row in xyz:
                                            dx, dy, dz = (
                                                v - w for v, w in zip(point, row)
                                            )
                                            d2 = dx * dx + dy * dy + dz * dz
                                            near |= abs(d2 - threshold) <= 8 * math.ulp(
                                                d2
                                            )
                                        if not near:
                                            raise AssertionError(
                                                (
                                                    "difference away from boundary",
                                                    size,
                                                    point,
                                                    threshold,
                                                    actual,
                                                    expected,
                                                )
                                            )
                                        boundary_differences += 1
                                        if len(examples) < 3:
                                            examples.append(
                                                {
                                                    "sites": size,
                                                    "point": point,
                                                    "threshold": threshold,
                                                    "candidate": actual,
                                                    "reference": expected,
                                                }
                                            )
                            queries_digest.update(
                                json.dumps(
                                    [points, edges], separators=(",", ":")
                                ).encode()
                            )
    # Strict 4 A neighbor boundary in both site and probe solvation counts.
    neighbor_checks = 0
    for n in range(7, 15):
        for spacing in (math.nextafter(4.0, 0), 4.0, math.nextafter(4.0, math.inf)):
            q = [1.0] + [0.0] * (n - 1)
            xyz = [[0.0, 0.0, 0.0]] + [[spacing, 0.0, 0.0]] * (n - 1)
            a = candidate.CoulombDielectric(q, xyz)
            b = reference.CoulombDielectric(q, xyz)
            for point in ((0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (4.0, 0.0, 0.0)):
                actual, expected = a(*point, 1.0), b(*point, 1.0)
                if not math.isclose(actual, expected, rel_tol=2e-12, abs_tol=2e-9):
                    raise AssertionError(
                        ("4 A boundary", n, spacing, point, actual, expected)
                    )
                neighbor_checks += 1
    return {
        "schema_version": 1,
        "seed": 960102,
        "scope": "copied explicit-charge local-solvation Coulomb model; synthetic scalar evaluations",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "rdkit": rdkit.__version__,
        "constructors": len(cases),
        "ordinary_comparisons": comparisons,
        "neighbor_boundary_checks": neighbor_checks,
        "max_absolute_error_kj_mol": max_absolute,
        "max_error_scaled_by_max_1_reference": max_scaled,
        "tolerance": {"relative": 2e-12, "absolute_kj_mol": 2e-9},
        "case_parameters_sha256": hashlib.sha256(
            json.dumps(cases, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest(),
        "query_stream_sha256": queries_digest.hexdigest(),
        "threshold_boundary_diagnostics": {
            "comparisons": boundary_count,
            "differences_above_tolerance": boundary_differences,
            "examples": examples,
            "interpretation": "diagnostics are not passing parity cases; differences must be within 8 ULP of an energy threshold boundary",
        },
        "exclusions": [
            "molecule overload",
            "charge generation",
            "arbitrary dielectric solver",
            "original threshold-boundary bit parity",
            "native build",
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
