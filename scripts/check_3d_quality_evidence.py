#!/usr/bin/env python3
"""Validate the bounded local P5 quality-evidence bundle.

This gate checks the already-recorded machine-readable evidence only. It does
not promote local results to force-field parity or independent-oracle
coverage; those boundaries remain explicit in the source records.
"""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS = ROOT / "benchmarks"


def load(name: str) -> dict:
    return json.loads((BENCHMARKS / name).read_text(encoding="utf-8"))


def main() -> int:
    failures: list[str] = []
    classes = load("2026-09-09-3d-class-failure-rates-v1.0.10.json")
    energy = load("2026-09-09-3d-energy-sanity-v1.0.10.json")
    diversity = load("2026-09-09-ensemble-diversity-v1.0.10.json")
    rmsd = load("2026-09-09-symmetric-rmsd-oracle-v1.0.10.json")
    torsion = load("2026-09-09-symmetric-torsion-distance-v1.0.10.json")

    totals = classes["totals"]
    if totals != {"checked": 58, "new_ok": 58, "new_failures": 0, "legacy_dg_ok": 50, "legacy_dg_failures": 8}:
        failures.append(f"unexpected class totals: {totals}")
    for arm, result in energy["arms"].items():
        if result != {"attempted": 63, "successful": 61, "finite": 61, "non_increasing": 61}:
            failures.append(f"unexpected energy arm {arm}: {result}")
    repro = diversity["reproducibility"]
    if repro != {"corpus_molecules_reproduced": 58, "corpus_molecules_checked": 58, "coordinate_mismatches": 0, "different_seed_non_aliasing": 58, "different_seed_checked": 58}:
        failures.append(f"unexpected diversity reproduction: {repro}")
    if rmsd["summary"] != {"total": 6, "within_tolerance": 5, "known_gap": 1, "unexplained_mismatches": 0}:
        failures.append(f"unexpected RMSD summary: {rmsd['summary']}")
    if torsion["checks"]["changed_torsion_corpus"] != {"fixtures": 10, "measured": 10, "tests": 1, "passed": 1}:
        failures.append(f"unexpected torsion corpus: {torsion['checks']['changed_torsion_corpus']}")

    if failures:
        print("3D quality evidence check failed:")
        print("\n".join(failures))
        return 1
    print("3D quality evidence OK: bounded diversity, class failure, energy, RMSD, and torsion records match the pinned results")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
