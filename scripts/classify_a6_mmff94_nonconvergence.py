#!/usr/bin/env python3
"""Pair pinned MMFF94 200/400-iteration rows and classify non-convergence."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

BASELINE_SHA = "c2525cf7a12f7fd77b1881ee868b49785b7b62b44e006eb708d2c630b6e0d850"
EXTENDED_SHA = "f5a12d84230ffe078dce14e9571d9d9f2b91b7db7bd1ada960212ceea0c67155"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_rows(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as stream:
        rows = [json.loads(line) for line in stream if line.strip()]
    identities = [(row["tier"], row["row_index"], row["name"], row["smiles"]) for row in rows]
    if len(rows) != 265 or len(set(identities)) != len(rows):
        raise ValueError(f"expected 265 unique rows in {path}")
    return rows


def force_band(force: float) -> str:
    if force < 0.01:
        return "lt_0_01"
    if force < 0.1:
        return "0_01_to_0_1"
    if force < 1:
        return "0_1_to_1"
    if force < 10:
        return "1_to_10"
    return "ge_10"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--extended", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--gate-v1031", action="store_true")
    args = parser.parse_args()

    hashes = {"baseline": sha256(args.baseline), "extended": sha256(args.extended)}
    if args.gate_v1031 and hashes != {"baseline": BASELINE_SHA, "extended": EXTENDED_SHA}:
        raise ValueError(f"raw 200/400-iteration files differ from pinned v1.0.31 inputs: {hashes}")
    baseline = read_rows(args.baseline)
    extended = read_rows(args.extended)
    counts: Counter[str] = Counter()
    force_bands: Counter[str] = Counter()
    nonconverged: list[dict] = []
    for before, after in zip(baseline, extended, strict=True):
        identity = ("tier", "row_index", "name", "smiles", "primary_category", "arm")
        if any(before[key] != after[key] for key in identity):
            raise ValueError(f"row identity or arm changed at {before['row_index']}")
        if before["status"] != "success" or after["status"] != "success":
            raise ValueError(f"non-success row at {before['row_index']}")
        old, new = before["force_field"], after["force_field"]
        if not old or not new:
            raise ValueError(f"missing force-field outcome at {before['row_index']}")
        old_ok, new_ok = bool(old["converged"]), bool(new["converged"])
        counts["converged_200"] += old_ok
        counts["converged_400"] += new_ok
        counts["gained_convergence"] += not old_ok and new_ok
        counts["lost_convergence"] += old_ok and not new_ok
        if new_ok:
            continue
        counts["not_converged_400"] += 1
        iterations = int(new["iterations"])
        if iterations == 400:
            stop_class = "iteration_cap"
        elif 0 <= iterations < 400:
            stop_class = "early_stop_unexplained"
        else:
            raise ValueError(f"unexpected iteration count at {before['row_index']}: {iterations}")
        counts[stop_class] += 1
        force = float(new["max_residual_force"])
        if not math.isfinite(force) or force < 0:
            raise ValueError(f"invalid residual force at {before['row_index']}")
        if stop_class == "iteration_cap":
            force_bands[force_band(force)] += 1
        else:
            counts["early_stop_same_iteration_as_200"] += iterations == int(old["iterations"])
        missing = new.get("missing_parameter_classes", [])
        counts["nonconverged_with_missing_parameters"] += bool(missing)
        energy_before = new["energy_before"]["total"]
        energy_after = new["energy_after"]["total"]
        if not (math.isfinite(energy_before) and math.isfinite(energy_after)):
            raise ValueError(f"nonfinite energy at {before['row_index']}")
        counts["nonconverged_with_energy_increase"] += energy_after > energy_before
        nonconverged.append({
            "row_index": before["row_index"], "name": before["name"],
            "category": before["primary_category"], "stop_class": stop_class,
            "iterations_200": old["iterations"], "iterations_400": iterations,
            "max_residual_force_400": force,
            "missing_parameter_kinds": sorted({item["kind"] for item in missing}),
        })
    report = {
        "schema": "a6-mmff94-v1031-nonconvergence-paired/v1",
        "scope": "published v1.0.31 wheel; same 265-row cohort; optimizer outcomes only",
        "raw_sha256": hashes,
        "counts": dict(sorted(counts.items())),
        "force_bands_at_400_cap": dict(sorted(force_bands.items())),
        "nonconverged_rows": nonconverged,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.gate_v1031 and (
        counts["converged_200"] != 100 or counts["converged_400"] != 164
        or counts["gained_convergence"] != 64 or counts["lost_convergence"] != 0
        or counts["iteration_cap"] != 97 or counts["early_stop_unexplained"] != 4
        or counts["early_stop_same_iteration_as_200"] != 4
        or counts["not_converged_400"] != 101
        or [r["row_index"] for r in nonconverged if r["stop_class"] == "early_stop_unexplained"]
        != [170, 191, 227, 234]
    ):
        raise ValueError("v1.0.31 200/400-iteration non-convergence gate failed; inspect report")
    print(json.dumps(report["counts"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
