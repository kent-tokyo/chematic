#!/usr/bin/env python3
"""
Benchmark import startup time: chematic vs RDKit.

Measures clean-process import latency using subprocess isolation so that
Python's module cache does not skew results.

Usage:
    python scripts/bench_startup.py
    python scripts/bench_startup.py --rdkit    # also time RDKit
    python scripts/bench_startup.py --json     # machine-readable output
    python scripts/bench_startup.py --runs 10  # more samples (default: 5)
"""

import argparse
import json
import subprocess
import sys
import time


SNIPPETS = {
    "chematic_import":     "import chematic",
    "chematic_first_parse": "import chematic; chematic.from_smiles('c1ccccc1')",
    "rdkit_import":        "from rdkit import Chem",
    "rdkit_first_parse":   "from rdkit import Chem; Chem.MolFromSmiles('c1ccccc1')",
}


def time_snippet(snippet: str, runs: int) -> tuple[float | None, list[float]]:
    times = []
    for _ in range(runs):
        t0 = time.perf_counter()
        r = subprocess.run([sys.executable, "-c", snippet],
                           capture_output=True, timeout=60)
        elapsed = time.perf_counter() - t0
        if r.returncode == 0:
            times.append(elapsed * 1000)
    if not times:
        return None, times
    sorted_times = sorted(times)
    return sorted_times[len(sorted_times) // 2], times   # median, raw samples


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rdkit", action="store_true", help="Also benchmark RDKit")
    ap.add_argument("--json",  action="store_true", help="Output JSON")
    ap.add_argument("--runs",  type=int, default=5,  help="Subprocess runs per measurement")
    args = ap.parse_args()

    keys = ["chematic_import", "chematic_first_parse"]
    if args.rdkit:
        keys += ["rdkit_import", "rdkit_first_parse"]

    results: dict[str, object] = {}
    rows: list[tuple[str, str]] = []

    for key in keys:
        t, samples = time_snippet(SNIPPETS[key], args.runs)
        results[key + "_ms"] = round(t, 1) if t is not None else None
        results[key + "_ms_samples"] = [round(s, 1) for s in samples]
        label = SNIPPETS[key][:55].ljust(55)
        val   = f"{t:>6.0f} ms" if t is not None else "  not found"
        rows.append((label, val))

    if args.rdkit:
        ci = results.get("chematic_import_ms")
        ri = results.get("rdkit_import_ms")
        if ci and ri:
            ratio = ri / ci
            results["import_speedup_x"] = round(ratio, 1)
            rows.append(("", ""))
            rows.append(("chematic import speedup vs RDKit", f"{ratio:.1f}×"))

    if args.json:
        print(json.dumps(results, indent=2))
    else:
        print(f"Startup benchmark  (median of {args.runs} subprocess runs)\n")
        for label, val in rows:
            if label:
                print(f"  {label}  {val}")
            else:
                print()


if __name__ == "__main__":
    main()
