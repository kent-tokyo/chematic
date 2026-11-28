#!/usr/bin/env python3
"""Measure independent generation stages on previously exposed development inputs."""
import argparse
import hashlib
import json
import platform
import statistics
import time
from pathlib import Path

import chematic
from chematic import chematic as extension


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("repeats must be positive")
    rows = json.loads(args.input.read_text(encoding="utf-8"))["results"]
    results = []
    for row in rows:
        molecule = chematic.from_smiles(row["smiles"]).add_hydrogens()
        timings = {name: [] for name in ("embed", "checked_0", "checked_200")}
        counts = {}
        result_hashes = {name: [] for name in timings}
        for rep in range(args.repeats):
            operations = [
                ("embed", lambda: molecule.etkdg_embed_corrected(random_seed=42)),
                ("checked_0", lambda: molecule.etkdg_embed_multiple_confs_mixed_checked_diagnostics(1, random_seed=42, max_mmff_iters=0)),
                ("checked_200", lambda: molecule.etkdg_embed_multiple_confs_mixed_checked_diagnostics(1, random_seed=42, max_mmff_iters=200)),
            ]
            operations = operations[rep % 3:] + operations[:rep % 3]
            for name, operation in operations:
                started = time.perf_counter()
                try:
                    value = operation()
                    counts[name] = len(value if name == "embed" else value[0])
                    result_hashes[name].append(hashlib.sha256(
                        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest())
                except Exception as error:
                    counts[name] = type(error).__name__ + ": " + str(error)
                    result_hashes[name].append("error: " + counts[name])
                timings[name].append(time.perf_counter() - started)
        results.append({"index": row["index"], "name": row["name"],
                        "seconds": timings, "return_count_or_error": counts,
                        "result_sha256": result_hashes})
        print(f"{len(results)}/{len(rows)}", flush=True)
    summary = {name: statistics.median(statistics.median(r["seconds"][name]) for r in results)
               for name in ("embed", "checked_0", "checked_200")}
    output = {"input_sha256": digest(args.input), "extension_sha256": digest(Path(extension.__file__)),
              "python": platform.python_version(), "platform": platform.platform(),
              "repeats": args.repeats, "summary_median_seconds": summary, "results": results}
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
