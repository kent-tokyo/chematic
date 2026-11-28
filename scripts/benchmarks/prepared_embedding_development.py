#!/usr/bin/env python3
"""Compare repeated single-attempt setup with reusable embedding setup."""
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


def result_hash(value):
    return hashlib.sha256(json.dumps(value, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.input.read_text(encoding="utf-8"))
    results = []
    for index, record in enumerate(source["results"]):
        molecule = chematic.from_smiles(record["smiles"]).add_hydrogens()
        start = time.perf_counter()
        prepared = molecule._timed_prepare_corrected(max_iterations=0)
        prepare_seconds = time.perf_counter() - start
        attempts = []
        for attempt, seed in enumerate((42, 84, 126)):
            operations = [("individual", lambda: molecule._timed_embed_corrected(random_seed=seed)),
                          ("prepared", lambda: prepared.embed(seed))]
            if (index + attempt) % 2:
                operations.reverse()
            entry = {"seed": seed}
            for name, operation in operations:
                start = time.perf_counter()
                coords = operation()
                entry[name] = {"seconds": time.perf_counter() - start,
                               "sha256": result_hash(coords), "embedded": coords is not None}
            entry["matches"] = entry["individual"]["sha256"] == entry["prepared"]["sha256"]
            attempts.append(entry)
        results.append({"index": record["index"], "name": record["name"],
                        "prepare_seconds": prepare_seconds, "attempts": attempts,
                        "individual_total_seconds": sum(a["individual"]["seconds"] for a in attempts),
                        "prepared_total_seconds": prepare_seconds + sum(a["prepared"]["seconds"] for a in attempts)})
        print(f"{len(results)}/{len(source['results'])}", flush=True)
    output = {"input_sha256": digest(args.input), "extension_sha256": digest(Path(extension.__file__)),
              "script_sha256": digest(Path(__file__)), "python": platform.python_version(),
              "platform": platform.platform(), "seeds": [42, 84, 126], "results": results,
              "median_individual_total_seconds": statistics.median(r["individual_total_seconds"] for r in results),
              "median_prepared_total_seconds": statistics.median(r["prepared_total_seconds"] for r in results)}
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
