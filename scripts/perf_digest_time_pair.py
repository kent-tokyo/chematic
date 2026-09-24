#!/usr/bin/env python3
"""Paired Rust cold-timing of two `tools/perf_digest` builds.

Runs `chematic-perf-digest time CORPUS` for a base and a head binary, one
operation at a time and alternating base/head for every round, so both
builds see the same machine state. Each call reports the best of REPS cold
passes (fresh molecule clones, empty derived caches) in microseconds per
molecule; the record keeps every round and the median of the rounds.

Output correctness is checked separately (`scripts/perf_digest_diff.sh`);
this script only times.

usage: perf_digest_time_pair.py --base BIN --head BIN --corpus FILE
           --ops a,b,c --rounds N --reps R --base-revision SHA
           --head-revision SHA --output-json OUT
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import statistics
import subprocess
import sys
from datetime import datetime, timezone

LINE = re.compile(r"^(\S+)\s+([0-9.]+) µs/mol \(cold\)")


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def run(binary: str, corpus: str, op: str, reps: int) -> float:
    env = dict(os.environ, ONLY=op, REPS=str(reps))
    out = subprocess.run([binary, "time", corpus], env=env, check=True,
                         capture_output=True, text=True).stdout
    for line in out.splitlines():
        m = LINE.match(line)
        if m and m.group(1) == op:
            return float(m.group(2))
    raise RuntimeError(f"no timing line for {op}: {out!r}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", required=True)
    ap.add_argument("--head", required=True)
    ap.add_argument("--corpus", required=True)
    ap.add_argument("--ops", required=True)
    ap.add_argument("--rounds", type=int, default=3)
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--base-revision", required=True)
    ap.add_argument("--head-revision", required=True)
    ap.add_argument("--output-json", required=True)
    args = ap.parse_args()

    results = {}
    for op in args.ops.split(","):
        base, head = [], []
        for _ in range(args.rounds):
            base.append(run(args.base, args.corpus, op, args.reps))
            head.append(run(args.head, args.corpus, op, args.reps))
        b, h = statistics.median(base), statistics.median(head)
        results[op] = {"base_us": base, "head_us": head, "base_median_us": b,
                       "head_median_us": h, "speedup": b / h if h else None}
        print(f"{op:22s} base {b:9.2f}  head {h:9.2f}  x{b / h:5.2f}", flush=True)

    record = {
        "schema": "chematic-perf-digest-time-pair/v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "base_revision": args.base_revision,
        "head_revision": args.head_revision,
        "corpus": {"path": os.path.relpath(args.corpus), "sha256": sha256(args.corpus)},
        "protocol": {"rounds": args.rounds, "reps_per_call": args.reps,
                     "per_call": "best of REPS cold passes (fresh clones), µs/molecule",
                     "aggregate": "median over rounds; base/head alternate every round"},
        "platform": platform.platform(),
        "cpu_count": os.cpu_count(),
        "ops": results,
    }
    with open(args.output_json, "w", encoding="utf-8") as fh:
        json.dump(record, fh, indent=1)
        fh.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
