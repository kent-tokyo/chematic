#!/usr/bin/env python3
"""Common steady-state native batches; retain actual counts and raw elapsed time."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys

WARMUP_NS = 500_000_000
TIMED_NS = 100_000_000


def normalize(row: dict, mode: str, minimum_repeats: int, records: int) -> dict:
    if type(minimum_repeats) is not int or minimum_repeats < 1 or type(records) is not int or records < 1:
        raise ValueError("positive integer workload required")
    def positive(key: str) -> int:
        value = row.get(key)
        if type(value) is not int or value <= 0:
            raise ValueError(f"invalid native {key}")
        return value
    actual = positive("repeats")
    elapsed = positive("elapsed_ns")
    if actual < minimum_repeats or row.get("mode") != mode:
        raise ValueError("native mode or repetition count differs")
    if type(row.get("records")) is not int or row["records"] != actual * records:
        raise ValueError("native record accounting differs")
    if row.get("minimum_warmup_ns") != WARMUP_NS or row.get("minimum_timed_ns") != TIMED_NS:
        raise ValueError("native time budget differs")
    if positive("warmup_repeats") < 1 or positive("warmup_elapsed_ns") < WARMUP_NS or elapsed < TIMED_NS:
        raise ValueError("native time budget was not met")
    if elapsed // actual <= 0:
        raise ValueError("normalized batch duration is not positive")
    normalized = dict(row)
    normalized.update(elapsed_ns_raw=elapsed, elapsed_ns=elapsed // actual,
                      elapsed_ns_per_batch=elapsed // actual, records_per_batch=records,
                      minimum_repeats=minimum_repeats,
                      elapsed_unit="nanoseconds per full batch; raw elapsed and actual counts retained")
    return normalized


def timing_row(binary: Path, corpus: Path, mode: str, repeats: int, records: int) -> dict:
    process = subprocess.run(["/usr/bin/time", "-l", str(binary), str(corpus), mode, str(repeats), str(records)],
                             check=True, capture_output=True, timeout=180)
    row = normalize(json.loads(process.stdout), mode, repeats, records)
    match = re.search(rb"(\d+)\s+maximum resident set size", process.stderr)
    if match is None:
        raise ValueError("Darwin peak RSS unavailable")
    row["peak_process_rss_bytes"] = int(match.group(1))
    return row


def main() -> None:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--kind", choices=["sdf", "v3000"], required=True)
    parser.add_argument("--repeats", type=int, default=3)
    opts, remainder = parser.parse_known_args()
    if opts.repeats < 1:
        raise ValueError("positive minimum repeats required")
    os.environ.update(CHEMATIC_BENCH_WARMUP_NS=str(WARMUP_NS), CHEMATIC_BENCH_TIMED_NS=str(TIMED_NS),
                      CHEMATIC_BENCH_WARM_DUMP="1", CHEMATIC_BENCH_MIN_REPEATS=str(opts.repeats))
    sys.argv = [sys.argv[0], *remainder, "--repeats", str(opts.repeats)]
    import bench_openbabel_sdf_throughput as sdf
    sdf.timing_row = timing_row
    if opts.kind == "sdf":
        runner = sdf
    else:
        import bench_openbabel_v3000_throughput as runner
        runner.timing_row = timing_row
    runner.main()
    position = sys.argv.index("--output")
    output = Path(sys.argv[position + 1])
    report = json.loads(output.read_text())
    assert report["complete"]
    report["schema"] = "openbabel-steady-throughput/v1"
    report["protocol"].update(minimum_repeats=report["protocol"].pop("repeats"),
        minimum_warmup_ns=WARMUP_NS, minimum_timed_ns=TIMED_NS,
        elapsed_unit="nanoseconds per full batch; raw native elapsed, actual repeats and records retained",
        warmup="at least 500ms in the measured mode; parse/roundtrip create fresh objects; prepared write warmed",
        rewrite_preflight="dump after the same write warmup/time budget, minimum repetitions from protocol")
    output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
