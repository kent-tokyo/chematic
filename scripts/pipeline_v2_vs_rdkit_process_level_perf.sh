#!/usr/bin/env bash
# Wave 1 follow-up: process-level (whole-corpus, separate-process) wall-clock
# performance, to sit alongside (not replace) the in-process per-molecule
# p50/p95/p99 already in the Wave 1 report. Each run is a fresh OS process
# reading the full committed corpus manifests; row-level output is discarded
# (already committed from the original Wave 1 run) -- only wall-clock time
# per whole-corpus run is recorded here.
#
# Run sequentially, never concurrently with another run of itself or with a
# build -- concurrent CPU contention would corrupt exactly the number this
# script exists to measure.
#
# Usage: bash scripts/pipeline_v2_vs_rdkit_process_level_perf.sh
# Output: validation/results/pipeline_v2_vs_rdkit_process_level_perf.json

set -euo pipefail
cd "$(dirname "$0")/.."

RUNS=5
OUT=validation/results/pipeline_v2_vs_rdkit_process_level_perf.json

echo "=== machine info ===" >&2
uname -a >&2
sysctl -n machdep.cpu.brand_string 2>/dev/null >&2 || true
sysctl -n hw.memsize 2>/dev/null | awk '{print $1/1024/1024/1024 " GB"}' >&2 || true
echo "load average before run: $(uptime)" >&2

CHEMATIC_BIN=./target/release/examples/pipeline_v2_vs_rdkit_dump
if [ ! -x "$CHEMATIC_BIN" ]; then
  echo "building chematic dump (release)..." >&2
  cargo build --release -p chematic-3d --example pipeline_v2_vs_rdkit_dump >&2
fi

chematic_times=()
for i in $(seq 1 "$RUNS"); do
  echo "chematic run $i/$RUNS..." >&2
  start=$(python3 -c "import time; print(time.monotonic())")
  "$CHEMATIC_BIN" > /dev/null 2>/dev/null
  end=$(python3 -c "import time; print(time.monotonic())")
  elapsed=$(python3 -c "print($end - $start)")
  echo "  ${elapsed}s" >&2
  chematic_times+=("$elapsed")
done

rdkit_times=()
for i in $(seq 1 "$RUNS"); do
  echo "rdkit run $i/$RUNS..." >&2
  start=$(python3 -c "import time; print(time.monotonic())")
  .venv/bin/python scripts/pipeline_v2_vs_rdkit_oracle.py > /dev/null 2>/dev/null
  end=$(python3 -c "import time; print(time.monotonic())")
  elapsed=$(python3 -c "print($end - $start)")
  echo "  ${elapsed}s" >&2
  rdkit_times+=("$elapsed")
done

echo "load average after run: $(uptime)" >&2

python3 - "$OUT" "${RUNS}" "${chematic_times[@]}" -- "${rdkit_times[@]}" << 'PYEOF'
import json
import statistics
import sys

out_path = sys.argv[1]
runs = int(sys.argv[2])
rest = sys.argv[3:]
sep = rest.index("--")
chematic_times = [float(x) for x in rest[:sep]]
rdkit_times = [float(x) for x in rest[sep + 1:]]

def summarize(times):
    return {
        "runs": len(times),
        "all_seconds": times,
        "median_seconds": statistics.median(times),
        "min_seconds": min(times),
        "max_seconds": max(times),
        "stdev_seconds": statistics.stdev(times) if len(times) > 1 else 0.0,
        "coefficient_of_variation": (statistics.stdev(times) / statistics.mean(times)) if len(times) > 1 else 0.0,
    }

result = {
    "methodology": "Whole-corpus (265 molecules), separate-process wall-clock "
    "per run, sequential (never concurrent with another run or a build). "
    "Includes process startup (Rust binary startup / Python+RDKit import) -- "
    "not steady-state-only.",
    "chematic": summarize(chematic_times),
    "rdkit": summarize(rdkit_times),
}

with open(out_path, "w") as f:
    json.dump(result, f, indent=2)
print(json.dumps(result, indent=2))
PYEOF
