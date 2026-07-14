#!/usr/bin/env bash
# Minimal self-check for criterion_gate.sh's block-pairing math (geomean +
# negate), without needing a real Criterion binary. Full end-to-end
# validation (does this actually avoid pseudo-replication on real timing
# noise, does it actually detect a real regression) was done manually against
# real compiled benchmarks on this machine -- see docs/cip_accurate_rfc.md's
# sibling issue #70 thread for those results; this script only guards the
# arithmetic against a future accidental edit.
set -euo pipefail
cd "$(dirname "$0")/.."

# Stub out run_point_estimate with fixed values keyed by binary path, so
# cmd_run_blocks's block-assembly logic runs unmodified against known inputs.
# (Extract to a real temp file rather than `source <(...)` -- macOS's bash
# 3.2 doesn't reliably persist function defs sourced via process substitution.)
extracted_fn=$(mktemp)
sed -n '/^cmd_run_blocks/,/^}/p' scripts/criterion_gate.sh > "$extracted_fn"
source "$extracted_fn"
rm -f "$extracted_fn"
run_point_estimate() {
  case "$1" in
    fast) echo 100 ;;
    slow) echo 200 ;;
  esac
}

out=$(mktemp)
cmd_run_blocks fast slow bench_name 1 0 0 1 abba "$out"
record=$(cat "$out")
echo "abba record: $record"

baseline=$(echo "$record" | jq '.baseline')
candidate=$(echo "$record" | jq '.candidate')
# fast=100 both runs -> geomean 100 -> negated -100; slow=200 both runs -> geomean 200 -> negated -200
[ "$baseline" = "-100" ] || { echo "FAIL: expected baseline -100, got $baseline"; exit 1; }
[ "$candidate" = "-200" ] || { echo "FAIL: expected candidate -200, got $candidate"; exit 1; }

cmd_run_blocks fast slow bench_name 1 0 0 1 baab "$out"
record_baab=$(cat "$out")
echo "baab record: $record_baab"
# order must not change which side is labeled baseline/candidate
[ "$(echo "$record_baab" | jq '.baseline')" = "-100" ] || { echo "FAIL: baab changed baseline label"; exit 1; }
[ "$(echo "$record_baab" | jq '.candidate')" = "-200" ] || { echo "FAIL: baab changed candidate label"; exit 1; }

rm -f "$out"
echo "OK"
