#!/usr/bin/env bash
# Process-level Criterion A/B runner for the CI regression gate (issue #70).
#
# Replaces the old approach of pairing one Criterion process's ~100 internal
# batch samples as if they were 100 independent trials (pseudo-replication --
# they all share one process's runner-environment state, so a single
# environment difference gets amplified into an extreme-looking uniform win
# rate). Instead: one Criterion *process run* -> one point estimate (median
# per-iteration time). N independent process-pair "ABBA" blocks (order
# A,B,B,A; geometric mean of each side's 2 runs within a block cancels slow
# thermal/frequency drift) -> N independent {baseline, candidate} records,
# fed to `veridict compare --metric sign-test`.
#
# Usage:
#   criterion_gate.sh run-blocks <bin_a> <bin_b> <bench_name> <n_blocks> \
#       <warm_up_secs> <measurement_secs> <sample_size> <order> <out_jsonl>
#
# <bin_a>/<bin_b> are compiled Criterion bench executable paths (from
# `cargo bench --no-run --message-format=json`, not a glob over
# target/release/deps -- multiple hash-suffixed binaries for the same bench
# name can coexist in a restored cache). Pass the SAME path for both to get a
# same-binary null control (locally validated: a same-binary run through this
# exact pipeline comes out "inconclusive", not "fail" -- see the criterion-gate
# job's null-control step).
#
# <order> is `abba` (block = a,b,b,a) or `baab` (block = b,a,a,b) -- output
# `baseline`/`candidate` fields always mean bin_a/bin_b regardless of which
# physically ran first, only the *order* of execution within the block
# flips. Used for two-stage confirmation: Stage 2 re-runs a Stage-1 fail
# candidate with reversed order, so a real regression must reproduce under
# genuinely different execution order, not just the same measurement twice.
set -euo pipefail

run_point_estimate() {
  local bin="$1" bench_name="$2" warm_up="$3" measurement="$4" sample_size="$5"
  local sample_json="target/criterion/${bench_name}/new/sample.json"
  "$bin" --bench "$bench_name" \
    --warm-up-time "$warm_up" --measurement-time "$measurement" \
    --sample-size "$sample_size" --noplot >/dev/null 2>&1
  jq '
    [ .times, .iters ] as $ti
    | [range(0; ($ti[0] | length))] | map($ti[0][.] / $ti[1][.]) | sort
    | .[length / 2 | floor]
  ' "$sample_json"
}

cmd_run_blocks() {
  local bin_a="$1" bin_b="$2" bench_name="$3" n_blocks="$4"
  local warm_up="$5" measurement="$6" sample_size="$7" order="$8" out="$9"
  : > "$out"
  for i in $(seq 1 "$n_blocks"); do
    local a1 b1 b2 a2
    if [ "$order" = "baab" ]; then
      b1=$(run_point_estimate "$bin_b" "$bench_name" "$warm_up" "$measurement" "$sample_size")
      a1=$(run_point_estimate "$bin_a" "$bench_name" "$warm_up" "$measurement" "$sample_size")
      a2=$(run_point_estimate "$bin_a" "$bench_name" "$warm_up" "$measurement" "$sample_size")
      b2=$(run_point_estimate "$bin_b" "$bench_name" "$warm_up" "$measurement" "$sample_size")
    else
      a1=$(run_point_estimate "$bin_a" "$bench_name" "$warm_up" "$measurement" "$sample_size")
      b1=$(run_point_estimate "$bin_b" "$bench_name" "$warm_up" "$measurement" "$sample_size")
      b2=$(run_point_estimate "$bin_b" "$bench_name" "$warm_up" "$measurement" "$sample_size")
      a2=$(run_point_estimate "$bin_a" "$bench_name" "$warm_up" "$measurement" "$sample_size")
    fi
    # Negate: veridict's mean-diff/sign-test treat a larger candidate-baseline
    # as an improvement, but for latency lower is better.
    jq -c -n --argjson a1 "$a1" --argjson a2 "$a2" --argjson b1 "$b1" --argjson b2 "$b2" --arg id "$i" '
      {
        id: $id,
        baseline: (-(($a1 * $a2) | sqrt)),
        candidate: (-(($b1 * $b2) | sqrt))
      }
    ' >> "$out"
  done
}

cmd_bin_path() {
  local crate_dir="$1" crate="$2" benchfile="$3"
  (cd "$crate_dir" && cargo bench -p "$crate" --bench "$benchfile" --no-run --message-format=json 2>/dev/null) \
    | jq -r 'select(.reason == "compiler-artifact" and .executable != null) | .executable' \
    | tail -1
}

case "${1:-}" in
  run-blocks)
    shift
    cmd_run_blocks "$@"
    ;;
  bin-path)
    shift
    cmd_bin_path "$@"
    ;;
  *)
    echo "usage: $0 run-blocks <bin_a> <bin_b> <bench_name> <n_blocks> <warm_up_secs> <measurement_secs> <sample_size> <order:abba|baab> <out_jsonl>" >&2
    echo "       $0 bin-path <crate_dir> <crate> <benchfile>" >&2
    exit 64
    ;;
esac
