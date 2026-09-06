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
#
#   criterion_gate.sh route-check <blocks_jsonl> <expected_count> <threshold>
#   criterion_gate.sh ratio-summary <blocks_jsonl> <expected_count>
#   criterion_gate.sh check-threshold <threshold>
#
# Stage-1 routing screen -- see cmd_route_check below for why Stage 1 can't
# use veridict's sign-test verdict directly. ratio-summary/check-threshold
# are the shared validation+computation helpers both route-check and the
# Stage 2 magnitude gate (in the workflow) call -- see cmd_ratio_summary.
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
    local started_at finished_at loadavg cpu_model steal_time
    started_at=$(date -u '+%Y-%m-%dT%H:%M:%SZ')
    if [ -r /proc/loadavg ]; then
      loadavg=$(awk '{print $1}' /proc/loadavg)
    elif command -v sysctl >/dev/null 2>&1; then
      loadavg=$(sysctl -n vm.loadavg 2>/dev/null | awk '{gsub(/[{}]/, ""); print $1}')
    else
      loadavg="unavailable"
    fi
    if [ -r /proc/cpuinfo ]; then
      cpu_model=$(awk -F: '/model name|Hardware|chip type/ {gsub(/^ +/, "", $2); print $2; exit}' /proc/cpuinfo)
    elif command -v sysctl >/dev/null 2>&1; then
      cpu_model=$(sysctl -n machdep.cpu.brand_string 2>/dev/null || true)
    fi
    cpu_model=${cpu_model:-unavailable}
    if [ -r /proc/stat ]; then
      steal_time=$(awk '/^cpu / {print $9; exit}' /proc/stat)
    else
      steal_time="unavailable"
    fi
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
    finished_at=$(date -u '+%Y-%m-%dT%H:%M:%SZ')
    jq -c -n --argjson a1 "$a1" --argjson a2 "$a2" --argjson b1 "$b1" --argjson b2 "$b2" \
      --arg id "$i" --arg order "$order" --arg started_at "$started_at" \
      --arg finished_at "$finished_at" --arg loadavg "$loadavg" \
      --arg cpu_model "$cpu_model" --arg steal_time "$steal_time" '
      {
        id: $id,
        baseline: (-(($a1 * $a2) | sqrt)),
        candidate: (-(($b1 * $b2) | sqrt)),
        execution_order: $order,
        started_at: $started_at,
        finished_at: $finished_at,
        environment: {
          loadavg: $loadavg,
          cpu_model: $cpu_model,
          proc_stat_steal_ticks: $steal_time
        }
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

# Shared helper (issue #70 follow-up): strict validation + median-ratio
# computation for a run-blocks jsonl file, used by both Stage 1's
# route-check and Stage 2's magnitude gate (in the workflow). A malformed or
# short/long block file must never silently degrade to a default decision
# ("no-route", "pass") -- that would hide a real infrastructure problem
# (a crashed run-blocks invocation, a truncated artifact) behind a verdict
# that looks like ordinary negative data. Every failure mode below exits
# non-zero with a message on stderr instead.
cmd_ratio_summary() {
  local jsonl="$1" expected_count="$2"
  if [ ! -s "$jsonl" ]; then
    echo "ratio-summary: $jsonl is missing or empty" >&2
    return 1
  fi
  jq -rs --argjson expected "$expected_count" '
    if length != $expected then
      error("expected exactly \($expected) records, got \(length)")
    else . end
    | (map(.id) | unique | length) as $uniq_ids
    | if $uniq_ids != length then
        error("duplicate id field -- ids must be unique")
      else . end
    | map(
        if (.baseline | type) != "number" or (.candidate | type) != "number" then
          error("baseline/candidate must be numbers")
        elif (.baseline | isinfinite or isnan) or (.candidate | isinfinite or isnan) then
          error("baseline/candidate must be finite")
        elif .baseline == 0 then
          error("baseline must not be zero")
        else . end
      )
    | map((.candidate | fabs) / (.baseline | fabs)) as $ratios
    | if ($ratios | map(isinfinite or isnan or . <= 0) | any) then
        error("computed ratio must be finite and positive")
      else . end
    | ($ratios | sort) as $sorted
    | ($sorted | length) as $n
    | if ($n % 2) == 1 then
        $sorted[($n - 1) / 2 | floor]
      else
        (($sorted[($n / 2 | floor) - 1] + $sorted[$n / 2 | floor]) / 2)
      end
  ' "$jsonl"
}

# Threshold values below 1 would mean "route/fail even when the candidate is
# faster or equal," which is never a sensible regression-detection bound --
# reject those (and non-numeric/non-finite input) before they can silently
# produce a nonsense comparison result downstream.
cmd_check_threshold() {
  local threshold="$1"
  case "$threshold" in
    ''|*[!0-9.]*)
      echo "threshold '$threshold' is not a valid number" >&2
      return 1
      ;;
  esac
  jq -n --argjson t "$threshold" '
    if ($t | isinfinite or isnan) then error("threshold must be finite, got \($t)")
    elif $t < 1 then error("threshold must be >= 1, got \($t)")
    else empty end
  '
}

# Stage-1 routing screen (issue #70 follow-up). Stage 1 only ever ran 3
# blocks through `veridict compare --metric sign-test`, but a sign test's
# strongest possible signal at n=3 (a unanimous 3-0 split) has a two-sided
# exact p-value of 2*(1/2)^3=0.25 -- it can NEVER cross a 95%-confidence fail
# bar, so Stage 1 could never produce a fail verdict for a regression of any
# size, making Stage 2 (the only place that sets any_fail) unreachable.
# Verified empirically: synthetic +5%/+10% regressions both came back
# "inconclusive" with byte-identical confidence bounds (issue #70 comment).
#
# Fix: Stage 1 stops asserting statistical significance and becomes a cheap
# routing screen instead -- does this benchmark look suspicious enough to
# spend Stage 2's 10-block budget on? Routes purely on the block-level
# latency ratio's median crossing <threshold>.
#
# A "route if all 3 blocks agree on direction, even below <threshold>" leg
# was tried first, to catch a small-but-consistent regression a magnitude
# rule alone would miss -- and was the actual shipped behavior in an earlier
# revision of this function. It's deliberately NOT here: an offline
# evaluation against 28 historical no-op runs (issue #70) showed pure
# direction-agreement is dominated by sampling noise at n=3 (21-22%
# false-routing, vs 4% for the magnitude threshold alone) -- with 16
# benchmarks routed independently per run, unanimity fires on noise often
# enough to be worse than not screening at all. Magnitude-only was the
# selected candidate; keep it that way unless a re-evaluation says otherwise.
# Stage 1 alone never sets any_fail; only Stage 2's real sign-test AND
# magnitude gate (with a null-control-clean check) does that -- see the
# workflow's Stage 2 loop, not this function.
cmd_route_check() {
  local jsonl="$1" expected_count="$2" threshold="$3"
  cmd_check_threshold "$threshold" || return 1
  local median
  median=$(cmd_ratio_summary "$jsonl" "$expected_count") || return 1
  jq -rn --argjson m "$median" --argjson t "$threshold" \
    'if $m >= $t then "route" else "no-route" end'
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
  route-check)
    shift
    cmd_route_check "$@"
    ;;
  ratio-summary)
    shift
    cmd_ratio_summary "$@"
    ;;
  check-threshold)
    shift
    cmd_check_threshold "$@"
    ;;
  *)
    echo "usage: $0 run-blocks <bin_a> <bin_b> <bench_name> <n_blocks> <warm_up_secs> <measurement_secs> <sample_size> <order:abba|baab> <out_jsonl>" >&2
    echo "       $0 bin-path <crate_dir> <crate> <benchfile>" >&2
    echo "       $0 route-check <blocks_jsonl> <expected_count> <threshold>" >&2
    echo "       $0 ratio-summary <blocks_jsonl> <expected_count>" >&2
    echo "       $0 check-threshold <threshold>" >&2
    exit 64
    ;;
esac
