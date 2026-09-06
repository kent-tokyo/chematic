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

# Keep the checked-in calibration contract parseable. The actual hosted
# experiments remain a separate evidence gate; this only prevents fixture
# drift from making the local arithmetic tests document-only.
calibration_manifest="validation/criterion-gate-calibration.json"
jq -e '.schema_version == 1 and (.cases | length) == 7' "$calibration_manifest" >/dev/null

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
[ "$(echo "$record" | jq -r '.execution_order')" = "abba" ] \
  || { echo "FAIL: block execution order metadata was not persisted"; exit 1; }
jq -e '.started_at | strings and length > 0' <<<"$record" >/dev/null \
  || { echo "FAIL: block start timestamp metadata was not persisted"; exit 1; }
jq -e '.finished_at | strings and length > 0' <<<"$record" >/dev/null \
  || { echo "FAIL: block finish timestamp metadata was not persisted"; exit 1; }
jq -e '.environment | has("loadavg") and has("cpu_model") and has("proc_stat_steal_ticks")' \
  <<<"$record" >/dev/null \
  || { echo "FAIL: environment metadata was not persisted"; exit 1; }

cmd_run_blocks fast slow bench_name 1 0 0 1 baab "$out"
record_baab=$(cat "$out")
echo "baab record: $record_baab"
# order must not change which side is labeled baseline/candidate
[ "$(echo "$record_baab" | jq '.baseline')" = "-100" ] || { echo "FAIL: baab changed baseline label"; exit 1; }
[ "$(echo "$record_baab" | jq '.candidate')" = "-200" ] || { echo "FAIL: baab changed candidate label"; exit 1; }

rm -f "$out"

# --- Shared ratio-summary/check-threshold/route-check (issue #70 follow-up) ---
# route-check calls check-threshold and ratio-summary by name, so extract all
# three together.
extracted_shared=$(mktemp)
sed -n '/^cmd_ratio_summary/,/^}/p;/^cmd_check_threshold/,/^}/p;/^cmd_route_check/,/^}/p' \
  scripts/criterion_gate.sh > "$extracted_shared"
source "$extracted_shared"
rm -f "$extracted_shared"

fixture_plus5=$(mktemp)
cat > "$fixture_plus5" << 'EOF'
{"id":"1","baseline":-11723.852948081518,"candidate":-12419.292912444453}
{"id":"2","baseline":-11710.817870733676,"candidate":-12516.890630251484}
{"id":"3","baseline":-11776.024760778471,"candidate":-13441.09866163192}
EOF
[ "$(cmd_route_check "$fixture_plus5" 3 1.04)" = "route" ] || { echo "FAIL: +5% fixture (median +6.9%) should route"; exit 1; }
rm -f "$fixture_plus5"

fixture_plus10=$(mktemp)
cat > "$fixture_plus10" << 'EOF'
{"id":"1","baseline":-11073.407252668874,"candidate":-12563.450800114455}
{"id":"2","baseline":-11155.27769228296,"candidate":-12628.101260372696}
{"id":"3","baseline":-11110.051949478835,"candidate":-12760.029773887674}
EOF
[ "$(cmd_route_check "$fixture_plus10" 3 1.04)" = "route" ] || { echo "FAIL: +10% fixture (median +13.5%) should route"; exit 1; }
rm -f "$fixture_plus10"

# +6% fixture (this PR's own required minimum test): clearly above the 1.04
# routing threshold, must route.
fixture_plus6=$(mktemp)
cat > "$fixture_plus6" << 'EOF'
{"id":"1","baseline":-10000.0,"candidate":-10600.0}
{"id":"2","baseline":-10000.0,"candidate":-10620.0}
{"id":"3","baseline":-10000.0,"candidate":-10580.0}
EOF
[ "$(cmd_route_check "$fixture_plus6" 3 1.04)" = "route" ] || { echo "FAIL: +6% fixture should route"; exit 1; }
rm -f "$fixture_plus6"

fixture_noop=$(mktemp)
cat > "$fixture_noop" << 'EOF'
{"id":"1","baseline":-11700.0,"candidate":-11705.0}
{"id":"2","baseline":-11710.0,"candidate":-11690.0}
{"id":"3","baseline":-11705.0,"candidate":-11720.0}
EOF
[ "$(cmd_route_check "$fixture_noop" 3 1.04)" = "no-route" ] || { echo "FAIL: clean no-op (mixed direction, sub-percent) should not route"; exit 1; }
rm -f "$fixture_noop"

# Regression test for a real incident: PR #117's own first CI run shipped a
# route-check with an "OR all 3 blocks agree on direction" leg (the offline
# eval had already disqualified that leg -- 21-22% no-op false-routing, see
# the comment above cmd_route_check -- but the shipped code kept it anyway).
# It routed `ecfp4_10mol` to Stage 2 purely on unanimous direction (median
# 1.0149, under the 1.04 threshold) on a PR touching zero chematic-fp files,
# and Stage 2 then confirmed a "fail" from real but spurious ~1% build-to-
# build variance -- a genuine false positive on an unrelated benchmark. This
# is that exact Stage-1 block data (issue #70): must NOT route under a pure
# magnitude threshold, even though all 3 blocks agree on direction.
fixture_ecfp4_incident=$(mktemp)
cat > "$fixture_ecfp4_incident" << 'EOF'
{"id":"1","baseline":-165481.56427608195,"candidate":-167948.2795937714}
{"id":"2","baseline":-165755.31364862694,"candidate":-171745.03897197713}
{"id":"3","baseline":-166218.16842646126,"candidate":-167825.3596853625}
EOF
[ "$(cmd_route_check "$fixture_ecfp4_incident" 3 1.04)" = "no-route" ] || { echo "FAIL: regression -- unanimous-but-under-threshold data must not route (issue #70 ecfp4_10mol false positive)"; exit 1; }
rm -f "$fixture_ecfp4_incident"

# --- Strict validation (issue #70 follow-up): malformed/insufficient input
# must fail loudly, never silently degrade to a default "no-route" ---

empty_jsonl=$(mktemp)
: > "$empty_jsonl"
if cmd_route_check "$empty_jsonl" 3 1.04 >/dev/null 2>&1; then
  echo "FAIL: empty Stage-1 jsonl should error, not return a verdict"; exit 1
fi
rm -f "$empty_jsonl"

too_few=$(mktemp)
cat > "$too_few" << 'EOF'
{"id":"1","baseline":-100.0,"candidate":-100.5}
{"id":"2","baseline":-101.0,"candidate":-99.0}
EOF
if cmd_route_check "$too_few" 3 1.04 >/dev/null 2>&1; then
  echo "FAIL: 2-record file (expected 3) should error, not return a verdict"; exit 1
fi
rm -f "$too_few"

too_many=$(mktemp)
cat > "$too_many" << 'EOF'
{"id":"1","baseline":-100.0,"candidate":-100.5}
{"id":"2","baseline":-101.0,"candidate":-99.0}
{"id":"3","baseline":-99.5,"candidate":-100.2}
{"id":"4","baseline":-99.0,"candidate":-100.0}
EOF
if cmd_route_check "$too_many" 3 1.04 >/dev/null 2>&1; then
  echo "FAIL: 4-record file (expected 3) should error, not return a verdict"; exit 1
fi
rm -f "$too_many"

zero_baseline=$(mktemp)
cat > "$zero_baseline" << 'EOF'
{"id":"1","baseline":0.0,"candidate":-100.5}
{"id":"2","baseline":-101.0,"candidate":-99.0}
{"id":"3","baseline":-99.5,"candidate":-100.2}
EOF
if cmd_route_check "$zero_baseline" 3 1.04 >/dev/null 2>&1; then
  echo "FAIL: baseline=0 should error (division by zero), not return a verdict"; exit 1
fi
rm -f "$zero_baseline"

if cmd_check_threshold "0.5" >/dev/null 2>&1; then
  echo "FAIL: threshold below 1 should be rejected"; exit 1
fi

# --- Stage 2 magnitude gate (issue #70 follow-up): the real ecfp4_10mol fix
# stopped Stage 1 from routing on noise, but the very next CI run confirmed a
# DIFFERENT benchmark (parse_smiles_10mol) via a legitimate route + a
# sign-test fail driven by real-but-spurious ~2.2% build/codegen variance
# between separately-compiled binaries -- Stage 2's sign-test alone is
# magnitude-blind. These fixtures pin the practical-effect gate that
# combines the sign-test verdict with a minimum median-ratio requirement. ---

# Permanent regression fixture: the exact parse_smiles_10mol Stage 2 data
# from the real incident (issue #70). Median ~1.022, all 10 blocks candidate
# slower (a real sign-test "fail" from veridict) -- must NOT meet the 1.04
# practical-effect threshold, so the final gate must not fail on this data.
fixture_stage2_incident=$(mktemp)
cat > "$fixture_stage2_incident" << 'EOF'
{"id":"1","baseline":-9076.104372053887,"candidate":-9201.91136904256}
{"id":"2","baseline":-8932.942526291896,"candidate":-9219.596574885541}
{"id":"3","baseline":-8910.356356885584,"candidate":-9143.826173380521}
{"id":"4","baseline":-8944.482964007622,"candidate":-9235.42058864827}
{"id":"5","baseline":-8952.404397055447,"candidate":-9106.500793595658}
{"id":"6","baseline":-8971.499186556917,"candidate":-9324.067453203195}
{"id":"7","baseline":-8918.297767483793,"candidate":-9081.170864737667}
{"id":"8","baseline":-9064.68757885133,"candidate":-9099.614237070064}
{"id":"9","baseline":-9012.72747438855,"candidate":-9246.215896065525}
{"id":"10","baseline":-9013.023872561964,"candidate":-9108.926275796783}
EOF
incident_median=$(cmd_ratio_summary "$fixture_stage2_incident" 10)
# Loose float comparison via awk -- bash has no native float arithmetic.
awk -v m="$incident_median" 'BEGIN { if (m < 1.021 || m > 1.023) exit 1; exit 0 }' \
  || { echo "FAIL: expected median ~1.022 for the pinned incident fixture, got $incident_median"; exit 1; }
awk -v m="$incident_median" 'BEGIN { exit !(m >= 1.04) }' \
  && { echo "FAIL: regression -- pinned incident median ($incident_median) must stay below STAGE2_FAIL_THRESHOLD=1.04, or the final gate will start blocking spurious ~2.2% build noise again"; exit 1; }
rm -f "$fixture_stage2_incident"

# Sibling fixture: a real ~6% Stage 2 difference (well above 1.04) DOES meet
# the practical-effect threshold -- confirms the gate isn't just permanently
# disabled, only insensitive to sub-threshold noise.
fixture_stage2_real=$(mktemp)
cat > "$fixture_stage2_real" << 'EOF'
{"id":"1","baseline":-10000.0,"candidate":-10600.0}
{"id":"2","baseline":-10000.0,"candidate":-10610.0}
{"id":"3","baseline":-10000.0,"candidate":-10590.0}
{"id":"4","baseline":-10000.0,"candidate":-10620.0}
{"id":"5","baseline":-10000.0,"candidate":-10580.0}
{"id":"6","baseline":-10000.0,"candidate":-10605.0}
{"id":"7","baseline":-10000.0,"candidate":-10595.0}
{"id":"8","baseline":-10000.0,"candidate":-10615.0}
{"id":"9","baseline":-10000.0,"candidate":-10585.0}
{"id":"10","baseline":-10000.0,"candidate":-10600.0}
EOF
real_median=$(cmd_ratio_summary "$fixture_stage2_real" 10)
awk -v m="$real_median" 'BEGIN { exit !(m >= 1.04) }' \
  || { echo "FAIL: expected median ~1.06 (>= 1.04 threshold) for the real-effect fixture, got $real_median"; exit 1; }
rm -f "$fixture_stage2_real"

# Stage-2 strict validation, same shape as Stage 1's: exactly 10 records
# required, not 9 or 11.
stage2_too_few=$(mktemp)
for i in 1 2 3 4 5 6 7 8 9; do
  printf '{"id":"%d","baseline":-10000.0,"candidate":-10600.0}\n' "$i" >> "$stage2_too_few"
done
if cmd_ratio_summary "$stage2_too_few" 10 >/dev/null 2>&1; then
  echo "FAIL: 9-record Stage-2 file (expected 10) should error"; exit 1
fi
rm -f "$stage2_too_few"

# The workflow-level null control must compare independently compiled main
# checkouts, not pass the same executable on both sides. Keep this contract
# checked locally so a future workflow edit cannot silently restore the old
# same-binary control that missed build/codegen variance.
workflow=.github/workflows/bench-pr-gate.yml
grep -Fq 'path: main-null' "$workflow" \
  || { echo "FAIL: Criterion workflow has no independent main-null checkout"; exit 1; }
grep -Fq 'null_control_bin_b=$(bin_path main-null' "$workflow" \
  || { echo "FAIL: Criterion workflow does not resolve the independent null binary"; exit 1; }
grep -Fq '"$null_control_bin_a" "$null_control_bin_b" "parse_smiles_10mol"' "$workflow" \
  || { echo "FAIL: Criterion workflow null control is not a two-build comparison"; exit 1; }

# Note: this script unit-tests cmd_ratio_summary/cmd_check_threshold/
# cmd_route_check's pure logic, including the exact incident data. The
# workflow-level combination logic in bench-pr-gate.yml --
#   any_fail=1 requires ALL of: Stage 1 routed (median >= STAGE1_ROUTE_THRESHOLD)
#   AND Stage 2 sign-test verdict == fail AND Stage 2 median >= STAGE2_FAIL_THRESHOLD
#   AND environment_contaminated == 0
#   -- a sign-test fail below STAGE2_FAIL_THRESHOLD is small-effect-inconclusive, not blocking
#   -- a sign-test fail with environment_contaminated == 1 is environment-inconclusive, not blocking
#   -- Stage 1 alone never sets any_fail (no code path in the Stage 1 loop touches it)
# -- lives in bench-pr-gate.yml's shell, not in this script, and is verified
# by code inspection plus real CI run data (issue #70) rather than by this
# harness -- same scope limitation the original block-pairing tests above
# already had.

echo "OK"
