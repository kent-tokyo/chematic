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

# --- Stage-1 routing screen (issue #70 follow-up) ---
# cmd_route_check replaced a sign test that could never fire at n=3 blocks
# (see criterion_gate.sh's comment above cmd_route_check for the math). These
# fixtures are the real block data from the +5%/+10% synthetic-regression CI
# runs and a synthetic clean no-op, so a future accidental edit that breaks
# routing shows up here without needing a real Criterion binary.
extracted_route=$(mktemp)
sed -n '/^cmd_route_check/,/^}/p' scripts/criterion_gate.sh > "$extracted_route"
source "$extracted_route"
rm -f "$extracted_route"

fixture_plus5=$(mktemp)
cat > "$fixture_plus5" << 'EOF'
{"id":"1","baseline":-11723.852948081518,"candidate":-12419.292912444453}
{"id":"2","baseline":-11710.817870733676,"candidate":-12516.890630251484}
{"id":"3","baseline":-11776.024760778471,"candidate":-13441.09866163192}
EOF
[ "$(cmd_route_check "$fixture_plus5" 1.04)" = "route" ] || { echo "FAIL: +5% fixture (median +6.9%) should route"; exit 1; }
rm -f "$fixture_plus5"

fixture_plus10=$(mktemp)
cat > "$fixture_plus10" << 'EOF'
{"id":"1","baseline":-11073.407252668874,"candidate":-12563.450800114455}
{"id":"2","baseline":-11155.27769228296,"candidate":-12628.101260372696}
{"id":"3","baseline":-11110.051949478835,"candidate":-12760.029773887674}
EOF
[ "$(cmd_route_check "$fixture_plus10" 1.04)" = "route" ] || { echo "FAIL: +10% fixture (median +13.5%) should route"; exit 1; }
rm -f "$fixture_plus10"

fixture_noop=$(mktemp)
cat > "$fixture_noop" << 'EOF'
{"id":"1","baseline":-11700.0,"candidate":-11705.0}
{"id":"2","baseline":-11710.0,"candidate":-11690.0}
{"id":"3","baseline":-11705.0,"candidate":-11720.0}
EOF
[ "$(cmd_route_check "$fixture_noop" 1.04)" = "no-route" ] || { echo "FAIL: clean no-op (mixed direction, sub-percent) should not route"; exit 1; }
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
[ "$(cmd_route_check "$fixture_ecfp4_incident" 1.04)" = "no-route" ] || { echo "FAIL: regression -- unanimous-but-under-threshold data must not route (issue #70 ecfp4_10mol false positive)"; exit 1; }
rm -f "$fixture_ecfp4_incident"

# Note: this script only unit-tests cmd_route_check's pure logic. The
# workflow-level invariants (Stage 1 alone never sets any_fail; any_fail=1
# only on a confirmed Stage 2 fail with a clean null control; a contaminated
# null control suppresses any_fail even on a Stage 2 fail) live in
# bench-pr-gate.yml's shell, not in this script, and are verified by code
# inspection plus real CI run data (issue #70) rather than by this harness --
# same scope limitation the original block-pairing tests above already had.

echo "OK"
