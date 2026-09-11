#!/usr/bin/env bash
# Run the same checks as CI locally. Usage: bash scripts/check.sh
set -e
# The workspace version is the default for current v1.0.12 evidence. A small
# number of independent oracle/binding reports remain intentionally pinned to
# v1.0.10; those checks opt into the historical version explicitly below.
HISTORICAL_VERSION=1.0.10
echo "=== fmt ===" && cargo fmt --all -- --check
echo "=== unsafe surface ===" && python3 scripts/check_unsafe_surface.py
echo "=== shared cross-binding manifest ===" && python3 scripts/check_cross_binding_manifest.py
echo "=== static binding surface ===" && python3 scripts/check_binding_surface.py >/dev/null
echo "=== WASM artifact boundary ===" && SCHEMATIC_BENCHMARK_VERSION=1.0.12 python3 scripts/check_wasm_artifact_boundary.py
echo "=== benchmark record index ===" && python3 scripts/check_benchmark_index.py
echo "=== cross-engine streaming matrix ===" && python3 scripts/validate_streaming_cross_engine_matrix.py
echo "=== same-process contract bundle ===" && python3 scripts/check_same_process_contracts.py
echo "=== 3D quality evidence bundle ===" && SCHEMATIC_BENCHMARK_VERSION="$HISTORICAL_VERSION" python3 scripts/check_3d_quality_evidence.py
echo "=== MMFF94 RDKit availability oracle ===" && SCHEMATIC_BENCHMARK_VERSION="$HISTORICAL_VERSION" python3 scripts/check_mmff94_rdkit_availability_oracle.py
echo "=== reaction SMARTS contract evidence ===" && SCHEMATIC_BENCHMARK_VERSION="$HISTORICAL_VERSION" python3 scripts/check_reaction_smarts_contract_evidence.py
echo "=== reaction application evidence ===" && SCHEMATIC_BENCHMARK_VERSION="$HISTORICAL_VERSION" python3 scripts/check_reaction_application_evidence.py
echo "=== Python binding contract evidence ===" && SCHEMATIC_BENCHMARK_VERSION="$HISTORICAL_VERSION" python3 scripts/check_python_binding_contract_evidence.py
echo "=== streaming safety manifest ===" && python3 scripts/check_streaming_format_limits.py --validate-only
echo "=== streaming safety evidence ===" && python3 scripts/check_streaming_safety_evidence.py
echo "=== streaming parser-entry evidence ===" && python3 scripts/check_streaming_parser_entry_evidence.py
echo "=== triclinic neighbor evidence ===" && SCHEMATIC_BENCHMARK_VERSION="$HISTORICAL_VERSION" python3 scripts/check_triclinic_neighbor_evidence.py
echo "=== MMFF94 issue #337 determinism evidence ===" && SCHEMATIC_BENCHMARK_VERSION="$HISTORICAL_VERSION" python3 scripts/check_mmff94_issue337_determinism_evidence.py
echo "=== roadmap disposition ===" && SCHEMATIC_BENCHMARK_VERSION="$HISTORICAL_VERSION" python3 scripts/check_roadmap_disposition.py
echo "=== workflow action pins ===" && python3 scripts/check_workflow_pins.py
echo "=== clippy ===" && cargo clippy --workspace --all-targets -- -D warnings
echo "=== test ===" && cargo test --workspace --lib --quiet
echo "=== test (integration) ===" && cargo test --workspace --tests --quiet
# chematic-py's pytest suite only runs here if the venv already has an
# editable `chematic` build installed (`.venv/bin/maturin develop --release
# -m crates/chematic-py/Cargo.toml`) AND pytest (pyproject.toml declares no
# test deps, so it's not guaranteed present) -- the maturin build is slow
# (full release compile), so it isn't triggered automatically by this fast
# local pre-commit gate. CI runs it unconditionally in its own job
# (.github/workflows/ci.yml).
if .venv/bin/python3 -c "import chematic, pytest" &>/dev/null; then
    echo "=== test (chematic-py pytest) ===" && .venv/bin/python3 -m pytest crates/chematic-py/tests/ -q
else
    echo "=== test (chematic-py pytest) === (skipped: .venv/bin/python3 lacks chematic and/or pytest -- run .venv/bin/pip install pytest && .venv/bin/maturin develop --release -m crates/chematic-py/Cargo.toml first to include it locally)"
fi
if command -v cargo-deny &>/dev/null || cargo deny --version &>/dev/null 2>&1; then
    echo "=== deny ==="
    deny_output=""
    if deny_output=$(cargo deny --all-features check 2>&1); then
        printf '%s\n' "$deny_output"
    elif grep -q "failed to acquire advisory database lock" <<<"$deny_output"; then
        printf '%s\n' "$deny_output"
        echo "=== deny (writable advisory-db fallback) ==="
        bash scripts/check_cargo_deny_local.sh
    else
        printf '%s\n' "$deny_output" >&2
        exit 1
    fi
else
    echo "=== deny === (skipped: cargo-deny not installed)"
fi
echo "=== publish graph ===" && python3 scripts/check_publish_graph.py
echo "=== comparison corpus manifest ===" && python3 validation/cosmolkit_comparison/validate_manifest.py
echo "=== version ==="
VER=$(grep '^version = ' Cargo.toml | head -1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+')
fail=0
grep -q "### v$VER release boundary" README.md || { echo "MISMATCH: README.md (expect v$VER release boundary)"; fail=1; }
grep -q "### v$VER の対応範囲" README_ja.md || { echo "MISMATCH: README_ja.md (expect v$VER の対応範囲)"; fail=1; }
grep -q "version: $VER"  CITATION.cff || { echo "MISMATCH: CITATION.cff (expect version: $VER)"; fail=1; }
grep -q "v$VER | Yes"    SECURITY.md  || { echo "MISMATCH: SECURITY.md (expect v$VER | Yes)"; fail=1; }
grep -q "\"version\": \"$VER\"" demo/pkg/package.json || { echo "MISMATCH: demo/pkg/package.json (expect version $VER)"; fail=1; }
for f in crates/*/Cargo.toml; do
    if grep -q 'path = "\.\./chematic-' "$f" && ! grep -q "version = \"$VER\"" "$f"; then
        echo "MISMATCH: $f (path-dependency versions not bumped to $VER)"; fail=1
    fi
done
[ $fail -eq 0 ] && echo "Version consistent: $VER" || { echo "Run: python scripts/bump_version.py"; exit 1; }
# Soft staleness check (warning only, doesn't fail the build): a version bump can't
# auto-write a new "Recent Development" prose entry, so this can't be a hard MISMATCH --
# but silent drift here is exactly how the section went 10 versions stale unnoticed.
TOP_DEV_VER=$(grep -oE '\*\*v[0-9]+\.[0-9]+\.[0-9]+\*\*' README.md | head -1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' || true)
if [ -n "$TOP_DEV_VER" ] && [ "$TOP_DEV_VER" != "$VER" ]; then
    echo "WARNING: README.md 'Recent Development' section's newest entry is v$TOP_DEV_VER, workspace is v$VER -- consider adding an entry."
fi
echo "All checks passed."
