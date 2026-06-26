#!/usr/bin/env bash
# Run the same checks as CI locally. Usage: bash scripts/check.sh
set -e
echo "=== fmt ===" && cargo fmt --all -- --check
echo "=== clippy ===" && cargo clippy --workspace --all-targets -- -D warnings
echo "=== test ===" && cargo test --workspace --lib --quiet
if command -v cargo-deny &>/dev/null || cargo deny --version &>/dev/null 2>&1; then
    echo "=== deny ===" && cargo deny check --all-features
else
    echo "=== deny === (skipped: cargo-deny not installed)"
fi
echo "All checks passed."
