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
echo "=== version ==="
VER=$(grep '^version = ' Cargo.toml | head -1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+')
fail=0
grep -q "chematic v$VER" README.md    || { echo "MISMATCH: README.md (expect chematic v$VER)"; fail=1; }
grep -q "chematic v$VER" README_ja.md || { echo "MISMATCH: README_ja.md (expect chematic v$VER)"; fail=1; }
grep -q "version: $VER"  CITATION.cff || { echo "MISMATCH: CITATION.cff (expect version: $VER)"; fail=1; }
grep -q "v$VER | Yes"    SECURITY.md  || { echo "MISMATCH: SECURITY.md (expect v$VER | Yes)"; fail=1; }
grep -q "\"version\": \"$VER\"" demo/pkg/package.json || { echo "MISMATCH: demo/pkg/package.json (expect version $VER)"; fail=1; }
for f in crates/*/Cargo.toml; do
    if grep -q 'path = "\.\./chematic-' "$f" && ! grep -q "version = \"$VER\"" "$f"; then
        echo "MISMATCH: $f (path-dependency versions not bumped to $VER)"; fail=1
    fi
done
[ $fail -eq 0 ] && echo "Version consistent: $VER" || { echo "Run: python scripts/bump_version.py"; exit 1; }
echo "All checks passed."
