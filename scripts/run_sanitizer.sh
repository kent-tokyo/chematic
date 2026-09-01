#!/usr/bin/env bash
# Run the sanitizer-backed Rust test suite locally or in CI.
# Usage: bash scripts/run_sanitizer.sh [address|leak|thread]
set -euo pipefail

sanitizer="${1:-address}"
case "$sanitizer" in
    address|leak|thread) ;;
    *)
        echo "usage: $0 [address|leak|thread]" >&2
        exit 2
        ;;
esac

if [[ "$(uname -s)" != "Linux" ]]; then
    echo "sanitizer runs require the Linux nightly toolchain; current host: $(uname -s)" >&2
    exit 2
fi

target="x86_64-unknown-linux-gnu"
export RUST_BACKTRACE=1
export RUSTFLAGS="-Zsanitizer=${sanitizer}"

cargo +nightly test \
    -Zbuild-std \
    --target "$target" \
    -p chematic-core \
    -p chematic-smiles \
    -p chematic-mol \
    --lib
