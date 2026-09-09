#!/usr/bin/env bash
set -euo pipefail

# cargo-deny expands deny.toml's `~/.cargo/advisory-db` independently from
# CARGO_HOME. Keep the user's existing registry/cache and Rust toolchain, but
# give the advisory database a writable temporary HOME so the lock is local.
ORIGINAL_HOME=${HOME:?HOME must be set}
ORIGINAL_CARGO_HOME=${CARGO_HOME:-"$ORIGINAL_HOME/.cargo"}
ORIGINAL_RUSTUP_HOME=${RUSTUP_HOME:-"$ORIGINAL_HOME/.rustup"}
SOURCE_DB="$ORIGINAL_HOME/.cargo/advisory-db"

if [[ ! -d "$SOURCE_DB" ]]; then
  echo "advisory database not found: $SOURCE_DB" >&2
  exit 2
fi

TEMP_ROOT=$(mktemp -d /private/tmp/chematic-cargo-deny.XXXXXX)
trap 'rm -rf "$TEMP_ROOT"' EXIT
mkdir -p "$TEMP_ROOT/.cargo"
cp -R "$SOURCE_DB" "$TEMP_ROOT/.cargo/advisory-db"

HOME="$TEMP_ROOT" \
CARGO_HOME="$ORIGINAL_CARGO_HOME" \
RUSTUP_HOME="$ORIGINAL_RUSTUP_HOME" \
cargo deny check --disable-fetch
