#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 4 ]; then
  echo "usage: $0 <manifest> <signature> <public-key> <expected-sha256-fingerprint>" >&2
  exit 2
fi

manifest=$1
signature=$2
public_key=$3
expected_fingerprint=$4

actual_fingerprint=$(
  openssl pkey -pubin -in "$public_key" -outform DER \
    | sha256sum | awk '{print tolower($1)}'
)
expected_fingerprint=$(printf '%s' "$expected_fingerprint" | tr '[:upper:]' '[:lower:]')
if [ "$actual_fingerprint" != "$expected_fingerprint" ]; then
  echo "release public-key fingerprint mismatch" >&2
  echo "expected: $expected_fingerprint" >&2
  echo "actual:   $actual_fingerprint" >&2
  exit 1
fi

openssl dgst -sha256 -verify "$public_key" -signature "$signature" "$manifest"
echo "release key fingerprint verified: $actual_fingerprint"
