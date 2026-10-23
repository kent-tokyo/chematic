#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 3 ]; then
  echo "usage: $0 <manifest> <private-key> <signature>" >&2
  exit 2
fi

manifest=$1
private_key=$2
signature=$3

# The private key is supplied by the release environment and is never copied
# into the repository. CI should publish the matching public key separately.
openssl dgst -sha256 -sign "$private_key" -out "$signature" "$manifest"
echo "provenance signature written: $signature"
