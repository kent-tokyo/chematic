#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 3 ]; then
  echo "usage: $0 <manifest> <signature> <public-key>" >&2
  exit 2
fi

manifest=$1
signature=$2
public_key=$3

# The public key is supplied by the release verifier. This checks the exact
# manifest bytes; it does not establish key custody or registry publication.
openssl dgst -sha256 -verify "$public_key" -signature "$signature" "$manifest"
