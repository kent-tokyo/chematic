#!/usr/bin/env bash
set -euo pipefail

# Inputs are reconstructed before scoring and checked against the preparation run.
os="$1"
dataset="$2"
source_commit="29e37e169668903d6b01490c565b2f4a2cf63e5e"
case "$os" in
  ubuntu-latest) expected_wheel_sha256="da93ee39d81727e3566a65b352657f745ed3ed91aa7fffd8a4dcbd784d8e8edb" ;;
  macos-latest) expected_wheel_sha256="416a8fe0ffce006259c74106f2f4d90849f472bacb94377dd779103c60c64693" ;;
  *) echo "unknown platform: $os" >&2; exit 1 ;;
esac
mkdir -p evaluation/wheel evaluation/reference
gh run download 45024603212 --repo "$GITHUB_REPOSITORY" \
  --name "fresh-evaluation-wheel-${os}-${source_commit}" --dir evaluation/wheel
test "$(cat evaluation/wheel/candidate/source-commit.txt)" = "$source_commit"
test "$(cat evaluation/wheel/candidate/wheel-sha256.txt)" = "$expected_wheel_sha256"
wheel=$(find evaluation/wheel -name '*.whl' -type f)
test -n "$wheel"
python3 - "$wheel" evaluation/wheel/candidate/wheel-sha256.txt <<'PY'
import hashlib
import pathlib
import sys
wheel, expected = map(pathlib.Path, sys.argv[1:])
assert hashlib.sha256(wheel.read_bytes()).hexdigest() == expected.read_text().strip()
PY
case "$dataset" in
  cod)
    gh run download 45024603212 --repo "$GITHUB_REPOSITORY" \
      --name "fresh-cod-reference-${source_commit}" --dir evaluation/reference
    echo '6e2f93fb05c8ff49aff24577aa62d3f04f6531b505e12a9f711c9b79cf3ba008  evaluation/reference/reference.sdf' | shasum -a 256 --check
    echo 'c9d185e075f83dacc48620d1b843727cbaf6729406aee0568fada9e30bfbff9d  evaluation/reference/split.json' | shasum -a 256 --check
    echo '03111f28aa56e26a0e895d9f55128c2e29fce6060502fbaa123f6f4f9c9e3f63  evaluation/reference/audit.json' | shasum -a 256 --check
    ;;
  macrocycle)
    gh run download 44955399937 --repo "$GITHUB_REPOSITORY" \
      --name macrocycle-pdb-prepared-3b3a5f827788c39195ccc58608e854e6ab543557 --dir evaluation/reference
    mv evaluation/reference/macrocycle-prepared.sdf evaluation/reference/reference.sdf
    mv evaluation/reference/macrocycle-prepared-split.json evaluation/reference/split.json
    mv evaluation/reference/macrocycle-prepared-invalid-audit.json evaluation/reference/audit.json
    echo '71bf05a54cd402146e8892acd88969ec52fd45d2a02c4194b7aa905fe900222d  evaluation/reference/reference.sdf' | shasum -a 256 --check
    echo 'ffb256b0a9dbc779a5c6b3dd2ac96deec6174dc25408dc1d0691efcef8c4c4ff  evaluation/reference/split.json' | shasum -a 256 --check
    echo '53cf0e7ff0c1ef42dbe8ee043e7fda762c8d0c1bbecd477794cd55a9a44326f2  evaluation/reference/audit.json' | shasum -a 256 --check
    ;;
  *) echo "unknown dataset: $dataset" >&2; exit 1 ;;
esac
printf '%s\n' "$wheel" > evaluation/wheel-path.txt
