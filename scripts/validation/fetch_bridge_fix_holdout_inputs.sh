#!/usr/bin/env bash
set -euo pipefail

os="$1"
candidate_commit="261e520ab4d0fce681b23c168bee0888c2103b0a"
case "$os" in
  ubuntu-latest)
    run_id=45032082181
    artifact="mixed-development-wheel-${candidate_commit}"
    expected_wheel="cb9d6b6df61db6b6f4abb3326f315f122ceb80afdc9ef2de909800c689dbeba6"
    ;;
  macos-latest)
    run_id=45033920398
    artifact="bridge-fix-wheel-macos-${candidate_commit}"
    expected_wheel="83eaf6994439b98647c65d34322ccb50489d6d537f0b396cba4674be6b14d20c"
    ;;
  *) echo "platform is not yet pinned: $os" >&2; exit 1 ;;
esac
mkdir -p bridge-input/wheel bridge-input/reference
gh run download "$run_id" --repo "$GITHUB_REPOSITORY" --name "$artifact" --dir bridge-input/wheel
test "$(cat bridge-input/wheel/source-commit.txt)" = "$candidate_commit"
wheel=$(find bridge-input/wheel -name '*.whl' -type f)
test -n "$wheel"
python3 - "$wheel" "$expected_wheel" <<'PY'
import hashlib
import pathlib
import sys
assert hashlib.sha256(pathlib.Path(sys.argv[1]).read_bytes()).hexdigest() == sys.argv[2]
PY
gh run download 45024603212 --repo "$GITHUB_REPOSITORY" \
  --name fresh-cod-reference-29e37e169668903d6b01490c565b2f4a2cf63e5e --dir bridge-input/reference
echo '6e2f93fb05c8ff49aff24577aa62d3f04f6531b505e12a9f711c9b79cf3ba008  bridge-input/reference/reference.sdf' | shasum -a 256 --check
gzip -dc validation/results/cod-unused-force-field-holdout-split-2026-11-22.json.gz > bridge-input/split.json
gzip -dc validation/results/cod-unused-force-field-holdout-audit-2026-11-22.json.gz > bridge-input/audit.json
echo '9cbf6b107745be663dab2f0edf9dc2e48e372d3f90319c5fd2bcd348b97764d5  bridge-input/split.json' | shasum -a 256 --check
echo 'c748851eb0bae87d61158ec59c12e41ece76c5b6951a97869e917c4d6a385cce  bridge-input/audit.json' | shasum -a 256 --check
printf '%s\n' "$wheel" > bridge-input/wheel-path.txt
printf '%s\n' "$expected_wheel" > bridge-input/wheel-sha256.txt
