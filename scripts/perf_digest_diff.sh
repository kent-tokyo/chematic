#!/usr/bin/env bash
# Differential output check for performance changes.
#
# Builds tools/perf_digest twice — against BASE_REV (a temporary git worktree)
# and against the current working tree — runs `digest` over the given corpora
# and compares every (corpus, row, op) output byte-for-byte. Writes a compact
# JSON summary (per-op row/diff counts, corpus and digest SHA-256, revisions);
# the raw digests (up to gigabytes) are deleted unless KEEP_DIGESTS=1.
#
# usage: scripts/perf_digest_diff.sh BASE_REV OUT.json CORPUS.smi [CORPUS.smi...]
#        ONLY=op1,op2 restricts the op set (default: all ops present in both).
#        KEEP_DIGESTS=1 keeps the raw base/head digests in the scratch directory.
set -euo pipefail

if [ "$#" -lt 3 ]; then
  echo "usage: $0 BASE_REV OUT.json CORPUS.smi [CORPUS.smi...]" >&2
  exit 2
fi
BASE_REV=$1; OUT=$2; shift 2
ROOT=$(git rev-parse --show-toplevel)
SCRATCH=$(mktemp -d "${TMPDIR:-/tmp}/perf-digest.XXXXXX")
BASE_SHA=$(git -C "$ROOT" rev-parse "$BASE_REV^{commit}")
HEAD_SHA=$(git -C "$ROOT" rev-parse HEAD)
DIRTY=$(git -C "$ROOT" status --porcelain --untracked-files=no -- crates | wc -l | tr -d ' ')

cleanup() {
  git -C "$ROOT" worktree remove --force "$SCRATCH/base-tree" >/dev/null 2>&1 || true
  # Raw digests can be gigabytes; keep them only on request.
  [ -n "${KEEP_DIGESTS:-}" ] || rm -f "$SCRATCH/base.tsv" "$SCRATCH/head.tsv"
}
trap cleanup EXIT

git -C "$ROOT" worktree add --detach -q "$SCRATCH/base-tree" "$BASE_SHA"

make_harness() { # $1 = harness dir, $2 = crates root, $3 = extra cargo flags
  mkdir -p "$1"
  cp -R "$ROOT/tools/perf_digest/src" "$1/"
  sed "s#\.\./\.\./crates/#$2/#g" "$ROOT/tools/perf_digest/Cargo.toml" > "$1/Cargo.toml"
  (cd "$1" && cargo build --release -q $3)
}
make_harness "$SCRATCH/h-base" "$SCRATCH/base-tree/crates" ""
make_harness "$SCRATCH/h-head" "$ROOT/crates" ""

# Absolute corpus paths so both binaries read the same files.
CORPORA=()
for c in "$@"; do CORPORA+=("$(cd "$(dirname "$c")" && pwd)/$(basename "$c")"); done

"$SCRATCH/h-base/target/release/chematic-perf-digest" digest "$SCRATCH/base.tsv" "${CORPORA[@]}" &
"$SCRATCH/h-head/target/release/chematic-perf-digest" digest "$SCRATCH/head.tsv" "${CORPORA[@]}" &
wait

python3 - "$SCRATCH" "$OUT" "$BASE_SHA" "$HEAD_SHA" "$DIRTY" "$ROOT" "${CORPORA[@]}" <<'PY'
import collections, hashlib, json, os, sys
scratch, out, base_sha, head_sha, dirty, root, *corpora = sys.argv[1:]
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()
rows = collections.Counter(); diffs = collections.Counter(); first = {}
# Stream both digests: they can be gigabytes (full SMARTS match maps).
with open(os.path.join(scratch, "base.tsv"), encoding="utf-8") as fa, \
     open(os.path.join(scratch, "head.tsv"), encoding="utf-8") as fb:
    for a, b in zip(fa, fb, strict=True):
        op = a.split("\t", 3)[2]
        rows[op] += 1
        if a != b:
            diffs[op] += 1
            first.setdefault(op, {"base": a[:400], "head": b[:400]})
summary = {
    "schema": "chematic-perf-digest-diff/v1",
    "base_revision": base_sha,
    "head_revision": head_sha,
    "head_crates_dirty_files": int(dirty),
    "only": os.environ.get("ONLY"),
    "corpora": [{"path": os.path.relpath(c, root), "sha256": sha(c)} for c in corpora],
    "digest_sha256": {"base": sha(os.path.join(scratch, "base.tsv")), "head": sha(os.path.join(scratch, "head.tsv"))},
    "ops": {op: {"rows": rows[op], "diffs": diffs[op]} for op in sorted(rows)},
    "total_rows": sum(rows.values()),
    "total_diffs": sum(diffs.values()),
    "first_diff_examples": first,
}
with open(out, "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=1)
    f.write("\n")
for op in sorted(rows):
    print(f"{op:20s} rows {rows[op]:7d} diffs {diffs[op]}")
print(f"total rows {summary['total_rows']} diffs {summary['total_diffs']}")
sys.exit(1 if summary["total_diffs"] else 0)
PY
