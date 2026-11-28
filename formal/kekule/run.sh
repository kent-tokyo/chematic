#!/usr/bin/env bash
# Differential test: Rust `kekulize` vs. the Lean 4 specification.
# Needs `lake` (Lean 4, version in lean/lean-toolchain) and cargo on PATH.
# Exit status is non-zero if Lean flags any Rust result as a bug.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="$(cd "$here/../.." && pwd)"
export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$root/target}"
random_cases="${KEKULE_RANDOM_CASES:-400000}"

(cd "$here/lean" && lake build)
(cd "$here/harness" && cargo build --release --quiet)
harness="$CARGO_TARGET_DIR/release/kekule-diff-harness"
check="$here/lean/.lake/build/bin/kekule-check"

status=0
for seed in 1 2 3 4 5; do
  echo "== random seed $seed ($random_cases cases)"
  "$harness" random "$seed" "$random_cases" | "$check" || status=1
done
for corpus in \
  scripts/fixtures/chembl_accuracy_corpus_4999.smi \
  scripts/fixtures/nci_first_5k_smiles_only.smi \
  scripts/fixtures/descriptor_census_corpus.smi \
  validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  validation/cosmolkit_comparison/corpora/stress-v1.smi \
  validation/rdkit_accuracy_split_v2/development.smi \
  validation/rdkit_accuracy_split_v2/exposed_holdout.smi; do
  echo "== $corpus"
  "$harness" smiles "$root/$corpus" | "$check" || status=1
done
exit "$status"
