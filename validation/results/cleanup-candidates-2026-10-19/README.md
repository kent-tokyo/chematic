# Sparse cleanup candidates — 2026-10-19

Motivation: [RDKit #9475](https://github.com/rdkit/rdkit/issues/9475), confirmed
open. Base `b7004777daf39a5f6cef03f56c4dcb68c12433f3`; candidate source/test/fixture
hashes, compiler/platform and settings are in `timing.json`. Its `checkout_head`
identifies the base before committing the measured working-tree change.

The production change replaces dense candidate flags plus a second atom scan
with ordered inline indices (capacity 8), and borrows unchanged adjacency in
both nitrogen passes. There is no persistent cache or public API change.
First-pass eligibility remains fixed even after charge/bond-order updates.

## Correctness

- Existing 1,131-row canonical fixture: **1,130 raw models** agree completely
  with the frozen former `clean_up`; one malformed ring-closure input is
  explicitly recorded as refused. This is raw projection, not a claim that all
  1,130 inputs subsequently pass full sanitization.
- Each model is compared after **two successive calls**, including returned
  `Result` and the whole internal `Mol` debug state: **2,260 fixture state
  comparisons**. The original expected canonical-output fixture is unchanged.
- 22 authored controls, also tested with toggled aromatic flags/cleared explicit
  valences, plus empty/inline-capacity/spill/large sparse generator controls.
  `N(=O)#N` explicitly checks that its second-phase rewrite survives the
  first-phase charge change. See `controls.json` for provenance, license,
  minimization, oracle and normalized hashes. No upstream sample is imported.
- Complete crate library/integration tests: **572 passed, 0 failed, 2 ignored**
  across 28 binaries. Three new functional tests are in the existing library
  target, so ordinary CI runs them; timing is explicitly opt-in.
- Release all-target Clippy, warnings denied: passed. Workspace Rust formatting,
  Python Black, product names, refactor budget and pinned MkDocs strict build:
  passed.

## Stage-only timing

Same release test binary contains the old and new functions. Each call receives
identical freshly cloned raw models, with preparation/cloning and destruction
outside the timer. Both paths are warmed; nine paired rounds alternate order,
300 passes per round.

| Packet | Models | Median old | Median new | Time reduction |
| --- | ---: | ---: | ---: | ---: |
| Existing fixture | 1,130 | 11.138377 ms | 4.814780 ms | 56.7731% |
| Authored controls | 22 | 0.405923 ms | 0.175412 ms | 56.7869% |

These are one-host cleanup-stage totals. No end-to-end parsing/sanitization
gain, performance against RDKit, representative Enamine REAL packet, physical
chemistry validation or published-package result is claimed.

## Reproduction

```sh
python3 scripts/measure_cleanup_candidates.py \
  --target-dir /private/tmp/chematic-conjugation-target \
  --output /private/tmp/chematic-cleanup-reproduction
```

For the full tests and Clippy, use the same `cargo_settings` environment from
`timing.json`:

```sh
cargo test -p chematic-smiles --offline --release --lib --tests
cargo clippy -p chematic-smiles --offline --release --all-targets -- -D warnings
python3.13 -m mkdocs build --strict --site-dir /private/tmp/chematic-cleanup-docs
```

Logs retain differential refusals, every timing round, all test summaries and
Clippy output. Timing source hashes pin the exact candidate; hosted CI, main
integration and package publication remain separate gates.
