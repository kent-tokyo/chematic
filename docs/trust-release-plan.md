# chematic 1.x Trust Release rules

Updated 2026-10-10 for the **v1.0.41** release line.

This document defines release rules. Priorities are in the
[roadmap](https://github.com/kent-tokyo/chematic/blob/main/ROADMAP.md), open
dependencies in the [open-work ledger](roadmap-open-work.md), and measurements
in [validation](validation.md).

## Objective

A Trust Release must make five things inspectable:

1. which APIs are stable, experimental, or unsupported;
2. which versions, artifacts, corpora, options, and failure policy were tested;
3. whether Rust, Python, Node, and WASM preserve result accounting;
4. whether malformed or uncertain chemistry fails with a typed outcome;
5. whether evidence comes from source, a built package, or a public channel.

## Required contracts

### Compatibility

Each comparison record must contain the chematic and comparator versions or
hashes, corpus identity, operation and options, supported domain, comparison
type, and counts for success, failure, refusal, and skipped inputs.

### Bindings

For batch APIs:

`input_count = success + failed + refused + skipped`

Results retain original indices, stage, and typed reason. Cancellation must
identify unprocessed input. `failed == 0` does not mean every row succeeded.

### Parser and runtime safety

Malformed-input tests run with explicit time, memory, depth, and size bounds.
The required outcome is zero panic/crash, terminal accounting for every input,
and typed errors. Rust alone is not safety evidence.

### Stereo and interchange

Regression suites cover atom-order and spelling permutations, file round trips,
CIP/E/Z, and declared representation limits. Lossy writers expose report or
strict modes; unsupported chemistry is never silently promoted to parity.

## Release-candidate gate

- Workspace tests, clippy, binding tests, documentation checks, dependency and
  license checks pass in the declared environments.
- Corrected behavior has focused regression tests and any required full-corpus
  rerun.
- README files, CHANGELOG, validation summary, benchmark index, package
  versions, and release metadata agree.
- Source-only measurements are not described as package or public-channel
  results.
- Tags, registry publication, GitHub Release, docs, and Pages are verified
  independently.

Open issues may remain when they are explicitly outside the release boundary.
A release note must not imply they are complete.

## Long-running gates

If a gate is stopped, record its command, source SHA, completed and unfinished
scope, effect on the conclusion, and restart command. An interrupted run is not
a pass.

## Records

- [Roadmap](https://github.com/kent-tokyo/chematic/blob/main/ROADMAP.md)
- [Open work](roadmap-open-work.md)
- [Compatibility scope](compatibility-scope.md)
- [Validation summary](validation.md)
- [Benchmark methodology](benchmark.md)
- [Dated benchmark records](https://github.com/kent-tokyo/chematic/tree/main/benchmarks)
- [Release history](https://github.com/kent-tokyo/chematic/blob/main/CHANGELOG.md)

Detailed historical plans remain in Git history and are not repeated here.
