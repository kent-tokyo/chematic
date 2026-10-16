# Benchmark method

chematic benchmarks are specific to a package or source revision, comparator,
corpus, host, runtime, options, and timed operation. A speed result counts only
when the stated work and outputs are equivalent.

The current release is **v1.1.6**. The v1.1.5 published-wheel chemistry
workflow completed all 11 hosted jobs; the latest durable checked-in
cross-platform packet remains v1.1.3. The broadest published-package RDKit
performance packet is v1.0.42. The v1.1.5 hosted run remains correctness and
failure-accounting evidence, not a v1.1.6 performance packet. See the
[record index](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/README.md) for current entry points
and the [complete inventory](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/archive-index.md) for all history.

## Required provenance

Every record must identify:

- chematic package hash or source revision;
- comparator package/build and version;
- corpus source, filtering, and hash;
- operation, options, warm-up, repetitions, and timeout;
- hardware, OS, runtime, thread count, and relevant compiler flags;
- counts for success, failure, refusal, unsupported, and skipped rows;
- exact, numeric, or semantic output-equivalence rule;
- raw paired observations and uncertainty calculation.

Registry availability, a source build, and a public package are separate facts.
Do not reuse a result across those boundaries.

## Timed boundaries

Report these separately when applicable:

1. cold process or browser startup;
2. package/WASM initialization;
3. parse-inclusive pipeline;
4. first operation after parsing;
5. prepared or repeated operation;
6. serialization or host-language conversion;
7. peak process RSS and, where measurable, library allocation.

Moving perception or caching work into parsing changes the question; it does
not make the prepared operation a fair replacement for a parse-inclusive one.

## Output gate

Before ranking speed:

- compare every successful output under the declared contract;
- include comparator errors and chematic refusals in the denominator;
- reject asymmetric work or hidden skipped rows;
- classify representation-only differences separately from graph or numeric
  differences;
- keep typed refusals visible rather than counting them as matches.

For chemistry with several acceptable answers, define the semantic observer
before running the benchmark. Canonical-SMILES string equality alone is not a
general identity test.

## Statistical gate

An equivalent-output speed win requires:

- at least 21 alternating-order blocks for a promoted claim;
- a win in every block;
- a paired 95% confidence lower bound above `1.0`;
- no output or failure-accounting regression;
- two recorded hosts and exact public artifacts for a broad release claim.

Short exploratory runs may guide development, but their record must say that
they are exploratory and must not use superiority language.

## 3D and force fields

Separate atom typing, parameter coverage, same-coordinate energy, gradient,
convergence, stereo, clash, conformer quality, time, and memory. Faster energy
evaluation is not faster optimization, and a better median cannot hide a
large quality tail. Inspected corpora are development data and cannot become
sealed holdouts again.

## File formats and streaming

Compare semantic graph, charge, isotope, bond order, stereo, coordinates,
properties, and format metadata before throughput. For batch readers, require:

`input = success + refusal + unsupported + invalid + internal_error`

Keep cold start, first record, steady-state throughput, round trip, large-file
behavior, and peak memory separate.

## Repository checks

```bash
python3 scripts/check_benchmark_index.py
python3 scripts/check_release_docs_consistency.py
python3 scripts/check_roadmap_disposition.py
```

Each dated record contains its own reproduction command and validation script.
Use that record rather than a moving command copied into an overview page.
