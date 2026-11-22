# Benchmark method

chematic benchmarks are specific to a package or source revision, comparator,
corpus, host, runtime, options, and timed operation. A speed result counts only
when the stated work and outputs are equivalent.

Start with the [benchmark index](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/README.md)
or the [historical inventory](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/archive-index.md).
For chemistry and release qualification, see [validation](validation.md).

| Question | Record | Scope |
|---|---|---|
| v1.2.2 published Python timing | [Descriptor and MCS comparison](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-11-21-v122-published-comparison.md) | One macOS arm64 host and five pinned wheels |
| Labute ASA timing | [v1.1.8 comparison](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-19-v118-published-labute-asa-cosmolkit.md) | One descriptor and one host |
| Published Linux and macOS timing | [v1.0.42 comparison](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-11-v1042-published-rdkit-packet.md) | Three chemistry lanes on two hosts |
| Open Babel file I/O | [Format comparison](openbabel-file-io.md) | File interchange |

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

### Open Babel comparison scope

Open Babel is included as the principal file-interchange comparator. It is
most useful for SDF/MOL, V3000, MOL2, CML, CDXML, PDB, mmCIF, and related
streaming or round-trip work. It is not inserted into descriptor, fingerprint,
reaction, or 3D tables unless the two libraries expose equivalent operations
and their outputs pass a declared equality gate.

The current Open Babel evidence has two distinct boundaries:

- the v1.0.42 records are historical, immutable package/source-era baselines;
- newer V2000, V3000, SDF, and stereo records are one-host source comparisons
  against Open Babel 3.2.1 and do not inherit published-package status.

Only rows with equivalent semantic graph, charge, isotope, bond order, stereo,
coordinates, properties, and required format metadata receive a speed rank.
Rows with output differences remain visible and unranked. Start with the
[Open Babel comparison program](openbabel-file-io.md) and the
[record index](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/README.md);
do not combine ratios from different
revisions, corpora, formats, or timing boundaries.

## Repository checks

```bash
python3 scripts/checks/check_benchmark_index.py
python3 scripts/checks/check_release_docs_consistency.py
python3 scripts/checks/check_roadmap_disposition.py
```

Each dated record contains its own reproduction command and validation script.
Use that record rather than a moving command copied into an overview page.
