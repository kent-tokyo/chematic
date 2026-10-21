# Benchmark method

chematic benchmarks are specific to a package or source revision, comparator,
corpus, host, runtime, options, and timed operation. A speed result counts only
when the stated work and outputs are equivalent.

The current release is **v1.2.0**. Its three-platform runtime smoke and
11-job public-wheel chemistry workflow passed, but its performance changes
remain source evidence until the exact public artifacts are timed. The latest
exact public-artifact performance record is the one-host macOS arm64 v1.1.8
Labute ASA comparison. The latest broad, checked-in cross-platform
public-wheel chemistry packet remains v1.1.3, and the broadest
published-package RDKit performance packet remains v1.0.42.
Newer source records do not inherit either package claim. See the
[record index](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/README.md)
for current entry points and the
[complete inventory](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/archive-index.md)
for retained history.

## Current evidence map

| Question | Best current entry | Evidence boundary |
|---|---|---|
| Does a current public wheel reproduce a measured speed result? | [v1.1.8 Labute ASA](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-19-v118-published-labute-asa-cosmolkit.md) | One macOS arm64 wheel, one host, one descriptor |
| What is the latest hosted public-wheel chemistry result? | [v1.2.0 chemistry run](https://github.com/kent-tokyo/chematic/actions/runs/39511354905) | Eleven jobs passed; this is pass/fail chemistry evidence, not performance evidence or a checked-in row packet |
| What is the broadest checked-in public-wheel chemistry packet? | [v1.1.3 artifact packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-14-v113-published-artifact-packet.md) | Four Python wheels with retained per-row evidence |
| What is the broadest published RDKit performance packet? | [v1.0.42 comparison](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-11-v1042-published-rdkit-packet.md) | Historical exact artifacts; do not transfer results to newer releases |
| What is the Open Babel comparison boundary? | [file-I/O program](openbabel-file-io.md) | File interchange only; source and historical package records are separate |
| Where are current source experiments? | [benchmark record index](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/README.md) | Development evidence unless an exact public artifact is named |

`current release`, `latest source`, and `broadest public packet` therefore
refer to different artifacts. Overview pages must not collapse them into one
claim.

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
python3 scripts/check_benchmark_index.py
python3 scripts/check_release_docs_consistency.py
python3 scripts/check_roadmap_disposition.py
```

Each dated record contains its own reproduction command and validation script.
Use that record rather than a moving command copied into an overview page.
