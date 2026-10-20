# SDF reader correctness and throughput follow-up

The native source candidate beats Open Babel 3.2.1 in all three lanes on the
additional fixed, unfiltered 1,000-record V2000 2D SDF slice. Internal canonical
graphs and every declared rewrite field agree with RDKit 2026.09.1. Every lane
wins all 21 rotating blocks and has a paired 95% lower bound above 1.0.
The original slice is retained: its candidate reader now agrees on all records,
while Open Babel still differs on two internal graphs and cannot receive a
parser speed ranking. This is bounded source evidence on one host.

## Final measurements

Times are median milliseconds per 1,000 records, from five full-batch repeats
per native process. Speedups are medians of paired Open Babel / candidate ratios,
with deterministic 10,000-resample paired bootstrap intervals. They need not
equal ratios of independent median times.

| Slice | Lane | chematic ms | Open Babel ms | Paired speedup [95% CI] | Faster blocks |
|---|---|---:|---:|---|---:|
| 0 | parse | 29.570 | 92.906 | Not ranked: Open Babel internal graphs differ | — |
| 0 | write | 3.822 | 23.436 | 6.171x [6.101, 6.201] | 21/21 |
| 0 | roundtrip | 33.571 | 116.489 | 3.468x [3.453, 3.492] | 21/21 |
| 5000 | parse | 7.598 | 41.110 | 5.399x [5.377, 5.442] | 21/21 |
| 5000 | write | 1.643 | 13.539 | 8.216x [8.122, 8.292] | 21/21 |
| 5000 | roundtrip | 9.391 | 55.410 | 5.916x [5.859, 5.949] | 21/21 |

All processes run sequentially; baseline, candidate and Open Babel positions
rotate evenly. Timers exclude file loading, startup, plugin discovery and
preparation. Both engines retain a full parsed batch. The Rust reader includes
its default stereo and geometry diagnostics. Prepared writing follows one
untimed warm write and uses the normal in-memory writer; roundtrip reads a new
batch and writes it within the timer. This is native library throughput, without
Python call overhead or filesystem output costs.

| Slice | Lane | chematic peak RSS median MiB | Open Babel peak RSS median MiB |
|---|---|---:|---:|
| 0 | parse | 15.77 | 52.78 |
| 0 | write | 10.88 | 39.69 |
| 0 | roundtrip | 15.78 | 52.89 |
| 5000 | parse | 11.38 | 41.36 |
| 5000 | write | 7.97 | 33.31 |
| 5000 | roundtrip | 11.38 | 41.41 |

Darwin `time -l` reports whole-process high-water RSS, including input, prepared
batch, warmup, plugins and allocator state. These figures do not establish
operation-only allocation costs or constant-memory streaming behavior.

## Source changes and agreement

The follow-up changes only local CTAB stereo perception in production, on top
of the [previous writer candidate](2026-10-14-openbabel-sdf-speed.md).
Degree and wedge eligibility are checked before collecting neighbors, and
eligible centers use a four-element inline buffer. Ordinary atoms avoid a
temporary neighbor allocation.

The reader also rejects an ambiguous CTAB center when two flat neighbor bonds
point in nearly the same direction. It compares squared distances between unit
directions against `1e-3`, following the overlap criterion in
[RDKit's stereo source](https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/Chirality.cpp).
Only owned wedges leave the drawing plane; wedge/flat XY overlap remains valid.
Generic non-CTAB local parity keeps its previous behavior. This is a narrow
geometry guard, not an implementation of RDKit's complete stereo algorithm.
The source link is upstream master; the recorded external oracle is the
installed RDKit 2026.09.1 package.

Original zero-based row 745 previously acquired an unwanted tetrahedral tag.
It now reports `DegenerateGeometry`, leaves that center unassigned, and retains
the other stereo centers. The exact original MOL block is a regression fixture.
Across both timing slices, baseline and candidate internal snapshots are byte
identical, including atoms, bonds, neighbor orders, coordinate bits, diagnostics,
metadata and properties. Both corrected builds match all 2,000 reference graphs.
Open Babel matches 998/1,000 original and 1,000/1,000 additional internal graphs;
the original differing rows remain 90 and 745. Nothing is filtered.

RDKit independently observes canonical isomeric graph, indexed atoms and bonds,
finite coordinates within `1e-4`, name, raw SD properties and explicit unknown
double-bond markers. Each engine's rewrites pass every declared field on all
2,000 records. Additional safe `either` markers are accepted, while explicit
input markers must survive. Output agreement does not override the separate
internal parser gate. Retained snapshots from 18 additional synthetic cases
also match the external oracle: 1-degree overlaps remain unassigned and
3-degree separations remain assigned across three scales and three rotations.

## Incremental optimization result

The measured baseline already includes the same overlap correctness fix and
previous writer changes. It retains the old unconditional neighbor `Vec`, so
this comparison isolates the allocation change without mixing different
chemical outputs.

| Slice | Lane | Baseline / candidate median [95% CI] | Faster blocks | Strict win |
|---|---|---|---:|---|
| 0 | parse | 1.009x [1.007, 1.013] | 19/21 | No |
| 0 | write | 0.988x [0.975, 0.993] | 3/21 | No |
| 0 | roundtrip | 1.007x [1.004, 1.013] | 18/21 | No |
| 5000 | parse | 1.026x [1.023, 1.027] | 21/21 | Yes |
| 5000 | write | 1.007x [0.999, 1.016] | 14/21 | No |
| 5000 | roundtrip | 1.016x [1.010, 1.020] | 18/21 | No |

The incremental optimization is not a universal speed win. Original-slice
prepared writing is approximately 1.2% slower in this measurement despite no
writer source change. The candidate retains strict Open Babel wins in all
eligible lanes. The earlier immutable packet and its rejected coordinate-parser
experiment are preserved unchanged.

## Provenance, validation and reproduction

Base commit: `ffe682c32cb4c76de0d309649427202a7d8e956d` (source version 1.1.3).
Apply the previous writer candidate patch, then either the follow-up baseline
or candidate patch in separate isolated checkouts. Both alternatives were
reconstructed and their source hashes verified. All 1,098 tracked build files
under `crates`, plus the workspace manifest and lockfile, were checked against
the base: only the MOL writer and local stereo perception differ. The added
harness, fixtures and tests are recorded separately. Other dirty workspace
changes, including Python caches and 3D work, are excluded.

Apply the small runner-dependency patch when reproducing the recorded runner
file hashes. It retains the already-present COSMolKit API spelling correction
in the imported statistics module; that adapter is not executed by this SDF
runner. The paired bootstrap function is unchanged. All 11 source-manifest paths
were reconstructed for both alternatives, including this dependency.

Build the unchanged native example in release mode with the recorded Rust
1.97.0 toolchain, opt-level 3, LTO and one codegen unit. The retained Open Babel
C++ harness uses the same installed 3.2.1 library and plugin as the earlier
packet; executable and library hashes are checked. Rebuild with fresh source
mtimes or a fresh target directory. Run the existing runner for each decompressed
corpus with `--blocks 21 --repeats 5 --records 1000` and this follow-up source
manifest. Source, binary and corpus hashes are checked before and after timing.

The exact corpora are retained:

- [Original slice](../validation/benchmark_corpora/openbabel-sdf-speed-1000-v1.sdf.gz):
  source rows 0–999, 27,435 atoms, 2,419,258 bytes, SHA256
  `8e7ad63db96f6fb87b87333ab1f7710a18bc3be5cb4676c15c345c09d4ce718b`.
- [Additional slice](../validation/benchmark_corpora/openbabel-sdf-speed-slice5000-1000-v1.sdf.gz):
  source rows 5000–5999, 15,291 atoms, 1,390,916 bytes, SHA256
  `38d1d169f40de9fe4c95ec475035658ffcb27a2c84c448e3c619c32aefe8fd60`.

Both are fixed slices of the existing browser 10k SMILES corpus, generated by
RDKit 2026.09.1 with identical coordinate settings. The generator now accepts
`--offset`; no validity or success filtering occurs. Both slices were fixed
before the follow-up timing. They are exposed development data, not sealed
holdouts. A subsequent attempt to generate all 10,000 V2000 records was refused
at source row 8341 because RDKit emits V3000 for that metal/dative-bond input.
The [failed generation record](2026-10-14-openbabel-sdf-read-followup-full10k.json)
retains the input, emitted block and failure without filtering or running a
full-corpus benchmark. This is a generator format boundary, not a reader result.

Isolated release tests for `chematic-perception` and `chematic-mol`: 1,225 passed,
zero failed, two ignored. Four added regression tests pass. Release Clippy for
both crates, all targets and warnings denied, passes. Five runner integrity
tests pass; changed Rust files pass rustfmt and the workspace diff has no
whitespace errors. Rewrites, reproduced hash-identical snapshots, oracle probes
and build/test logs are retained as lossless gzip artifacts. Archival additions
to the reports do not change raw timing rows.

- [Original-slice timing blocks](2026-10-14-openbabel-sdf-read-followup-slice0-paired21.json)
- [Additional-slice timing blocks](2026-10-14-openbabel-sdf-read-followup-slice5000-paired21.json)
- [Source/compiler/library/corpus manifest](2026-10-14-openbabel-sdf-read-followup-source.json)
- [Artifact hashes, reconstruction and oracle verification](2026-10-14-openbabel-sdf-read-followup-verification.json)
- [Corrected baseline patch](2026-10-14-openbabel-sdf-read-followup-baseline.patch)
- [Optimized candidate patch](2026-10-14-openbabel-sdf-read-followup-candidate.patch)
- [Exact runner dependency patch](2026-10-14-openbabel-sdf-read-followup-runner-dependency.patch)
- [Prior writer patch](2026-10-14-openbabel-sdf-speed-candidate.patch)
- [Native Rust harness](../crates/chematic-mol/examples/openbabel_sdf_throughput.rs),
  [native Open Babel harness](../tools/openbabel_sdf_throughput.cpp),
  [runner](../scripts/bench_openbabel_sdf_throughput.py), and
  [corpus generator](../scripts/generators/build_openbabel_sdf_speed_corpus.py)

This packet establishes the listed native 2D V2000 SDF wins. Other formats,
larger and fresh corpora, malformed inputs, cold startup, Python/CLI overhead,
streaming memory, 3D workflows, additional hosts, and published package artifacts
still need their own equivalent-output comparisons. General Open Babel speed
superiority is not established.
