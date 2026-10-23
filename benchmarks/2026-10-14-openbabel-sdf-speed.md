# Diverse SDF batch throughput vs Open Babel 3.2.1

The isolated source candidate passes the final-output semantic gate on all
1,000 unfiltered V2000 2D SDF records. Prepared SDF rewriting and native
read-plus-rewrite both beat the installed Open Babel library in every one of
21 rotating blocks, with paired 95% lower bounds above 1.0. Internal parse
stereo differences remain, so parser timings have no equivalent-graph win claim.
This is one-host source evidence, not a released-package or toolkit-wide result.

## Final measurement

Times are medians in milliseconds per 1,000 records. Speedups are medians of
paired Open Babel / candidate ratios, with deterministic 10,000-resample paired
bootstrap intervals. A ratio of independent medians need not equal the paired
median ratio. Each native process times five repetitions of the full batch.

| Lane | chematic ms | Open Babel ms | Paired speedup [95% CI] |
|---|---:|---:|---|
| parse | 29.837 | 92.931 | Not ranked: internal parse graph differs |
| write | 3.788 | 23.378 | 6.151x [6.140, 6.232], 21/21 faster |
| roundtrip | 33.792 | 116.697 | 3.442x [3.435, 3.457], 21/21 faster |

Both engines retain the full parsed batch. Input file loading, process startup,
plugin discovery and preparation are outside the native timers. The default
Rust reader includes stereo/geometry diagnostics; it is not the fast supplier
path. Prepared rewriting uses each engine's parsed batch, one untimed warm
write and its normal in-memory SDF writer (reusable Rust String versus native
OBConversion WriteString). Round-trip timing reads a fresh batch and writes it.
The three engine orders rotate evenly; all processes run sequentially. This
compares native APIs, without Python call overhead or filesystem output costs.

Whole-process high-water RSS includes input loading, the prepared batch, warmup,
plugins, allocator state and timed work. It is measured by Darwin time -l and
is not an operation-only or constant-memory streaming measurement.

| Lane | chematic peak RSS median MiB | Open Babel peak RSS median MiB |
|---|---:|---:|
| parse | 15.78 | 52.78 |
| write | 10.92 | 39.72 |
| roundtrip | 15.81 | 52.89 |

## Meaning and source changes

RDKit 2026.09.1 independently observes indexed atoms (element, isotope, formal
charge, radical count and H count), indexed bonds and bond orders, canonical
isomeric graph, finite coordinates within 1e-4, name, raw SD properties and
explicit unknown-double-bond markers. Additional safe `either` markers are
accepted; an explicit source marker cannot disappear. Both candidate and Open
Babel rewritten SDFs pass all fields on 1,000/1,000 records. No input is excluded.

The baseline writer incorrectly makes the unspecified azo N=N on zero-based
row 766 definite E. The candidate marks undeclared stereogenic double bonds
`either` when supplied coordinates could otherwise invent E/Z. Declared stereo
continues through the existing depiction path. A small atom-row fast path also
writes the common zero mass/charge fields directly. The production change is
limited to the MOL2000 writer; unrelated dirty workspace changes were excluded.
This fixes the single rewritten-graph defect, rather than removing its input.
The baseline writer therefore performs semantically different work and has no
ranked write or round-trip speed comparison with the candidate.

Independent internal graph observations are intentionally separate from output
SDF observations. Rust canonical-SMILES observations differ from the reference
on row 745, where RDKit rejects overlapping-neighbor tetrahedral geometry.
Open Babel canonical-SMILES observations differ on rows 90 and 745. These are
observational parse-graph differences, not filtered inputs. Rewriting can hide
such a difference when a consumer reinterprets the drawing; the runner must
not use output-only agreement to declare a parser win. Rust baseline/candidate
snapshots match byte-for-byte, including atoms, bonds, neighbor orders, exact
coordinate bits, diagnostics, metadata and properties. The parser is unchanged.

A fixed-decimal coordinate parser and property-token allocation experiment was
removed after its 21-block development interval crossed parity (1.010x,
95% CI [0.997, 1.013], 14/21 faster). Its [raw development result](../validation/results/openbabel-sdf-speed-development-2026-10-14.json.gz)
and [source patch](../validation/results/openbabel-sdf-speed-development-2026-10-14.patch.gz)
are retained, with no parser optimization claim.

## Provenance and validation

The base is `ffe682c32cb4c76de0d309649427202a7d8e956d`. Baseline and candidate use
the identical Rust example instrumentation in isolated source trees. The
candidate applies only the writer change and its regression tests. Release
builds use Rust 1.97.0, opt-level 3, LTO and one codegen unit. The C++ harness
uses Apple clang 21.0.0, -O3 and the installed Homebrew Open Babel 3.2.1 library.
The host is macOS 27.0.1 arm64 with 10 logical CPUs. Library, plugin, executable,
corpus, source and patch hashes are recorded. This is not a published Rust crate,
Python wheel, npm package, or hosted CI result.

The corpus is the first 1,000 rows of the existing RDKit browser 10k SMILES set,
without validity or success filtering, converted by the retained RDKit 2026.09.1
generator to 2D V2000 SDF. It contains 27,435 atoms, 7–137 per record, names and
one `source_index` SD field per record. Uncompressed size is 2,419,258 bytes and
SHA256 is `8e7ad63db96f6fb87b87333ab1f7710a18bc3be5cb4676c15c345c09d4ce718b`.
The [exact compressed corpus](../validation/benchmark_corpora/openbabel-sdf-speed-1000-v1.sdf.gz)
is retained. This inspected development corpus is not a fresh holdout.

Final isolated `cargo test --release -p chematic-mol`: 970 passed, zero failures
or ignored tests. Regression tests check unspecified azo preservation, declared
E/Z preservation and supplied-coordinate preservation. Release Clippy with all
targets and warnings denied passes. Five benchmark-gate tests pass, including
rejection of invalid timing/accounting, missing rows/properties, nonfinite or
malformed coordinates and parse-only stereo differences. Changed Rust files
pass rustfmt. Raw outputs, parse snapshots, mismatch details and logs are kept
as lossless gzip artifacts listed in the raw report.

- [Final raw timing blocks and artifact hashes](2026-10-14-openbabel-sdf-speed-paired21.json)
- [Source, compiler, library and corpus manifest](2026-10-14-openbabel-sdf-speed-source.json)
- [Reconstructable source patch](https://github.com/kent-tokyo/chematic/blob/0bc130f24a7d163535e0f4ad59eb5a155847cf07/benchmarks/2026-10-14-openbabel-sdf-speed-candidate.patch)
- [Rust native harness](../crates/chematic-mol/examples/openbabel_sdf_throughput.rs)
- [Open Babel native harness](../tools/openbabel_sdf_throughput.cpp)
- [Runner](../scripts/bench_openbabel_sdf_throughput.py) and [corpus generator](../scripts/generators/build_openbabel_sdf_speed_corpus.py)

## Reproduce

Use the base revision plus the candidate patch in an isolated checkout. Preserve
a baseline executable built from the base with only the same example harness;
then build the candidate with a fresh target directory or source-mtime-safe
rebuild. Pin the recorded toolchain and Open Babel library. Decompress the exact
corpus, or regenerate it using the recorded RDKit version and verify its hash.

```sh
cargo build --release --locked -p chematic-mol --example openbabel_sdf_throughput
c++ -std=c++17 -O3 $(pkg-config --cflags openbabel-3) \
  tools/openbabel_sdf_throughput.cpp $(pkg-config --libs openbabel-3) \
  -o /tmp/openbabel-sdf
python scripts/bench_openbabel_sdf_throughput.py \
  --baseline /path/to/baseline-sdf --candidate /path/to/candidate-sdf \
  --openbabel /tmp/openbabel-sdf --corpus /path/to/corpus-1000.sdf \
  --source-manifest benchmarks/2026-10-14-openbabel-sdf-speed-source.json \
  --output /path/to/paired21.json --blocks 21 --repeats 5
```

Darwin RSS collection needs the normal read-only sysctl access used by time -l.
Keep source files identical to the manifest; the runner checks their hashes
before and after measurement. It also checks binary and corpus hashes, exact
record accounting and positive native timings, retaining incomplete status on
failure. The final report's retained-artifact list is lossless archival
postprocessing; timing rows are unchanged.

## Remaining boundaries

Full Open Babel superiority is not established. Internal parse stereo parity,
other Tier-A formats on diverse corpora, larger and fresh holdouts, malformed
input, streaming memory, cold startup, Python/CLI overhead, 3D workflows,
additional hosts and exact published-package comparisons remain open. Existing
tiny-fixture format results retain their own scope. The current result supports
only prepared 2D V2000 SDF rewriting and native read-plus-rewrite throughput on
this fixed corpus and host.
