# V3000 native block throughput and semantic boundaries

The reader improves 1.055x and 1.087x over the corrected baseline on the
two diverse slices, with 21/21 faster blocks. In the separate controlled
chain-size profile it improves 2.116x; Open Babel ratios there are 86.675x
for parsing, 2.522x for writing and 49.318x for roundtrip. Those controlled
results do not establish diverse-corpus or toolkit-wide Open Babel superiority.

This packet extends the Open Babel comparison from V2000 SDF to native V3000
MOL block APIs. The optimized reader preserves the corrected baseline's exact
internal snapshots. Both source builds match every declared read/rewrite field
on two fixed, unfiltered 1,000-record slices. Open Babel has internal graph and
rewrite failures on those slices, so they cannot establish an equivalent-output
Open Babel speed win. A separate controlled chain-size profile exercises the
dense-index fast path with equal outputs in all engines.

## Implementation

The reader borrows ordinary V30 logical lines and key/value values. Only a
continued line needs an owned joined string. Source-id maps track whether every
id appended is exactly `1..=len`; only with that proof do they use direct
indexing. Sparse, permuted, zero and duplicate ids keep the old first-match
lookup. Duplicate bond ids and missing endpoints still fail. Existing physical
line checks, tokenization, limits, stereo and geometry diagnostics remain.

Two correctness fixes precede the performance comparison in both source arms:

- Supported non-aromatic `VAL` overrides are applied after bonds are known,
  retaining the source attribute. This prevents an added implicit H on the
  sulfur/iodine radical cases at source rows 5374, 5572 and 5645. `VAL=0` keeps
  inference, while the zero sentinels suppress H. Aromatic/query and unsupported
  valence representations retain their previous boundary. See
  [RDKit's parser and property processing](https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/FileParsers/MolFileParser.cpp).
- Explicit double-bond `CFG=2` remains in normalized bond metadata and survives
  writing, without duplication when depiction regenerates it. Original row 766
  otherwise becomes definite E solely from its coordinates. The generic writer
  still permits geometry-defined E/Z when no explicit unknown marker exists.
  The full regression suite checks that established behavior.

The measured baseline already has both fixes, with the old owned logical lines
and linear id maps. Only borrowing and lookup optimizations distinguish the
candidate, so a chemical behavior change is not presented as a speed gain.

## Corpora and protocol

The diverse profiles are the same exposed SMILES slices at offsets 0 and 5000,
1,000 rows each, SHA256
`f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`.
RDKit 2026.09.1 uses the same 2D coordinate settings as the earlier SDF packet,
but `forceV3000=True`. Every source row is retained. Names carry absolute source
indices; SD properties are absent because these are individual MOL API calls.
SGROUP/COLLECTION/query semantics are not established by these ordinary blocks.

The separate scalability profile contains all 64 linear alkanes of lengths
16, 32, …, 1024: 33,280 atoms total, including blocks too large for V2000.
Its [explicit SMILES list](../validation/benchmark_corpora/openbabel-v3000-linear-chains-v1.smi)
is retained. This controlled profile was declared after observing the diverse
Open Babel semantic failures, before its formal timing. It measures contiguous
id scaling on saturated chains; it is not diverse chemistry, a sealed holdout,
or a replacement for the original failed gates. No source row or timing block
is filtered from any of the three profiles.

Both native harnesses load and split the framed block collection before timing.
The `$$$$` separators frame the API batch, rather than measuring an SDF supplier.
Each parser receives borrowed individual MOL blocks and retains the entire
parsed batch. Rust runs the normal diagnostic V3000 reader. Both writers return
new native MOL strings; prepared writing has one untimed warm write. Roundtrip
parses a new batch and writes it inside the timer. File I/O, framing, process
startup and initial plugin discovery are excluded symmetrically. These numbers
are not interchangeable with the earlier SDF supplier timings.

Twenty-one sequential fresh-process blocks rotate baseline, candidate and
Open Babel evenly. Each process times five full-batch repetitions. Ratios are
paired comparator/candidate medians with deterministic 10,000-resample paired
bootstrap 95% intervals. A strict win needs matching outputs, a lower bound
above 1, and 21/21 faster blocks. Whole-process high-water RSS includes input,
framing, prepared records, warmup, plugins and timed work; it is not streaming
memory or operation-only allocation evidence.

## Final measurements

Times are medians in milliseconds per full batch (1,000 diverse records or 64 chains). Ratios are paired medians, not ratios of independent median times.

| Profile | Lane | Candidate ms | Open Babel ms | Open Babel / candidate [95% CI] | Faster blocks |
|---|---|---:|---:|---|---:|
| 0 | parse | 34.901 | 113.031 | Not ranked: semantic gate fails | — |
| 0 | write | 7.522 | 17.519 | Not ranked: semantic gate fails | — |
| 0 | roundtrip | 42.778 | 130.528 | Not ranked: semantic gate fails | — |
| 5000 | parse | 10.832 | 53.041 | Not ranked: semantic gate fails | — |
| 5000 | write | 3.448 | 10.655 | Not ranked: semantic gate fails | — |
| 5000 | roundtrip | 14.311 | 63.887 | Not ranked: semantic gate fails | — |
| chains | parse | 9.442 | 819.566 | 86.675x [86.453, 87.197] | 21/21 |
| chains | write | 7.527 | 18.983 | 2.522x [2.509, 2.538] | 21/21 |
| chains | roundtrip | 16.976 | 841.159 | 49.318x [49.161, 49.787] | 21/21 |

All three controlled-chain lanes pass the strict Open Babel gate. The six diverse-profile Open Babel lanes remain unranked; their raw times do not imply an equivalent-output win.

### Incremental optimization

| Profile | Lane | Corrected baseline / candidate [95% CI] | Faster blocks | Strict win |
|---|---|---|---:|---|
| 0 | parse | 1.055x [1.051, 1.058] | 21/21 | Yes |
| 0 | write | 1.007x [0.998, 1.008] | 14/21 | No |
| 0 | roundtrip | 1.037x [1.032, 1.043] | 21/21 | Yes |
| 5000 | parse | 1.087x [1.082, 1.091] | 21/21 | Yes |
| 5000 | write | 0.991x [0.988, 0.994] | 5/21 | No |
| 5000 | roundtrip | 1.060x [1.056, 1.067] | 21/21 | Yes |
| chains | parse | 2.116x [2.102, 2.127] | 21/21 | Yes |
| chains | write | 0.995x [0.986, 1.007] | 9/21 | No |
| chains | roundtrip | 1.616x [1.606, 1.619] | 21/21 | Yes |

Reader and roundtrip improvements pass all 21 blocks on all three profiles. Prepared writing has no incremental strict win. Additional-slice writing is 0.9% slower [0.988, 0.994] in paired ratios despite identical writer source in both arms; this regression is retained. The complete candidate is not a universal incremental speed win.

### Whole-process memory

| Profile | Lane | Candidate median peak MiB | Open Babel median peak MiB |
|---|---|---:|---:|
| 0 | parse | 16.47 | 51.94 |
| 0 | write | 11.12 | 38.97 |
| 0 | roundtrip | 16.50 | 52.06 |
| 5000 | parse | 11.39 | 40.92 |
| 5000 | write | 7.88 | 32.92 |
| 5000 | roundtrip | 11.38 | 40.98 |
| chains | parse | 14.91 | 50.78 |
| chains | write | 10.33 | 39.38 |
| chains | roundtrip | 14.91 | 51.52 |


## Independent semantics

RDKit 2026.09.1 observes canonical isomeric graph, indexed atom/bond properties,
finite coordinates within `1e-4`, names, public properties and explicit unknown
double-bond markers. Additional safe `either` markers are accepted; explicit
input markers must survive. Internal Rust snapshots also retain exact coordinate
bits, metadata, neighbor orders and diagnostics, and are compared byte-for-byte.
Output agreement never substitutes for internal parser agreement.

The Open Babel 3.2.1 V3000 writer emits `CFG=6` in 153 original-slice records,
which the strict RDKit observer rejects. The accepted remainder contains one
graph/unknown-marker difference. On the additional slice, all outputs parse,
but 137 graph/unknown-marker differences remain. Internal canonical graphs
differ on 3 original and 137 additional rows. The complete row indices and raw
outputs are retained. A rejected output keeps its original row position and
fails the gate; it is never silently dropped or counted as agreement.

The `CFG=6` output also reproduces with a fresh ordinary Open Babel CLI call on
retained source row 73 using `obabel -imol input.mol -omol -x3 -O output.mol`.
This observation does not depend on the custom harness's warm write. No output
marker is repaired to obtain a comparison win, and no upstream message or issue
is posted as part of this packet.

## Provenance and validation

Base commit is `ffe682c32cb4c76de0d309649427202a7d8e956d`, source version 1.1.3.
Apply the recorded prior SDF writer, reader follow-up and runner-dependency
patches, then either new baseline or candidate patch in isolated checkouts.
The source/compiler/library manifest records each patch and file hash. Other
dirty Python, chemistry, 3D and file-format work is excluded.

Exploratory binaries from a reused build target were discarded after dependency
drift from the preceding Python build was detected. Every final source arm was
rebuilt in a dedicated initially empty target, with fresh source mtimes. The
final clean-build logs and binary hashes distinguish those artifacts. A
provisional row-745 difference from the discarded build is not a new parser
defect or a claimed fix. No exploratory timing becomes a formal observation.

Rust 1.97.0 release settings are opt-level 3, LTO and one codegen unit. The C++
harness uses `-std=c++17 -O3` and the installed Homebrew Open Babel 3.2.1 library.
Native library and MOL/canonical-SMILES plugin hashes are checked before and
after timing. The wrapper also checks the frozen build sources. All 1,098
tracked Rust/workspace files were checked against the pinned base: only the
previous MOL writer/stereo changes and these two V3000 source files differ.
This is one macOS arm64 host and native source builds, not a published crate,
wheel, hosted CI, release, commit or merge.

Final isolated `cargo test --release --locked -p chematic-mol`: 979 passed,
zero failed or ignored. Eight added tests cover id mapping, historical duplicate
atom handling, duplicate/missing references, continued properties, total-valence
overrides, radical fixtures and unknown azo stereo retention. Release Clippy
for all MOL targets passes with warnings denied. Ten benchmark-gate tests pass,
including rejected-row accounting, original row indices and independent parser
versus writer gates. Changed Rust files pass rustfmt and the diff has no
whitespace errors.

## Retained artifacts

- [Source/compiler/library/corpus manifest](2026-10-14-openbabel-v3000-speed-source.json)
- [Reconstruction, full mismatch rows, artifact hashes and validation](2026-10-14-openbabel-v3000-speed-verification.json)
- [0 raw 21-block record](2026-10-14-openbabel-v3000-speed-0-paired21.json)
- [5000 raw 21-block record](2026-10-14-openbabel-v3000-speed-5000-paired21.json)
- [chains raw 21-block record](2026-10-14-openbabel-v3000-speed-chains-paired21.json)
- [baseline reconstructable patch](2026-10-14-openbabel-v3000-speed-baseline.patch)
- [candidate reconstructable patch](2026-10-14-openbabel-v3000-speed-candidate.patch)
- [0 exact corpus](../validation/benchmark_corpora/openbabel-v3000-speed-slice0-1000-v1.molblocks.gz)
- [5000 exact corpus](../validation/benchmark_corpora/openbabel-v3000-speed-slice5000-1000-v1.molblocks.gz)
- [chains exact corpus](../validation/benchmark_corpora/openbabel-v3000-speed-chains64-v1.molblocks.gz)
- [Rust native harness](../crates/chematic-mol/examples/openbabel_v3000_throughput.rs), [Open Babel native harness](../tools/openbabel_v3000_throughput.cpp), [runner](../scripts/bench_openbabel_v3000_throughput.py), and [generator](../scripts/build_openbabel_v3000_speed_corpus.py)

Reproduce from the pinned base plus the three prerequisite patches and one alternative patch. Build the unchanged example with the recorded toolchain, compile the C++ harness using the installed `pkg-config openbabel-3` flags, and decompress the exact profile. Run the V3000 runner with `--blocks 21 --repeats 5`, `--records 1000` or `64`, the exact three executables and the source manifest. The retained wrapper adds frozen-build checks around the same runner. Raw snapshots/rewrites, full mismatch indices, CLI reproduction and build/test logs are losslessly archived. Archival fields do not change timing rows.

## Remaining work

General Open Babel speed superiority remains unproven. The controlled V3000
profile cannot replace the diverse semantic gates. Other production formats,
general 3D and query/SGROUP/COLLECTION inputs, malformed-input and termination
work, fresh corpora, cold startup, Python/CLI overhead, constant-memory streaming,
additional hosts and exact published packages need their own comparisons.
