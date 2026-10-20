# Native conjugation assignment — 2026-10-19

Related tracking issue: [RDKit #9475](https://github.com/rdkit/rdkit/issues/9475),
verified OPEN on 2026-10-19. This implements a local invariant-state cache and
a single/dative-only early return in chematic's native RDKit-compatible
sanitizer. Public APIs and the existing conjugation rules remain unchanged.

## Correctness and build boundary

- Fresh Rust source build on macOS arm64, Rust 1.97.0; no published Python/native
  wheel is used as the implementation fixture.
- `cargo test -p chematic-smiles --offline --release --lib --tests`: **573 passed,
  0 failed, 2 ignored**, across 28 targets. Full output is in `tests.log`.
  Four new functional tests cover the frozen-loop differential, predicate
  independence, cache-size/element controls and the single-bond fast path.
- One ignored test is the opt-in timing added here; it was separately run
  successfully. The other ignored test is the pre-existing pipeline profiler.
- The fixed 1,131-input fixture prepares 1,125 internal models; the six
  parser/sanitizer refusals and their reasons are retained in `differential.log`.
  They are not counted as chemical matches. Refusals arise before the changed
  assignment stage.
- The test retains the pre-change loop independently. For all-false and
  all-true initial conjugated flags, all 2,250 fixture comparisons agree on
  the complete internal model (atoms, bonds, adjacency and ring state). All
  2,250 immediate repeat calls are idempotent.
- Further controls include empty/64/65/256-atom models, large polyenes,
  radicals, charges, H removal, silicon, hypervalent sulfur, coordination and
  aromatic flags independent of bond types.
- Release Clippy with `--all-targets --no-deps -- -D warnings` passed; output
  is retained in `clippy.log`. Rustfmt, Black on the measurement script,
  `git diff --check` and the MkDocs site build passed. Recorded source hashes
  match the final measured implementation and frozen-loop test source.

The cache reads only fields unchanged by conjugation assignment. Single/dative
bonds contribute at most one valence unit per endpoint, while the old loop
requires at least 1.5 for its initiating bond; its aromatic-flag reset is still
performed before the early return. No persistent molecule cache is introduced.

## Timings

`timing.log` and `timing.json` contain the raw 27 paired rounds, compiler/profile
settings, source/fixture/lockfile hashes and base checkout identity. Each group
has nine rounds, alternating which implementation runs first. Each timed lane
performs 300 passes; preparation, cloning and Debug comparisons are untimed.
The frozen reference and candidate run in the same release test binary.

| Group | Models | Reference median (s) | Candidate median (s) | Time reduction |
|---|---:|---:|---:|---:|
| Full prepared fixture | 1,125 | 0.052714542 | 0.023106583 | 56.1666% |
| Single/dative types only | 445 | 0.005055875 | 0.000691625 | 86.3204% |
| Other bond types present | 680 | 0.047212500 | 0.022639292 | 52.0481% |

These figures measure only `set_conjugation` over prepared models. They are
not end-to-end parser/sanitizer performance, application throughput, current
published release measurements or evidence of superiority over RDKit. The
upstream Enamine REAL packet is not used. The corpus includes small compatibility
and boundary examples and is not a representative drug-screening benchmark.
LTO is disabled and codegen units set to 16; other profiles/hardware can differ.

## Reproduce

```sh
python3 scripts/benchmarks/measure_conjugation_cache.py \
  --target-dir /tmp/chematic-conjugation-target \
  --output /tmp/chematic-conjugation-evidence
```

The script uses an offline release build with two jobs, no debug information,
LTO disabled, 16 codegen units and incremental compilation disabled. It first
runs functional checks, then the ignored timing test by itself. The target
directory is regenerable build output. The committed raw logs concern the
candidate working-tree source identified by hashes, based on main
`444ee7ab4e68d4b6a10aa30fed50eb829e8d4114`; they are not a release artifact.
