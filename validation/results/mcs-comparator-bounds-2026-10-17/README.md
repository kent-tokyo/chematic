# MCS comparator enums and admissible search bounds

Final base: `abcc3569e7f49d15e911945509a5083ffff8a63d` (integration PR #1016).
Initially reproduced on `a07a1294c550d5fa55e5fdac66c5755556d76ff2`; the MCS
engine/binding source is unchanged between those bases. Tests and the native
extension were rebuilt after rebasing onto the final base.
[RDKit #9601](https://github.com/rdkit/rdkit/issues/9601) remains OPEN and tracks
typing enum-or-subclass MCS comparators. This candidate implements the native
enum choices, alongside legacy literal strings, and fixes comparator-dependent
MCS size underestimation. It does not implement RDKit callback/subclass adapters.

## Reproduction and results

The previous native fixture returned one MCS atom for `CCC`/`NNN` in both `any`
and `any_heavy_atom` modes, versus RDKit's three; `CCCO`/`NCCN` returned three
versus four. `baseline.json` records nine conditions (six mismatches). The MCS
engine and Python binding source were byte-identical between that fixture's
`cc72a803` source and this candidate's base. Both old seed/frontier bounds used
element buckets even when the comparator ignored element identity.

The fix chooses bound buckets according to the active comparator. Element mode
retains its old bound; any mode combines all atoms; native heavy-only mode combines
non-H atoms. Additional compatibility constraints can only lower achievable size,
so ignoring them in the bound remains conservative. Helper extraction keeps the
existing production-file refactor allowance from growing.

- Five new Rust regressions: four failed before the fix, all passed afterward.
- Full `chematic-smarts` lib: **222 passed**, no ignored tests or failures.
- Full `chematic-chem` lib: **962 passed, 1 existing ignored**, no failures.
- Python: **264 passed**, including **34 new** comparator tests, 68 module-function
  tests and 162 RDKit compatibility tests. Uses a fresh locally built candidate
  native-InChI extension, not a published wheel.
- Strict mypy accepts the enum module/public consumer; all **six** negative
  consumer calls fail with expected argument-type diagnostics.
- Clippy, format, refactor-budget, stub, unsafe-surface, cross-binding manifest,
  binding-surface and product-spelling gates; strict MkDocs; one executed guide
  block and three new rendered public API anchors.

## RDKit comparison and scope

RDKit **2026.09.1**, deterministic seed 9601, 35 short acyclic molecular inputs,
165 pair/triple groups including fixed regressions, both input orders, three atom
and two bond modes: **1,980 atom/bond-count comparisons**. No observed timeouts.
The raw RDKit comparisons match **1,976/1,980**; all four residuals are retained
in `oracle.json` and the complete rows are in `oracle.jsonl`.

The residuals are the existing explicit-H semantic difference: chematic's native
heavy mode excludes H altogether; RDKit's similarly named `CompareAnyHeavyAtom`
allows H-to-H matching. Native `ANY` still includes H. An additional independent
RDKit reference on the heavy induced graphs verifies chematic's documented mode
definition: **1,980/1,980** match this native profile. The script reports both the
raw comparison and the profile comparison; the latter is not full RDKit parity.

Only counts are compared. The existing Python result is a first-input-side
concrete projection, not generalized SMARTS encoding all accepted alternatives.
No generalized-query, atom-mapping, canonical-winner, performance or general FMCS
equivalence claim is made. Existing charge/ring/isotope/chirality option behavior
and native limits are unchanged. Python's checked boolean reports timeouts, not
Rust's separate resource-limited outcome. Dynamic runtime string compatibility
also permits a wrong enum family if its string value is valid; static typing
rejects that family mismatch.

## Commands

Set `CANDIDATE_PY` to a directory containing this candidate's Python package and
fresh native extension, and `RDKIT_PY` to the RDKit 2026.09.1 package directory.
The run used Python 3.13.6 / pytest 9.1.1; local builds use macOS arm64 linker
flags. The source/native hashes and build command are in `provenance.json`.

The host's default Apple linker produced a library rejected by dyld with a
misaligned LINKEDIT string-pool error (`default-loader.log`). Final validation
uses Rust's bundled LLD 22.1.6 and the installed macOS 26.5 SDK; the library is
copied directly from that link without binary patching/stripping. The final log
also records older deployment-target warnings for cached native-InChI C objects;
execution is verified on this macOS 27 host only, with no older-OS compatibility
or packaging claim. Repository linker/SDK configuration is unchanged.

```sh
MCS_MACHO_LD="$(rustc --print sysroot)/lib/rustlib/aarch64-apple-darwin/bin/gcc-ld/ld64.lld"
SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX26.5.sdk \
  PYO3_PYTHON=/opt/homebrew/bin/python3.13 CARGO_PROFILE_DEV_DEBUG=0 \
  CARGO_INCREMENTAL=0 cargo rustc -p chematic-py --features native-inchi -- \
  -C link-arg=-undefined -C link-arg=dynamic_lookup \
  -C "link-arg=-fuse-ld=$MCS_MACHO_LD"
cargo test -p chematic-smarts --lib
cargo test -p chematic-chem --lib
cargo clippy -p chematic-smarts -- -D warnings
PYTHONPATH="$CANDIDATE_PY" python3.13 -m pytest \
  crates/chematic-py/tests/test_mcs_comparators.py \
  crates/chematic-py/tests/test_module_functions.py \
  crates/chematic-py/tests/test_rdkit_compat.py -q -p no:cacheprovider
PYTHONPATH="$CANDIDATE_PY:$RDKIT_PY" python3.13 \
  validation/results/mcs-comparator-bounds-2026-10-17/oracle.py
MYPYPATH="$CANDIDATE_PY" python3.13 -m mypy --no-incremental \
  --config-file validation/results/mcs-comparator-bounds-2026-10-17/mypy.ini \
  crates/chematic-py/python/chematic/_mcs_comparators.py \
  validation/results/mcs-comparator-bounds-2026-10-17/typing_valid.py
# typing_invalid.py: expect exit 1 and six errors.
python3.13 -m mkdocs build --strict --site-dir /tmp/chematic-mcs-comparator-site
```

The worktree omits historical validation output and benchmarks via sparse
checkout to avoid exhausting local storage; the binding inventory required by
the gate was restored byte-for-byte from the base. Existing datasets, evidence
and other working trees were not removed.
