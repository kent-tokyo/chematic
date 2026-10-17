# Checked explicit core removal

Base: `abcc3569e7f49d15e911945509a5083ffff8a63d`.
RDKit [#9601](https://github.com/rdkit/rdkit/issues/9601) remains OPEN and lists
typed `ReplaceCore` match/matches inputs. No `ReplaceCore` / `replace_core` API
was found in the base's crates/docs. This candidate adds a native checked subset
with an explicit core-index-to-source-index sequence and index-labelled dummies.
It does not fix or replace the whole upstream wrapper checklist.

## Results

- New native module reuses checked bridge-single-bond fragmentation, then removes
  mapped core atoms and their own-side dummies with order-preserving atom deletion.
  The existing fragment graph validator becomes crate-visible; its behavior stays
  unchanged. Return values include source-to-output atom indices and original /
  output attachment provenance.
- Four new native tests; full chem lib **966 passed, 1 existing ignored**, no
  failures. Covers retained atom fields/H counts/caller tags, shifted numbering,
  ring-core and whole-source removal, invalid maps and unsupported edits.
- Fresh local native-InChI extension: Python **259 passed** (42 new core tests,
  55 fragment regressions, 162 existing RDKit compatibility tests).
- Strict mypy accepts the public helper and positive consumer; all **six**
  negative calls/assignments fail as expected. Runtime additionally excludes bool
  indices, which the Python type system otherwise treats as integers.
- Clippy for both changed Rust crates; format, stub, cross-binding, binding,
  unsafe, refactor-budget and product spelling gates; strict MkDocs; two executed
  guide blocks and three rendered API anchors.

## Independent reference

RDKit **2026.09.1**, 32 source graphs and 14 core graphs. An independent bounded
RDKit match enumeration selects concrete atom-compatible maps with single,
non-ring boundary bonds. For **415 accepted mappings**, canonical SMILES graph
identity and labelled attachment/source-map metadata match
`Chem.ReplaceCore(..., labelByIndex=True)`; four outputs are empty. For the
411 nonempty outputs, independent V3000 reading also checks indexed retained
atom numbers/isotopes/charges/maps/H counts and each dummy's label/anchor.
Full conditions are in `oracle.jsonl`, summary in `oracle.json`.

Another **12** RDKit-perceived aromatic-core mappings into a native Kekulé aspirin
graph are verified as explicit refusals, not counted as successful graph parity.
The matching contract compares native aromatic flags/order and does not perceive
or normalize aromaticity. The lowercase aromatic aspirin representation is also
present in the corpus. No native refusal inside the declared comparison profile
is silently discarded by the oracle.

## Boundaries

- Source <=4096 atoms / 8192 bonds; core 1..1024 atoms / <=8192 bonds; <=1024 cuts.
- Both inputs must be concrete and stereo-free. Only non-ring plain single
  boundary bonds; complete cyclic cores can still be removed.
- Explicit mapping validates element, native aromaticity, charge/isotope and core
  bonds/orders, ignoring H counts and atom maps during matching. Retained source
  values survive. Extra source edges among mapped atoms are allowed.
- No query cores, automatic match search, stereo transfer, alternate dummy bond /
  label policies, replaceDummies/requireDummyMatch or batch API.
- Existing V2000/MolJSON isotope-dummy export residuals remain; use SMILES. Only
  the named V3000 indexed fields are validated, not dummy query metadata or 3D.
- New graph and independent Python property store; all source atoms not in the
  core remain in order, including unrelated components. No performance claim,
  registry/wheel claim or broad ReplaceCore parity claim.

## Reproduction

Use Python 3.13.6 / pytest 9.1.1, a candidate package containing this source and
fresh native extension (`CANDIDATE_PY`), and RDKit 2026.09.1 (`RDKIT_PY`). Exact
source/native SHA-256 values, platform and build command are in `provenance.json`.
The local macOS 27 run uses installed SDK 26.5 and Rust bundled LLD 22.1.6 because
the default Apple linker previously emitted a dyld-rejected LINKEDIT layout.
The final build log retains cached native-InChI deployment-target warnings;
execution is checked only on this host, not older macOS. No repository linker
settings were changed and no binary patching/stripping was applied.

```sh
CORE_MACHO_LD="$(rustc --print sysroot)/lib/rustlib/aarch64-apple-darwin/bin/gcc-ld/ld64.lld"
SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX26.5.sdk \
  PYO3_PYTHON=/opt/homebrew/bin/python3.13 CARGO_PROFILE_DEV_DEBUG=0 \
  CARGO_INCREMENTAL=0 cargo rustc -p chematic-py --features native-inchi -- \
  -C link-arg=-undefined -C link-arg=dynamic_lookup \
  -C "link-arg=-fuse-ld=$CORE_MACHO_LD"
cargo test -p chematic-chem --lib
cargo clippy -p chematic-chem -p chematic-py -- -D warnings
PYTHONPATH="$CANDIDATE_PY" python3.13 -m pytest \
  crates/chematic-py/tests/test_core_replacement.py \
  crates/chematic-py/tests/test_bond_fragment.py \
  crates/chematic-py/tests/test_rdkit_compat.py -q -p no:cacheprovider
PYTHONPATH="$CANDIDATE_PY:$RDKIT_PY" python3.13 \
  validation/results/checked-core-removal-2026-10-17/oracle.py
MYPYPATH="$CANDIDATE_PY" python3.13 -m mypy --no-incremental \
  --config-file validation/results/checked-core-removal-2026-10-17/mypy.ini \
  crates/chematic-py/python/chematic/_core_replacement.py \
  validation/results/checked-core-removal-2026-10-17/typing_valid.py
# typing_invalid.py: expect exit 1 and six errors.
python3.13 -m mkdocs build --strict --site-dir /tmp/chematic-core-site
```

The worktree omits historical results and benchmarks via sparse checkout. Original
working-tree document changes and all other working trees/evidence are preserved.
The directory date follows the UTC session date; local testing crossed into
2026-10-18 JST.
