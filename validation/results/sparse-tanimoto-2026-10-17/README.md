# Sparse Tanimoto candidate — 2026-10-17

Base: `de41dadd8f123a4b1368937c7dd86df823434588` (main).
Native API: `sparse_tanimoto_matrix_checked`. The base exposes dense fingerprint
matrices, but no supplied sparse on-bit condensed matrix API. No baseline native
binary was retained or rebuilt; this is a source-level capability addition.

The upstream snapshot confirms RDKit #9601 is open and lists typed explicit /
sparse bit-vector inputs as remaining wrapper work. RDKit already calculates
sparse matrices. This candidate covers the native bounded-input slice only.

## Numeric evidence

`summary.json` records 240 deterministic datasets, 2,880 native calls and 960
live RDKit oracle calls. Five native input forms (lists, nested generators,
NumPy uint32 arrays, ExplicitBitVect and SparseBitVect sources), both metrics
and reversed row orders were checked. Ordinary nonempty-union comparisons had
zero failures and maximum absolute error zero. Counts exclude the separate high
ID and empty/singleton diagnostics. `corpus-results.json.gz` retains inputs and
all ordinary outputs; this is an original set corpus, not chemical fingerprint
accuracy evidence. No timing or speed claim is made.

The separate diagnostics disclose:

- Native empty/empty similarity is 1, consistent with existing BitVec2048;
  RDKit 2026.09.1 returns 0 on this host.
- Native singleton input returns an empty list; RDKit refuses it with an invariant
  exception. The printed upstream exception in `probe.log` is expected/caught.
- Sparse direct-source IDs near 2**31 match RDKit without dense expansion.
  This RDKit build exposes bit 2**32-2 as signed -2 through GetOnBits(); native
  unsigned extraction refuses that view. Full u32 IDs supplied directly as
  unsigned integer rows match the independent set reference.

The candidate native binary hash, Python, NumPy and RDKit versions are recorded
in `summary.json`. The binary is local and is not a wheel/publication artifact.
Source hashes are recorded in `source-hashes.json`.

## Local checks

Rust: 346 existing library tests and 7 new integration tests passed. Python:
40 new API tests and 25 existing screening/registration tests passed. Strict
Clippy on both changed crates (all targets, no deps) passed with warnings denied.
Source/documentation gates passed; logs are retained alongside this record.

## Reproduction

Use Python 3.13 with pytest, NumPy and RDKit 2026.09.1. Build the extension with
Cargo using `CARGO_INCREMENTAL=0`, `CARGO_BUILD_JOBS=1`, dev/test debug=0,
opt-level=1, codegen-units=1, dev strip=none, and a writable CARGO_TARGET_DIR.
On this macOS host the command was:

```sh
cargo rustc -p chematic-py --offline -- \
  -C strip=none -C link-arg=-undefined -C link-arg=dynamic_lookup
```

Copy the resulting dylib as `chematic/chematic.so` beside the Python package
files in a temporary directory. Set PYTHONPATH to that directory and the isolated
RDKit package target. No existing installed chematic package was overwritten.

```sh
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.13 -m pytest \
  crates/chematic-py/tests/test_sparse_tanimoto.py \
  crates/chematic-py/tests/test_screening_enrichment.py -q -p no:cacheprovider
python3.13 validation/rdkit_issues/sparse-tanimoto/probe.py <output-directory>
cargo test -p chematic-fp --offline --test sparse_matrix --lib
cargo clippy -p chematic-fp -p chematic-py --all-targets --offline --no-deps -- -D warnings
```

External pytest plugins were disabled because the host's old pytest-asyncio
plugin fails collection with pytest 9. Native Python tests do not require it.
Rust tests exercise exact raw/work limits, next-item refusal, u32::MAX, ordering,
duplicate idempotence, source preservation and existing dense native convention.
Local source gates cover formatting, refactor budget, stubs, binding inventory,
cross-binding manifest, unsafe surface, canonical product spelling and strict
MkDocs generation. Hosted CI and main integration are separate from these checks.
