# GAMESS explicit internal-coordinate templates — 2026-10-17

The bounded native/Python candidate emits `COORD=UNIQUE` Cartesian `$DATA`
with C1 and explicit `$ZMAT IZMAT` definitions. It implements the `$ZMAT`
conversion capability requested by Open Babel #2523. It does not emit
`COORD=ZMT` internal-coordinate `$DATA`, select a method/basis, or execute a
GAMESS job. The source-frame profile keeps atom indices tied to IZMAT.

Evidence:

- Rust: 805 existing library, 9 new GAMESS and 13 existing Z-matrix tests pass.
- Python: 133 selected tests pass, including 20 new GAMESS tests, using the
  freshly built extension.
- 640 original synthetic records (1/2/3/4/8 atoms, varying earlier reference
  chains, 8/17 significant digits) match the direct plain-XYZ adapter exactly.
  Seventeen-digit Cartesian coordinates round-trip exactly.
- 3,584 internal values measured using RDKit 2026.09.1 agree with the native
  checked converter; maximum difference 3.979039320256561e-13 in Å/degrees.
  IZMAT codes and atom tuples are separately decoded and checked against the
  caller's reference chains. These are geometry checks, not molecule validity.
- 80 GAMESS-input → XYZ cross-reads with Open Babel 3.2.1 preserve element
  order and Cartesian coordinates within 1e-4 Å. Open Babel does not validate
  the execution of the `$ZMAT` optimizer definitions.
- The same Open Babel version's default XYZ → GAMESS output has no `$ZMAT`;
  the original water probe and full output are in `openbabel-baseline.json`.
- Strict all-target Clippy for both changed crates, formatting, refactor,
  source stubs, binding inventory, cross-binding manifest, unsafe surface,
  product naming, workflow pins and strict MkDocs pass.

`water.xyz` / `water.inp` are an original small reviewable example, not a
third-party calculation result. Source/extension/compiler/reference-manual
hashes are in `manifest.json`. The native inventory adds exactly two Python
functions over the base main. No package publication is claimed.

Reproduction (from the candidate worktree):

```sh
export CARGO_INCREMENTAL=0 CARGO_BUILD_JOBS=1
export CARGO_PROFILE_DEV_DEBUG=0 CARGO_PROFILE_TEST_DEBUG=0
export CARGO_PROFILE_DEV_STRIP=none
export CARGO_PROFILE_DEV_OPT_LEVEL=1 CARGO_PROFILE_TEST_OPT_LEVEL=1
export CARGO_PROFILE_DEV_CODEGEN_UNITS=1 CARGO_PROFILE_TEST_CODEGEN_UNITS=1
export CARGO_TARGET_DIR=/private/tmp/chematic-sparse-tanimoto-target
export PYO3_PYTHON=/Library/Frameworks/Python.framework/Versions/3.13/bin/python3
cargo test -p chematic-mol --offline --lib --test gamess_zmat_contract --test gzmat_precision
cargo clippy -p chematic-mol -p chematic-py --all-targets --offline --no-deps -- -D warnings
cargo rustc -p chematic-py --offline -- -C strip=none -C link-arg=-undefined -C link-arg=dynamic_lookup
# Copy built dylib as chematic.so beside candidate Python package files.
export PYTHONPATH=/private/tmp/chematic-gamess-zmatrix-python:/private/tmp/chematic-sparse-tanimoto-rdkit
export PYTHONDONTWRITEBYTECODE=1 PYTEST_DISABLE_PLUGIN_AUTOLOAD=1
python3.13 -m pytest crates/chematic-py/tests/test_gamess_zmat.py crates/chematic-py/tests/test_gzmat_precision.py crates/chematic-py/tests/test_xyz_symbol_case.py crates/chematic-py/tests/test_new_formats.py -q -p no:cacheprovider
python3.13 validation/validate_gamess_zmat.py
```

The plain XYZ adapter refuses extra columns and key/value comments; it does
not silently discard extxyz cells or per-atom properties. Reference selection
remains the caller's responsibility. GAMESS runtime acceptance and optimizer
behavior remain unverified on this host.
