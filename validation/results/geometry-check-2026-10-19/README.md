# Supplied-geometry screen — 2026-10-19

Source base: `b7004777`, Open Babel issue
[#2651](https://github.com/openbabel/openbabel/issues/2651), confirmed open.
Candidate source hashes and the freshly built Python extension hash are in
`environment.json`. The reported native version is the source manifest's 1.1.9,
not evidence of a published artifact.

## Results

- Rust `chematic-3d --lib --tests`: **764 passed, 0 failed, 4 ignored**, across
  12 test binaries, including all 13 new `geometry_check` tests.
- Python geometry screen + existing checked-SASA suite: **78 passed**
  (35 new geometry-screen cases, 43 existing SASA cases).
- Release `cargo clippy` for `chematic-3d` and `chematic-py`, all targets,
  `-D warnings`: passed.
- Whole-workspace format check, Python Black on both new Python files, source
  stub gate, refactor budget, product names, binding inventory and the existing
  Open Babel file-I/O contract: passed.
- MkDocs strict build with the repository-pinned top-level documentation
  package versions (1.6.1 / 9.7.7 / 1.0.6 / 12.0.1): passed.
- Independent script: **360 indexed-distance reports** match Python `math.dist`
  and RDKit 2026.09.1's independently parsed indexed atom/bond graphs. All
  nonbonded pair sets and distances match over ten original molecular controls,
  twelve coordinate packets each and three cutoffs. Returned bond-outlier
  distances are also checked independently.
- **192 nondegenerate tetrahedral reports** match RDKit 2026.09.1's coordinate
  stereoperception: two explicit-ligand and two implicit-H declarations, each
  mirrored and transformed by 24 proper signed-axis rotations plus translation.
  Planar/unknown and E/Z controls are separately exercised in the unit tests.

The independent script does not validate the physical quality of the random
coordinates, fit the reference bond lengths, or establish conformer quality.
No universal validity, energy, performance or Open Babel-wide parity claim.

## Reproduction

All new molecular controls are small original authored examples. No external
corpus or coordinates are redistributed; no held-out molecules are imported.
The random-distance packet uses Python's fixed seed 2651.

Rust commands used the same environment:

```sh
export CARGO_TARGET_DIR=/private/tmp/chematic-conjugation-target
export CARGO_BUILD_JOBS=2 CARGO_INCREMENTAL=0
export CARGO_PROFILE_RELEASE_DEBUG=0 CARGO_PROFILE_RELEASE_LTO=false
export CARGO_PROFILE_RELEASE_CODEGEN_UNITS=16
export PYO3_PYTHON=python3.13
cargo test -p chematic-3d --offline --release --lib --tests
cargo clippy -p chematic-3d -p chematic-py --release --offline --all-targets -- -D warnings
cargo rustc -p chematic-py --release --offline -- -C strip=none -C link-arg=-undefined -C link-arg=dynamic_lookup
```

The last command builds the actual candidate extension. The final-crate
`strip=none` override was necessary for this host's macOS 27 LINKEDIT loader
error; the unstripped build imports successfully. The host's removed
`-ld_classic` option and stripping the resulting file did not fix it. This is
consistent with [Rust #157750](https://github.com/rust-lang/rust/issues/157750).
No repository-wide linker or release settings were changed.

Stage **all current Python sources**, excluding pre-existing `.so` files and
`__pycache__`, into a temporary `chematic` package, then link the fresh cdylib as
`chematic/chematic.cpython-313-darwin.so`. The full package initializer is used;
there is no published-wheel bootstrap or mocked native implementation.

```sh
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 PYTHONPATH=/private/tmp/chematic-geometry-check-python python3.13 -m pytest -p no:cacheprovider -q crates/chematic-py/tests/test_geometry_check.py crates/chematic-py/tests/test_sasa_checked.py
PYTHONPATH=/private/tmp/chematic-geometry-check-python:/private/tmp/chematic-rdkit-2026091/lib/python3.13/site-packages python3.13 scripts/validation/validate_geometry_check.py
python3.13 -m mkdocs build --strict --site-dir /private/tmp/chematic-geometry-check-docs
```

`rust-tests.log`, `python-tests.log`, `clippy.log` and
`independent-checks.json` retain the successful native/test outputs. These are
local source checks; hosted CI and main integration are separate gates.
