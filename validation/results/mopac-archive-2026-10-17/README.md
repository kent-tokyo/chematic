# MOPAC archive reader validation — 2026-10-17

The candidate adds a field-preserving reader for one archive restart deck.
It does not convert internal coordinates, infer bonds, evaluate ancillary
sections, or establish calculation convergence. Source/artifact hashes and
compiler versions are in `manifest.json`; upstream issue state is recorded
separately. No upstream implementation code or attachment contents are copied
into the candidate.

- Rust: 805 existing library tests and 7 new archive integration tests pass.
  Existing Cartesian integration results are in `rust-cartesian-tests.log`.
- Python: 277 selected archive/Cartesian/format/cross-binding tests pass,
  including 9 new archive tests, using the freshly built extension.
- Deterministic corpus: 512 original Cartesian/internal archives preserve
  parsed values, references, flags, announced charges and raw sections.
  Sixteen extracted Cartesian decks cross-read with Open Babel 3.2.1 and
  agree in coordinates within 1e-4 Å. This checks coordinate interpretation,
  not energy, topology, symmetry execution or internal-coordinate conversion.
- The separately downloaded issue attachment parses as four atoms, with the
  expected final reference tuple and internal-coordinate values, retaining
  its ancillary SYMMETRY section. Open Babel 3.2.1 `-imopout -oxyz` emits
  no molecule for that archive. Archive/zip hashes are in `corpus.json`.
- Strict Clippy for both changed crates/all targets, formatting, refactor
  budget, stubs, binding inventory, cross-binding manifest, unsafe surface,
  naming, workflow pins and strict MkDocs pass.

Commands (from the candidate worktree):

```sh
# Local low-disk build settings, used for Rust tests and the extension:
export CARGO_INCREMENTAL=0 CARGO_BUILD_JOBS=1
export CARGO_PROFILE_DEV_DEBUG=0 CARGO_PROFILE_TEST_DEBUG=0
export CARGO_PROFILE_DEV_STRIP=none
export CARGO_PROFILE_DEV_OPT_LEVEL=1 CARGO_PROFILE_TEST_OPT_LEVEL=1
export CARGO_PROFILE_DEV_CODEGEN_UNITS=1 CARGO_PROFILE_TEST_CODEGEN_UNITS=1
export CARGO_TARGET_DIR=/private/tmp/chematic-sparse-tanimoto-target
export PYO3_PYTHON=/Library/Frameworks/Python.framework/Versions/3.13/bin/python3
cargo test -p chematic-mol --offline --lib --test mopac_archive_contract
cargo test -p chematic-mol --offline --test mopac_cartesian_contract
cargo clippy -p chematic-mol -p chematic-py --all-targets --offline --no-deps -- -D warnings
cargo rustc -p chematic-py --offline -- -C strip=none -C link-arg=-undefined -C link-arg=dynamic_lookup
# Copy the dylib as chematic.so beside candidate Python package source files.
export PYTHONPATH=/private/tmp/chematic-mopac-archive-python
export PYTHONDONTWRITEBYTECODE=1 PYTEST_DISABLE_PLUGIN_AUTOLOAD=1
python3.13 -m pytest crates/chematic-py/tests/test_mopac_archive.py crates/chematic-py/tests/test_mopac_cartesian.py crates/chematic-py/tests/test_new_formats.py crates/chematic-py/tests/test_format_convert.py crates/chematic-py/tests/test_cross_binding_contract.py -q -p no:cacheprovider
python3.13 validation/validate_mopac_archive.py --upstream-zip /private/tmp/chematic-mopac-example.zip
```

The corpus command also works without the optional upstream zip argument.
The baseline source has no archive reader; the existing Cartesian reader is
retained. Local results are distinct from hosted CI and a published artifact.
