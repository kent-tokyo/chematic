# Caller-supplied MMFF partial charges, 2026-10-17

Base: `bef0d008d49f3d432fc35a900ccbf6acab2cd957`.
`source-sha256.json` records the modified implementation, tests, probe,
manifests, all MMFF implementation source and reused coordinate fixture.
`summary.json` records the compiled native hash, host and measured results.

## Scope and upstream observation

[Open Babel #2496](https://github.com/openbabel/openbabel/issues/2496) was
OPEN when inspected; `upstream-issue.json` is the dated API snapshot. The
upstream issue concerns GAFF charge-model and charge-file inputs. This change
adds a transferable construction-time input to chematic's prepared
MMFF94/MMFF94s fields, not GAFF or an external charge-file parser.

The authored six-atom methanol C++ probe creates a fresh Open Babel 3.2.1
molecule/field for each charge scale, preassigns atom charges, and marks them
perceived before setup. `openbabel-probe.json` retains the output:

- GAFF charge scales 0/1/2: electrostatics always `8.5165206621963172`, total
  energy always `71.301943158404939` in the backend's native energy units.
- MMFF94 charge scales 0/1/2: electrostatics always zero, total energy always
  `17.652556788138931` in the backend's native energy units.

This observes preassigned charges not controlling those prepared calculations;
it does not reproduce the issue's separate charge-file atom-mapping error.
No upstream source implementation or attached example was copied.

## Local validation

- Rust: **293 passed, zero failed/ignored**: 283 library unit tests, 5 existing
  dielectric integration tests and 5 new custom-charge integration tests.
- Python: **45 passed**: 19 new custom-charge cases and 26 existing dielectric
  cases, against the compiled native extension.
- Independent analytic matrix: **504 passed** across 7 geometries, 2 MMFF
  variants, 2 cutoffs, 2 interfragment policies, 3 charge arrays and 3 dielectric
  values. The probe independently computes graph-BFS exclusions and the
  buffered Coulomb energy/gradient. It subtracts a zero-charge field to isolate
  electrostatics while keeping all bonded/vdW terms present.
- Maximum absolute energy error: `7.105427357601002e-14` kcal/mol; maximum
  gradient-component error: `2.708944180085382e-14` kcal/mol/Å. Both thresholds
  are `1e-8`, and every input/measurement is in `analytic.json.gz`.
- Full-model finite differences: **27 axes passed**, step `1e-5` Å, maximum
  gradient difference `7.336671181690235e-8`, threshold `1e-5`.
- Strict all-target Clippy for both changed crates, formatting, source gates
  and strict MkDocs passed. Compressed logs retain their outputs.
- Absence baseline: the retained PR #973 native, built from the same main base
  plus unrelated SASA work, rejects the new keyword with `TypeError`. This
  expected single test failure is retained and is not a fresh main build.
  Base source also has no custom-charge constructor/keyword. Baseline native
  SHA256: `5b3d449d0923a9e9e05a98316ceabf75df4e6d9b632277956594939ebaca91cc`.

The supplied charge arrays are functional controls, not calibrated physical
models. These checks establish the explicit input contract and algebra on the
recorded cases, not GAFF parity, scientific accuracy or general replacement
of standard MMFF charges. No net-charge normalization is applied. The old
default constructor still assigns MMFF charges and keeps its near-zero cutoff;
the custom path retains nonzero supplied charges and refuses nonrepresentable
included-pair coefficients. Energy/gradient finiteness at arbitrary geometries
is outside that coefficient check.

## Reproduction

Run from this checkout using Python 3.13 on macOS arm64:

```sh
export CARGO_TARGET_DIR=/private/tmp/chematic-mmff-charges-target
export CARGO_INCREMENTAL=0 CARGO_BUILD_JOBS=1
export CARGO_PROFILE_DEV_DEBUG=0 CARGO_PROFILE_TEST_DEBUG=0
export CARGO_PROFILE_DEV_STRIP=none
export CARGO_PROFILE_DEV_OPT_LEVEL=1 CARGO_PROFILE_TEST_OPT_LEVEL=1
export CARGO_PROFILE_DEV_CODEGEN_UNITS=1 CARGO_PROFILE_TEST_CODEGEN_UNITS=1
export PYO3_PYTHON=/Library/Frameworks/Python.framework/Versions/3.13/bin/python3
cargo test -p chematic-ff --offline --lib --test mmff_dielectric --test mmff_partial_charges
cargo clippy -p chematic-ff -p chematic-py --all-targets --offline --no-deps -- -D warnings
cargo rustc -p chematic-py --offline -- -C strip=none -C link-arg=-undefined -C link-arg=dynamic_lookup
mkdir -p /private/tmp/chematic-mmff-charges-python/chematic
cp "$CARGO_TARGET_DIR/debug/libchematic.dylib" /private/tmp/chematic-mmff-charges-python/chematic/chematic.so
cp crates/chematic-py/python/chematic/{__init__.py,__init__.pyi,rdkit_compat.py,py.typed} /private/tmp/chematic-mmff-charges-python/chematic/
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 PYTHONPATH=/private/tmp/chematic-mmff-charges-python python3.13 -m pytest -q crates/chematic-py/tests/test_mmff_partial_charges.py crates/chematic-py/tests/test_mmff_dielectric.py
PYTHONPATH=/private/tmp/chematic-mmff-charges-python python3.13 validation/openbabel_issues/mmff-partial-charges/probe.py --output /private/tmp/mmff-charge-analytic.json
clang++ -std=c++17 validation/openbabel_issues/mmff-partial-charges/openbabel_probe.cpp -I/opt/homebrew/Cellar/open-babel/3.2.1/include/openbabel3 -L/opt/homebrew/Cellar/open-babel/3.2.1/lib -lopenbabel -o /private/tmp/mmff-charge-ob-probe
/private/tmp/mmff-charge-ob-probe
```

The pinned coordinate corpus was reused from the earlier dielectric controls;
its source/hash and generation provenance remain in
`validation/openbabel_issues/mmff-dielectric/`. No new embedding randomness is
used by the custom-charge probe.
