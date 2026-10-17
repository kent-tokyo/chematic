# Checked bridge fragmentation / typed caller counters

Candidate source based on `c6d10dc0218f190f88fb01f1e43fee0631575c22`.
Upstream [RDKit #9601](https://github.com/rdkit/rdkit/issues/9601) was OPEN at
inspection; it lists `FragmentOnBonds.cutsPerAtom` among remaining typed mutable
outputs. `upstream-issue.json` records the inspection. This does not assert that
fragmentation itself is missing from RDKit or that chematic covers all its modes.

## Results

- Fresh macOS arm64 native PyO3 extension built from candidate source, native
  InChI enabled; not a wheel or registry/release verification.
- Rust chem library: 962 passed, 1 existing ignored test, 0 failures; four new
  bridge fragmentation tests are included. Focused run: 4 passed.
- Python: 217 passed (55 new fragmentation + 162 existing RDKit compatibility).
- RDKit 2026.09.1 oracle: 30 source molecules, 780 nonempty cut/label combinations
  match counts, canonical SMILES graph identity and V3000 indexed atom/bond fields.
  Another 90 no-op conditions check chematic's defined clone/zero behavior against
  a reference clone; **RDKit FragmentOnBonds rejects empty cuts**.
- Source bond indices are mapped by endpoints into RDKit's independent parse;
  label orientation is adjusted for any reversed source bond. Zero, high isotope
  65535, custom and default labels; multiple/adjacent/reordered cuts; aromatic
  substituents, charges, atom maps, isotopes and explicit hydrogens are covered.
- Query metadata of V3000 `R` dummies, coordinates, stereo transfer, ring cuts,
  alternative bond orders and chemistry outside the bounded profile are untested.
- Strict mypy accepts helper/public consumer and rejects all six negative consumer
  lines. Runtime rejects bool even though Python's type system treats bool as int.
- Source stub, cross-binding manifest, unsafe surface, binding surface, refactor
  budget, product spelling and formatting gates; Clippy for both changed crates;
  strict MkDocs; three executed guide blocks and three rendered API anchors.

**Known export residuals:** pre-existing V2000 writer corrupts some isotope dummy
mass offsets, and MolJSON serializes dummy placeholder carbon. `export_residuals`
records an actual reproducer. Neither serializer is used as fragment graph oracle.
No writer changes are included. Use SMILES; validated V3000 fields are listed above.

## Reproduction

Set `CANDIDATE_PY` to a package directory containing candidate Python source/stubs
and the freshly built native extension, and `RDKIT_PY` to RDKit 2026.09.1's package
directory. The run used Python 3.13.6 and pytest 9.1.1. Native artifact and source
hashes are in `provenance.json`.

```sh
PYO3_PYTHON=/opt/homebrew/bin/python3.13 CARGO_PROFILE_DEV_DEBUG=0 \
  CARGO_INCREMENTAL=0 cargo rustc -p chematic-py --features native-inchi -- \
  -C link-arg=-undefined -C link-arg=dynamic_lookup
cargo test -p chematic-chem --lib
cargo clippy -p chematic-chem -p chematic-py -- -D warnings
PYTHONPATH="$CANDIDATE_PY" python3.13 -m pytest \
  crates/chematic-py/tests/test_bond_fragment.py \
  crates/chematic-py/tests/test_rdkit_compat.py -q -p no:cacheprovider
PYTHONPATH="$CANDIDATE_PY:$RDKIT_PY" python3.13 \
  validation/results/bond-fragment-counters-2026-10-17/oracle.py
PYTHONPATH="$CANDIDATE_PY:$RDKIT_PY" python3.13 \
  validation/results/bond-fragment-counters-2026-10-17/export_residuals.py
MYPYPATH="$CANDIDATE_PY" python3.13 -m mypy --strict --follow-imports=silent \
  --follow-imports-for-stubs --ignore-missing-imports --no-incremental \
  crates/chematic-py/python/chematic/_bond_fragment.py \
  validation/results/bond-fragment-counters-2026-10-17/typing_valid.py
# Replace typing_valid.py by typing_invalid.py: expect exit 1 and six errors.
python3.13 -m mkdocs build --strict --site-dir /tmp/chematic-bond-fragment-site
```

On macOS the explicit dynamic lookup link flags are necessary for this local
`cargo rustc` extension build. CI/wheel builds use their platform packaging flow.
No performance comparison is claimed.
