# Supplied-candidate tautomer selection — 2026-10-17

Source implementation motivated by open RDKit
[#9601](https://github.com/rdkit/rdkit/issues/9601). Base main is
`a4ca27e63bdf689f5f52e36bbddc02d6862f65ab`.

## Observations

- Main source had no public supplied-candidate selector. Three tests against
  the retained PR #968 native module failed with `AttributeError`. That module
  has an unrelated GRO/LAMMPS addition on the same base; it is not a fresh
  rebuild of unmodified main. `baseline.log.gz` retains the output.
- Both RDKit 2025.09.3 and 2026.09.1 terminated with SIGSEGV (`returncode=-11`)
  for `TautomerEnumerator().PickCanonical([])`. Each check ran in a separate
  process with a 15-second timeout and core dumps disabled. No large-input
  crash or application process was exercised.
- The new native function returns a typed empty-input error; Python raises
  `ValueError` with no winner. This adds a missing capability safely rather
  than fixing a pre-existing chematic null-pointer crash.

## Passed validation

- Five new Rust integration tests.
- 34 new Python tests and nine existing tautomer tests: 43 passed, no skips.
- Two live reference versions, each with 175 candidate sets and three input
  forms (list, tuple, generator): 1,050 matching selections in total. Score
  comparisons were also checked for all supplied rows: 882 matches.
- Corpus includes keto/enol and heteroaromatic forms, permutations, duplicate
  ties, tetrahedral/geometric stereo, isotope and atom-map labels, ions and an
  empty molecule. It also includes mixed unrelated candidates deliberately:
  no tautomer-family membership claim is made. Selected source-row properties
  verify ownership and deterministic duplicate handling.
- Clippy for both changed crates, all targets, with warnings denied. Format,
  Python stub syntax, binding inventory, cross-binding manifest, product names,
  unsafe-surface policy and strict MkDocs build passed.

The separate mixed specified/unspecified stereo tie records show the boundary
of the claim: the existing lexical policy chooses `CC(F)Cl`, while RDKit
2026.09.1 prefers `C[C@H](F)Cl`. It is excluded from the matching-corpus count.
The function does not apply RDKit's final stereo-cleanup transformation or
support custom scoring. Selection preserves the chosen source graph and typed
properties; it does not estimate physical stability.

## Reproduction

The retained native package is `/private/tmp/chematic-pick-tautomer-python`;
its extension SHA-256 is in `summary.json` and both compressed reference
reports (`rdkit-*.json.gz`). Exact
Rust/Python source-file hashes are in `source-sha256.json.gz`.

Build/test used one job, no incremental cache or debug information, opt-level 1
and one codegen unit. The CPython 3.13 build used `PYO3_PYTHON` pointing to
`/Library/Frameworks/Python.framework/Versions/3.13/bin/python3` and:

```sh
cargo test -p chematic-smiles --offline --test tautomer_pick
cargo rustc -p chematic-py --offline -- \
  -C strip=none -C link-arg=-undefined -C link-arg=dynamic_lookup
cargo clippy -p chematic-smiles -p chematic-py --all-targets --offline --no-deps -- -D warnings
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 \
PYTHONPATH=/private/tmp/chematic-pick-tautomer-python:/private/tmp/chematic-path-rdkit \
python3.13 -m pytest -q crates/chematic-py/tests/test_tautomer_pick.py \
  crates/chematic-py/tests/test_rdkit_tautomer.py
PYTHONPATH=/private/tmp/chematic-pick-tautomer-python:/private/tmp/chematic-path-rdkit \
python3.13 validation/rdkit_issues/tautomer_pick/probe.py --output /private/tmp/tautomer-reference.json
```

`/private/tmp/chematic-path-rdkit` supplies RDKit 2026.09.1; omitting it selects
this host's installed 2025.09.3. Logs are compressed without modifying their
contents. The candidate build initially hit low disk space; the final build
log records the successful resumed build. Only checked unused build caches
were removed; the native module, source and evidence were retained.

These are local source checks. No wheel publication, upstream change, hosted
CI completion or main integration is established by this packet.
