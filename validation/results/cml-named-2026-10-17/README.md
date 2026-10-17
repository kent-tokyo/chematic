# Named CML atoms — 2026-10-17

Local source implementation for the atom-name preservation request in open
[Open Babel #2120](https://github.com/openbabel/openbabel/issues/2120).
The exact base commit, native binary hash, source hashes and tool versions are
recorded in `summary.json` and `source-sha256.json.gz`.

## Baseline and scope

The retained PR #971 native module has no named-CML writer: the new fixture
test fails with `AttributeError`. This retained build contains unrelated
tautomer selection and is not a fresh rebuild of current main. Current main
source also lacks this API. The unchanged default writer emits `a1`, `a2`,
`a3`, confirmed again by the candidate probe.

Open Babel 3.2.1 converts an original synthetic sequential-ID water MOL2 to
CML with the same generated IDs, replacing `O_site`, `H_one`, `H-two`. Both
bonds are present in this diagnostic. The separate native fixture uses source
IDs 9, 72 and 200 to exercise correspondence; it is not used to infer an
Open Babel bond-parser defect. No upstream attachment or external molecule
corpus is copied.

This adds a graph writer and a documented MOL2 recipe, not a complete MOL2
metadata converter. Source partial charges, force-field types, substructures,
opaque data and source-only identities are not emitted. Unsupported stereo,
maps, tags, wildcard atoms and special bonds are refused. The existing CML
reader's graph/2D projection does not expose IDs or 3D coordinates.

## Independent validation

- Six native regression tests and 64 Python tests (35 new, 29 existing H
  accounting) pass with no failures or skips.
- 129 original cases validate IDs, matching bond references, all coordinate
  components, graph formula and canonical graph identity. These cover atom
  permutations, charged/isotopic/aromatic and ordinary graphs, and coordinate
  scales 1, 1e-200 and 1e200. One additional named-water fixture checks MOL2
  source names, arbitrary source IDs and exact coordinates.
- Every generated case passes the independently downloaded official CML
  Schema 3, checked by lxml 6.0.2. The schema URL and SHA-256 are retained;
  the external schema itself is not copied into the repository.
- `schema-dot-boundary.json` reproduces a specification discrepancy: the
  molecular convention permits dots in atom IDs, while Schema 3's atomRefs2
  list pattern rejects them. The writer uses their common subset and refuses
  dotted names even on isolated atoms, without automatic renaming.
- Clippy with warnings denied for both changed crates/all targets, formatting,
  stubs, binding inventory/manifest, product naming, unsafe-surface checks and
  strict MkDocs pass.

`comparison.json.gz` retains raw Open Babel CML/stderr, candidate CML, native
and schema hashes, and all case rows. Logs are compressed without changing
content. Source validation does not establish published-package chemistry,
WASM parity, hosted CI completion or main integration.

## Reproduction

Builds use one job, no incremental cache/debug data, opt-level 1 and one
codegen unit. `PYO3_PYTHON` points at the host's CPython 3.13.6. Commands:

```sh
cargo test -p chematic-mol --offline --test cml_named
cargo rustc -p chematic-py --offline -- \
  -C strip=none -C link-arg=-undefined -C link-arg=dynamic_lookup
cargo clippy -p chematic-mol -p chematic-py --all-targets --offline --no-deps -- -D warnings
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 PYTHONPATH=/private/tmp/chematic-cml-named-python \
python3.13 -m pytest -q crates/chematic-py/tests/test_cml_named.py \
  crates/chematic-py/tests/test_cml_hydrogen_contract.py
curl -L --fail https://www.xml-cml.org/schema/schema3/schema.xsd -o /tmp/cml-schema3.xsd
PYTHONPATH=/private/tmp/chematic-cml-named-python \
python3.13 validation/openbabel_issues/cml-named/probe.py \
  --schema /tmp/cml-schema3.xsd --output /tmp/cml-named-comparison.json
```

The retained source-built native module is
`/private/tmp/chematic-cml-named-python/chematic/chematic.so`.
The isolated worktree excludes large historic results and benchmark corpora
because the host is low on disk; the small ethanol V3000 fixture required by
existing all-target checks is included. Existing worktrees, source/evidence and
concurrent builds were preserved.
