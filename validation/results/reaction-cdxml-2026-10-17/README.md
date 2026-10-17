# Checked reaction CDXML candidate — 2026-10-17

Base: `ad0283a8101dfc6f7d5ea6b6bdea94a35f8404f0` (main).
The base exposes molecular `write_cdxml`, but no reaction CDXML writer.
This is a source-level capability addition; a main-baseline native extension
was not rebuilt. The candidate runtime is isolated from installed packages.

Related open issue: Open Babel #2703. The saved upstream snapshot records its
state and URL. Three original reaction-SMILES probes reproduce the missing
scheme/step/arrow on Open Babel 3.2.1. The resulting CDXML and stderr are retained.
The new API preserves caller roles and layout with an explicit limited graph
profile; it does not implement automatic depiction or every reaction annotation.

## External-reader evidence

The original `validation/openbabel_issues/reaction-cdxml/probe.py` builds 256
documents from eight symbolic recipes, all role-order reversals and four layout
transforms (translation, rotation, reflection/scaling included). These examples
exercise data transport; they do not establish chemical feasibility or balance.

- XML checks cover ordered role references, duplicate components, globally
  unique IDs, fragment-local bonds, exact coordinates and arrow endpoints.
- Indigo 1.46.0 independently reads reactants, agents and products correctly
  for all 256 documents, compared with its canonicalization of original inputs.
- RDKit 2026.09.1 reads reactants/products correctly for all 256 documents.
  It omits agents in all 160 documents containing agents on this host; these
  differences are explicitly retained in `summary.json` and raw records.
  Zero scoped failures does not mean full RDKit role parity.
- Four separate atom-state diagnostics use the unchanged legacy molecular
  writer. Indigo changes `[O]`, `[N]` and `[H]` on read; `[CH3]` is preserved.
  The new API conservatively refuses all four via its explicit-H/isolated-H
  profile, rather than silently changing their state. Other molecular writer
  behavior is not changed by this candidate.

`corpus-results.json.gz` contains complete original inputs, emitted documents and
external-reader role outputs. `summary.json` records versions, candidate binary
hash, Open Babel results and diagnostics. `source-hashes.json` records source and
build settings. The candidate binary is not a published wheel.

## Local validation

802 existing MOL library tests plus 10 new integration tests passed. Python:
46 new API tests plus 83 existing CDX/CDXML/RXN tests passed. Strict Clippy on
both changed crates, all targets, passed with warnings denied. Formatting,
refactor budget, source stubs, binding inventory, cross-binding manifest,
unsafe surface, product spelling, workflow pins and strict MkDocs passed.
Logs are retained with trailing whitespace normalized. Hosted CI and main
integration are separate from these local checks.

The new Rust target is listed explicitly in CI's integration matrix; Python's
existing full-suite job discovers the new test file.

## Reproduction

Use Python 3.13 with NumPy, pytest, RDKit 2026.09.1 and Indigo 1.46.0; Open Babel
CLI 3.2.1 is used for the issue diagnostic. The probe currently locates obabel
at `/opt/homebrew/bin/obabel` on this macOS arm64 host.

Build with a writable CARGO_TARGET_DIR and CARGO_INCREMENTAL=0,
CARGO_BUILD_JOBS=1, dev/test debug=0, opt-level=1, codegen-units=1,
and dev strip=none. On this host:

```sh
cargo rustc -p chematic-py --offline -- \
  -C strip=none -C link-arg=-undefined -C link-arg=dynamic_lookup
```

Copy the dylib as `chematic/chematic.so` beside package Python files in a temporary
package directory and put that directory plus the isolated RDKit package on
PYTHONPATH. Existing installed chematic packages are not overwritten.

```sh
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.13 -m pytest \
  crates/chematic-py/tests/test_reaction_cdxml.py \
  crates/chematic-py/tests/test_chemdraw_block.py \
  crates/chematic-py/tests/test_cdxml_stereo_options.py \
  crates/chematic-py/tests/test_cross_binding_contract.py::test_python_binding_matches_shared_rxn_document_contract \
  crates/chematic-py/tests/test_cross_binding_contract.py::test_python_rxn_document_contract_is_loss_aware \
  -q -p no:cacheprovider
cargo test -p chematic-mol --offline --test reaction_cdxml --lib
cargo clippy -p chematic-mol -p chematic-py --all-targets --offline --no-deps -- -D warnings
python3.13 validation/openbabel_issues/reaction-cdxml/probe.py <output-directory>
```

External pytest plugins are disabled because the host's old pytest-asyncio
plugin fails collection with pytest 9. These tests do not require that plugin.
