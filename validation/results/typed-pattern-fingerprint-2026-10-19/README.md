# Typed Pattern fingerprint facade

Base: `00596fb1ea5397524f98e3789c41f17ee8616fdc`. Checked 2026-10-19.
RDKit [#9601](https://github.com/rdkit/rdkit/issues/9601) remains open and
includes Pattern fingerprint mutable-output typing. Chematic had a native
default calculator but no `rdkit_compat.PatternFingerprint` entry point.

## Behavior and evidence

The new facade validates its molecule, width, boolean mode, structural bit mask
and plain-list/NumPy counter inputs before native calculation. The pinned RDKit
Pattern implementation does not use `atomCounts` or `setOnlyBits` after size
checks; the facade preserves this explicitly, without inventing atom counters or
mask filtering. Official source URL, Git blob ID and SHA256 are recorded in
`upstream-source.json` for tag `Release_2026_09_1`. Options with different
semantics from Layered fingerprints must not share Layered's output mutation.

The supported widths are powers of two in 1..2048. Folding is exact since each
width divides 2048. RDKit accepts other widths and tautomer mode; this facade
rejects them. Query/bundle overloads are absent. Graph caps (1024 atoms/4096
bonds) are not a new matcher visit budget or deadline.

- **377 tests passed**, including 54 new adapter cases and existing native
  fingerprint, Layered counter and RDKit-compatibility tests.
- **2736/2736** complete bit-vector comparisons match RDKit **2026.09.1**:
  76 molecules × 12 widths × three input profiles. The samples are the first
  32 rows of each hash-pinned descriptor/ChEMBL corpus plus twelve authored
  simple probes. These are development fixtures, not outcome-blind holdouts.
- **1824** supplied counter/mask profiles preserve seed values and match bits
  despite an empty mask. Profiles are plain lists and reversed strided uint32
  views; default calls are counted only in bit-vector comparisons.
- Scoped strict mypy passes the new helper and valid public consumer. The
  invalid consumer must report four argument errors. Legacy imports are silent;
  whole-package strict typing is not claimed.

Tests also cover result ownership, explicit unsupported widths/modes, empty
graphs/destinations, mask validation/iteration failure, error identity, no native
call on invalid options, unchanged counters on native failure, and layout/type
limits. An unrelated source's own side effects are not rolled back.

## Local runtime boundary and reproduction

Tests load current Python sources from `/private/tmp/chematic-pattern-python`
with the existing v1.1.7 native extension from
`/private/tmp/chematic-blog-v117-venv`. Its hash is in `oracle-summary.json`.
The Pattern core and fingerprint binding files have no diff between tag v1.1.7
and the base. This is not a fresh native build or published-wheel validation.
No Rust code changed. The existing chemistry/model compatibility limits of the
native calculator remain; this sample's zero mismatches is not universal parity.
Unit tests use NumPy 2.5.3; the reference environment supplies NumPy 2.5.4.

```sh
PYTHONDONTWRITEBYTECODE=1 PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 \
PYTHONPATH=/private/tmp/chematic-pattern-python \
python3.13 -m pytest crates/chematic-py/tests/test_pattern_fingerprint_adapter.py \
  crates/chematic-py/tests/test_fp.py crates/chematic-py/tests/test_layered_counts.py \
  crates/chematic-py/tests/test_rdkit_compat.py -q -p no:cacheprovider

PYTHONDONTWRITEBYTECODE=1 \
PYTHONPATH=/private/tmp/chematic-pattern-python:/private/tmp/chematic-rdkit-2026091/lib/python3.13/site-packages \
python3.13 validation/results/typed-pattern-fingerprint-2026-10-19/oracle.py

MYPYPATH=/private/tmp/chematic-pattern-python \
python3.13 -m mypy \
  --config-file validation/results/typed-pattern-fingerprint-2026-10-19/mypy.ini \
  crates/chematic-py/python/chematic/_pattern_inputs.py \
  validation/results/typed-pattern-fingerprint-2026-10-19/typing_valid.py
```
