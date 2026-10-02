# Reaction compatibility: checked source profile (2026-10-02)

This is an **unpublished source-candidate** result, not a rerun of published
v1.0.30. It builds on the published-artifact packet and retains the same 57
historical plus 26 stratified exposed cases. The older `run_smirks` API is
unchanged; `run_smirks_checked(..., rdkit_compat=True)` is opt-in and returns
typed outcomes and diagnostics.

The local CPython 3.13 arm64 development wheel has SHA-256
`37a854ef2c8fc534d96a53c2e1f2644fb95d79354bf60e6e291c34a31a345269`.
The [83-row report](../validation/results/v1.0.30-source-checked-reaction-83-strata.json)
has SHA-256
`e676f63364d1375b8ce8222c80c83dd2acf02369ae02d3a52c62cbafd2035078`.
The report pins RDKit 2026.03.6 and both fixture hashes; it records each
reaction, oracle product set, candidate product set or refusal reason, and
stratum. These are product-graph/multiplicity comparisons, not reaction yield
or selectivity measurements.

| Outcome | Cases | Meaning |
|---|---:|---|
| Semantic match | 76 | Same canonicalized product-graph sets and multiplicity |
| Typed unsupported | 3 | Reactant tetrahedral `@`/`@@` template semantics differ from RDKit 2026.03.6; opt-in profile declines them |
| Typed refusal | 1 | Isotope-labeled split candidate fails CheMatic product-valence validation; filtered product is reported, not hidden |
| Jointly invalid input | 3 | Both systems reject the malformed or missing-reactant case |
| Wrong confident | 0 | No confident product mismatch in this bounded corpus/profile |

The three stereo cases include one that the prior source candidate matched,
but the checked profile deliberately refuses all reactant-side tetrahedral
templates until their semantics are adjudicated. Product-side stereo and E/Z
templates are not blanket-refused. `partial_products` and truncated matches
are explicitly non-success statuses; the 83-case corpus does not exercise
truncation. The binding reports accepted, applied, valence-rejected and
truncated match counts. The runner's default mode still measures the legacy
published API and must not be switched retroactively.

Reproduce with a source wheel built from this branch and an isolated Python
environment containing RDKit 2026.03.6:

```sh
python scripts/run_reaction_compatibility_v2.py \
  --artifact /path/to/source-candidate.whl \
  --output /path/to/reaction-checked.json \
  --expected-rdkit 2026.03.6 --checked-rdkit-compat
```

Verification on the candidate: 231 Rust reaction tests passed; five targeted
Python checked-API tests passed; the complete Python suite passed 1,000/1,000.
The reported wheel was rebuilt after adding fail-closed truncation handling.
Product atom-map provenance was **not** compared with RDKit, and no broad
SMIRKS compatibility or published-package claim follows from these 83 cases.
