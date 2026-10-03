# Reaction 83: current-source graph, origin and template-map audit

This is **source-candidate evidence**, not a v1.0.31 published-package rerun.
The pinned baseline is the [published v1.0.30 83-row audit](2026-10-02-v1.0.30-published-artifact-gates.md).
Source commit `ab4eaf9d63a6ffb075d9eba6a9dd0240ffa18856`, Rust 1.97.0,
macOS arm64, and RDKit 2026.03.6 (wheel SHA-256
`e16c467cb254a223e59a0cf81358c6b39da15a99d2909170d693e95778fddb41`)
were used with the unchanged 57+26 exposed fixtures. The pinned input hashes,
row classifications and output hashes are in the [machine-readable summary](../validation/results/v1.0.31-source-reaction-83-checked-summary.json).

| Outcome | Published v1.0.30 Rust | Current source, native | Current source, opt-in RDKit-checked |
|---|---:|---:|---:|
| Graph + atom origin + template map agree | 73 | 77 | 76 |
| Confident graph/origin/map difference | 6 | 2 | 0 |
| Typed unsupported | 0 | 0 | 3 |
| Diagnosed valence refusal | 1 | 1 | 1 |
| Jointly invalid input | 3 | 3 | 3 |

In native mode the two remaining confident differences are the alanine
reactant-tetrahedral cases. The opt-in RDKit 2026.03.6 profile refuses *three*
reactant-tetrahedral templates with the shared typed reason
`chiral_reactant_template_semantics`: one previously matching case and both
native differences. It does not blanket-refuse product-side stereo or E/Z.
The existing `v2_isotope_methanol_split` valence refusal remains separate; it
is not counted as parity. The three jointly invalid inputs are not successes.
Among the six published confident residuals, source matching now preserves
both ether-cleavage graph/origin assignments, the styrene template-map
orientation and the isotope no-match decision. No ordinary reaction output
was changed by the checked-profile addition.

The [source runner](../tools/source_reaction_gate/Cargo.toml) reuses the
published runner logic but links the local crate; its lockfile is committed.
From repository root, build both binaries with the same label and profile:

```sh
CHEMATIC_REACTION_GATE_LABEL='chematic 1.0.31 source ab4eaf9d checked' CHEMATIC_REACTION_GATE_PROFILE=rdkit-2026.03.6 cargo run --manifest-path tools/source_reaction_gate/Cargo.toml --locked --bin reactions -- validation/reaction_product_parity_cases.json validation/reaction_product_parity_strata_v2.json OUT/rows.json OUT/summary.json
CHEMATIC_REACTION_GATE_LABEL='chematic 1.0.31 source ab4eaf9d checked' CHEMATIC_REACTION_GATE_PROFILE=rdkit-2026.03.6 cargo run --manifest-path tools/source_reaction_gate/Cargo.toml --locked --bin reaction_template_maps -- validation/reaction_product_parity_cases.json validation/reaction_product_parity_strata_v2.json OUT/map-rows.json OUT/map-summary.json
python3 scripts/reaction_template_map_gate.py --rust-rows OUT/map-rows.json --rust-summary OUT/map-summary.json --baseline-rows OUT/rows.json --output OUT/report.json --expected-rdkit 2026.03.6 --expected-crate 'chematic 1.0.31 source ab4eaf9d checked' --expected-version 1.0.31 --expected-profile rdkit-2026.03.6
```

`OUT` must be an existing directory; run the Python command in an isolated
environment containing pinned RDKit 2026.03.6. Full raw rows are regenerable
but not checked in for this source diagnostic. Their SHA-256 values are in the
summary. The published-v1.0.30 comparator default was rerun unchanged and
still reports 73 exact, six confident differences, one refusal and three
jointly invalid. Rust reaction tests passed 234/234, Python binding type-check
passed, and four comparison-script tests passed. No new Python wheel or npm
artifact was produced; cross-binding and published-package acceptance remain
open. This corpus does not establish general SMIRKS parity, yield or
selectivity accuracy.
