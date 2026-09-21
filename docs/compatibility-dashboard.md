# Compatibility dashboard

Generated from checked-in manifests by `python3 scripts/generate_compatibility_dashboard.py`.
This is a compatibility-contract dashboard, not a universal RDKit parity or speed claim.

- Target version: `1.0.19`
- Regeneration: deterministic, offline, clean-checkout compatible
- Contract manifest: `validation/cross_binding_contract.json` (SHA-256 `2c08569608a26d40239c2b331ae9890257406afd672ca3fca35334610ba31966`)
- Streaming matrix: `validation/results/cross-engine-matrix-v1.0.13.json` (SHA-256 `84870b6ee1b402327fec33df7fc76cb0f36301f595758d575467066818405330`)

## Shared binding contract

Operation inventory: `58` currently shared operations; validate with `python3 scripts/check_cross_binding_manifest.py`.

| Area | Checked-in assertions | Status |
|---|---:|---|
| `fixtures` | 4 | covered by versioned fixture contract |
| `descriptor_contract` | 4 | covered by versioned fixture contract |
| `standardization_contract` | 10 | covered by versioned fixture contract |
| `fingerprint_contract` | 6 | covered by versioned fixture contract |
| `fingerprint_detail_contract` | 4 | covered by versioned fixture contract |
| `reaction_application_contract` | 21 | covered by versioned fixture contract |
| `adversarial` | 8 | covered by versioned fixture contract |

## Streaming record/failure contract

Pinned matrix: 20 repetitions across 10 formats.
- Matrix status: `historical for target`; matrix target is `1.0.13` and was not remeasured for the current release.
The engine/process boundaries remain explicit; records and failures are the only cross-engine claims.

| Format | Expected records | Engines with zero failures |
|---|---:|---:|
| `sdf` | 40 | 3/3 |
| `mol` | 40 | 3/3 |
| `xyz` | 40 | 3/3 |
| `extxyz` | 40 | 2/2 |
| `v3000` | 20 | 3/3 |
| `mol2` | 20 | 3/3 |
| `cml` | 20 | 2/2 |
| `cdxml` | 20 | 2/2 |
| `mmcif` | 20 | 2/2 |
| `pdb` | 20 | 2/2 |

## Boundaries

- A `covered` row means the checked-in assertion and its local consumer tests exist; it does not mean every edge case is supported.
- Missing optional engines, external review, publication, and broader corpus quality remain separate release gates.
- Regenerate after changing either source manifest and review the resulting digest changes before committing.

## API profile Compatibility Contract

Each row separates support status, API profile, comparator lane, and measurement coverage. `not_measured` is a visible gap, not zero coverage or a passing result.

| Operation | Support | Profile | Comparator | Coverage |
|---|---|---|---|---|
| `smiles_parse_write` | `stable` | `native` | — | not_measured |
| `canonical_identity` | `experimental` | `native` | RDKit `2025.09.3` | partial |
| `aromaticity` | `experimental` | `rdkit_compatibility` | RDKit `2025.09.3` | partial |
| `cip_labels` | `experimental` | `rdkit_compatibility` | RDKit `2025.09.3` | partial |
| `smarts_substructure` | `experimental` | `rdkit_compatibility` | RDKit `2025.09.3` | 2025.09.3: 155651 cells, 145588 matches, 21 residuals, 10042 oracle parse errors; 2026.03.6 separate lane: 155651 cells, 155629 matches, 22 residuals |
| `morgan_ecfp` | `experimental` | `rdkit_compatibility` | RDKit `2025.09.3` | partial |
| `mol_sdf_v2000` | `stable` | `native` | — | not_measured |
| `mol_sdf_v3000` | `experimental` | `native` | — | partial |

Evidence paths, API names, settings, and exact/numeric/semantic comparison rules live in `validation/compatibility_profiles.json`; validate them with `python3 scripts/check_compatibility_profiles.py`.

## RDKit accuracy profiles

This section is generated from `validation/manifests/rdkit_accuracy_v2.json`. It records declared evidence lanes; it does not convert development or exposed data into a sealed evaluation.

- Manifest status: `development_not_sealed`
- Comparator: `RDKit 2025.09.3`

| Operation | Scope | Split | Expected rows | Profile |
|---|---|---|---:|---|
| `descriptor_eight_field` | `rdkit_descriptor_semantics_v1` | `development` | 5000 | `rdkit_compatibility` |
| `descriptor_unused_holdout_workflow` | `rdkit_descriptor_semantics_v1` | `reserved_unused_workflow_probe` | 4 | `rdkit_compatibility` |
| `descriptor_eight_field_binding_exposed` | `rdkit_descriptor_binding_consistency_v1` | `exposed_holdout_not_sealed` | 7737 | `binding_consistency` |

A missing profile is an explicit contract gap: native, RDKit-compatible, and binding-consistency lanes must not be conflated in a scorecard.
