# Compatibility dashboard

Generated from checked-in manifests by `python3 scripts/generate_compatibility_dashboard.py`.
This is a compatibility-contract dashboard, not a universal RDKit parity or speed claim.

- Target version: `1.0.14`
- Regeneration: deterministic, offline, clean-checkout compatible
- Contract manifest: `validation/cross_binding_contract.json` (SHA-256 `633e7a8ecfe2d0595d94d697742e3cbb2cbfaa6d95ad9de81a19a90434cc7930`)
- Streaming matrix: `validation/results/cross-engine-matrix-v1.0.13.json` (SHA-256 `84870b6ee1b402327fec33df7fc76cb0f36301f595758d575467066818405330`)

## Shared binding contract

Operation inventory: `57` currently shared operations; validate with `python3 scripts/check_cross_binding_manifest.py`.

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
