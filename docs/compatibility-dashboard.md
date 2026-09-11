# Compatibility dashboard

Generated from checked-in manifests by `python3 scripts/generate_compatibility_dashboard.py`.
This is a compatibility-contract dashboard, not a universal RDKit parity or speed claim.

- Target version: `1.0.12`
- Regeneration: deterministic, offline, clean-checkout compatible
- Contract manifest: `validation/cross_binding_contract.json` (SHA-256 `223d37f4af3ef2a0b64141e75d8661fbeb72acb686603fb7d5fe0a5c813f404e`)
- Streaming matrix: `validation/results/cross-engine-matrix-v1.0.12.json` (SHA-256 `fc8ad57ebfad30145436d98776363b0cf6852c5430f65c3f210074a15d580526`)

## Shared binding contract

Operation inventory: `56` currently shared operations; validate with `python3 scripts/check_cross_binding_manifest.py`.

| Area | Checked-in assertions | Status |
|---|---:|---|
| `fixtures` | 4 | covered by versioned fixture contract |
| `descriptor_contract` | 4 | covered by versioned fixture contract |
| `standardization_contract` | 10 | covered by versioned fixture contract |
| `fingerprint_contract` | 6 | covered by versioned fixture contract |
| `fingerprint_detail_contract` | 4 | covered by versioned fixture contract |
| `reaction_application_contract` | 18 | covered by versioned fixture contract |
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
