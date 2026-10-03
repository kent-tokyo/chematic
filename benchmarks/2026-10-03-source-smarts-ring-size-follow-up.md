# Source-only SMARTS ring-size follow-up

This follows the [310,000-cell source profile](2026-10-03-source-smarts-310k-profiles.md)
against pinned RDKit 2026.03.6. It is **not** a published package or a
general SMARTS parity result. The native matcher and its 200 match-set / 43
Boolean v1.0.30 residuals are unchanged. Only the opt-in RDKit-style matcher
now evaluates `[kN]` using a bounded symmetrized ring set, with a typed
resource-limit error rather than a partial answer.

| Opt-in outcome across 10,000 molecules × 31 queries | Before | After |
|---|---:|---:|
| Published/default match-set residuals corrected | 154 | 164 |
| Remaining wrong-confident match sets | 29 | 19 |
| Published/default Boolean residuals corrected | 37 | 37 |
| Remaining wrong-confident Booleans | 6 | 6 |
| Typed refusal cells | 18 | 18 |
| Newly wrong-confident cells | 0 | 0 |

The ten corrected match sets are `[k6]` (9) and `[k5]` (1); they already
had correct hit/no-hit Booleans. One `[k6]` match-set residual remains (row
3498). Remaining match-set errors are concentrated in `[R1]` (3), `[R2]`
(6), `[R3]` (5), plus four organometallic/ring-topology cells. The six
Boolean errors are unchanged. The 18 refusals are the prior ambiguous
ring-count cells; 17 had wrong match sets in the default lane and one had a
correct match set. All 18 previously had correct Booleans. This is a smaller
wrong-confident domain, **not** an unconditional 310,000-cell parity claim.

The targeted regression is a bridged boronate molecule where atom 46 is in
a six-membered ring selected by RDKit but absent from chematic's plain SSSR
basis. The opt-in lane finds it; native matching remains unchanged. A second
test requires typed refusal when ring candidate generation exceeds its cap.

To reproduce, run `rdkit_parity_dump` on
`validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi`, then
`scripts/check_smarts_310k_profiles.py --dump <jsonl> --output <summary.json>`
as in the parent profile. The complete external source dump had SHA-256
`58636b38b443d858f362581e554fe39945452986c32f53bf330f49e74f9f808a`;
the checker summary had SHA-256
`85fa18328c51e3e6e107857fcf796d86d056800c0b28cc876c9926ce549449d2`.
These intermediates are not checked in. The checker verified corpus/query
hashes, exact source-default baseline cell identities, row completeness, and
zero newly wrong-confident cells. No Python/npm binding or release artifact
was measured here.
