# Canonical identity focused gate — v1.0.10

The focused local canonical identity suite passed all 36 tests:

- canonical E/Z residuals: 6/6
- canonical idempotency corpus: 2/2
- canonical robustness: 24/24
- Kekulé S0 stereo invariance: 4/4

Reproduce with:

```sh
CARGO_TARGET_DIR=/private/tmp/chematic-target-canonical-audit cargo test -p chematic-smiles --offline --test canonical_robustness --test canonical_ez_residual --test canonical_idempotency_corpus --test kekule_s0_stereo_invariance
```

This is focused local regression evidence, not completion of the full atom-order/E/Z roadmap item. The three aromatic-stash residuals remain fail-closed, and the broader corpus gate remains open.
