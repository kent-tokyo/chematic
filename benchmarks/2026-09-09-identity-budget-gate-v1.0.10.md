# Identity budget gate — v1.0.10

The focused local identity gate passed all 263 tests:

- `chematic-smarts` library: 184/184
- `chematic-smiles` canonical module: 79/79

Reproduce with:

```sh
CARGO_TARGET_DIR=/private/tmp/chematic-target-identity-budget cargo test -p chematic-smarts --offline --lib
CARGO_TARGET_DIR=/private/tmp/chematic-target-identity-budget cargo test -p chematic-smiles --offline --lib canonical::tests
```

The SMARTS suite covers exact-output matching, visit budgets, MCS timeouts and
resource limits, and ring-model exhaustion. The canonical module suite covers
the current atom-order, E/Z, stereo, idempotency, and fail-closed regressions.

This is focused local evidence, not completion of the broader canonical
atom-order/E/Z roadmap item. The three aromatic-stash residuals remain
intentionally fail-closed, and the broader stable-operation manifest remains
open.
