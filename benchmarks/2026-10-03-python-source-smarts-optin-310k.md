# Opt-in Python SMARTS profile: 310,000-cell source-wheel gate

This is a **locally built source wheel**, not the v1.0.31 PyPI wheel or any
other published artifact. The base is post-PR #733 main `896abb8b`, plus the
opt-in Python binding in this change. RDKit is pinned to **2026.03.6**. The
10,000-SMILES corpus, 31 queries, oracle and published v1.0.30 baseline are
SHA-256-pinned in the [machine-readable result](../validation/results/v1.0.31-source-python-smarts-parity-310k.json).

| Outcome | Cells |
|---|---:|
| Exact match sets vs RDKit | 309,982 / 310,000 |
| Typed unsupported (`ring_model_ambiguous`) | 18 |
| Wrong confident match sets | 0 |
| Original 200 match-set residuals corrected | 183 |
| Original 200 residuals refused | 17 |
| Previously correct cells refused | 1 |
| Original 43 hit/no-hit residuals corrected | 43 |

The 18 refusals are `[R1]`, `[R2]`, `[R3]` on six charged polycycles. They
are **not** exact matches. The default Python `find_matches` remains unchanged.
The new `find_matches_rdkit_parity` returns `status`, `matches`, and a stable
`reason`: `matches=[]` means a completed no-match, while `matches=None` means
typed unsupported/refusal. This API is a bounded profile, not general RDKit
SMARTS compatibility. The [ring-family diagnosis](2026-10-03-smarts-ring-family-diagnostic.md)
explains why the 18 cells remain guarded.

Reproduce with a source-built wheel installed into an isolated environment
that also has RDKit 2026.03.6:

```sh
python scripts/check_python_smarts_parity_310k.py \
  --wheel dist/hba-wheel/chematic-*.whl \
  --module-root .venv \
  --output validation/results/smarts-checked-source-wheel-310k-local.json
```

The checker validates import provenance, wheel and input hashes, all 310,000
match sets, the exact 18 refusal cells, and the predeclared outcome counts.
Post-merge PR [#735](https://github.com/kent-tokyo/chematic/pull/735),
run `37101314517`, built release-profile source wheels on Linux and macOS
arm64. Both retained reports (`rdkit-hba-source-wheel-5k` and
`rdkit-hba-source-wheel-5k-macos`) were downloaded and independently read:
each reports 309,982 exact, 18 typed unsupported, 183 corrected original
match-set residuals, 43 corrected Boolean residuals and zero unexpected
cells. All 76 run checks passed (two intentionally skipped). These are
**CI source wheels**, not PyPI artifacts. A published wheel, npm/WASM profile
and independent held-out corpus remain unverified.
