# Opt-in SMARTS `[k6]` mixed-size replacement

Source-only comparison against pinned RDKit **2026.03.6** on the exposed
10,000-molecule × 31-query corpus. This is not a published package or a claim
of complete SMARTS equivalence. The candidate starts from main `6185ea96`.

For corpus row 3498, RDKit's `[k6]` set includes atom 22 in the six-ring
1–2–20–21–22–23. The shared symmetrized ring set retained only the four-ring
through atom 22, but the existing bounded SMARTS-specific replacement selector
already found the six-ring. The opt-in `[kN]` path now considers the accepted
replacement rings as well. Native SMARTS and the shared perception API are
unchanged.

| Pinned 310,000-cell lane | Before | Candidate |
|---|---:|---:|
| Published/default match-set residual baseline | 200 | 200 |
| Opt-in corrected match sets | 170 | **171** |
| Opt-in wrong-confident match sets | 13 | **12** |
| Opt-in corrected hit/no-hit results | 38 | 38 |
| Opt-in wrong-confident hit/no-hit results | 5 | 5 |
| Opt-in typed-refusal cells | 18 | 18 |

No newly wrong-confident cell appeared. The remaining 12 match-set cells are
the three fullerene-like rows 3837/3976/4004 (`[R2]` and `[R3]`) and the
organometallic row 8341 (six ring and ring-bond predicates). Five of these
also differ in hit/no-hit. The 18 refusals are not exact matches.

Reproduce from this source candidate:

```sh
cargo run -p chematic-smarts --example rdkit_parity_dump -- \
  validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  > /tmp/chematic-smarts-k6-dump.jsonl
python3 scripts/check_smarts_310k_profiles.py \
  --dump /tmp/chematic-smarts-k6-dump.jsonl \
  --output /tmp/chematic-smarts-k6-profile.json
```

The completed local dump SHA-256 was
`6a978e86be4d6b5706e39daf552f067c3fd7522dedb78bbc7df9ff5956e02d0c`;
the checked JSON profile SHA-256 was
`e83957d13ad2afde7c739b4c38982dd46d4a7a2e36c0bfe8c4ebde176417d50f`.
The checker also verified unchanged 200/43 source-default cell identities,
input and oracle hashes, row order, and the complete 10,000-row footer.
