# Opt-in SMARTS dense-cage ring-count follow-up

Source-only comparison against pinned RDKit **2026.03.6** on the exposed
10,000-molecule × 31-query corpus. This is not a published package or a
general SMARTS-parity claim. The candidate starts from main `da70d2c6`.

Rows 3837, 3976 and 4004 contain closely related fullerene-like fused cages.
In row 3837, RDKit records 33 rings while CheMatic's minimal and prior
symmetrized sets recorded 32. A missing six-membered ring raises six atoms
from `[R2]` to `[R3]`. The graph-wide D2 search started from peripheral
substituents, not the cage. The opt-in ring-count model now peels degree-<3
vertices from the basis-ring bond graph and, only when its earlier passes found
no extra ring, blocks each ring edge at surviving 3-core roots for a bounded
second-shortest-ring search. Candidate rings still must pass the existing
same-size, shared-bond and unique-bond-preservation rules. Native SMARTS and
general ring perception are unchanged.

| Pinned 310,000-cell lane | Main before | Candidate |
|---|---:|---:|
| Published/default differing match sets | 200 | 200 |
| Opt-in corrected match sets | 171 | **177** |
| Opt-in wrong-confident match sets | 12 | **6** |
| Opt-in corrected hit/no-hit results | 38 | **41** |
| Opt-in wrong-confident hit/no-hit results | 5 | **2** |
| Opt-in typed-refusal cells | 18 | 18 |

No newly wrong-confident cell appeared. The six remaining match-set cells
(`*@*`, `*!@*`, `[R1]`, `[R2]`, `[x2]`, `[x3]`) and two Boolean cells
(`[R1]`, `[x2]`) are all row 8341, an Fe-containing organometallic graph.
The 18 charged-polycycle refusals are not exact matches.

Reproduce from this source candidate:

```sh
cargo run -p chematic-smarts --example rdkit_parity_dump -- \
  validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  > /tmp/chematic-smarts-cage-dump.jsonl
python3 scripts/check_smarts_310k_profiles.py \
  --dump /tmp/chematic-smarts-cage-dump.jsonl \
  --output /tmp/chematic-smarts-cage-profile.json
```

The completed local dump SHA-256 was
`2c5b4858afe71082faf889e9652bec8c17ec50172fb461bc5367c89d442fe563`;
the checked JSON profile SHA-256 was
`2dd1994460ca17f916350c7ceffaa6d492be0c2becfe8bf9b222e9a509107404`.
The checker verified pinned input/oracle hashes, complete rows, unchanged
200/43 source-default cell identities and no opt-in regression. It does not
validate 3D force-field typing or published binding behavior.
