# Opt-in SMARTS organometallic ring-view follow-up

Source-only comparison on the pinned RDKit **2026.03.6** oracle and the exposed
10,000-molecule × 31-query corpus. The candidate is based on the dense-cage
ring fix merged as PR #726. This is not a registry-package or binding result.

The last wrong-confident row (8341) has a carbon anion single-bonded to Fe in
the input SMILES. RDKit sanitization changes that hypervalent nonmetal-to-metal
single bond to dative before ring perception ([pinned RDKit implementation](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/Code/GraphMol/MolOps.cpp),
[sanitization contract](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/Docs/Book/RDKit_Book.rst)). The source-only opt-in matcher
now makes a private bond-order view for the unambiguous `[C-]`--Fe valence-four
case. Ring-bond and `[xN]` queries exclude that dative bond, even if both
endpoints also occur in another ring. The caller's molecule, native SMARTS,
and general SMILES parsing are unchanged. This deliberately does not claim
RDKit's full organometallic cleanup, including multiple-metal rank selection.

| Pinned 310,000-cell lane | PR #726 source | Candidate source |
|---|---:|---:|
| Published/default differing match sets | 200 | 200 |
| Opt-in corrected match sets | 177 | **183** |
| Opt-in wrong-confident match sets | 6 | **0** |
| Opt-in corrected hit/no-hit results | 41 | **43** |
| Opt-in wrong-confident hit/no-hit results | 2 | **0** |
| Opt-in typed-refusal cells | 18 | 18 |

Of the 200 default-wrong match sets, 17 are now a typed refusal rather than
an exact result. One further refused cell was default-correct. No new
wrong-confident cell appeared. The 18 charged-polycycle refusals still need
chemical adjudication; they do **not** count as matches. No published Rust,
Python, npm or WASM artifact has been rerun for this source candidate.

Reproduce from this source candidate:

```sh
cargo run -p chematic-smarts --release --example rdkit_parity_dump -- \
  validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  > /tmp/chematic-smarts-organometallic-dump.jsonl
python3 scripts/check_smarts_310k_profiles.py \
  --dump /tmp/chematic-smarts-organometallic-dump.jsonl \
  --output /tmp/chematic-smarts-organometallic-profile.json
```

The completed local dump SHA-256 is
`8713c456adaec21f0ff60441152a593453612e5f24210d21a5ccd2ca06ae9b91`;
the checked profile SHA-256 is
`7627beefc27068c34e985e1e8b7ec8eb5dc1764929350149af072719a0bafb87`.
The checker verified fixed input/oracle hashes, all rows, unchanged 200/43
default cell identities, and no opt-in regression. The SMARTS crate test suite
passed. This packet does not address the reaction published-artifact rerun or
the separate A6 3D quality exits.
