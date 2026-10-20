# RDKit Cleanup: bounded acid/base and metal follow-up

Date: 2026-10-19. The source candidate starts from
`00596fb1ea5397524f98e3789c41f17ee8616fdc`; the
[row-level result](2026-10-19-rdkit-cleanup-acid-metal-followup.json) records
the changed source hashes. The comparator is RDKit `2026.09.1` on the pinned
4,072-row stress corpus.

The earlier 25 coordination-tag path differences were already closed by the
[non-tetrahedral follow-up](2026-10-18-rdkit-cleanup-nontetra-followup.md).
This follow-up improves exact Cleanup agreement from **3,923 to 3,933/4,069
(96.66%)**. Ten rows improve and zero regress against the immediately previous
source result.

After integration with the separate acid-anion change, the
[current-source record](../validation/results/rdkit-standardization-current-source-2026-10-20.json)
reports **3,931/4,069** exact rows. It fixes three metal rows and regresses
none against the integrated October 19 source record (3,928/4,069). The
earlier 3,933 figure describes this branch before integration.

| Class | Fixed rows | Remaining rows |
|---|---:|---:|
| Existing-anion acid/base transfer | 7 | 26 |
| Metal charge/valence representation | 3 | 101 across disconnection, charge, bond-order and stereo classes |
| Other representation/normalization | 0 | 9 |

The acid/base rule is intentionally narrower than the public heuristic
`reionize()`: it preserves total charge, requires an existing non-aromatic C/N
anion with an available proton, and moves that charge only to an unambiguous
neutral oxygen acid. Applying the generic Reionize stage directly after
Cleanup was rejected because it fixed none of this lane and regressed 845
previously exact rows. The metal changes are three RDKit-versioned
representations: permanganate, copper cyanide and all-dative neutral magnesium.

The 136 residual rows remain classified in the JSON: 131 are not implemented,
four are representation-only and one is a declared policy difference. This is
agreement with a pinned RDKit output contract, not an independent chemistry
truth or a toolkit-wide accuracy result.

Reproduce with RDKit `2026.09.1`:

```sh
cargo build -p chematic-chem --example cleanup_rdkit_dump --offline
target/debug/examples/cleanup_rdkit_dump \
  < validation/cosmolkit_comparison/corpora/stress-v1.smi \
  > /tmp/chematic-cleanup-candidate.jsonl
python scripts/rdkit_cleanup_source_gate.py \
  --dump /tmp/chematic-cleanup-candidate.jsonl \
  --previous-source benchmarks/2026-10-18-rdkit-cleanup-alkali-oxygen-followup.json \
  --output /tmp/chematic-cleanup-candidate-gate.json
```

Focused Cleanup tests pass. The full changed-crate and repository checks are
reported with the change rather than inferred from this comparison.
