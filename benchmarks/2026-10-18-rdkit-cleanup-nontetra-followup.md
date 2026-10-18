# RDKit Cleanup: isolated coordination tags

Date: 2026-10-18. Candidate source built from `d8c4475a` with the changes in this branch. The comparator is RDKit `2026.09.1` on the pinned 4,072-row stress corpus. The [row-level result](2026-10-18-rdkit-cleanup-nontetra-followup.json) includes source hashes and exact mismatches.

The two parsers agree on all 4,072 inputs; RDKit Cleanup accepts 4,069. CheMatic's exact Cleanup output rises from **3,891 to 3,916/4,069 (96.24%)** against the prior source result. **25 rows are fixed and zero regress** on the identical scored set. **153** output differences remain.

The fixed rows are fully occupied non-tetrahedral metal complexes. Cleanup already retains their coordination tags in the graph after metal disconnection, but the general RDKit-compatible SMILES writer sanitizes an isolated metal as though its SMILES had been parsed anew, dropping the tag. The Cleanup-specific writer now retains the tag for isolated centers with no vacant coordination sites. The ordinary writer keeps its existing parse-and-sanitize contract.

One three-ligand square-planar case remains different: RDKit writes `[Pt@SP2+3]` after disconnection, while this writer emits `[Pt+3]`. Reusing its original `@SP1` tag would assert a different configuration; the vacant-site permutation needs a separate remap. The case stays in the mismatch list.

Reproduce with a Python environment containing RDKit `2026.09.1`:

```sh
cargo build -p chematic-chem --example cleanup_rdkit_dump --offline
target/debug/examples/cleanup_rdkit_dump \
  < validation/cosmolkit_comparison/corpora/stress-v1.smi \
  > /tmp/chematic-cleanup-nontetra.jsonl
python scripts/rdkit_cleanup_source_gate.py \
  --dump /tmp/chematic-cleanup-nontetra.jsonl \
  --previous-source benchmarks/2026-10-18-rdkit-cleanup-n-oxide-followup.json \
  --output /tmp/chematic-cleanup-nontetra-gate.json
```

This is RDKit agreement for a selected source corpus, not an independent chemistry truth or a toolkit-wide accuracy comparison. Reaction and R-group source code are unchanged in this follow-up; their separate results are in the [workflow packet](2026-10-18-rdkit-workflows-followup.md).
