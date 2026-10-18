# RDKit Cleanup: pyridinium phenoxide resonance forms

Date: 2026-10-18. Candidate source builds on `2ab1e555` with this branch's changes. The comparator is RDKit `2026.09.1` on the pinned 4,072-row stress corpus. The [row-level result](2026-10-18-rdkit-cleanup-pyridone-followup.json) records source hashes and every remaining difference.

The two parsers agree on all 4,072 inputs; RDKit Cleanup accepts 4,069. Exact Cleanup output agreement rises from **3,916 to 3,920/4,069 (96.34%)** against the previous source result. **Four rows improve and zero regress** on the identical scored set; **149** output differences remain. The fixed indices are 157, 178, 180, and 929.

These four inputs represent ortho or para pyridinium phenoxide resonance forms. Cleanup now converts the exocyclic O− single bond to carbonyl and neutralizes the ring N+ and O− charges. The rule recognizes both aromatic and Kekulé input. It is limited to a simple six-membered C5N aromatic ring with one qualifying exocyclic O−. Meta substitution, fused systems, and ambiguous multiple sites are left alone. A separate 15-input check against RDKit `2026.09.1` covered both positive and negative boundaries with 15 exact results.

Reproduce with a Python environment containing RDKit `2026.09.1`:

```sh
cargo build -p chematic-chem --example cleanup_rdkit_dump --offline
target/debug/examples/cleanup_rdkit_dump \
  < validation/cosmolkit_comparison/corpora/stress-v1.smi \
  > /tmp/chematic-cleanup-pyridone.jsonl
python scripts/rdkit_cleanup_source_gate.py \
  --dump /tmp/chematic-cleanup-pyridone.jsonl \
  --previous-source benchmarks/2026-10-18-rdkit-cleanup-nontetra-followup.json \
  --output /tmp/chematic-cleanup-pyridone-gate.json
```

Local checks: `chematic-chem` unit tests 981 passed and 1 ignored; Clippy with warnings denied and rustfmt passed. The worktree is sparse, so repository-wide benchmark-index validation cannot resolve unrelated omitted documentation and validation files; both new record links are checked directly. This result measures RDKit output agreement for a selected source corpus, not independent chemical truth or toolkit-wide accuracy superiority.
