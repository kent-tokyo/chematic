# RDKit Cleanup: isolated alkali with one hydroxyl site

Date: 2026-10-18. Candidate source: `5fec6a07`.
The comparator is RDKit `2026.09.1` on the pinned 4,072-row stress corpus.
The [row-level result](2026-10-18-rdkit-cleanup-alkali-oxygen-followup.json)
records source hashes, the corpus hash, and every remaining difference.

Both parsers agree on all 4,072 inputs. RDKit Cleanup accepts 4,069 of them.
Exact Cleanup output agreement increases from **3,922 to 3,923/4,069
(96.41%)** against the previous source result. Row 2796 improves, zero rows
regress, and 146 output differences remain.

The previous isolated-metal rule charged a bare K atom in the same mixture as
a phenol but left the phenol neutral. The output gained a net positive formal
charge. RDKit pairs the new K cation with a phenoxide anion in that row. This
follow-up applies that charge pairing only to one originally isolated bare
Li/Na/K atom and one unambiguous neutral carbon-bound OH site in an otherwise
uncharged mixture. Multiple metal atoms, competing OH sites, explicit OH
hydrogen counts, and already charged mixtures remain outside this rule.
Focused Rust tests check alcohol, phenol, carboxylic acid, water, an existing
Na cation, competing OH sites, and Cleanup idempotence.

Reproduce with a Python environment containing RDKit `2026.09.1`:

```sh
cargo build -p chematic-chem --example cleanup_rdkit_dump --offline
target/debug/examples/cleanup_rdkit_dump \
  < validation/cosmolkit_comparison/corpora/stress-v1.smi \
  > /tmp/chematic-cleanup-alkali-oxygen.jsonl
python scripts/rdkit_cleanup_source_gate.py \
  --dump /tmp/chematic-cleanup-alkali-oxygen.jsonl \
  --previous-source benchmarks/2026-10-18-rdkit-cleanup-isolated-metals-followup.json \
  --output /tmp/chematic-cleanup-alkali-oxygen-gate.json
```

Local checks: `chematic-chem` unit tests 984 passed and 1 ignored; Clippy,
rustfmt, and refactor budget passed. This is output agreement on a selected
corpus, not independent chemical truth or broad accuracy superiority over
RDKit. Reaction and R-group source behavior is unchanged by this follow-up.
