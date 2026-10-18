# RDKit Cleanup: originally isolated metals

Date: 2026-10-18. Candidate source: `1a0bfab5fcc0cd786f68563b64e303a0208d8a32`.
The comparator is RDKit `2026.09.1` on the pinned 4,072-row stress corpus.
The [row-level result](2026-10-18-rdkit-cleanup-isolated-metals-followup.json)
records source hashes, the corpus hash, and every remaining difference.

Both parsers agree on all 4,072 inputs. RDKit Cleanup accepts 4,069 of them.
Exact Cleanup output agreement increases from **3,920 to 3,922/4,069
(96.39%)** against the previous source result. Two rows improve (indices 437
and 737), zero regress, and 147 output differences remain.

RDKit Cleanup charges a bare, originally isolated Li/Na/K atom to +1 and a
bare Mg/Ca atom to +2. The new rule follows that boundary. It does not change
explicit hydrides, bonded metals, Rb/Cs, or metals newly isolated by bond
disconnection. The focused Rust test also checks isotope and atom-map labels,
a separate stereocenter, and Cleanup idempotence against RDKit outputs.

Reproduce with a Python environment containing RDKit `2026.09.1`:

```sh
cargo build -p chematic-chem --example cleanup_rdkit_dump --offline
target/debug/examples/cleanup_rdkit_dump \
  < validation/cosmolkit_comparison/corpora/stress-v1.smi \
  > /tmp/chematic-cleanup-isolated-metals.jsonl
python scripts/rdkit_cleanup_source_gate.py \
  --dump /tmp/chematic-cleanup-isolated-metals.jsonl \
  --previous-source benchmarks/2026-10-18-rdkit-cleanup-pyridone-followup.json \
  --output /tmp/chematic-cleanup-isolated-metals-gate.json
```

Local checks: `chematic-chem` unit tests 982 passed and 1 ignored; Clippy,
rustfmt, refactor budget, and the benchmark index passed. This is output
agreement on a selected corpus, not independent chemical truth or broad
accuracy superiority over RDKit. Reaction and R-group source behavior is
unchanged by this follow-up.
