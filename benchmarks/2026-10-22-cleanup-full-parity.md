# RDKit Cleanup parity on the pinned stress corpus

The Rust source at `8e4898c7d2db3a74a516d182a2259d87b06eaa28` and a
CPython 3.13 wheel built from that commit both match RDKit 2026.09.1 Cleanup
on **4,069/4,069 comparable rows**. The frozen `stress-v1.smi` corpus has
4,072 rows. All 4,072 inputs produce the same parsed canonical SMILES in both
toolkits; RDKit Cleanup raises on three rows (indices 427, 975, 1326), so
those rows have no comparison result.

The [Rust source report](../validation/results/rdkit-cleanup-source-full-parity-2026-10-22.json)
records the corpus and source hashes, every excluded row, and the comparison
with the previous 4,053/4,069 source result. All 16 previous mismatches are
fixed, with zero regressions on this corpus. The [Python wheel report](../validation/results/python-cleanup-full-parity-2026-10-22.json)
records the native extension hash and independently confirms zero Cleanup
mismatches. The wheel SHA-256 is
`e9dc30abef4aa856c49f8bda212d7eb61e5cb81b821a19cdce9ec10aac4c2d85`.

The remaining differences covered acid/base handling (three rows), canonical
representation (two), and metal charge, disconnection, or stereo handling
(eleven). The changes align Cleanup-specific normalization, metal
disconnection, acid/base matching, and the SMILES writer with RDKit's
in-memory Cleanup result. Standalone Normalize remains **4,069/4,069** on
the comparable rows; standalone Reionize is **4,064/4,069**.

Three RDKit Cleanup output SMILES (indices 1497, 1630, 2677) cannot be parsed
back by RDKit. CheMatic matches those output strings exactly. This is output
parity with this pinned RDKit version and corpus, not independent chemical
validation or a claim about unseen molecules.

Validation: `cargo test -p chematic-chem --lib --locked --offline` passed
(1,031 passed, one ignored); `cargo test -p chematic-smiles --lib --locked
--offline` passed (435 passed, three ignored); Clippy passed with `-D
warnings`; and the final wheel passed 74 Python molecule and pickle tests.
The source and wheel comparison commands were:

```sh
cargo run -p chematic-chem --release --offline --locked --example cleanup_rdkit_dump < validation/cosmolkit_comparison/corpora/stress-v1.smi > /private/tmp/cleanup.jsonl
PYTHONPATH=/private/tmp/chematic-rdkit-2026091/lib/python3.13/site-packages /usr/local/bin/python3 scripts/rdkit_cleanup_source_gate.py --dump /private/tmp/cleanup.jsonl --output /private/tmp/source-report.json
PYTHONPATH=/private/tmp/chematic-wheel-install-20261022:/private/tmp/chematic-rdkit-2026091/lib/python3.13/site-packages /usr/local/bin/python3 scripts/rdkit_cleanup_parity.py --sample-size 4072 --output /private/tmp/python-report.json
```
