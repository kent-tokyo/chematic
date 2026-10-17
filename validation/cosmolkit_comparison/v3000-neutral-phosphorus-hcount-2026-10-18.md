# Neutral phosphorus V3000 `HCOUNT` interoperability, 2026-10-18

CheMatic's V3000 writer emitted `HCOUNT=1` for neutral `[PH]` even when
ordinary phosphorus valence already determines that hydrogen. COSMolKit
0.5.0 treats this property as a query feature and refuses to return a concrete
molecule. For `NC(N)=NCCCC(N)[PH](=O)O`, deleting that one redundant property
preserves the exact RDKit 2026.09.1 read-back identity and lets COSMolKit
read the same identity. The writer now omits it only for neutral, nonaromatic
phosphorus with ordinary bonds when CheMatic's inferred hydrogen count equals
the stored count. Non-inferred cases, including `[PH2](=O)O`, retain it.

## Same-input comparison

The baseline CLI was built from main `9d599051b6d3b73dcaa6b6f0272551b41d453857`
(SHA-256 `6b78a74ff20655e8a5d63f7d1b6d472b3d8681b9c98e626c452e5b9bcd4c278a`).
The candidate CLI SHA-256 is
`9cd35038d0637ce5f5321eba30b0b555b9ba16a60f85d9e3574648e27e812957`.
Both runs used `scripts/compare_rdkit_cosmolkit_v3000.py` (SHA-256
`fe323c8bd00a775a1e0e640814db90c5e9bcc401afce75d0e7cc3cafac64ca9c`),
RDKit 2026.09.1 and COSMolKit 0.5.0. The script gives all writers the same
RDKit-normalized SMILES. RDKit rereads the resulting V3000; the separate
interoperability column uses COSMolKit to read CheMatic's block, then RDKit
to normalize its SMILES. Exact canonical isomeric identity is required.

| Corpus | Comparable | CheMatic block → COSMolKit, before → after | CheMatic block → RDKit | RDKit writer → RDKit | COSMolKit writer → RDKit |
| --- | ---: | ---: | ---: | ---: | ---: |
| Browser 10k | 10,000 | 9,957 → **9,965** | 10,000 | 9,991 | 9,992 |
| `stress-v1` | 4,070 | 3,523 → **3,530** | 4,025 | 3,995 | 3,990 |
| ChEMBL 5k, candidate only | 5,000 | **5,000** | 5,000 | 4,981 | 4,993 |

Browser corpus SHA-256:
`f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`.
Stress corpus SHA-256:
`879166ab5481cac805ecd2d77f37d628383ec124085b90e2ee1659d55ce76277`.
ChEMBL corpus SHA-256:
`1c47371dcbe37f4e0a141bf545b72bf238de2761fa3894fa251a552d84728d3e`.
The baseline and candidate full reports were compared by source row. Every
expected identity and all three writer/read-back identities stayed the same.
Only the CheMatic-block/COSMolKit-reader identity changed: eight gains on
Browser lines 301, 324, 340, 7605, 7606, 7609, 7611 and 9610; seven gains
on stress lines 1499, 2411, 2856, 2871, 3286, 3325 and 3653. There were
no regressions on either full corpus. The Browser and stress reports can be
regenerated with the command below and their `rows` compared by `id`.

```sh
cargo build --offline -p chematic-cli --bin chematic
python scripts/compare_rdkit_cosmolkit_v3000.py \
  --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  --cli target/debug/chematic --report /tmp/browser-v3000.json
python scripts/compare_rdkit_cosmolkit_v3000.py \
  --corpus validation/cosmolkit_comparison/corpora/stress-v1.smi \
  --cli target/debug/chematic --report /tmp/stress-v3000.json
```

The writer regression test checks inferred and non-inferred phosphorus
hydrogens, including a charged control. `cargo test --offline -p chematic-mol
--lib --quiet` passed 814 tests; `cargo clippy --offline -p chematic-mol
--lib -- -D warnings`, `cargo fmt --all -- --check` and `git diff --check`
passed. These source-built, exposed-corpus results measure a V3000
interchange path with RDKit as input normalizer and read-back judge. They do
not establish independent chemical truth or toolkit-wide accuracy superiority.
