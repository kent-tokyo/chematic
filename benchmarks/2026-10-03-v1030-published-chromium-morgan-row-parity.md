# Published v1.0.30 browser Morgan bit parity

The [20-pair speed record](2026-10-03-v1030-rdkitjs-published-chromium-prepared-split20.md)
checks an aggregate fingerprint digest. This separate correctness gate compares
all 2,048 Morgan bits **for each molecule in Chromium**, using the same
published `@kent-tokyo/chematic` 1.0.30 and official `@rdkit/rdkit` 2026.03.6
JS/WASM files. The npm tarball SHA-256 is
`fd427a161c11b14e6cca51e78e0c35ff05bbc52faf0539a2ddc2b7090e534201`;
`scripts/check_v1030_published_browser_morgan_rows.py` verifies it along with
both packages' JS/WASM hashes and the input hashes. The comparison runs both
modules in one page and does not measure speed.

| Corpus and operation | Exact | Typed refusal | Wrong bit vector |
|---|---:|---:|---:|
| First 250 descriptor-census rows, direct | 250/250 | 0 | 0 |
| Same 250 rows, prepared | 250/250 | 0 | 0 |
| Browser 10k corpus, direct | 9,999/10,000 | 1 | 0 |
| Browser 10k corpus, prepared | 9,999/10,000 | 1 | 0 |

The 10k refusal is row 8341, a high-degree Fe coordination structure.
CheMatic returns `unsupported RDKit coordination sanitization`; RDKit.js
returns a fingerprint. It is an explicit non-parity case, not a correct
fingerprint. Both CheMatic entry points account for the same row. This is an
exposed corpus; it is not sealed validation or general fingerprint parity.

Raw records: [250 direct](2026-10-03-v1030-published-chromium-morgan-row-direct.json),
[250 prepared](2026-10-03-v1030-published-chromium-morgan-row-prepared.json),
[10k direct](2026-10-03-v1030-published-chromium-morgan-10k-direct.json), and
[10k prepared](2026-10-03-v1030-published-chromium-morgan-10k-prepared.json).
Their SHA-256 values and all counts are pinned by
`scripts/check_v1030_artifact_packet.py`. Re-run one lane with
`scripts/check_browser_rdkit_ecfp4_parity.py`, the archived corpus,
`--schematic-dir` pointing at the extracted pinned npm tarball,
`--schematic-tarball` naming that tarball, and `--schematic-operation direct`
or `prepared`; use `--allow-unsupported` only for the 10k corpus.

This closes the browser row-level output check for these two corpora and
configured Morgan operation. It does not prove canonical SMILES equality,
other operations, cross-browser speed, independent-host replication, or
library-level memory use. P0.2 remains open on those dimensions.
