# Open Babel MOL stereo preservation: paired source check (2026-10-15)

CheMatic source revision `06da70001b47c22ecb4803cbc0957666f230ce83`
(v1.1.4 source), Open Babel 3.2.1, and RDKit 2025.09.3 were run on one
macOS arm64 host. The [raw record](2026-10-15-openbabel-stereo-paired.json)
retains all 310 paired rows, input hashes, output SMILES, tool versions and
binary hashes.

## Method

The fixed 155-row [CIP label corpus](../validation/cip_label_corpus.jsonl)
has SHA-256 `14ce78b96a24bc66b65fd963ffc7e236e1153f72c006933555335d0bdd61730e`.
RDKit writes each molecule as a 2D V2000 and V3000 MOL block. A row is compared
only if RDKit rereads that exact block to the source's isomeric semantic
identity. CheMatic and Open Babel then read the same file bytes and emit SMILES;
RDKit checks each output against that source identity. Canonical SMILES spelling
differences alone do not count as mismatches. No source rows were excluded.

| Profile | Both match source | CheMatic only matches | Open Babel only matches | Neither matches |
| --- | ---: | ---: | ---: | ---: |
| MOL V2000 | 143 | 12 | 0 | 0 |
| MOL V3000 | 146 | 9 | 0 | 0 |

The 155 corpus rows contain 122 distinct RDKit semantic identities. The
CheMatic-only rows represent four distinct identities in V2000 and three in
V3000; several corpus lines repeat a structure. Every disagreement is
classified as a stereo descriptor assignment mismatch after the same RDKit
source/reopen contract. This is a preservation result under that contract,
not an independent adjudication of absolute configuration.

Reproduce with the pinned tool versions and a CheMatic CLI built from the
recorded source revision:

```bash
python3 scripts/compare_openbabel_stereo_formats.py \
  --cli /path/to/chematic \
  --obabel /path/to/obabel \
  --output /tmp/openbabel-stereo-paired.json
```

The result is limited to one host, the stated 155-row corpus, RDKit-generated
2D MOL input, and CLI MOL-to-SMILES conversion. It is not a published-package,
cross-platform, broad-format, performance, or toolkit-wide superiority claim.
