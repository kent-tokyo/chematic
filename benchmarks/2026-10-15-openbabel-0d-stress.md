# Open Babel 0D V2000 molecular identity, stress corpus (2026-10-15)

This is a source and CLI comparison on one macOS arm64 host. Open Babel 3.2.1
writes a zero-coordinate V2000 MOL from each eligible source SMILES. RDKit
2025.09.3 canonicalizes the source, Open Babel's own reread SMILES, and the
SMILES emitted by two CheMatic CLI binaries. Both CheMatic binaries read the
same MOL bytes for each case. A case is scored only when Open Babel's own
reread retains the source isomeric identity.

The input is the checked-in [4,072-row stress corpus](../validation/cosmolkit_comparison/corpora/stress-v1.smi)
(SHA-256 `879166ab5481cac805ecd2d77f37d628383ec124085b90e2ee1659d55ce76277`).
The script selects the 1,313 unique SMILES containing `@`, `/`, or `\`;
three overlap the earlier 122-source CIP comparison. It keeps V2000 records
with atom parity or a directional bond marker. Open Babel MOL generation is
retried up to three times when it fails to produce a complete header. The
[raw JSON](2026-10-15-openbabel-0d-stress.json) records every scored input
hash, source and output identity, marker count, and excluded source with its
reason.

## Result

Of the 1,313 selected sources, 110 produced no complete Open Babel MOL after
three attempts, and 123 MOL records had neither parity nor directional bond
markers. Of the remaining 1,080 records, 229 failed the Open Babel self-reread
source gate. The 851 scored records include 809 with atom parity and 349 with
directional bonds; these categories overlap.

| CheMatic V2000 reader | Source identity retained |
| --- | ---: |
| Before the 0D reader changes | 63 / 851 |
| After the 0D reader changes | 850 / 851 |

The modified reader retains explicit-H parity, three-neighbor sulfur parity,
shared E/Z markers in conjugated systems, and CTfile default-valence hydrogen
counts for Na, K, Si, As, Se, and Sn. The single remaining mismatch concerns
the location of an aromatic-ring cation. Open Babel's MOL reread returns the
source aromatic-cation identity; CheMatic and a direct RDKit MOL read give the
same Kekulé cation location. The RDKit MOL read drops other stereo on this
zero-coordinate record, so this observation does not independently adjudicate
the complete molecular identity.

The paired results, input hashes, and every excluded source and reason are in
the raw JSON. A mismatch means the output SMILES does not canonicalize to the
source identity under RDKit; it does not by itself establish which toolkit has
the chemically correct interpretation of an ambiguous MOL record.

## Reproduction

Build CheMatic CLI binaries at the before and after source revisions, then
run with Open Babel 3.2.1 and RDKit 2025.09.3:

```bash
python3 scripts/compare_openbabel_0d_stereo.py \
  --corpus validation/cosmolkit_comparison/corpora/stress-v1.smi \
  --stereo-only \
  --before /path/to/chematic-before \
  --after /path/to/chematic-after \
  --report /tmp/openbabel-0d-stress.json
```

This is an exposed, existing stress corpus and a format-specific source
comparison. It does not prove independent absolute stereochemical truth,
published-package behavior, or toolkit-wide superiority over Open Babel.
