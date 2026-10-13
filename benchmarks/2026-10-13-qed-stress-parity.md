# QED unusual-SMILES parity follow-up

This record measures the source candidate at commit `3f9b5cd5` against
RDKit 2026.03.6 on the frozen 4,072-row unusual-SMILES corpus. It is a
reference-compatibility result, not an experimental-property-accuracy or
published-package claim.

## Result

| Metric | Public v1.1.2 | Source candidate |
|---|---:|---:|
| HBD agreement | 4,072 / 4,072 | 4,072 / 4,072 |
| QED agreement on comparable rows | 4,021 / 4,071 | **4,071 / 4,071** |
| QED disagreements | 50 | **0** |
| RDKit QED errors | 1 | 1 |

RDKit's non-comparable row is corpus row 1,327. Its QED calculation raises a
`KekulizeException`; it is not counted as a CheMatic disagreement.

## Root causes closed

1. QED component properties were split across the original graph and an
   RDKit-parity aromatic view. All eight properties now use one reconstructed
   RDKit-sanitized graph.
2. RDKit suppresses neutral divalent O/S two-electron donation when the first
   selected ring containing the atom has at least nine atoms. The parity
   aromaticity path now applies the same ring-local rule.
3. RDKit assigns zero mass to a dummy atom but still counts its bracket
   hydrogens. `rdkit_molecular_weight` now handles `[*H]` and `[*H2]`
   accordingly.

## Reproduction and provenance

- corpus: `validation/cosmolkit_comparison/corpora/stress-v1.smi`
- corpus SHA-256:
  `879166ab5481cac805ecd2d77f37d628383ec124085b90e2ee1659d55ce76277`
- local source wheel SHA-256:
  `01d752e4157d60fc14a13f2678e8cc9460ca9519f186474eaa370d243f4332a4`
- Python 3.13.6, macOS arm64
- detailed result:
  [`validation/results/qed-stress-rdkit-source-2026-10-13.json`](../validation/results/qed-stress-rdkit-source-2026-10-13.json)

The dump and comparator were run with:

```bash
python scripts/chematic_chemistry_dump.py \
  --corpus validation/cosmolkit_comparison/corpora/stress-v1.smi \
  --output /tmp/qed-source-final-dump.jsonl \
  --descriptors-only
python scripts/compare_chemistry_dump_rdkit.py \
  --dump /tmp/qed-source-final-dump.jsonl \
  --output /tmp/qed-source-final-vs-rdkit.json \
  --descriptors-only
```

Rust, Python and WASM regression tests cover the three residual families.
The source candidate also passed the complete `chematic-perception` library
tests, the `chematic-wasm` native tests, focused Python descriptor tests and
Clippy for all changed crates.
