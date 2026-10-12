# RDKit workflow comparison: reactions, standardization, and R-groups

Date: 2026-10-12  
CheMatic source revision: `77d138adb7829e8c9f3dca128ee792ccb4e08a20`  
Comparator: RDKit `2026.03.6`  
Host: macOS arm64, Python 3.13.6

This is bounded source-candidate evidence. It does not establish complete
RDKit reaction, MolStandardize, or R-group parity. The raw record includes the
source revision, host, corpus hashes, every timing sample, mismatch examples,
and the exact timed-row selection.

## Method

- Accuracy is evaluated before timing. Only exact-output rows are timed.
- The engines alternate first position over 21 blocks. Intervals are
  deterministic bootstrap intervals over paired RDKit/CheMatic ratios.
- Parse-inclusive and preparsed/prepared work are separate lanes.
- The 256-row standardization slice is selected evenly across the frozen
  4,072-row unusual-SMILES corpus. Two rows whose parsed input semantics differ
  between engines are excluded before standardization is assessed.
- R-group output retains `R1`/`R2` labels. Unordered substituent sets do not
  count as agreement.

## Results

| Workflow | Accuracy boundary | Parse-inclusive RDKit/CheMatic | Prepared RDKit/CheMatic |
|---|---:|---:|---:|
| Reaction | 80/80 exact product sets; 3/3 jointly invalid | 2.17x (2.13–2.20) | 1.44x (1.43–1.47) |
| Largest-fragment fixture | CheMatic 39/39; RDKit 38/39 | 1.59x (1.58–1.65) | — |
| Labelled R-group | both 16/16 | 58.20x (57.02–59.68) | 60.04x (56.60–61.50) |

Every timing lane above favored CheMatic in all 21 blocks. The prepared
R-group lane uses the new reusable `PreparedRGroupCore`; the free functions
share the same implementation and output contract.

## Standardization profile

The table below reports operation agreement after restricting the sample to
254 rows with identical parsed RDKit-compatible input semantics. Speed is
measured only on those exact rows, so mismatches cannot become performance
wins.

| Operation | Exact output | Parse-inclusive speedup | Preparsed speedup |
|---|---:|---:|---:|
| Normalize | 254/254 | 6.23x | 5.49x |
| Reionize | 203/254 | 4.46x | 3.39x |
| FragmentParent | 233/254 | 4.60x | 5.00x |
| ChargeParent | 227/254 | 7.41x | 8.14x |
| IsotopeParent | 233/254 | 5.40x | 6.29x |
| StereoParent | 233/254 | 5.99x | 7.21x |
| Cleanup | 185/254 | 1.26x | 1.13x |

All 14 standardization timing lanes favored CheMatic in 21/21 blocks, and
every paired 95% interval stayed above 1.0. Normalize reached 254/254 after
neutral three-coordinate sulfoxides were changed from `S=O` to RDKit's
charge-separated `[S+]-[O-]` form.

The remaining differences are not presented as CheMatic accuracy wins. They
include metal coordination, unusual valence, proton-placement, and broader
Cleanup policy differences. Each must be independently classified before any
chemistry rule is changed.

## Validation

- `cargo test -p chematic-chem --lib`: 943 passed, 1 ignored.
- Python binding suite: 1,179 passed, 1 skipped.
- `cargo clippy -p chematic-chem -p chematic-py --all-targets -- -D warnings`:
  passed.
- Prepared and free-function R-group rows/columns are regression-tested as
  exactly equal.

## Reproduce

```bash
python scripts/bench_rdkit_workflows.py \
  --blocks 21 \
  --output benchmarks/2026-10-12-rdkit-workflows-expanded-source.json
```

See the [raw JSON](2026-10-12-rdkit-workflows-expanded-source.json) for full
samples and residual examples.
