# RDKit workflow comparison: reactions, standardization, and R-groups

Date: 2026-10-12  
CheMatic source revision: `8b529ce6`  
CheMatic wheel SHA-256: `3d945852bc57eac2dde418febd98c7c629db21446a8d5cd962106a454c4ed1ea`  
Comparator: RDKit `2026.03.6`  
Host: macOS arm64, Python 3.13.6

This is a bounded source-candidate comparison. It does not establish complete
RDKit reaction, standardization, or `RGroupDecomposition` parity.

## Method

- Accuracy is evaluated before timing.
- Timing includes only rows where both engines return the same declared
  semantic result.
- Each lane uses 31 alternating-order blocks. The reported interval is a
  deterministic bootstrap interval over paired RDKit/CheMatic ratios.
- Reaction timing is split into parse-inclusive and prepared-input lanes.
- The R-group gate compares labelled `R1`/`R2` payloads. Sorting payloads while
  discarding labels is not accepted as parity.

Corpora and their SHA-256 hashes are recorded in the accompanying
[raw JSON](2026-10-12-rdkit-workflows-source.json).

## Results

| Lane | Accuracy boundary | RDKit/CheMatic median | Paired 95% interval | Winning blocks |
|---|---:|---:|---:|---:|
| Reaction, parse-inclusive | 80/80 semantic matches; 3/3 jointly invalid | 2.19x | 2.15–2.22x | 31/31 |
| Reaction, prepared | same 80 valid rows | 1.44x | 1.43–1.46x | 31/31 |
| Largest-fragment standardization | CheMatic 39/39; RDKit 38/39 | 1.59x on 38 common-exact rows | 1.59–1.65x | 31/31 |
| Labelled R-group decomposition | both 16/16 | 57.27x | 56.31–57.92x | 31/31 |

The reaction speed result is not an artifact of charging RDKit for template
construction while CheMatic uses its template cache: the prepared lane keeps
reactants parsed and RDKit reactions compiled before timing, and still favors
CheMatic in every block.

## R-group correction

The prior gate compared each row as an unordered set of substituents. That
missed a real error on a symmetric two-site phenyl core: CheMatic returned the
right two fragments but assigned ethyl and hydroxyl to the opposite `R1`/`R2`
columns. The stricter labelled gate reproduces the old wheel at 15/16 and the
source candidate at 16/16.

The implementation now:

1. preserves exact previous-column matches as the strongest alignment signal;
2. uses ECFP4 similarity to established column fragments to break remaining
   symmetric ties, matching RDKit's bounded row alignment in the fixture; and
3. matches the original molecule graph directly instead of serializing and
   reparsing a canonical SMILES for every row.

Rust regression tests cover the corrected row and atom-order-equivalent input
spellings. The full `chematic-chem` suite, its corpus standardization tests,
clippy, and the Python labelled R-group tests pass.

## Boundaries and next gates

- R-group coverage remains one mapped core with terminal wildcard or mapped
  core-atom labels. Multi-core/MCS alignment, tautomer/core enumeration,
  non-terminal groups, and multiple groups per site remain unsupported.
- The 39-row standardization lane covers largest-fragment selection, not the
  complete cleanup/normalization/tautomer pipeline.
- Reaction accuracy is the existing 83-row checked gate, not a broad reaction
  database. Atom provenance and template-map gates remain separate.
- These are source-built macOS results. Published-wheel and second-host gates
  are required before release-level performance claims.

Reproduce with:

```bash
python scripts/bench_rdkit_workflows.py \
  --blocks 31 \
  --output benchmarks/2026-10-12-rdkit-workflows-source.json
```
