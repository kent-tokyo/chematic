# RDKit comparison — chematic descriptor accuracy

> Historical validation snapshot. These measurements were generated for
> chematic **v0.4.20** and are not evidence for the current v1.0.3 release.
> For current scope and evidence, see
> [`docs/compatibility-scope.md`](../../docs/compatibility-scope.md),
> [`docs/validation.md`](../../docs/validation.md), and
> [`docs/benchmark.md`](../../docs/benchmark.md).

Reference tool: **RDKit 2026.03.3**  
Historical chematic version: **v0.4.20**
Corpus: `scripts/rdkit_ref_properties.tsv` (175 drug-like molecules)  
Full accuracy history: [benchmarks/](../../benchmarks/)

---

## Per-descriptor results (175-mol corpus)

| Descriptor | chematic API | Agreement | Tolerance | Notes |
|------------|-------------|-----------|-----------|-------|
| Molecular weight | `mol.mw` | **100%** | exact | |
| Heavy atom count | `mol.heavy_atom_count` | **100%** | exact | |
| HBD (Lipinski) | `mol.hbd` | **100%** | exact | |
| HBA (Ertl) | `mol.hba` | **100%** | exact | SMARTS-based |
| Aromatic ring count | `mol.aromatic_ring_count` | **100%** | exact | |
| TPSA | `mol.tpsa` | **100%** | ±0.1 Å² | |
| LogP (Crippen) | `mol.logp` | ~99% | ±0.3 | |

## Known divergence: TPSA on large corpora

On the 4,999-mol ChEMBL subset, TPSA agreement is **93.3%** (±0.1 Å²).

The 6.7% divergence is concentrated in:
- N-oxide groups (e.g., pyridine N-oxide)
- Sulfonamide patterns with unusual connectivity
- Certain charged nitrogen heterocycles

Both chematic and RDKit implement the Ertl (2000) SMARTS-based TPSA algorithm.
The divergence arises from edge-case differences in how SMARTS matching handles
polar surface contributions of unusual functional groups.

Ongoing calibration work is tracked in the CHANGELOG under `chematic-chem`.

## Known divergence: aromaticity model

chematic applies Hückel 4n+2 per SSSR ring independently.  
RDKit uses fused-ring electron delocalization.

Visible differences in: pyridone, quinolone, indolizine, and similar N-heterocycles.  
Impact on practical drug-like compounds: negligible (100% agreement on 175-mol corpus).

## Data files

| File | Contents | Rows |
|------|----------|------|
| `../../scripts/rdkit_ref_properties.tsv` | MW, LogP, TPSA, HAC, HBD, HBA per molecule | 175 |
| `../../scripts/rdkit_ref_tanimoto.tsv` | Tanimoto similarity matrix (subset) | 51 |
| `../../scripts/chematic_vs_rdkit.tsv` | Side-by-side comparison (rdXXX vs chXXX columns) | 175 |

## Reproduce

```bash
# Regenerate RDKit reference data
python scripts/gen_rdkit_reference.py   # writes scripts/rdkit_ref_properties.tsv

# Run comparison and print per-descriptor agreement
python scripts/rdkit_benchmark.py

# Side-by-side diff with mismatch details
python scripts/bench5k.py ~/Downloads/SMILES.csv --detail
```
