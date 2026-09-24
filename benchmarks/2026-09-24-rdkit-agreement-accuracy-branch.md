# RDKit output agreement — post-v1.0.22 accuracy branch

This record measures how often chematic's output equals RDKit 2026.03.6 for
RDKit-defined operations, on v1.0.22 (`main` at `297fec4d`) and on the
accuracy branch at `ea967d47`, and re-runs the Python per-operation timing
matrix on the branch to check that the accuracy fixes did not cost speed. It
is a source-level record on a shared 2-vCPU x86_64 VM, not a published-package
or release claim.

## Scope

| Item | Value |
|---|---|
| Base | `297fec4d2e932d50dbe1ed51c7eff9cec36d1686` (v1.0.22) |
| Branch | `ea967d476fa787bde64d001e752eb5cca0705b14` |
| Wheels | `maturin build --release` of `crates/chematic-py` at each revision (SHA-256 in the JSON) |
| RDKit | 2026.03.6 (PyPI) |
| Python / host | 3.11.15, Linux x86_64 VM, 2 vCPU "Intel(R) Xeon(R) Processor @ 2.10GHz" |
| Corpora | `scripts/chembl_accuracy_corpus_4999.smi` (aromatic SMILES), `scripts/nci_first_5k_smiles_only.smi` (mostly Kekulé SMILES), `scripts/descriptor_census_corpus.smi`; SHA-256 in the JSON |
| Agreement script | [`scripts/check_rdkit_output_agreement.py`](../scripts/check_rdkit_output_agreement.py) |
| Agreement raw data | [v1.0.22](2026-09-24-rdkit-agreement-v1.0.22-297fec4d.json), [branch](2026-09-24-rdkit-agreement-branch-ea967d47.json) |
| Timing raw data | [branch matrix](2026-09-24-python-op-matrix-vs-rdkit-branch-ea967d47.json) (script [`bench_python_op_matrix_vs_rdkit.py`](../scripts/bench_python_op_matrix_vs_rdkit.py)) |

## Row-level agreement with RDKit (rows equal / rows compared)

Each corpus row that both libraries parse is compared once. "new" marks an
API the base does not have. Bold marks an improvement.

| Operation | Criterion | ChEMBL 5k | NCI 5k (Kekulé SMILES) | Descriptor census 5k |
|---|---|---:|---:|---:|
| `morgan_r2_2048(rdkit-compatible)` | bit-identical | 5000 / 5000 | 4990 / 4991 | 5000 / 5000 |
| `atom_pair(rdkit-compatible)` | bit-identical | 4249 → **5000** / 5000 | 4074 → **4990** / 4991 | 4429 → **5000** / 5000 |
| `torsion(rdkit-compatible)` | bit-identical | 4248 → **4999** / 5000 | 4114 → **4991** / 4991 | 4428 → **4999** / 5000 |
| `maccs` | all 166 keys identical | 4714 → **4995** / 5000 | 1572 → **4989** / 4991 | 4597 → **5000** / 5000 |
| `logp` | abs. diff <= 1e-6 | 5000 / 5000 | 4439 → **4986** / 4991 | 5000 / 5000 |
| `mr` | abs. diff <= 0.01 | 5000 / 5000 | 4439 → **4986** / 4991 | 5000 / 5000 |
| `rdkit_tpsa` | abs. diff <= 0.1; N/O-only TPSA | new: 5000 / 5000 | new: 4990 / 4991 | new: 5000 / 5000 |
| `tpsa(S/P included)` | abs. diff <= 0.1 vs includeSandP=True | 5000 / 5000 | 4697 → **4990** / 4991 | 5000 / 5000 |
| `qed` | abs. diff <= 1e-3 | 3201 → **4962** / 5000 | 2594 → **4842** / 4991 | 3469 → **4988** / 5000 |
| `rdkit_mw` | abs. diff <= 1e-3 | 5000 / 5000 | 4991 / 4991 | 5000 / 5000 |
| `hbd` | exact | 5000 / 5000 | 4984 / 4991 | 5000 / 5000 |
| `rotatable_bonds` | exact | 5000 / 5000 | 4894 / 4991 | 4999 / 5000 |
| `murcko_scaffold` | same RDKit canonical SMILES | 1591 → **4998** / 5000 | 3466 → **4923** / 4991 | 1742 → **4998** / 5000 |
| `pains_passes` | vs RDKit PAINS FilterCatalog | 5000 / 5000 | 4801 → **4802** / 4991 | 4990 / 5000 |
| `brenk_passes` | vs RDKit BRENK FilterCatalog (chematic uses its own alert list) | 2550 → **2565** / 5000 | 3196 → 3182 / 4991 | 2617 → **2626** / 5000 |

### What changed

- **Atom pair / torsion (RDKit-compatible)**: RDKit's `numPiElectrons` returns
  0 for atoms it perceives as SP3; sulfonyl/sulfinyl S, phosphoryl P and
  seleninyl Se carry formal double bonds but are SP3 in RDKit. The exactly
  decidable part of RDKit's hybridization rule is now ported, and the
  fingerprints run on the RDKit-perceived aromatic view so Kekulé spellings
  agree.
- **MACCS**: key 26 uses RDKit's verbatim `[#6]=;@[#6](@*)@*` (the parser now
  accepts compound bond expressions), key 125 counts all-aromatic-bond rings
  of the symmetrized SSSR, key 1 is unset as in RDKit (`'?'`), and matching
  runs on the perceived aromatic view.
- **QED**: every property follows `rdkit.Chem.QED.properties` (QED acceptor
  SMARTS, N/O-only TPSA, the aliphatic-ring-deletion AROM count, all 116
  alerts verbatim with isotope-aware isotope alerts). The unspecified-bond fix
  below closed the remaining alert over-matches on ChEMBL.
- **Murcko scaffold**: linkers of any length (previously only single-atom
  linkers), exocyclic double-bonded atoms kept, stereo preserved by editing
  in place, chiral tags cleared on atoms that lose a substituent (RDKit
  behaviour).
- **Kekulé input**: Crippen hydrogen typing and TPSA nitrogen typing use the
  perceived environment.
- **SMARTS semantics**: an unspecified bond matches single or aromatic bonds
  (Daylight/RDKit), not any bond. LogP/MR, MACCS and PAINS agreement is
  unchanged by this; QED improves; the chematic-specific Brenk list moves by
  +15 / −14 / +9 rows on the three corpora (it is not RDKit's Brenk catalogue,
  so its agreement is shown for information only).
- **New `rdkit_tpsa`**: RDKit-default N/O-only TPSA; `tpsa` still includes S/P.

## Python timing matrix on the branch

Same protocol as the
[2026-09-24 perf-branch record](2026-09-24-python-op-matrix-vs-rdkit-perf-branch.md)
(cold inputs, 3 repeats, median, ChEMBL 5k). The comparison column is that
record's `8cc01365` run, whose crates differ from v1.0.22 only in version
strings; both runs were made on the same VM type in separate sessions, so
compare ratios, not absolute microseconds.

| | `8cc01365` | Branch `ea967d47` |
|---|---:|---:|
| Checked operations | 43 | 44 |
| Faster than RDKit (median) | 39 | 38 |
| …with full output agreement | 21 | 22 |
| …with full agreement and every repeat faster | 18 | 21 |

Cost of the accuracy fixes: the RDKit-compatible atom-pair and torsion
fingerprints and MACCS now run on the RDKit-perceived aromatic view, so on a
cold molecule they pay that perception (it is memoized per molecule and
shared with the descriptors). Their prepared-input ratios therefore fall —
atom pair 8.44x → 2.72x, torsion 27.03x → 4.20x, MACCS 1.41x → 1.25x — while
staying faster than RDKit and becoming bit-identical far more often. TPSA and
the new `rdkit_tpsa` remain slower than RDKit on prepared inputs (RDKit's
aromaticity is computed inside `MolFromSmiles`); the `parse+tpsa` row
includes parsing on both sides.

| Operation | Class | v1.0.22-equivalent `8cc01365` ratio (agreement) | Branch `ea967d47` ratio (agreement) |
|---|---|---:|---:|
| `parse_smiles` | same-task | 33.01x (—) | 28.81x (—) |
| `canonical_smiles(prepared)` | same-task | 1.45x (—) | 1.39x (—) |
| `parse+canonical_smiles` | same-task | 3.81x (—) | 3.64x (—) |
| `mol_block_write` | not-equivalent | 209.40x (—) | 293.88x (—) |
| `mol_block_read` | same-task | 3.42x (—) | 3.80x (—) |
| `inchi` | not-equivalent | 1.85x (0/1000) | 2.01x (0/1000) |
| `inchikey` | not-equivalent | 2.16x (0/1000) | 2.11x (0/1000) |
| `mw` | checked | 1.49x (5000/5000) | 1.54x (5000/5000) |
| `exact_mass` | checked | 1.68x (5000/5000) | 1.45x (5000/5000) |
| `logp` | checked | 3.06x (5000/5000) | 3.23x (5000/5000) |
| `mr` | checked | 2.86x (5000/5000) | 3.08x (5000/5000) |
| `tpsa` | checked | 0.96x (3659/5000) | 0.88x (3659/5000) |
| `rdkit_tpsa` | checked | new | 0.81x (5000/5000) |
| `hbd` | checked | 7.87x (5000/5000) | 7.99x (5000/5000) |
| `hba` | checked | 3.17x (3641/5000) | 2.55x (3641/5000) |
| `rotatable_bonds` | checked | 4.19x (5000/5000) | 3.10x (5000/5000) |
| `fsp3` | checked | 2.00x (5000/5000) | 1.98x (5000/5000) |
| `ring_count` | checked | 0.02x (4946/5000) | 0.03x (4946/5000) |
| `aromatic_ring_count` | checked | 0.14x (4940/5000) | 0.13x (4940/5000) |
| `qed` | checked | 2.36x (3201/5000) | 2.16x (4962/5000) |
| `kappa2` | checked | 21.26x (7/5000) | 16.12x (7/5000) |
| `chi1v` | checked | 3.70x (3380/5000) | 3.10x (3380/5000) |
| `bertz_ct` | checked | 296.38x (0/5000) | 356.51x (0/5000) |
| `formula` | checked | 1.57x (4992/5000) | 1.53x (4992/5000) |
| `num_stereocenters` | checked | 5.35x (4932/5000) | 5.89x (4932/5000) |
| `lipinski_bundle(mw,logp,hbd,hba)` | checked | 3.17x (3641/5000) | 3.14x (3641/5000) |
| `parse+logp` | checked | 3.98x (5000/5000) | 4.35x (5000/5000) |
| `parse+tpsa` | checked | 4.14x (3659/5000) | 4.46x (3659/5000) |
| `parse+ring_count` | checked | 5.60x (4946/5000) | 5.97x (4946/5000) |
| `parse+aromatic_ring_count` | checked | 3.79x (4940/5000) | 4.58x (4940/5000) |
| `morgan_r2_2048(native ecfp4)` | not-equivalent | 1.67x (—) | 1.83x (—) |
| `morgan_r2_2048(rdkit-compatible)` | checked | 0.97x (5000/5000) | 0.98x (5000/5000) |
| `maccs` | checked | 1.41x (4714/5000) | 1.25x (4995/5000) |
| `atom_pair(rdkit-compatible)` | checked | 8.44x (4249/5000) | 2.72x (5000/5000) |
| `torsion(rdkit-compatible)` | checked | 27.03x (4248/5000) | 4.20x (4999/5000) |
| `rdkit_fp(rdkit-compatible)` | checked | 1.18x (2000/2000) | 1.12x (2000/2000) |
| `pattern_fp(rdkit-compatible)` | checked | 1.32x (2000/2000) | 1.12x (2000/2000) |
| `tanimoto_1xN(per target, compatible Morgan)` | checked | 1.33x (50/50) | 0.97x (50/50) |
| `has_substructure c1ccccc1` | checked | 2.33x (2000/2000) | 2.74x (2000/2000) |
| `has_substructure [OH]` | checked | 4.06x (2000/2000) | 3.40x (2000/2000) |
| `has_substructure C(=O)N` | checked | 2.74x (1997/2000) | 2.88x (2000/2000) |
| `has_substructure [#7;R]` | checked | 1.10x (2000/2000) | 1.05x (2000/2000) |
| `has_substructure c1ccc2ccccc2c1` | checked | 1.17x (2000/2000) | 1.11x (2000/2000) |
| `has_substructure [CX3](=O)[OX2H1]` | checked | 2.27x (2000/2000) | 2.68x (2000/2000) |
| `has_substructure [NX3;H2,H1;!$(NC=O)]` | checked | 3.35x (2000/2000) | 4.52x (2000/2000) |
| `has_substructure *~*~*~*~*~*` | checked | 2.29x (2000/2000) | 2.58x (2000/2000) |
| `find_matches [#6]~[#7]` | checked | 7.30x (2000/2000) | 5.50x (2000/2000) |
| `murcko_scaffold` | checked | 3.53x (1591/5000) | 4.98x (4998/5000) |
| `add_hydrogens` | same-task | 4.72x (—) | 4.36x (—) |
| `remove_hydrogens` | same-task | 13.80x (—) | 13.76x (—) |
| `largest_fragment` | checked | 23.29x (2000/2000) | 24.26x (2000/2000) |
| `neutralize` | checked | 11.59x (2000/2000) | 7.78x (2000/2000) |
| `standardize vs Cleanup` | not-equivalent | 0.37x (—) | 0.45x (—) |
| `standardize vs Cleanup+Uncharge+Canonicalize` | checked | 21.16x (241/300) | 18.36x (241/300) |
| `canonical_tautomer` | checked | 17.05x (241/300) | 18.00x (241/300) |
| `brics_fragments` | not-equivalent | 70.61x (—) | 85.59x (—) |
| `cip_labels` | not-equivalent | 46.98x (—) | 53.51x (—) |
| `sssr_rings` | not-equivalent | 1.35x (—) | 1.63x (—) |
| `svg_depiction` | not-equivalent | 15.94x (—) | 17.33x (—) |
| `2d_layout` | not-equivalent | 11.97x (—) | 11.53x (—) |
| `embed_3d` | not-equivalent | 189.27x (—) | 230.85x (—) |
| `embed+minimize_mmff94` | not-equivalent | 0.10x (—) | 0.09x (—) |
| `mcs_pair` | not-equivalent | 2.24x (—) | 1.78x (—) |

## Limits

- One shared 2-vCPU VM; no confidence intervals; ratios near 1 are ties.
- Agreement is measured on three 5k corpora; residuals remain (for example
  hypervalent halogen oxoacids for Crippen, phosphazene bracket hydrogens for
  Murcko, RDKit-specific `CalcNumHBA`, which this record does not cover).
- Brenk and TPSA-with-S/P have no RDKit-default counterpart and are
  informational.
