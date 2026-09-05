# Validation Report

Summary of descriptor accuracy against RDKit on a ChEMBL-derived corpus.

The figures below are a dated measurement snapshot, not a claim that every
current workspace revision has been re-measured. Reproduce them with the
version, commit, corpus, and commands shown here.

**Environment:** Python 3.13.6, Apple M4, chematic v0.18.0 (commit `24a9239`), RDKit 2026.03.4, measured 2026-08-23T00:20:36Z

---

## Descriptor Accuracy (4,999-molecule ChEMBL subset)

| Descriptor | Agreement | Tolerance | Notes |
|---|---|---|---|
| Molecular weight | **99.82%** (4990/4999) | ±0.01 Da | vs `Descriptors.MolWt` |
| Heavy atom count | **100%** (4999/4999) | exact | |
| H-bond donors (HBD) | **100%** (4999/4999) | exact | |
| H-bond acceptors (HBA) | **100%** (4999/4999) | exact | |
| TPSA | **100%** (4999/4999) | ±0.1 Å² | |
| LogP (Crippen) | **100%** (4999/4999) | exact* | max Δ = 1.10e-13 |
| MR (molar refractivity) | **100%** (4999/4999) | ±0.01 | |
| Fsp3 | **100%** (4999/4999) | ±0.001 | |
| Aromatic ring count | **100%** (4999/4999) | exact | |
| Aliphatic ring count | **100%** (4999/4999) | exact | |
| Saturated ring count | **100%** (4999/4999) | exact | |
| Rotatable bonds | **100%** (4999/4999) | exact | |
| Num heteroatoms | **100%** (4999/4999) | exact | |
| Num spiro atoms | **100%** (4999/4999) | exact | |
| Num bridgehead atoms | **100%** (4999/4999) | exact | bond-intersection algorithm |
| Num amide bonds | **100%** (4999/4999) | exact | |
| Arom./aliph. heterocycles | **100%** (4999/4999) | exact | |
| [nH] SMARTS match | **100%** (4999/4999) | precision & recall | TP=467 TN=4532 FP=0 FN=0 |
| Num stereocenters (legacy)  | **99.96%** (4997/4999) | exact† | vs `CalcNumAtomStereoCenters` |
| Num stereocenters (new CIP) | 98.6% (4929/4999) | exact† | vs `FindPotentialStereo` |

20 of 20 tested metrics reach ≥98.6% on the 4,999-molecule ChEMBL corpus.
chematic stereocenters is calibrated between legacy (99.96%) and new-CIP (98.6%) oracles.

---

## Stereocenters — Oracle Calibration

chematic's stereocenter count is calibrated between two RDKit oracles:

| Oracle | Agreement | Count | Notes |
|---|---|---|---|
| Legacy `CalcNumAtomStereoCenters` | **99.96%** (4997/4999) | 4997/4999 | 68 molecule where chematic is more accurate (legacy under-counts) |
| New CIP `FindPotentialStereo` | 98.6% (4929/4999) | 4929/4999 | 0 molecules where chematic correctly agrees with legacy (new CIP over-counts cage systems) |
| Consensus (all three agree) | 98.6% (4929/4999) | 4929/4999 | molecules where legacy, new CIP, and chematic all agree |

**Oracle disagreements:** 68 molecules where legacy ≠ new CIP.
- 68 where legacy under-counts a pseudoasymmetric polyester (chematic and new CIP both correctly return 4; legacy returns 2)
- 0 where new CIP over-counts cage/adamantane-like systems (chematic and legacy correctly agree on fewer stereocenters)

---

## CIP R/S/E/Z Label Agreement

A distinct metric from stereocenter *count* agreement above: given a stereocenter
both chematic and RDKit agree exists, does chematic assign the same R/S/E/Z label?
Measured via `chematic-cip`'s `corpus_snapshot` example
(`assign_cip_accurate_experimental`, the production-path engine) against the same
4,999-molecule corpus, cross-checked against a freshly regenerated `rdCIPLabeler`
oracle by `scripts/cip_accurate_full_corpus_report.py`:

| Oracle | Agreement | Count |
|---|---|---|
| Modern `rdCIPLabeler` | **99.64%** | 4171/4186; 15 phosphorus rows fail closed as representation-unstable |

A supplementary spot-check against RDKit's older `AssignStereochemistry`/`_CIPCode`
algorithm (R/S atom stereocenters only, not the E/Z bond stereocenters `rdCIPLabeler`
also covers above) gives a consistent **99.78%** (4150/4159), confirming the two
RDKit-side oracles agree with each other and with chematic to within a similar
margin. Both figures are a substantial improvement over a prior snapshot's
96.30%/96.83%, reflecting CIP-engine fixes landed in the interim releases.

```bash
cargo run -p chematic-cip --release --example corpus_snapshot -- \
    --candidate scripts/chembl_accuracy_corpus_4999.smi /tmp/candidate.tsv
.venv/bin/python scripts/cip_accurate_full_corpus_report.py \
    /tmp/candidate.tsv /tmp/candidate.tsv scripts/chembl_accuracy_corpus_4999.smi
```

(Passing the same file as both `baseline` and `candidate` reports zero
regressions trivially and surfaces the direct oracle-agreement percentage in the
`candidate_correct` line -- this script's real purpose is comparing two engine
variants, reused here for a single-engine accuracy figure.)

---

## Reproduce

```bash
# Requires RDKit; corpus is committed at scripts/chembl_accuracy_corpus_4999.smi
.venv/bin/python scripts/bench5k.py scripts/chembl_accuracy_corpus_4999.smi
.venv/bin/python scripts/bench5k.py scripts/chembl_accuracy_corpus_4999.smi --detail
.venv/bin/python scripts/bench5k.py scripts/chembl_accuracy_corpus_4999.smi --json validation/results/bench5k_latest.json
python3 scripts/gen_validation_report.py validation/results/bench5k_latest.json
```

Reference TSV files: `scripts/rdkit_reference_*.tsv` (generated by `scripts/gen_rdkit_reference.py`).

---

*\* LogP max |Δ| = 1.10e-13 — within float64 rounding error. bench5k.py uses ±0.01 as the test threshold.*

*† Stereocenters: see Oracle Calibration section above.*

---

## Known Limitations

- **Kekulization**: 1 of 5,000 tested molecules — `[H][H]` (no heavy atoms; IUPAC InChI library constraint). Returns `KekuleError` explicitly.
- **Aromaticity model**: Hückel 4n+2 per SSSR ring; RDKit uses fused-ring delocalization. Visible in pyridone, quinolone, indolizine.
- **InChI**: Pure-Rust implementation is approximate. Use `native-inchi` feature for standard-compliant InChI/InChIKey.

---

*Validation corpus: ChEMBL-derived 4,999-molecule SMILES set (`scripts/chembl_accuracy_corpus_4999.smi`, 5,000 raw lines; `csv.DictReader` treats the first line as a header, so 4,999 molecules are actually evaluated -- see the Reproduce section). Details: [`benchmark.md`](benchmark.md) · [`rdkit-comparison.md`](rdkit-comparison.md)*
