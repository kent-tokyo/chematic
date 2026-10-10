# Direct ring and descriptor fast paths — 2026-10-10

This record evaluates commit `49fcfff7` against RDKit 2026.03.1 and
COSMolKit 0.5.0rc15. It is source-candidate evidence, not a published-package
or universal speed claim.

## Scope

- Artifact: CPython 3.13 macOS arm64 wheel built from the clean commit above,
  SHA-256
  `86da1280e3547833fe6a22c86929318f96e315df2fd04a645bd56735547c647f`.
- Corpus: the first 1,000 rows of `scripts/chembl_accuracy_corpus_4999.smi`.
- Timing: 21 process-isolated, rotating-order blocks with paired bootstrap
  95% intervals. [Raw JSON](2026-10-10-direct-descriptor-fastpaths-candidate.json).
- Accuracy: direct RDKit comparison on the complete 5,000-row ChEMBL corpus
  and 4,991 rows parsed by both engines from the NCI corpus.

## Changes

- SSSR ring count reuses the exact cycle rank already encoded by SMILES ring
  closures; unsupported bond orders and edited molecules retain the general
  graph path.
- RDKit-default TPSA computes only the N/O atom types it returns, in the same
  atom and floating-point accumulation order.
- Whole-molecule Labute ASA accumulates its total directly instead of
  allocating a second per-atom output vector. The per-atom API is unchanged.

## Result

The parse-inclusive pipelines remain faster than both comparators in all 21
blocks. The direct descriptor calls are smaller but still slower, so this
change does not close those prepared-operation roadmap gates.

| Operation | chematic pipeline ms | COSMolKit ms | RDKit ms | pipeline speedup vs COSMolKit | pipeline speedup vs RDKit |
|---|---:|---:|---:|---:|---:|
| TPSA | 6.888 | 122.865 | 69.171 | 17.99x [17.66, 18.06] | 10.09x [10.03, 10.17] |
| Labute ASA | 6.771 | 122.036 | 67.998 | 18.19x [17.98, 18.37] | 10.13x [10.07, 10.19] |

Direct-call paired speedup intervals remain below 1.0: TPSA is 0.32x
versus COSMolKit and 0.71x versus RDKit; Labute ASA is 0.19x and 0.43x.
These are retained as explicit deficits, not reported as wins.

## Accuracy gate

RDKit-compatible ring count, TPSA, non-chiral Morgan and chiral Morgan remain
exact on 5,000/5,000 ChEMBL rows and 4,991/4,991 comparable NCI rows. The full
row-level records are
[`direct-fastpaths-49fcfff7-chembl.json`](../validation/results/direct-fastpaths-49fcfff7-chembl.json)
and
[`direct-fastpaths-49fcfff7-nci.json`](../validation/results/direct-fastpaths-49fcfff7-nci.json).
The affected Rust library tests and clippy pass.

## Boundary

This patch reduces real work without changing outputs, but it does not prove
a complete operation-only speed win. Symmetric ring counts, descriptor
aromaticity, chiral Morgan and reaction application remain the next profiling
targets. Published artifacts and another host must be measured separately.
