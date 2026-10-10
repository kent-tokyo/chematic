# RDKit/COSMolKit performance and accuracy gate start

Date: 2026-10-10. This record defines the starting point for a bounded attempt
to beat the pinned distributed comparators without reducing correctness. It is
not a claim of universal superiority.

## Scope

- Candidate: clean commit `0972ca64d6781ebb96a61327763dce712599dc4a`.
  macOS arm64 CPython 3.13 wheel SHA-256:
  `afab8432c8ede45c96cbc739c86fd690218f27caa3c9256968a3f0d5b813eae4`.
- Comparators: RDKit 2026.09.1 and COSMolKit 0.5.0rc22.
- Host: Apple silicon macOS, CPython 3.13.6.
- Accuracy: ChEMBL 5,000 and NCI 5,000 input rows, using the repository's
  RDKit-agreement script.
- Performance: first 1,000 ChEMBL rows, 21 alternating-order blocks. Parsing,
  prepared operation time, and parse-inclusive pipeline time are separate.
- A speed ratio is invalid when either engine reports a row error.

## Accuracy preservation

The candidate retained exact RDKit agreement for Morgan radius 2 (plain and
chiral), atom pair, torsion, MACCS, LogP, MR, both TPSA definitions, QED,
molecular weight, HBD, rotatable bonds, and ring count on every jointly parsed
row: 5,000/5,000 ChEMBL and 4,991/4,991 NCI. Murcko, PAINS, and Brenk retain
their existing documented residuals; this work does not claim to close them.

## 21-block performance result

All six parse-inclusive pipelines beat both comparators in all 21 blocks. The
prepared-operation result is intentionally less favorable:

| Prepared operation | vs COSMolKit median (95% CI) | vs RDKit median (95% CI) | Result |
|---|---:|---:|---|
| parse | 18.120x (17.921–18.167) | 18.348x (18.239–18.525) | gate passed |
| TPSA | 0.453x (0.446–0.470) | 0.999x (0.977–1.020) | open |
| Labute ASA | 0.268x (0.263–0.275) | 0.524x (0.521–0.531) | open |
| ring count | 0.050x (0.049–0.051) | 0.035x (0.033–0.038) | open |
| chiral Morgan | 0.636x (0.633–0.644) | 0.967x (0.953–0.974) | open |
| amide/amine reaction | 0.765x (0.762–0.769) | 1.189x (1.182–1.198) | RDKit only |

“Passed” requires all 21 blocks to win and the paired 95% lower bound to be
above 1.0. The current candidate therefore does **not** satisfy the complete
prepared-operation target.

Machine-readable evidence:

- [`2026-10-10-rdkit-cosmolkit-championship-accuracy.json`](2026-10-10-rdkit-cosmolkit-championship-accuracy.json)
- [`2026-10-10-rdkit-cosmolkit-championship-performance.json`](2026-10-10-rdkit-cosmolkit-championship-performance.json)

## Changes measured

- Compute and memoize the shared RDKit-model safety predicate once during
  validated SMILES parsing.
- Avoid a second cleanup scan when that predicate proves an ordinary molecule.
- Preserve RDKit bond accumulation order for Labute ASA without allocating a
  complete bond-order vector, and keep common per-atom scratch data inline.
- Update the COSMolKit rc22 adapter. The previous method name caused every
  chiral-Morgan call to fail; the benchmark now suppresses speed ratios for any
  operation with engine errors.

## Next order

1. cache an exact, mutation-safe prepared ring model rather than recomputing it;
2. reduce Labute/TPSA Python-call and RDKit-view overhead without moving costly
   cyclic perception into parsing;
3. profile chiral Morgan and reaction allocation against rc22;
4. rerun a clean source commit, public artifacts, a second host, and peak memory.
