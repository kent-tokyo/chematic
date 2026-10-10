# RDKit/COSMolKit performance and accuracy gate start

Date: 2026-10-10. This record defines the starting point for a bounded attempt
to beat the pinned distributed comparators without reducing correctness. It is
not a claim of universal superiority.

## Scope

- Candidate: local source candidate based on `origin/main`; final clean commit
  and wheel hashes are recorded in the adjacent JSON files.
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
| parse | 18.082x (17.984–18.156) | 18.300x (18.144–18.483) | gate passed |
| TPSA | 0.469x (0.458–0.479) | 1.014x (1.003–1.036) | RDKit only |
| Labute ASA | 0.266x (0.262–0.270) | 0.522x (0.511–0.532) | open |
| ring count | 0.049x (0.048–0.050) | 0.033x (0.033–0.035) | open |
| chiral Morgan | 0.638x (0.633–0.640) | 0.966x (0.958–0.971) | open |
| amide/amine reaction | 0.768x (0.762–0.770) | 1.194x (1.189–1.198) | RDKit only |

“Passed” requires all 21 blocks to win and the paired 95% lower bound to be
above 1.0. The current candidate therefore does **not** satisfy the complete
prepared-operation target.

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
