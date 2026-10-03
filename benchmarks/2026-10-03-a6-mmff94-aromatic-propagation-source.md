# MMFF94 aromatic-state propagation: source candidate (2026-10-03)

This is a source-only A6 typing check, not a published-package, geometry,
energy, or convergence result. Baseline: v1.0.31-source before this change;
oracle: RDKit 2026.03.6. The pinned 10,000-SMILES corpus has SHA-256
`f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`.

RDKit's MMFF aromaticity pass reads a neighbor's current aromatic flag after
earlier accepted rings have updated it. CheMatic read the frozen Kekule input
flag. The correction propagates the accepted-ring state for compact ring
systems. Cage systems with more than 20 symmetrized rings and macrocycles with
a ring of at least 16 atoms retain the prior behavior: their ring-model parity
has not been established. This is a conservative support boundary, not an
assertion that their existing types are all correct.

| Paired source comparison | Before | Candidate |
| --- | ---: | ---: |
| Atoms compared | 405,727 | 405,727 |
| Different numeric types vs RDKit | 2,964 | 459 |
| RDKit type 37 / CheMatic type 2 | 1,322 | 16 |
| Rows typed | 9,774 | 9,774 |
| CheMatic typing errors/refusals | 22 | 22 |

At exact input-row and atom-index granularity, 2,505 previously wrong types
became equal, zero previously equal types became unequal, and one atom changed
from one wrong type to another (row 8173, atom 8: RDKit 51, CheMatic 7→6).
The 459 residuals therefore remain open. The four source rows that regressed
under unrestricted propagation were three fullerene-like cages (3837, 3976,
4004) and one 16-membered macrocycle (7573); the conservative boundary
prevents those regressions on this corpus.

Reproduce in an environment with RDKit 2026.03.6:

```sh
python scripts/mmff94_atom_type_oracle_tsv.py \
  --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  --output /tmp/chematic-mmff94-oracle.tsv
SCHEMATIC_MMFF94_ORACLE_TSV=/tmp/chematic-mmff94-oracle.tsv \
  cargo test -p chematic-ff --test mmff_type_census -- --ignored --nocapture
```

The generated oracle TSV has SHA-256
`d77f7b5e8a0911da60674a0d1c63904176fe953fae06347384d04ae1c8144c10`.
Re-run the published wheel and independent A6 quality gates before any release
claim. No improvement in optimization or geometry is inferred from type parity.
