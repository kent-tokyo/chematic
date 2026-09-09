# Symmetry-aware RMSD oracle — v1.0.10 candidate

The automorphism-aware conformer RMSD implementation was compared with an
independent RDKit `GetBestRMS` calculation on six deterministic conformer
pairs. The second conformer applies a fixed rigid transform and, where
applicable, swaps automorphism-equivalent atoms.

| Result | Count |
|---|---:|
| Total cases | 6 |
| Agreement within 0.001 Å | 5 |
| Documented symmetrization gap | 1 |
| Unexplained mismatches | 0 |

The one expected discrepancy is acetate: RDKit additionally applies
`symmetrizeConjugatedTerminalGroups`, which treats resonance-equivalent
terminal atoms as interchangeable. chematic does not yet apply that
representation-level preprocessing and therefore reports the explicit
topological result. It is retained as a known gap, not counted as parity.

## Reproduction

```text
cargo run --release -p chematic-3d --example rmsd_symmetric_oracle_dump --offline \
  > /private/tmp/rmsd_symmetric_oracle_dump-v1010.jsonl
python3 scripts/rmsd_symmetric_oracle_check.py \
  /private/tmp/rmsd_symmetric_oracle_dump-v1010.jsonl
```

This closes only the bounded RMSD oracle slice. Full conformer-quality,
force-field, and symmetry-aware TFD corpus gates remain open.
