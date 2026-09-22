# MMFF94 current-source same-coordinate energy — 2026-09-23

This diagnostic reruns the fixed 265-molecule Tier A+B manifests after the
#337 typing and charge corrections. RDKit 2026.03.6 generates one explicit-H
coordinate set per row; pinned RDKit and source-built CheMatic 1.0.19 then
evaluate those identical coordinates. It measures total MMFF94 energy
agreement only—not conformer quality, minimization, stereo, or speed.

## Before and after

The comparison uses the same seed, manifests, 262 comparable rows, two RDKit
embedding failures, and one declared unsupported probe as the latest v1.0.13
diagnostic.

| Metric | Earlier v1.0.13 diagnostic | Current source | Change |
|---|---:|---:|---:|
| Median absolute delta | 0.250084 kcal/mol | 0.217545 kcal/mol | 13.0% lower |
| p90 absolute delta | 9.239645 kcal/mol | 1.144679 kcal/mol | 87.6% lower |
| Maximum absolute delta | 32.670488 kcal/mol | 9.874153 kcal/mol | 69.8% lower |
| Rows within 1 kcal/mol | 193/262 | 232/262 | +39 |
| Rows within 5 kcal/mol | 213/262 | 260/262 | +47 |
| Rows above 5 kcal/mol | 49/262 | 2/262 | -47 |

The two remaining total-energy residuals above 5 kcal/mol are rows 166
(`chembl_tier_b_0101`, 9.874153 kcal/mol) and 231
(`chembl_tier_b_0166`, 6.505362 kcal/mol). Each successful row preserves the
CheMatic bond, angle, stretch-bend, torsion, out-of-plane, van der Waals, and
electrostatic breakdown. RDKit's Python force-field API does not expose an
equivalent per-term oracle in this packet, so these two rows remain open for
term-level adjudication rather than being assigned to a specific term.

## Provenance and boundary

- CheMatic source revision: `8ca05c8b2636e5e1e52ef2c37c3c6c2ce5cff41a`
- RDKit wheel: `rdkit==2026.3.6`, SHA-256
  `e16c467cb254a223e59a0cf81358c6b39da15a99d2909170d693e95778fddb41`
- CheMatic extension SHA-256:
  `e340782b114b23bcfa7365f85182d504fbd9a1e152edcd574de74914614cd6dd`
- Row accounting: 265/265 terminal; 262 comparable, two embed failures, one
  declared unsupported; no rows dropped.

This is current-source evidence, not a published-package result or an A6 exit.
Convergence/timeout, gradient, stereo, conformer quality, and speed remain
separate gates in issue #637.

Evidence:

- `validation/results/mmff94-same-explicit-h-energy-current-main-v1.0.19-2026-09-23.json`
- `validation/results/mmff94-same-explicit-h-energy-current-main-v1.0.19-2026-09-23.jsonl`
