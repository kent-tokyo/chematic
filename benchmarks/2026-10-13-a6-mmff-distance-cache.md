# A6 MMFF distance-cache development record

Date: 2026-10-13

This is source-candidate development evidence against RDKit 2026.09.1. The
265-row A6 set and holdouts v1–v3 have been inspected and are not final or
sealed evaluation data.

## Change

RDKit-compatible MMFF94 now reuses each interatomic distance from an energy
evaluation in the gradient evaluated immediately afterwards. The cache belongs
to one minimization, compares coordinates bit for bit, and uses generation
markers so changing coordinates does not clear an atom-count-squared array.
Public standalone energy and gradient APIs are unchanged. RDKit's special
stale-distance termination remains explicit.

The diagnostic `onlyInitialEmbedding` boundary also showed that seeded metric
embedding agrees with RDKit 2026.09.1 to floating-point noise. The observed
conformer divergence begins in the ETK minimization stage, not bounds, random
distance selection, or eigendecomposition.

## Fixed development comparison

- Host: macOS arm64, Python 3.13.6
- CheMatic: source candidate based on merge `b78e33bf`
- Oracle: RDKit 2026.09.1
- Corpus: frozen 265-row A6 development set
- Work: 20 ETKDGv3 candidates, seed 42, 200 MMFF94 iterations
- Quality profile: if no candidate converged, retry the lowest-energy
  non-converged candidate for 100 iterations
- Timing: three alternating blocks

The 100-row result is identical before and after the cache for every recorded
per-molecule CheMatic and RDKit field. The 265-row quality profile produced:

| Metric | CheMatic | RDKit |
|---|---:|---:|
| Seconds, block 1 | 55.523 | 56.443 |
| Seconds, block 2 | 55.085 | 56.035 |
| Seconds, block 3 | 54.661 | 56.325 |
| Median seconds | 55.085 | 56.325 |
| Generated candidates | 5,260 | 5,258 |
| Converged candidates | 2,855 | 2,859 |

CheMatic was faster in all three blocks; the median ratio was 1.023x. Its best
energy was at or below RDKit on 192/265 rows, within 1 kcal/mol on 261/265,
within 2 kcal/mol on 263/265, and the maximum positive delta was 2.969
kcal/mol. The default no-retry profile remains 257/265 within 1 kcal/mol.

This record does not establish complete or universal superiority. Acceptance
still requires a candidate commit, a once-only outcome-blind v4 set excluding
all previous development/evaluation rows, peak memory, a second host, at least
21 paired blocks for the final performance claim, and exact published
artifacts.
