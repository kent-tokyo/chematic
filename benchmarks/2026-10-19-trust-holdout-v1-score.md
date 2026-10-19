# Frozen trust holdout v1: one-shot score

Date: 2026-10-19. The candidate is the immutable v1.1.8 main commit
`00596fb1ea5397524f98e3789c41f17ee8616fdc`; the Python artifact reports
version `1.1.8`. RDKit `2026.09.1` is the pinned compatibility oracle. The
frozen manifest SHA-256 is
`18ac28b8674fc202c9c8efebc1102e4de6d2d8cab7ff1bb2fa1e7b5bba6262b6`.
No chemistry implementation was changed after inspecting these results.

Every observable is reported independently as success, mismatch, refusal,
invalid input or internal error. There were no refusals, invalid inputs or
internal errors in these three lanes.

| Lane / observable | Success | Mismatch | Total |
|---|---:|---:|---:|
| MCS parse | 64 | 0 | 64 |
| MCS completion | 64 | 0 | 64 |
| MCS atom count | 31 | 33 | 64 |
| MCS bond count | 31 | 33 | 64 |
| SDF parse | 64 | 0 | 64 |
| SDF graph | 64 | 0 | 64 |
| SDF title | 64 | 0 | 64 |
| SDF properties | 64 | 0 | 64 |
| SDF stereo | 64 | 0 | 64 |
| 3D parse | 96 | 0 | 96 |
| 3D stereo graph | 96 | 0 | 96 |
| 3D embed status | 96 | 0 | 96 |
| 3D conformer count | 95 | 1 | 96 |
| 3D coordinates | 4 | 92 | 96 |

The coordinate column compares all pairwise distances with a maximum absolute
tolerance of `1e-4` Å, so it is invariant to rigid rotation and translation.
It is intentionally much stricter than a conformer-quality or energy test.
The score shows three different conclusions: SDF attribute/stereo preservation
generalizes on this holdout, MCS size parity does not, and multi-conformer
ETKDG status/stereo parity does not imply coordinate parity.

The first draft score was discarded before publication because the scorer
mapped the MCS `order` profile incorrectly, counted SDF `_Name` twice and used
raw rather than rigid-invariant coordinates. Only comparator instrumentation
was corrected; the candidate and frozen rows were unchanged. The complete
row-level result is
[`validation/results/2026-10-19-trust-holdout-v1-v118.json`](../validation/results/2026-10-19-trust-holdout-v1-v118.json).

This is a one-shot compatibility evaluation, not training data and not an
independent physical-quality benchmark. Future work may fix the implementation
on new development data, but this frozen score remains historical and must not
be rewritten.

