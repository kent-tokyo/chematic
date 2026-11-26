# Cooperative minimization development check

This local check compares the unchanged `minimize_` API on the parent and
candidate source wheels, and the candidate's `minimize_with_timeout_` with a
60-second budget. The input is the first 20 phosphorus/sulfur records of the
previously exposed [force-field development set](2026-11-27-forcefield-hetero-development.md),
using its saved explicit-hydrogen coordinates. This is a selected development
sample, not an independent holdout. Each of MMFF94 and UFF was run three times
per molecule in each block. Blocks ran in parent, candidate, candidate, parent
order; the candidate block also interleaved 0.1 ms and 1 ms deadline probes.

All 120 status, energy and coordinate digests per block agreed exactly across
both APIs and both source builds. There were no exceptions. Across the two
candidate blocks, the 0.1 ms probes timed out in all 80 cases. At 1 ms, 76
timed out and four converged first. For every timed-out case with an available
energy, reevaluation at the returned coordinates agreed exactly. The largest
recorded ratio of method elapsed time to requested deadline was 1.152. This is
only a small local sample; it is not a guarantee of a maximum overrun.

The median call times per block were 7.51–7.53 ms for parent/candidate
`minimize_`, and 7.30–7.63 ms for the candidate's cooperative method,
depending on force field and block. Other tests were running on the same host,
and there was no controlled CPU isolation. These observations are not a general
speed comparison. Field preparation was outside the timed method.

The host was macOS 27.0.1 arm64 with Python 3.13.6. Both local wheels were
built with release LTO disabled and 16 codegen units; neither is a published
wheel. The parent was commit `b994f5e`; its wheel SHA-256 was
`913134369c6e406459366d0f7973087d7b7e25dd55323a84cfc28ad6781a139b`.
The candidate wheel used for the four timing blocks, before an internal
function/test-file split, had SHA-256
`504dca735c1d612301d6ea8afb177fc3b10d88e1c897e53dc9ef361d24dac87e`.
The final source wheel, after that split, had SHA-256
`a29de59e5d58e63e09e4f315147948bb3d01ff34f686ffdf81f9a4007aea9c51`.
A separate final-wheel block again matched all 120 parent digests for each API;
it is retained as `coop-final.json` and is not included in the timings above.
Input and coordinate SHA-256 values were
`86218b8c8dac403fc025e0cf9d739c0f1cbcc680cbe8aa139be9bdf7fd1d14f4` and
`6b642f40b8007bbe217c653ee115f3078ed624c0f8c914a6d0299eb2c08257ef`.

The [probe](../scripts/validation/measure_cooperative_minimization.py), SHA-256
`72006ffd8051136a31bf4ad726bf0e41254aca5f69fec81458da854846ddf8bf`, accepts
`old` or `new` and an output path. Set `PYTHONPATH` to the corresponding wheel
installation and run it from the repository root. The raw blocks are in
`validation/results/2026-11-27-cooperative-minimize/`. Each block was run
with `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, and `RAYON_NUM_THREADS` set to 1.
