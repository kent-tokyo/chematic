# Timed unoptimized previews: flexible development probe

The optional `return_preview` mode records a corrected-gradient embedding
before MMFF setup and minimization, after checking declared stereo and finite
coordinates. The eventual optimization consumes those same coordinates. On
timeout or worker failure a pending preview, if any, is returned separately
from completed candidates, with no energy or convergence claim.

This experiment used the 20 previously exposed flexible development molecules
in `validation/results/2026-11-26-flexible-development/input.json`, not an
independent holdout. Each molecule was run twice with preview disabled and twice
with it enabled, alternating order within each run, requesting 20 attempts and giving the
Python process 0.2 seconds including startup. MMFF was limited to 200
iterations, with no repair retries. These times are not directly comparable to
the earlier reserved benchmark, whose budget started after worker initialization.
Both modes used the same local candidate wheel SHA-256
`e9e094c87278cf83bd3a075794a74b83f0bbbf4ab1e45b095edf2c101bdd2a75`,
not a published release; it was built with LTO disabled and 16 codegen units.
The candidate extension SHA-256 was
`093950a628c99b2abd3dc73a8344d1146405d1d0f26b504e87844407c5a5ac2c`.
The host was macOS 27.0.1 arm64, Python 3.13.6 and RDKit 2026.09.2.

All 80 calls timed out. The two sequential runs produced:

| Run and mode | Inputs with completed candidates | Completed candidates | Inputs with a preview | Inputs with any candidate |
| --- | ---: | ---: | ---: | ---: |
| First, preview disabled | 13/20 | 23 | 0 | 13/20 |
| First, preview enabled | 13/20 | 22 | 7 | 15/20 |
| Second, preview disabled | 13/20 | 23 | 0 | 13/20 |
| Second, preview enabled | 14/20 | 23 | 6 | 15/20 |

All previews passed the independent strict declared-stereo and coarse UFF
bond/angle/clash screens; none of the screens had incomplete parameter
coverage. No preview or completed candidate in either run passed those screens
and was within 1 Å of the reference crystal after symmetry-aware heavy-atom
alignment. The second run additionally stored candidate coordinates and status
fields so its screens can be recomputed. The first run's scorer snapshot is
retained because it did not store coordinates.

The 52 stored candidate screens in the second run were recomputed from the
stored coordinates and matched the recorded values on this host. These
results describe two small, non-isolated development runs. Previews are
unoptimized, the geometry test is coarse, and declared-stereo coverage has
limits. The experiment does not establish a general accuracy or timing
improvement, and the modes can produce different numbers of completed
candidates because of time variation and preview overhead. No RDKit conformer
generation comparison was performed here. The input, per-record outputs, scorer snapshot and screening results are in
`validation/results/2026-11-27-timed-preview/`; the current runner is
`scripts/benchmarks/timed_preview_development.py`.
