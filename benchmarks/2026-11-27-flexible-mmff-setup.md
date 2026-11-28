# Flexible conformer development: MMFF setup reuse

The checked mixed fallback previously constructed the same MMFF field twice
from the same first embedded conformer when MMFF setup succeeded. It now uses
the first field for optimization. The field is still rebuilt at a converged
candidate before its final gradient check: its nonbonded pair list depends on
coordinates and cannot in general be reused across conformers.

A local exploratory comparison used the 20 previously exposed flexible
molecules in `validation/results/2026-11-26-flexible-development/input.json`.
The before run loaded a source-built macOS arm64 extension SHA-256
`af8b3a615fc163232d6d277f11bf5829d7334b94eef5adf986dfad5f37af3893`;
the after run loaded an installed development wheel SHA-256
`f9a04b8dfaf3d64cc173ee884cc4835191a545db1d014e80ab970c0677ddbdeb`,
containing extension SHA-256
`17f0d46b52b0eb6e8bb30edae084f0534404f01a0483549029502c5923a2403f`.
The wheel is not a published release and was built with LTO disabled and 16
codegen units. Python 3.13.6 ran on macOS 27.0.1 arm64 with Rayon, OpenMP and
OpenBLAS restricted to one thread. Each call was timed three times, rotating
the order of standalone corrected embedding, checked embedding with zero
MMFF iterations, and checked embedding with up to 200 MMFF iterations.
These are independent calls, not an exact decomposition of one call.

| Run order | Embed median | Checked, 0 iterations median | Checked, 200 iterations median |
| --- | ---: | ---: | ---: |
| Before 1 | 51.97 ms | 52.31 ms | 78.36 ms |
| After 1 | 52.94 ms | 53.55 ms | 80.47 ms |
| After 2 | 53.05 ms | 53.62 ms | 80.40 ms |
| Before 2 | 53.51 ms | 53.98 ms | 81.57 ms |

The medians are medians of the per-molecule timing medians. Each operation's
serialized outputs matched exactly across all three repeats and all four
blocks, including the error for one standalone embedding failure. This is
output equivalence for these inputs only. Timings changed with run order even
for the unchanged standalone embedding. This small, non-isolated run does not
establish a speedup or general numerical equivalence. The experiment input,
per-call times and output digests are stored under
`validation/results/2026-11-27-flexible-profile/`; the runner is
`scripts/benchmarks/flexible_profile_development.py`.
