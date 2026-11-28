# Optional regular-geometry screening of completed timed conformers

The input is the first 12 records from the preceding fixed 30-molecule exposed
COD development selection in `../2026-11-28-timed-preview-regular/results.json`
(SHA-256 `29280526dc0b3256d1c15104d8ef36fa31436549807fa9a5445028c87c265ca3`).
This small selection is not a molecular holdout or a representative benchmark.
Each molecule is run once in normal batch mode and once in `return_partial`
mode, with one requested conformer, seed 42, at most 30 MMFF iterations and a
20-second process deadline. These two modes have different established
embedding profiles and are compared only against the same mode.

On the same macOS arm64 host with CPython 3.13.6, all 24 prior-wheel default
runs completed and returned one conformer. The candidate wheel's 24 default
runs had equal responses after excluding elapsed time. The candidate's 24
screened runs also completed and returned the same conformers; none were
rejected by the optional screen. After excluding elapsed time and the new
screening-count fields, their responses matched the candidate defaults. This is
not a runtime comparison, a crystal-structure accuracy test, or evidence that
the screen never rejects a real conformer.

A separately retained distorted-coordinate fixture verifies a nonzero
rejection. Worker tests also exercise force-field-unavailable fallback rows,
unknown backend rejection counts, impossible/mixed checkpoint counts, and a
worker killed after a completed screened attempt. The counter is disjoint from
backend rejection counts; partial recovery reports only completed attempts.
These are behavior tests, not a general stereo or force-field validation.

The prior wheel was built locally from commit `c33ef42b` with SHA-256
`51ee2466457abbdfc350446c09abe714b70f25af7cc2dbdbf5aabef13b0b78c2`.
The candidate was built from that commit plus the accompanying source changes,
SHA-256 `7f0c988b4d1a537542e88df914291f742b005c252a719c1088b8e790c1b3625b`.
Both are local debug wheels, not published distributions. `before.json` and
`after.json` include raw responses and elapsed times for traceability; timings
were not compared. To reproduce with each wheel installed:

```sh
python validation/results/2026-11-28-timed-conformer-regular/run.py before > before.json
python validation/results/2026-11-28-timed-conformer-regular/run.py after > after.json
```

The preview and completed-conformer screens are separate options. A passing
screen does not prove force-field convergence, physical geometry, or coverage
of every stereo class, and the optional screen may reject plausible structures.
