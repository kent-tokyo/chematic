# Prepared corrected embedding: flexible development comparison

Timed preview attempts previously reconstructed identical distance bounds and
torsion knowledge for every seed. The preview worker now prepares this
immutable data once per call and reuses it for each corrected-gradient
embedding. Preparation is inside the existing deadline. The non-preview
partial mode retains its existing execution path, and the worker still checks
stereo, records a provisional geometry and runs MMFF as before.

The local development experiment used the same 20 previously exposed flexible
molecules in `validation/results/2026-11-26-flexible-development/input.json`.
For each molecule, a single prepared object and the existing individual-call
API embedded seeds 42, 84 and 126. The timed total for the prepared path
includes its one-time initialization. The order of individual and prepared
calls alternated; SMILES parsing and MMFF were outside both timings. All 60
serialized coordinate results, including one unsuccessful embedding, matched
exactly. This is output agreement on this input and host, not a general proof.

The median of the 20 per-molecule three-attempt totals was 172.58 ms for
individual setup and 166.20 ms including prepared setup. The median paired
per-molecule relative change was -2.87%. This was a single, non-isolated
macOS arm64 development run; it does not prove a general speedup or a
whole-workflow improvement. It includes no independent accuracy comparison.

The candidate was a local source-built wheel, not a published release. Wheel
SHA-256: `913134369c6e406459366d0f7973087d7b7e25dd55323a84cfc28ad6781a139b`;
extension SHA-256:
`36680f037e9510e7f245a0ebb2b6b5a75b5e5076fb71ae959b6a29817ad4685b`.
It used Python 3.13.6 on macOS 27.0.1 arm64 with Rayon, OpenMP and OpenBLAS
restricted to one thread; the wheel was built with LTO disabled and 16
codegen units. Per-attempt times and result hashes are saved in
`validation/results/2026-11-27-prepared-embedding/development.json`; the
runner is `scripts/benchmarks/prepared_embedding_development.py`.

## Whole timed-preview development probe

The old and prepared implementations were also run in old/new/new/old order
on the same 20 exposed development molecules. Each block requested 20
attempts per molecule, with a 0.2-second deadline including process startup;
the preview and non-preview modes alternated their call order. This is the
same host and screening protocol as
`benchmarks/2026-11-27-timed-preview-development.md`, not an independent
holdout. Both local source-built wheels were built with LTO disabled and 16
codegen units. The old wheel SHA-256 was
`4d23fadd0b8c1258b157447206bc6a211b066d1873a3ed390ca952ae156b41e3`;
the new wheel is identified above.

| Block | Inputs with any preview-mode candidate | Completed preview-mode candidates | Pending previews |
| --- | ---: | ---: | ---: |
| Old 1 | 15/20 | 26 | 4 |
| New 1 | 15/20 | 24 | 5 |
| New 2 | 15/20 | 25 | 4 |
| Old 2 | 15/20 | 23 | 6 |

Every call timed out. Across blocks, 139 pairs of corresponding completed
candidate rows and 22 pairs of pending previews from the same attempt had
identical coordinates and, for completed candidates, identical status and
energy. The small time-limited experiment shows no detectable end-to-end
improvement; its outcomes are affected by process startup and timing noise.
It does not establish general accuracy or speed. Full per-record candidates,
screens and provenance are in `timed-old1.json`, `timed-new1.json`,
`timed-new2.json`, and `timed-old2.json` alongside the embedding-only data.
