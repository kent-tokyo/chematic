# Partial conformer checkpoint development comparison

The opt-in timed Python wrapper previously rewrote its complete accumulated
ensemble after every attempt. It now appends one newline-terminated record
per completed attempt. The reader ignores a final unterminated record from an
interrupted write, rejects malformed or out-of-order completed records, and
reconstructs the same public result fields. Worker failures and unreadable
checkpoints remain distinct. This is an internal temporary format, not a
persistent file format or a change to the default non-partial API.

The checkpoint-only experiment used a synthetic worker returning 60 fixed
coordinate triples per attempt, with no embedding, minimization or chemistry.
On macOS 27.0.1 arm64, Python 3.13.6, two repetitions per size gave:

| Attempts | Before: median write | After: median write | Before: median read | After: median read |
| --- | ---: | ---: | ---: | ---: |
| 100 | 0.512 s | 0.0057 s | 0.0017 s | 0.0019 s |
| 400 | 8.067 s | 0.0224 s | 0.0069 s | 0.0077 s |
| 800 | 32.376 s | 0.0446 s | 0.0140 s | 0.0157 s |

Times include JSON serialization and filesystem operations; the two modes
ran sequentially on a non-isolated local host. This demonstrates checkpoint
scaling in the synthetic workload, not an end-to-end molecular speedup. The
before wrapper was from commit `a61648f`, SHA-256
`5094be14ae4d03510b45266f24de47c316e832afd6bdf81ca79abe9b3aa3c0dd`;
the measured after wrapper was SHA-256
`9f3215c93a9455fa845dc2a335ab71087176270d845610024ea98c13fd7aa2fb`.
The scripts and individual timings are stored in
`validation/results/2026-11-27-partial-checkpoint/` and
`scripts/benchmarks/partial_checkpoint_scaling.py`.

For a separate output check, both wrappers processed the 20 previously
exposed flexible development molecules, three corrected-gradient attempts per
molecule, same seed, with the source-built extension SHA-256
`af8b3a615fc163232d6d277f11bf5829d7334b94eef5adf986dfad5f37af3893`.
All 20 full result objects, including conformer coordinates and diagnostics,
matched exactly after JSON serialization, and both formats reconstructed
their own results. Only output hashes and counts, not coordinates, are kept
in `equivalence.json`. This is not independent chemical accuracy evidence.
