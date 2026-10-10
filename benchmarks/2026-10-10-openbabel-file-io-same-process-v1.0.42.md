# Open Babel 3.2.1 same-process file-I/O comparison

## Result

CheMatic is faster in every one of the 21 alternating blocks for all 12
checked lanes. The paired 95% speedup lower bound also exceeds 1.0 in every
lane.

| Format | Parse | Write | Round trip |
|---|---:|---:|---:|
| V3000 | 6.38x | 4.06x | 5.67x |
| MOL2 | 2.79x | 3.11x | 2.86x |
| CML | 4.24x | 3.27x | 4.00x |
| CDXML | 5.77x | 2.59x | 4.78x |

Values are median paired speedups. The weakest paired 95% lower bound is
2.56x for CDXML write. The raw block timings, fixture and executable hashes,
and every lower bound are in the [JSON record](2026-10-10-openbabel-file-io-same-process-v1.0.42.json).

## Method

- CheMatic source: clean commit `5f247956519a1650e1e4c9caa526859ed2e48667`
- Comparator: Open Babel 3.2.1
- Host: macOS 27.0.1, arm64
- Repetitions: 5,000 operations inside each process
- Blocks: 21, alternating which engine runs first
- Timer boundary: the in-process operation loop; executable startup and plugin
  discovery are excluded
- Semantic prerequisite: the same four fixtures pass the common-observer
  semantic gate before timing

## Boundary

This is source-candidate evidence for one small checked-in fixture per format.
It establishes bounded hot-loop latency for V3000, MOL2, CML, and CDXML. It
does not establish large-file throughput, peak-memory superiority, broad-corpus
format fidelity, the remaining Tier-A formats, or published-package performance.
