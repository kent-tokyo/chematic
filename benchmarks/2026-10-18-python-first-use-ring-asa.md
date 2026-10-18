# Python first-use Labute ASA and ring count

> Historical record. The later
> [10,000-row RDKit/COSMolKit follow-up](2026-10-19-python-labute-asa-cosmolkit.md)
> supersedes the current-candidate ASA conclusion. This 1,000-row result is
> retained to show the earlier development boundary.

This source-candidate record tests the two first-use Python operations that
still trailed RDKit in the 2026-10-18 performance audit. It is evidence for
one corpus and one macOS arm64 host, not a published-package or toolkit-wide
claim.

## Candidate and method

- chematic source: `6b0db101` (base `cfcbbdfc`)
- comparator: RDKit `2026.09.1`
- additional context: COSMolKit `0.5.0`
- Python: CPython `3.13.6`
- host: macOS `27.0.1`, Apple arm64, 10 logical CPUs
- corpus: rows 0–999 of `rdkit-js-browser-10k-v1.smi`
- corpus SHA-256:
  `f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`
- protocol: 20 rotating fresh-process blocks
- raw record: [JSON](2026-10-18-python-first-use-ring-asa.json)

The candidate prepares the RDKit-compatible aromatic-view decision and ring
count at the Python SMILES boundary for molecules with at most 128 atoms, 192
bonds, and 16 SMILES ring-closure bonds. Larger or unusually cyclic inputs
retain the prior lazy path. Graph mutation continues to invalidate every
derived cache. The Rust parser itself remains lazy.

The locally built wheel required a metadata-only Mach-O alignment repair to
load with this host linker. The repaired wheel SHA-256 is
`b5b42c77bebc8a7347e0b6e544d28d684acbd5ef3675e69fd0e154facac5dcb8`;
the loaded extension SHA-256 is
`118ca4519720914f0b572af3e0978421e957be713a9437c6674eea988bb9857f`.
No Rust code or chemistry data changed during that repair.

## Results

Times are medians in milliseconds per 1,000 molecules. Speedup is
`RDKit / chematic`.

| Operation | Boundary | chematic | RDKit | Speedup | 95% CI | Winning blocks |
|---|---|---:|---:|---:|---:|---:|
| Labute ASA | first use | 0.668 | 1.252 | **1.908x** | 1.829–1.928 | 20/20 |
| Labute ASA | repeated use | 0.044 | 0.117 | **2.665x** | 2.541–2.739 | 20/20 |
| Labute ASA | parse inclusive | 7.092 | 66.300 | **9.321x** | 9.252–9.497 | 20/20 |
| ring count | first use | 0.054 | 0.076 | **1.419x** | 1.316–1.529 | 20/20 |
| ring count | repeated use | 0.039 | 0.066 | **1.726x** | 1.673–1.781 | 20/20 |
| ring count | parse inclusive | 6.586 | 65.798 | **9.991x** | 9.886–10.111 | 20/20 |

Candidate and RDKit outputs agree on all 1,000 rows for both operations.
Median parse time is 6.695 ms for chematic and 63.422 ms for RDKit. Median
whole-process peak RSS is 50.6 MiB and 97.6 MiB respectively.

COSMolKit is faster on first-use Labute ASA in this lane (0.570 ms versus
0.668 ms), while chematic is faster after reuse and parse-inclusive. Ring
timings are not ranked against COSMolKit because its ring count differs from
RDKit on six rows.

## Interpretation and limits

This closes the measured first-use deficits against RDKit on this fixed
Python/macOS lane without changing numeric results. It does so by moving work
that RDKit already performs during `MolFromSmiles` into chematic's bounded
Python parse path. Therefore first-use and parse-inclusive timings must always
be reported together; this is not a claim that the underlying standalone Rust
descriptor became faster.

The JSON contains identical `baseline` and `candidate` arms because the exact
pre-change wheel was unavailable during the final run. Those arms must not be
used as before/after evidence. The RDKit comparisons are paired within this
run. Previous pre-change timings came from a separate audit and are not
treated as a paired comparison here.

Before a release claim, repeat the lane from an unmodified public wheel on
hosted Linux and another macOS run. Keep oversized/high-cycle inputs, peak
memory, and parser resource limits in the regression set.
