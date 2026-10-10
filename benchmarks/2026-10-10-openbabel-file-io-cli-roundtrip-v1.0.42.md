# Open Babel file-I/O CLI round-trip candidate — 2026-10-10

This is a source-candidate measurement, not a v1.0.42 registry-artifact claim.
It measures a fresh process that parses one checked-in molecule and writes the
same format. Process startup is included. The candidate source was clean at
`c69f7bac0befafc320f46a3bf6e431f2c923f5f5`; the measured Open Babel executable
reported `Open Babel 3.2.1 -- Jul 11 2026 -- 19:21:27`.

## Prerequisite

Both engines' rewritten files were first read by the same Rust semantic probe.
CheMatic preserved the observed fields for all four fixtures. Open Babel tied
on V3000, CML, and CDXML, and changed the MOL2 residue label (`LIG1` to
`LIG11`). See the [semantic record](2026-10-10-openbabel-file-io-semantics-v1.0.42.json).

## Method

- Host: macOS 27.0.1 arm64, Darwin 27.0.0.
- Toolchain: rustc 1.97.0, Python 3.13.6.
- Binary: release `file_io_semantic_probe` versus `/opt/homebrew/bin/obabel`.
- Work: read one V3000, MOL2, CML, or CDXML fixture and write the same format.
- Schedule: 21 blocks, 10 fresh-process invocations per engine per block,
  alternating which engine ran first.
- Win rule: every paired block faster and the paired speedup 95% lower bound
  above 1.0.

## Result

| Format | Median CheMatic per invocation | Median Open Babel per invocation | Median paired speedup | 95% lower bound | All 21 blocks faster |
|---|---:|---:|---:|---:|---:|
| V3000 | 3.53 ms | 220.65 ms | 65.38x | 60.48x | yes |
| MOL2 | 3.40 ms | 255.77 ms | 66.88x | 59.33x | yes |
| CML | 2.20 ms | 140.11 ms | 65.50x | 62.52x | yes |
| CDXML | 2.25 ms | 140.12 ms | 63.84x | 62.66x | yes |

The strict rule passes in all four lanes. The magnitude is dominated by Open
Babel CLI startup and plugin loading; it must not be described as a parser or
writer throughput ratio.

## Boundary

This record does not measure parser-only work, writer-only work, warm-process
batch throughput, large files, peak RSS, malformed inputs, PDB, mmCIF, total
format breadth, or published packages. The raw alternating samples are in
[the JSON record](2026-10-10-openbabel-file-io-cli-roundtrip-v1.0.42.json).
