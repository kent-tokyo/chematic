# A6 MMFF gradient-rescue development record

Date: 2026-10-13

This is source-candidate development evidence against RDKit 2026.09.1. The
265 rows were already inspected and are not a sealed evaluation set. This
record does not establish universal or published-package superiority.

## Candidate

- Base: `d376aabe7d743524867b298370bcc41c11e89a4d`
- Implementation commit: `fe14064b`
- Python wheel SHA-256:
  `c876ecc28a4e8bc6f0c7595259b8f6da717436774fa200a61eeb6b47efb98cc2`
- RDKit: 2026.09.1
- Host: macOS arm64, Python 3.13.6
- Tier A manifest SHA-256:
  `6a478ea0f5d4ef067a4d1739e77a7209e8f76ecaa837e3f487c723dd6f465d6b`
- Tier B manifest SHA-256:
  `ec93ac160aa91c3c4416d57d3000d6e344e0b75cdf5c122211364c7495dce239`

The measured wheel was built immediately before the final public-input guard
that rejects a non-finite gradient threshold. The measured finite-threshold
path and optimization code are identical to `fe14064b`; this is not an exact
clean-commit artifact measurement.

## Fixed profile

Both engines generate 20 ETKDGv3 candidates per molecule with seed 42 and run
200 MMFF94 iterations. CheMatic then retries the lowest-energy candidate for
100 iterations only when:

1. the initial ensemble has no converged candidate; and
2. that candidate's already-computed RDKit normalized gradient test is at
   least `0.05`.

The threshold and budget were fixed from this exposed development corpus. A
retry is retained only when it has finite, non-increasing energy. Existing
APIs and default behavior remain unchanged.

## Result

| Metric | CheMatic | RDKit |
|---|---:|---:|
| Generated candidates | 5,260 | 5,258 |
| Converged candidates | 2,854 | 2,859 |
| Best energy within 1 kcal/mol of RDKit | 261/265 | reference |
| Best energy within 2 kcal/mol of RDKit | 263/265 | reference |
| At or below RDKit's best energy | 186/265 | reference |
| Maximum positive delta | +2.969 kcal/mol | reference |

Three alternating timing blocks were:

| Block | CheMatic (s) | RDKit (s) | Faster |
|---:|---:|---:|---|
| 1 | 55.464 | 58.394 | CheMatic |
| 2 | 57.784 | 58.801 | CheMatic |
| 3 | 57.769 | 60.524 | CheMatic |

The median ratio is `1.0179x` RDKit/CheMatic. This is a positive development
result, but three blocks on one host are insufficient for the repository's
formal speed-win contract. A fresh outcome-blind set, exact clean artifacts,
21 alternating blocks, a second host, and peak-memory accounting remain open.

Raw result:
[`a6-mmff-gradient-rescue-development-r3-2026-10-13.json`](../validation/results/a6-mmff-gradient-rescue-development-r3-2026-10-13.json)
(SHA-256 `166609165164d903386a760fe06092689f60af5b50bbc4a90ad2aeff9c6c79b8`).
