# CDK 2.13 / Indigo 1.46.0 bounded baseline

This record establishes a reproducible development baseline against CDK 2.13
and Indigo 1.46.0. It is deliberately a smoke comparison, not a universal
superiority claim.

## Scope and provenance

- CheMatic source commit: `ca1efa847a007b7be30ee73b427d504da0aedbc7`
- CheMatic wheel SHA-256:
  `44178cd07d8fd343ff0f8ccc78a50f23b9a39028b1eb1408da49fc8ae6dbd34f`
- CDK: `org.openscience.cdk:cdk-bundle:2.13`; adapter jar SHA-256:
  `f5e82454a60f0a8d1a1303373d41350b5dc8d50f8eff9cfe08a48407df9df997`
- Indigo: `epam.indigo==1.46.0`
- Accuracy reference: RDKit `2026.03.6`
- Corpus: 24 checked-in stratified rows, SHA-256
  `6bd9b6c34d42e9645fd4a8a50e4dbd1565b42b813fcf5b21d69a4786878498bd`
- Host: Apple M4, arm64, macOS 27.0.1, Python 3.13.6, Temurin JDK
  21.0.12.1
- Timing: persistent processes, 21 rotating-order blocks, separate
  parse-inclusive and prepared lanes, paired bootstrap intervals.

The adapters call public APIs only. No CDK or Indigo implementation code was
copied. A speed row is rankable only when both engines produce the pinned
reference result on every corpus row.

## Accuracy

CheMatic matched the pinned RDKit result on all 24 rows for all 15 measured
operations. The comparator residuals are:

| Operation | CheMatic | CDK | Indigo |
|---|---:|---:|---:|
| canonical stability, formula, parse, six SMARTS queries | 24/24 | 24/24 | 24/24 |
| HBD | 24/24 | 24/24 | 16/24, 4 errors |
| HBA | 24/24 | 16/24 | 10/24 |
| molecular weight | 24/24 | 20/24 | 22/24 |
| TPSA | 24/24 | 22/24 | 20/24, 4 errors |
| MOL V2000 semantic round trip | 24/24 | 19/24 | 19/24 |
| ring count, diagnostic only | 24/24 | 23/24 | 23/24 |

The ring-count row is not an accuracy verdict: CDK and Indigo expose SSSR
counts while RDKit and CheMatic use the pinned RDKit ring model. Cubane differs
by basis definition.

The direct HBD implementation was also checked separately on the pinned
10,000-row browser corpus: 10,000/10,000 values matched RDKit 2026.03.6.

This gate found and fixed two CheMatic defects: V2000 `M  ISO` records were
lost, and coordinate-derived E/Z directions did not reach the RDKit-compatible
SMILES model. The fixed cases have Rust regression tests.

## Speed

The table reports competitor time divided by CheMatic time. Values above 1.0
favour CheMatic. All listed confidence intervals are entirely above 1.0.

| Equivalent operation | Lane | CDK / CheMatic | Indigo / CheMatic | CDK blocks won |
|---|---|---:|---:|---:|
| parse | pipeline | 1.65x | 7.23x | 21/21 |
| canonical stability | pipeline | 8.24x | 2.51x | 21/21 |
| formula | pipeline | 5.38x | 20.50x | 21/21 |
| HBD | pipeline | 2.31x | not rankable | 21/21 |
| HBD | prepared | 6.60x | not rankable | 21/21 |
| six SMARTS queries | pipeline | 1.25–1.34x | 11.86–15.66x | five 21/21; `[#7]` 20/21 |
| six SMARTS queries | prepared | 1.45–2.38x | 24.44–156.04x | 21/21 each |

A separate 1,000-iteration hot run confirmed prepared HBD at 6.82x and the
three atomic-number SMARTS queries at 1.43–1.62x in the prepared lane. A
5,000-iteration diagnostic retained one losing block for pipeline and prepared
`[#8]`, despite medians of 1.30x and 1.68x and 95% lower bounds above 1.0.
Under the predeclared every-block rule, that lane remains open.

## Boundary and next gate

The evidence proves only the named operations, versions, corpus, artifact and
host. It does not cover the full CDK or Indigo feature sets. A public claim
requires a frozen broad corpus, adversarial strata, peak RSS, a second host,
and exact published CheMatic artifacts. Descriptor rows on which a comparator
does not match the reference are accuracy evidence but are not ranked for
speed.

Raw evidence:

- [`2026-10-10-cdk-indigo-smoke-scorecard.json`](2026-10-10-cdk-indigo-smoke-scorecard.json)
- [`2026-10-10-cdk-indigo-smoke-benchmark.json`](2026-10-10-cdk-indigo-smoke-benchmark.json)
- [`2026-10-10-cdk-indigo-smoke-hot-benchmark.json`](2026-10-10-cdk-indigo-smoke-hot-benchmark.json)

