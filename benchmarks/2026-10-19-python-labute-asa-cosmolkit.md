# Python Labute ASA versus RDKit and COSMolKit

**Decision: adopt the source candidate.** On the fixed 10,000-molecule corpus,
the candidate returns the same finite Labute ASA values as its source baseline,
RDKit, and COSMolKit. It is faster than both comparators in all 20 paired
blocks for first use, repeated use, and parse-inclusive use.

This is one-host source-wheel evidence. It is not a published-package,
second-host, or toolkit-wide performance claim.

## Reader summary

This result is relevant when a Python workflow parses many SMILES strings and
calculates RDKit-compatible Labute ASA values, or repeatedly calculates ASA on
already parsed molecules.

- **Correctness:** all 10,000 values agree with RDKit and COSMolKit.
- **First calculation:** 4.102 ms for 10,000 parsed molecules, 2.358x faster
  than RDKit and 1.114x faster than COSMolKit in this run.
- **Parse plus calculation:** 66.107 ms for 10,000 SMILES, 6.390x faster than
  RDKit and 6.229x faster than COSMolKit in this run.
- **Availability:** these speed figures apply to source candidate `524a0210`.
  They do not yet describe the current registry package.

The Python API is intentionally small:

```python
import chematic

mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
area = mol.labute_asa
```

The same protocol has now been rerun against the
[published PyPI 1.1.8 wheel](2026-10-19-v118-published-labute-asa-cosmolkit.md).
That release artifact reproduces the RDKit wins. Its first-use comparison with
COSMolKit is statistically tied rather than a strict win, so use the published
record for package-selection claims and this document for source-adoption
provenance.

## Scope and artifacts

- candidate: `524a0210` (`perf(descriptors): accelerate first-use Labute ASA`)
- source baseline: `85f5000b`
- RDKit: `2026.09.1`
- COSMolKit: `0.5.0`
- Python: CPython `3.13.6`
- host: macOS `27.0.1`, Apple arm64, 10 logical CPUs
- corpus: all 10,000 rows of `rdkit-js-browser-10k-v1.smi`
- corpus SHA-256:
  `f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`
- protocol: 20 rotating, sequential, fresh-process blocks
- raw record: [JSON](2026-10-19-python-labute-asa-cosmolkit.json)
- raw-record SHA-256:
  `d02769a161f48541e5fc6b7fef739ce505bc219aa8220bd8181eb29306a31670`

The candidate wheel SHA-256 is
`fcfd8c20815d22c2efa4f8a83f4ca7bf52039950df7245f38ba5af3751958d3d`;
its loaded extension SHA-256 is
`67aa67aaa342d6a62e98f9d20a6b2733aeef0f29e2d1651bc6d58ccf6f42844d`.
The source baseline wheel SHA-256 is
`b5b42c77bebc8a7347e0b6e544d28d684acbd5ef3675e69fd0e154facac5dcb8`.
Both local wheels required the same metadata-only Mach-O string-table
alignment repair to load on this host. The repair did not alter Rust code or
chemistry data.

## Method

Each block runs the four arms in a rotating order. Output conversion and
comparison are outside the timed interval.

- **First use:** parse fresh molecules before timing, then call Labute ASA once.
- **Repeated use:** parse and prime fresh molecules before timing, then call it
  again.
- **Parse inclusive:** parse and call Labute ASA inside the same timed loop.

An equivalent-output speed win requires 10,000/10,000 comparable values, a
paired 95% confidence lower bound above `1.0`, and 20/20 faster blocks. The
speedup is `competitor time / candidate time`.

## Results

Median time is milliseconds per 10,000 molecules.

| Boundary | Candidate | Source baseline | RDKit | COSMolKit |
|---|---:|---:|---:|---:|
| First use | **4.102** | 24.736 | 9.614 | 4.579 |
| Repeated use | **0.661** | 0.828 | 0.949 | 0.756 |
| Parse inclusive | **66.107** | 67.380 | 422.047 | 412.234 |

| Boundary | Comparator | Speedup | Paired 95% CI | Faster blocks | Decision |
|---|---|---:|---:|---:|---|
| First use | source baseline | 6.030x | 6.000–6.064 | 20/20 | win |
| First use | RDKit | 2.358x | 2.322–2.414 | 20/20 | win |
| First use | COSMolKit | 1.114x | 1.106–1.136 | 20/20 | win |
| Repeated use | source baseline | 1.266x | 1.224–1.283 | 20/20 | win |
| Repeated use | RDKit | 1.422x | 1.388–1.461 | 20/20 | win |
| Repeated use | COSMolKit | 1.140x | 1.104–1.160 | 20/20 | win |
| Parse inclusive | source baseline | 1.020x | 1.012–1.022 | 20/20 | win |
| Parse inclusive | RDKit | 6.390x | 6.337–6.406 | 20/20 | win |
| Parse inclusive | COSMolKit | 6.229x | 6.200–6.245 | 20/20 | win |

Candidate, source baseline, and COSMolKit each match RDKit on all 10,000
Labute ASA values within the lane's `1e-9` tolerance. Candidate and baseline
values are exactly equal. Median parse-only time is 64.372 ms for the
candidate, 43.879 ms for the baseline, 402.675 ms for RDKit, and 418.002 ms
for COSMolKit. The candidate deliberately moves bounded model preparation to
the Python parse boundary; first-use results must therefore be reported with
parse-inclusive results.

Median whole-process peak RSS is 98.0 MiB for both chematic arms, 302.0 MiB
for RDKit, and 399.7 MiB for COSMolKit. This includes imports and all measured
operations, so it is descriptive rather than a per-operation memory ranking.

## Ring-count residuals are not a chematic speed deficit

Ring-count timing remains separate from this ASA result. On the same fixed
corpus, chematic and RDKit `2026.09.1` agree on all 10,000 rows. COSMolKit
differs on six rows:

| Zero-based row | Corpus line | chematic | RDKit | COSMolKit |
|---:|---:|---:|---:|---:|
| 9 | 10 | 10 | 10 | 8 |
| 23 | 24 | 22 | 22 | 10 |
| 28 | 29 | 22 | 22 | 9 |
| 29 | 30 | 22 | 22 | 10 |
| 30 | 31 | 22 | 22 | 10 |
| 34 | 35 | 11 | 11 | 9 |

The corpus hash and row numbers identify every structure without duplicating
long SMILES in this report. Because equivalent outputs are a prerequisite for
speed ranking, no chematic/COSMolKit ring-count speed result is reported.
Changing chematic to reproduce these six COSMolKit values would create six
differences from the pinned RDKit contract; no such change was made.

## What this record does not establish

- It does not show that every chematic descriptor is faster.
- It does not compare full application pipelines, database I/O, or parallel
  execution.
- It does not establish behavior on Linux, Windows, or a second macOS host.
- It does not turn the six unequal COSMolKit ring counts into a correctness or
  speed claim.
- It does not replace published-wheel chemistry and compatibility gates.

## Implementation and checks

The candidate:

- reuses a bounded RDKit-compatible scalar model at the Python parse boundary;
- retains lazy parsing in the low-level Rust API;
- removes duplicate bond lookups and unnecessary Labute distance clamps;
- scans small ring-closure sets with an ordered cursor; and
- invalidates all prepared state after graph mutation.

`cargo fmt --check`, `git diff --check`, Clippy for the modified crates, the
new model-preparation cache test, and 11 Labute/VSA-focused tests pass. A full
workspace test was not rerun in this low-disk session, so hosted CI and an
exact public-artifact rerun remain release gates.

## Reproduction

```text
python scripts/bench_python_speed_frontier.py \
  --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  --limit 10000 --blocks 20 --ops labute_asa \
  --base-revision 85f5000b \
  --baseline-wheel <baseline-wheel> --baseline-path <baseline-install> \
  --candidate-wheel <candidate-wheel> --candidate-path <candidate-install> \
  --output benchmarks/2026-10-19-python-labute-asa-cosmolkit.json
```

Before a release claim, repeat this protocol with exact public artifacts in
the available hosted CI environment. The performance claim remains explicitly
scoped to this single macOS host unless another comparable host is measured.
Keep first-use, repeated-use, parse-inclusive, parse-only, and memory
observations separate.
