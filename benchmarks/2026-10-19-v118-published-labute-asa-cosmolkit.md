# Published chematic 1.1.8 Labute ASA comparison

**Result:** the PyPI macOS arm64 wheel reproduces the 10,000-row Labute ASA
contract and the main release performance claim against RDKit. All values
agree within `1e-9`. The wheel is faster than RDKit in every one of the 20
paired blocks for first use, repeated use, and parse-inclusive use.

Against COSMolKit, the wheel wins all 20 blocks for repeated and
parse-inclusive use. First use is statistically tied in this run: the median
ratio is `1.009x`, the paired 95% interval crosses `1.0`, and chematic is
faster in 11/20 blocks. This lane is therefore **not** called a win.

This is one-host evidence for one published wheel and one operation. It does
not establish cross-platform or toolkit-wide superiority.

## Reader summary

- **Correctness:** chematic, the source baseline, and COSMolKit each agree
  with RDKit on 10,000/10,000 finite values.
- **First calculation:** chematic takes 4.555 ms per 10,000 already parsed
  molecules, 2.120x faster than RDKit. COSMolKit is statistically tied.
- **Repeated calculation:** chematic takes 0.608 ms, 1.581x faster than RDKit
  and 1.241x faster than COSMolKit.
- **Parse plus calculation:** chematic takes 65.848 ms, 6.429x faster than
  RDKit and 6.244x faster than COSMolKit.

The Python API measured here is:

```python
import chematic

mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
area = mol.labute_asa
```

## Artifacts and environment

- chematic: PyPI `1.1.8`, CPython 3.13 macOS 11+ arm64 wheel
- published filename:
  [`chematic-1.1.8-cp313-cp313-macosx_11_0_arm64.whl`](https://files.pythonhosted.org/packages/70/c6/a7469ae7e844cbd382807d89a96ce91ad08d333a3c25ecf3c61a7388eb33/chematic-1.1.8-cp313-cp313-macosx_11_0_arm64.whl)
- PyPI wheel SHA-256:
  `7dc15a10323ded4644291e530f559c4c84f2529f6fe6c05f0dce3b8fe09c2d5e`
- installed extension SHA-256:
  `ed5a04679987ee9b0a0a99daadca9f7a907b7af00e9e0c4180eaf49b9dbfcc55`
- PyPI upload time: `2026-10-18T22:21:48.417424Z`
- RDKit: `2026.09.1`
- COSMolKit: `0.5.0`
- comparison source baseline: chematic `1.1.7`, source revision `85f5000b`
- Python: CPython `3.13.6`
- host: macOS `27.0.1`, Apple arm64, 10 logical CPUs
- corpus: all 10,000 rows of `rdkit-js-browser-10k-v1.smi`
- corpus SHA-256:
  `f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`
- protocol: 20 rotating, sequential, fresh-process blocks
- [raw JSON](2026-10-19-python-labute-asa-cosmolkit-public-wheel.json)
- raw-record SHA-256:
  `710d7cb31a1f7b9784fbc5e81e04e9d3a33e7b70f317862b75bb837d0d800769`

The extension extracted from the downloaded wheel and the installed extension
have the same SHA-256. The wheel digest also matches PyPI's version metadata.
The benchmark workers recorded the public wheel's installed module path and
version on every block.

## Method

Every block runs chematic 1.1.7, public chematic 1.1.8, RDKit, and COSMolKit
in a rotating order. Each arm uses a fresh process. Native return objects are
retained during timing; output normalization and comparison happen afterward.

- **First use:** parse fresh molecules before timing, then call Labute ASA
  once.
- **Repeated use:** parse and prime fresh molecules before timing, then call
  Labute ASA again.
- **Parse inclusive:** parse and calculate Labute ASA inside the timed loop.

A speed win requires equivalent output, a paired 95% confidence lower bound
above `1.0`, and 20/20 faster blocks. Speedup is comparator time divided by
chematic time.

## Results

Median time is milliseconds per 10,000 molecules.

| Boundary | Public 1.1.8 | Source 1.1.7 | RDKit | COSMolKit |
|---|---:|---:|---:|---:|
| First use | **4.555** | 24.755 | 9.672 | 4.561 |
| Repeated use | **0.608** | 0.801 | 0.961 | 0.753 |
| Parse inclusive | **65.848** | 67.291 | 422.677 | 411.566 |

| Boundary | Comparator | Median ratio | Paired 95% CI | Faster blocks | Decision |
|---|---|---:|---:|---:|---|
| First use | source 1.1.7 | 5.363x | 5.331-5.468 | 20/20 | win |
| First use | RDKit | 2.120x | 2.072-2.179 | 20/20 | win |
| First use | COSMolKit | 1.009x | 0.980-1.056 | 11/20 | tie |
| Repeated use | source 1.1.7 | 1.322x | 1.238-1.390 | 20/20 | win |
| Repeated use | RDKit | 1.581x | 1.516-1.695 | 20/20 | win |
| Repeated use | COSMolKit | 1.241x | 1.162-1.324 | 20/20 | win |
| Parse inclusive | source 1.1.7 | 1.021x | 1.015-1.030 | 20/20 | win |
| Parse inclusive | RDKit | 6.429x | 6.373-6.471 | 20/20 | win |
| Parse inclusive | COSMolKit | 6.244x | 6.216-6.328 | 20/20 | win |

Median parse-only time is 63.733 ms for public chematic 1.1.8, 43.973 ms for
the source 1.1.7 baseline, 402.737 ms for RDKit, and 418.965 ms for COSMolKit.
The release moves bounded descriptor preparation into the Python parse
boundary, so the first-use result must always be read together with the
parse-inclusive result.

Median whole-process peak RSS is 97.7 MiB for public chematic 1.1.8,
97.8 MiB for the source baseline, 301.2 MiB for RDKit, and 399.4 MiB for
COSMolKit. Imports and every measured mode are included, so this is
descriptive process memory rather than isolated operation memory.

## Ring-count boundary

Ring count is not part of this timing result. On the same corpus, chematic and
RDKit agree on all 10,000 rows while COSMolKit differs on six rows. Equivalent
output is a prerequisite for speed ranking, so no chematic/COSMolKit
ring-count speed result is reported. The exact rows remain recorded in the
[source-candidate report](2026-10-19-python-labute-asa-cosmolkit.md).

## Reproduction

```text
python scripts/bench_python_speed_frontier.py \
  --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  --limit 10000 --blocks 20 --ops labute_asa \
  --base-revision 85f5000b \
  --baseline-wheel <source-1.1.7-wheel> \
  --baseline-path <source-1.1.7-install> \
  --candidate-wheel <downloaded-PyPI-1.1.8-wheel> \
  --candidate-path <public-1.1.8-install> \
  --rdkit-path <pinned-RDKit-install> \
  --cosmolkit-path <pinned-COSMolKit-install> \
  --output benchmarks/2026-10-19-python-labute-asa-cosmolkit-public-wheel.json
```

The raw record contains all 20 timings, operation outputs as per-engine
digests, module paths, extension hashes, process RSS, load averages, and the
full paired summaries. A second host and other published wheel platforms are
still needed before making a broader performance claim.
