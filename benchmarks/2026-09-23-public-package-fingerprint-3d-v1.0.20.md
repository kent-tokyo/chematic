# Public-package fingerprint and 3D comparison — v1.0.20

This record uses registry-installed `@kent-tokyo/chematic@1.0.20` and
`chematic==1.0.20`. Results are scoped to the named corpus, operation, host,
and comparator; they do not establish universal superiority over RDKit.

## Fixed environment

- Host: macOS 26.5.2 arm64, CPython 3.13.6, Google Chrome/Chromium.
- Comparators: `@rdkit/rdkit@2026.3.6` and `rdkit==2026.3.6`.
- Fingerprints: exposed 10,000-SMILES corpus, radius 2, 2,048 bits, identical
  packed 256-byte output, 20 fresh browser processes per arm.
- 3D: fixed 265-molecule Tier A+B corpus, seed `20260801`, at most eight
  embedding attempts, and the same external geometry/stereo scorer.
- PyPI wheel SHA-256:
  `7eacc40e9eb76e2253add9b318b2291c490c53cfbcae59070a1e409cd65f646b`.
- npm package integrity:
  `sha512-HhzFqZwKnCdDGbbyqiaGRfhYooSt6R7Cr06hC+qJD2lGYw0GiPRlJ4huWli6T5Y3ErkMqWbA89vkE5RwKT82Dw==`.

## Compatible Morgan result

| Operation | CheMatic p50 ms/mol | RDKit.js p50 ms/mol | Geometric speedup | 95% lower bound |
|---|---:|---:|---:|---:|
| Parse + prepare + Morgan | 0.071567 | 0.099155 | **1.398x** | **1.363x** |
| Prepared Morgan | 0.006876 | 0.024327 | **3.511x** | **3.407x** |

All 20 paired repetitions favored CheMatic. Minimum observed speedups were
1.210x and 2.898x. The prepared API matched RDKit.js bit-for-bit on all 9,999
supported rows. One high-coordinate Fe(II) row remains the declared typed
coordination refusal. The published WASM is 4,210,700 bytes raw and 1,542,710
bytes under gzip-9; its SHA-256 is
`0f477fdbe8a2e9378c47829dfbf19c63320522061c50e2bbf8b2528e7d3b4c15`.

## 3D result

`usable` means successful, independently sound, and stereo-clean under the
shared external scorer. Timing ratios use paired common successes.

| Lane | CheMatic usable | RDKit usable | Geometric speedup | 95% lower bound |
|---|---:|---:|---:|---:|
| Raw embedding | 212/265 | 264/265 | 3.968x | 3.599x |
| UFF | 263/265 | 264/265 | 3.292x | 3.023x |
| UFF best-of-10 | **265/265** | 264/265 | **2.359x** | **2.198x** |
| MMFF94, lightweight/ignore stereo | 211/265 | 264/265 | 1.514x | 1.411x |
| UFF stereo-safe | 262/265 | 264/265 | 2.589x | 2.358x |
| MMFF94 stereo-safe | **265/265** | 264/265 | 0.944x | 0.861x |

The published MMFF94 stereo-safe lane succeeds on all 265 inputs, is
independently sound and stereo-clean on all 265, and has zero gross-clash
rows. This confirms that the historical 12 typed failures and four clash rows
are closed in the public package. Its p50 is 32.613 ms versus RDKit's 28.732 ms,
so quality-equivalent MMFF94 is still slightly slower and is not a speed win.
The faster lightweight MMFF94 lane is not quality-equivalent.

## Evidence

- `validation/results/competitive-browser-rdkitjs-published-npm-10k-chromium-20rep-prepared-v1.0.20-2026-09-23.json`
- `validation/results/competitive-browser-rdkitjs-ecfp4-prepared-parity-10k-v1.0.20-2026-09-23.json`
- `validation/results/public-package-3d-chematic-v1.0.20-2026-09-23.{jsonl,meta.json}`
- `validation/results/public-package-3d-common-scored-v1.0.20-vs-rdkit-2026.3.6-2026-09-23.jsonl`
- `validation/results/public-package-3d-summary-v1.0.20-vs-rdkit-2026.3.6-2026-09-23.json`
- Comparator: `validation/results/public-package-3d-rdkit-2026.3.6-2026-09-22.{jsonl,meta.json}`

## Decision

The public package establishes scoped speed wins for parse-inclusive and
prepared compatible Morgan, plus UFF best-of-10 with one more usable row than
RDKit. It closes the named MMFF94 quality cohort, but does not establish a
quality-equivalent MMFF94 speed win or general 3D superiority. The two
same-coordinate MMFF94 energy residuals, term-level adjudication, timeout, and
broader conformer-quality gates remain open.
