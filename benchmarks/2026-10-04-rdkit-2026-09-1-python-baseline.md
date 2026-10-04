# RDKit 2026.09.1 rebaseline — published Python baseline (2026-10-04)

This records the **old-version baseline only**. Official `rdkit==2026.9.1`
was not on PyPI at the recorded check, so no Python 2026.09.1 or nanobind
result, old/new speed ratio, or full rebaseline is claimed. The separate
[npm/WASM lane](2026-10-04-rdkit-2026-09-1-npm-rebaseline.md) already measures
both RDKit.js versions; it cannot substitute for the Python binding.

## Pinned artifacts and conditions

Official PyPI wheels were downloaded with `pip download --only-binary=:all:`
and installed into a fresh CPython 3.13.6 virtual environment with
`pip install --no-index --find-links`. This was macOS arm64. The exposed,
non-sealed 10,000-input corpus is
`validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi` (SHA-256
`f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`).
The 31 queries are pinned by `validation/rdkit_rebaseline_smarts_queries.json`.

| Wheel | SHA-256 |
|---|---|
| `rdkit-2026.3.6-cp313-cp313-macosx_11_0_arm64.whl` | `e16c467cb254a223e59a0cf81358c6b39da15a99d2909170d693e95778fddb41` |
| `chematic-1.0.33-cp313-cp313-macosx_11_0_arm64.whl` | `c541d42f51302039d178d546f02bdbb235925749ed970fb3388f040ef41f0898` |
| `numpy-2.5.3-cp313-cp313-macosx_14_0_arm64.whl` | `f9a2353b37a1a9e78fd82b27ad7e2a32a2d036604d18f02b05e3136c62ca3b09` |
| `pillow-12.3.0-cp313-cp313-macosx_11_0_arm64.whl` | `d69141514cc30b774ceea5e3ed3a6635c8d8a96edf664689b890f4089111fb35` |

The executed RDKit runtime reports `2026.03.6`; the binding callable is a
`Boost.Python.function`. The CheMatic runtime reports `1.0.33`. Wheel paths
in the raw records are host-local metadata; the hashes, versions, corpus and
command catalog are the portable identity checks.

## Complete-row chemistry results

| Check | Published CheMatic 1.0.33 vs RDKit 2026.03.6 |
|---|---:|
| Parsed in both and retained as one row | 10,000 / 10,000 |
| Semantic SMILES round-trip | 10,000 / 10,000 |
| Exact canonical SMILES spelling | 64 / 10,000 |
| Graph-corresponded CIP labels exact | 9,989 / 10,000 |
| Typed CIP abstentions (not exact) | 5 |
| Imine E/Z mismatches (not exact) | 6 |
| Morgan radius 2, 2048-bit exact | 9,999 / 10,000 |
| Typed Fe coordination refusal (not exact) | 1 |
| SMARTS atom-set mismatches | 200 / 310,000 cells |

The five CIP abstentions are 2960, 4419, 4439, 4643 and 4644. The six
imine E/Z mismatches are 1206, 1213, 1214, 1287, 1370 and 1371: CheMatic
labels E where RDKit labels Z. The Morgan refusal is row 8341. This run
corrects a comparator error: `STEREOTRANS` is **not** interchangeable with
the CIP E label. The runner now reads RDKit bond `_CIPCode`; older records
using the cis/trans shortcut must not be reused as E/Z parity evidence.
No canonical-spelling difference is counted as a molecular
identity error: each row's semantic round-trip is checked separately. The
200 SMARTS cells are **Python binding** results; do not merge them with the
194/195 mismatches in the published WASM comparison. The raw JSONL retains
every row, match-set difference, and typed reason. This exposed corpus does
not prove general RDKit compatibility.

## Python binding boundary

The same wheel was probed for `str`, `bytes`, `PathLike`, keyword, invalid
type/SMILES, and list/tuple/generator Morgan inputs. Seven sequential
10,000-row samples were recorded for scalar parse, scalar Morgan and batch
Morgan. Their median whole-batch times on this host were 808.27, 293.12 and
287.65 ms respectively. These are **old-RDKit-only diagnostic timings**, not
a paired comparison or a claim that CheMatic is faster. The new RDKit wheel
must be run with the same input and contract probe before any version-change
or wrapper-overhead interpretation.

## Evidence and remaining gate

- [Wheel/runtime provenance](../validation/results/rdkit-rebaseline-python-provenance-v1.0.33-vs-2026.03.6-2026-10-04.json)
- [10k chemistry summary](../validation/results/rdkit-rebaseline-python-chemistry-v1.0.33-vs-2026.03.6-2026-10-04.json) and [all rows](../validation/results/rdkit-rebaseline-python-chemistry-v1.0.33-vs-2026.03.6-2026-10-04.jsonl.gz). The latter is `gzip -n -9` of the runner's JSONL; the summary records both uncompressed and compressed SHA-256.
- [Python contract and seven timing samples](../validation/results/rdkit-rebaseline-python-contract-vs-2026.03.6-2026-10-04.json)
- Reproduction argv: `validation/rdkit_rebaseline_execution.json`, lane
  `python-2026.03.6-distributed`. Integrity and row-accounting check:
  `python3 scripts/check_rdkit_2026_09_rebaseline_evidence.py`.

When a versioned 2026.09.1 Python wheel exists, hash and run it in a separate
environment, then compare complete rows and classify changes. Native C++
requires its own separately pinned binary or build packet. Neither pending
lane is marked complete by this baseline.
