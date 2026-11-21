# Published Python package comparison, 2026-11-21

## Current COSMolKit 0.6.3 rerun

After COSMolKit 0.6.3 became the latest stable PyPI release, all five arms were rerun in 20 rotating fresh-process blocks on the same macOS arm64/Python 3.13 host. The other arms are the published chematic 1.2.0, 1.2.1 and 1.2.2 wheels and RDKit 2026.9.2. The [registry record](../validation/results/published-benchmark-wheel-registry-cos063-2026-11-21.json) verifies all five wheel hashes against PyPI and records that RDKit 2026.9.2 and COSMolKit 0.6.3 were the latest stable versions when checked.

All five packages matched TPSA, ring count and Morgan-2/2048 on the same 1,000 pinned inputs under the recorded tolerance. The fresh-molecule TPSA column follows one warm-up call on a separate molecule. Peak memory covers imports and all worker operations, not TPSA alone.

| Package | TPSA with parse, ms | TPSA on fresh molecules, ms | TPSA after preparation, ms | Process peak, MiB |
|---|---:|---:|---:|---:|
| chematic 1.2.0 | 5.168 | 1.580 | 0.045 | 54.5 |
| chematic 1.2.1 | 5.230 | 1.575 | 0.045 | 54.6 |
| chematic 1.2.2 | 4.968 | 1.494 | 0.042 | 54.6 |
| RDKit 2026.9.2 | 76.240 | 2.553 | 0.478 | 97.6 |
| COSMolKit 0.6.3 | 58.720 | 0.927 | 0.082 | 122.5 |

The separate MCS run uses the same five published wheels and 20 rotating blocks. Every timed call returns SMARTS, including concrete-molecule-to-SMARTS conversion for chematic 1.2.0 and 1.2.1. The table reports the benzene/toluene case, whose query predicates and input matches agree for all five. First parse covers all five cases, and RDKit's parse includes lazy initialization.

| Package | Import, ms | First parse of five cases, ms | First MCS, ms | Repeated MCS, ms | Process peak, MiB |
|---|---:|---:|---:|---:|---:|
| chematic 1.2.0 | 85.91 | 0.090 | 0.1317 | 0.0345 | 44.2 |
| chematic 1.2.1 | 85.68 | 0.089 | 0.3982 | 0.2917 | 44.2 |
| chematic 1.2.2 | 85.75 | 0.087 | 0.0539 | 0.0100 | 43.6 |
| RDKit 2026.9.2 | 53.12 | 668.955 | 0.1429 | 0.0499 | 41.5 |
| COSMolKit 0.6.3 | 21.67 | 0.780 | 0.1125 | 0.0320 | 43.5 |

On the five cases, the query-predicate and input-match checks pass for cases 0 and 2 in 1.2.0, cases 0–3 in 1.2.1, and all five in 1.2.2, RDKit and COSMolKit. Earlier releases fail to serialize case 4 as SMARTS; invalid or unequal results are excluded from paired ratios. These numbers are medians for a small case set on a shared host; the earlier run below is a separate measurement and its timings should not be combined with this one. The full samples, load observations, wheel hashes, case results and paired ratios are in [descriptor](2026-11-21-v122-published-speed-cos063.json) and [uniform MCS](2026-11-21-v122-published-mcs-uniform-cos063.json).

## Earlier COSMolKit 0.6.2 run

The measured packages are chematic 1.2.0, 1.2.1 and 1.2.2, RDKit 2026.9.2, and COSMolKit 0.6.2 on macOS arm64 with Python 3.13. Their wheel hashes and individual measurements are retained in the adjacent JSON files. Each block ran all five packages in a rotating order in fresh processes. Medians below describe this host and the specified inputs.

The [registry verification](../validation/results/published-benchmark-wheel-registry-verification-2026-11-21.json) checks the exact filename and SHA-256 of all five measured wheels against their PyPI JSON entries and both timing records. It also rechecks the complete 26-file v1.2.2 PyPI hash set against the release evidence. The dated record retains the registry URLs and response hashes; CI checks the retained records for consistency without making live registry requests.

For 1,000 rows of the pinned RDKit.js corpus, all packages produced identical TPSA, ring-count and Morgan-2/2048 results under the comparison tolerance recorded in the raw data. Twenty blocks completed. A disk-full interruption followed block 15; package extension hashes and output values were checked again before the last five blocks. Load averages changed during the run.

| Package | TPSA with parse, ms | TPSA on fresh molecules, ms | TPSA after preparation, ms | Process peak, MiB |
|---|---:|---:|---:|---:|
| chematic 1.2.0 | 5.370 | 2.340 | 0.046 | 54.7 |
| chematic 1.2.1 | 5.741 | 2.270 | 0.047 | 54.7 |
| chematic 1.2.2 | 5.269 | 1.766 | 0.045 | 54.7 |
| RDKit 2026.9.2 | 70.979 | 2.492 | 0.425 | 97.5 |
| COSMolKit 0.6.2 | 57.208 | 1.069 | 0.082 | 122.5 |

The fresh-molecule TPSA column runs after a warm-up call on a separate molecule. The peak covers imports and all measured operations in the worker, including other descriptors and fingerprints. See `2026-11-21-v122-published-speed-frontier.json` for measurements, agreement checks, corpus hash and wheel identities.

The separate five-case MCS probe measures process import, first parsing of the entire case set, the first search on benzene/toluene, a repeated search on those prepared molecules, and whole-process peak RSS. These are 20-block medians.

| Package | Import, ms | First parse of five cases, ms | First MCS, ms | Repeated MCS, ms | Process peak, MiB |
|---|---:|---:|---:|---:|---:|
| chematic 1.2.0 | 56.23 | 0.060 | 0.0732 | 0.0185 | 41.6 |
| chematic 1.2.1 | 55.53 | 0.058 | 0.2390 | 0.1857 | 41.6 |
| chematic 1.2.2 | 56.26 | 0.064 | 0.0651 | 0.0126 | 41.5 |
| RDKit 2026.9.2 | 33.71 | 412.594 | 0.0972 | 0.0322 | 40.0 |
| COSMolKit 0.6.2 | 13.99 | 0.505 | 0.0777 | 0.0196 | 41.8 |

RDKit's first parse includes lazy initialization. The MCS case set is small. The benzene/toluene results all match both inputs and contain six atoms and six bonds. For the triangle/chain case, 1.2.0 returned two atoms while the other packages returned three. For the fused aromatic case, all three chematic versions returned a concrete SMILES that RDKit could not parse; RDKit and COSMolKit returned eight-atom SMARTS queries. Concrete molecule and SMARTS outputs have different expressive power. The raw `2026-11-21-v122-published-mcs.json` records per-case times and query-predicate validation after conversion to SMARTS.

An additional 20-block run uses 1.2.2's `find_mcs_smarts_checked` with the bond-first objective. The earlier releases use their available concrete-molecule API. RDKit and COSMolKit use SMARTS. For the five cases, the 1.2.2, RDKit and COSMolKit query graphs have matching atom and bond predicates and match every input. The earlier 1.2.0 result meets that check on cases 0 and 2, and 1.2.1 meets it on cases 0–3 after its concrete molecule is converted to SMARTS. The fused aromatic concrete result remains unparseable for both earlier releases. Times below are medians for the benzene/toluene case from this separate run.

| Package and MCS output | Import, ms | First parse of five cases, ms | First MCS, ms | Repeated MCS, ms | Process peak, MiB |
|---|---:|---:|---:|---:|---:|
| chematic 1.2.0, molecule | 51.86 | 0.059 | 0.0732 | 0.0181 | 43.7 |
| chematic 1.2.1, molecule | 52.65 | 0.060 | 0.2374 | 0.1819 | 43.7 |
| chematic 1.2.2, SMARTS | 52.27 | 0.060 | 0.0349 | 0.0063 | 43.3 |
| RDKit 2026.9.2, SMARTS | 32.73 | 408.340 | 0.0879 | 0.0317 | 41.2 |
| COSMolKit 0.6.2, SMARTS | 13.34 | 0.490 | 0.0755 | 0.0191 | 43.1 |

The raw `2026-11-21-v122-published-mcs-smarts.json` includes the five wheel hashes, validated installed files, per-case measurements and paired ratios. Timings include creation of each library's native result; conversion of a concrete molecule to SMARTS and query comparison happen after timing.

A third 20-block run uses SMARTS as the timed output for all five packages. For 1.2.0 and 1.2.1, the timed call includes converting the concrete MCS molecule with `rdkit_smarts()`. The 1.2.2 lane uses `find_mcs_smarts_checked` with the bond-first objective. RDKit and COSMolKit return SMARTS directly. Post-timing validation compares atom and bond query predicates and input matches. Successful query predicates agree for 1.2.0 on cases 0 and 2, for 1.2.1 on cases 0–3, and for 1.2.2, RDKit and COSMolKit on all five cases. The older versions return a serialization error on case 4. The timing table uses case 0, benzene/toluene, where all five queries pass.

| Package | Import, ms | First parse of five cases, ms | First MCS, ms | Repeated MCS, ms | Process peak, MiB |
|---|---:|---:|---:|---:|---:|
| chematic 1.2.0 | 52.30 | 0.057 | 0.0871 | 0.0217 | 44.0 |
| chematic 1.2.1 | 51.75 | 0.055 | 0.2560 | 0.1839 | 44.0 |
| chematic 1.2.2 | 52.22 | 0.060 | 0.0364 | 0.0062 | 43.5 |
| RDKit 2026.9.2 | 32.87 | 408.218 | 0.0885 | 0.0313 | 41.4 |
| COSMolKit 0.6.2 | 13.18 | 0.488 | 0.0735 | 0.0191 | 43.2 |

The full measurements, wheel verification, validation outcomes, paired ratios and excluded cases are in `2026-11-21-v122-published-mcs-uniform.json`. The five cases are a limited input set.
