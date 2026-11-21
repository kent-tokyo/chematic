# Published Python package comparison, 2026-11-21

The measured packages are chematic 1.2.0, 1.2.1 and 1.2.2, RDKit 2026.9.2, and COSMolKit 0.6.2 on macOS arm64 with Python 3.13. Their wheel hashes and individual measurements are retained in the adjacent JSON files. Each block ran all five packages in a rotating order in fresh processes. Medians below describe this host and the specified inputs.

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
