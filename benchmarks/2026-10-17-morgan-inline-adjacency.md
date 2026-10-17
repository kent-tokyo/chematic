# Private inline adjacency: failed performance gate

Baseline `c84c68f297f95562d1d29c75234f2190cd5f9d75`; measured candidate `61d91ccb7091c74d9823dda25035d232b492e31b`. The private RDKit model changes each atom's adjacency from `Vec<usize>` to `SmallVec<[usize;4]>`, spilling above four entries. Bond order, mutation rules, stereo perception and public `RdkitMolView` types are preserved. No new dependency or molecule-output cache.

**Performance decision: rejected.** Both first-use chiral Morgan intervals must have lower bounds above 1 and all 42 source medians must be at least 0.99. The second primary interval fails the improvement rule and eight lanes fail the median floor. Do not rerun or alter the gate after observing this epoch.

This evidence-only change excludes the runtime optimization and its representation-specific tests. The [source experiment PR #1014](https://github.com/kent-tokyo/chematic/pull/1014) was closed after the performance gate failed. The source commits, complete negative results and reproducibility evidence remain available there; none of the historical timing observations is reclassified as a passing result.

## Evidence and scope

[Numeric packet](2026-10-17-morgan-inline-adjacency.json.gz) ([auditor](../scripts/check_morgan_numeric_packet.py)) retains every batch duration/count/output hash, all 78 lanes and raw-journal SHA-256 commitments. Full row outputs, wheel/native binaries, source/runtime manifests, build logs and local audits remain local; this smaller export permits numeric consistency checks but cannot independently establish the chemistry of unexported outputs.

One Apple M4 arm64 host (10 CPUs, 16 GiB, macOS 27.0.1), CPython 3.13.6, RDKit 2026.09.1, COSMolKit 0.5.0. Source wheels only; no published-package or universal superiority claim. Main advanced after the baseline to `c6d10dc0218f190f88fb01f1e43fee0631575c22`; this rejected experiment is not a measurement against that newer main.

- Source wheels: old baseline retained with original empty-target build/source proof; candidate rebuilt from an empty target. Identical release opt-level 3, LTO=true, codegen-units=1, strip=none, one Cargo job, incremental off, Rust 1.97.0, native-inchi enabled, same Python interpreter. ZIP-installed native hashes reverified; runtime/source seals checked before/after both profiles. No earlier elapsed sample reused.
- Corpus: unfiltered 1,000-row windows at offsets 0 and 5,000 of `validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi`. All 24 permutations of four arms once in each window; seeds 20261017+offset; each arm occurs in each position six times.
- First use parses fresh molecules outside each timed full batch; repeated use parses and primes outside timing; parse inclusive includes fresh parsing. Each batch calls the public API once per row. Accumulate independently timed full batches for at least 50 ms (Morgan) or 20 ms (source controls). Preparation and normalization gaps are excluded, so this is not continuous steady-state or process cold-start timing.
- Native return conversion remains in the public API timing boundary; canonical on-bit normalization occurs outside timing. Validate actual outputs after every batch, including every repeat/mode. Full row outputs match across all four arms on both timed windows; source controls match source exactly.
- 10,000 paired bootstrap resamples per lane (fixed comparator-index seeds); ratios are comparator time / candidate time, so >1 favors candidate. Raw unrounded paired medians decide the gate. An external strict win additionally requires lower CI>1 and 24/24 winning blocks. No RSS ranking or cross-host extrapolation.
- Local independent audit: 192 journals, 2880 accounting cells, 187,961 timed full batches, 187,961,000 processed rows. No missing or partial journals. No own build, profiler, source edit or heavy audit during formal timing.
- Correctness: all 14,072 ordinary/stress rows × four R2/R3, chiral/nonchiral fingerprints = 56,288 cells exactly equal between baseline, candidate and RDKit. COSMolKit ordinary 40,000 cells equal; each stress operation retains six refusals and each chiral radius also retains seven bit mismatches. These are not counted as agreement.
- Regression checks: 547 Rust tests passed, 0 failed, 3 existing ignored; all-target smiles Clippy and format/refactor checks pass. Python binding subset 329 passed. Seven measurement-driver tests pass. Frozen original add/remove bodies are whitespace-normalized exact against baseline, with high-degree spill, self-loop, clone, single/bulk removal, remapping and empty-removal checks. Numeric export auditor also rejects corrupted batch accounting.

## Primary source lanes

| Offset | Source / candidate | Paired 95% CI | Wins / 24 |
|---:|---:|---|---:|
| 0 | 1.013527 | [1.009430, 1.017613] | 21 |
| 5000 | 0.992919 | [0.987785, 1.000875] | 8 |

## Every source lane

All eight floor failures and every negative interval are retained below.

| Offset | Operation | Mode | Source / candidate | Paired 95% CI | Wins / 24 |
|---:|---|---|---:|---|---:|
| 0 | morgan2_2048 | first_use | 0.991984 | [0.988092, 0.997124] | 5 |
| 0 | morgan2_2048 | repeated_use | 0.981502 | [0.973383, 0.987211] | 2 |
| 0 | morgan2_2048 | parse_inclusive | 1.001847 | [0.998356, 1.004417] | 15 |
| 0 | morgan2_chiral | first_use | 1.013527 | [1.009430, 1.017613] | 21 |
| 0 | morgan2_chiral | repeated_use | 0.980930 | [0.971601, 0.986384] | 3 |
| 0 | morgan2_chiral | parse_inclusive | 1.011423 | [1.009689, 1.014896] | 23 |
| 0 | morgan3_2048 | first_use | 0.993678 | [0.990554, 0.997706] | 7 |
| 0 | morgan3_2048 | repeated_use | 0.985553 | [0.980043, 0.991274] | 3 |
| 0 | morgan3_2048 | parse_inclusive | 1.000207 | [0.991615, 1.002652] | 12 |
| 0 | tpsa | first_use | 0.983139 | [0.967695, 0.992715] | 3 |
| 0 | tpsa | repeated_use | 1.002830 | [0.999163, 1.007066] | 16 |
| 0 | tpsa | parse_inclusive | 1.005659 | [0.990642, 1.012635] | 14 |
| 0 | labute_asa | first_use | 1.007035 | [0.999678, 1.012003] | 17 |
| 0 | labute_asa | repeated_use | 1.001408 | [0.998028, 1.009763] | 13 |
| 0 | labute_asa | parse_inclusive | 1.004244 | [0.998693, 1.009943] | 16 |
| 0 | num_rings | first_use | 0.987907 | [0.977944, 0.994948] | 5 |
| 0 | num_rings | repeated_use | 0.993301 | [0.987070, 0.999136] | 7 |
| 0 | num_rings | parse_inclusive | 0.995340 | [0.988993, 1.011780] | 10 |
| 0 | maccs | first_use | 1.001425 | [0.999022, 1.004239] | 14 |
| 0 | maccs | repeated_use | 1.001955 | [0.998331, 1.004740] | 16 |
| 0 | maccs | parse_inclusive | 1.002148 | [1.000449, 1.004682] | 17 |
| 5000 | morgan2_2048 | first_use | 0.997318 | [0.995418, 1.004936] | 11 |
| 5000 | morgan2_2048 | repeated_use | 0.988140 | [0.983898, 0.998367] | 6 |
| 5000 | morgan2_2048 | parse_inclusive | 0.995470 | [0.988940, 0.999856] | 7 |
| 5000 | morgan2_chiral | first_use | 0.992919 | [0.987785, 1.000875] | 8 |
| 5000 | morgan2_chiral | repeated_use | 0.986155 | [0.979416, 0.995725] | 5 |
| 5000 | morgan2_chiral | parse_inclusive | 1.000887 | [0.994121, 1.005059] | 12 |
| 5000 | morgan3_2048 | first_use | 0.996982 | [0.994060, 1.003737] | 10 |
| 5000 | morgan3_2048 | repeated_use | 0.989757 | [0.986953, 1.000415] | 8 |
| 5000 | morgan3_2048 | parse_inclusive | 1.003450 | [1.001391, 1.007212] | 18 |
| 5000 | tpsa | first_use | 1.004194 | [0.998131, 1.009007] | 16 |
| 5000 | tpsa | repeated_use | 1.004489 | [0.995291, 1.010085] | 14 |
| 5000 | tpsa | parse_inclusive | 1.002435 | [0.987493, 1.024395] | 13 |
| 5000 | labute_asa | first_use | 0.995922 | [0.990681, 1.005445] | 11 |
| 5000 | labute_asa | repeated_use | 1.006843 | [1.002098, 1.008968] | 18 |
| 5000 | labute_asa | parse_inclusive | 1.002215 | [0.993823, 1.010697] | 14 |
| 5000 | num_rings | first_use | 1.002894 | [0.996788, 1.010346] | 14 |
| 5000 | num_rings | repeated_use | 1.001372 | [0.997804, 1.009234] | 14 |
| 5000 | num_rings | parse_inclusive | 1.004267 | [0.992816, 1.011632] | 13 |
| 5000 | maccs | first_use | 1.000699 | [0.997398, 1.006282] | 13 |
| 5000 | maccs | repeated_use | 1.003076 | [0.997966, 1.010118] | 16 |
| 5000 | maccs | parse_inclusive | 1.001261 | [0.995826, 1.003798] | 13 |

## Every external Morgan lane

An external win does not rescue a failed source adoption gate. In particular first-use chiral Morgan at offset 0 remains slower than COSMolKit; offset 5,000 also misses the 24/24 strict-win condition.

| Offset | Comparator | Operation | Mode | Comparator / candidate | Paired 95% CI | Wins / 24 | Strict win |
|---:|---|---|---|---:|---|---:|---|
| 0 | rdkit | morgan2_2048 | first_use | 3.304433 | [3.284130, 3.327492] | 24 | yes |
| 0 | cosmolkit | morgan2_2048 | first_use | 2.083180 | [2.069346, 2.097416] | 24 | yes |
| 0 | rdkit | morgan2_2048 | repeated_use | 5.309868 | [5.282137, 5.333502] | 24 | yes |
| 0 | cosmolkit | morgan2_2048 | repeated_use | 3.216857 | [3.199774, 3.241691] | 24 | yes |
| 0 | rdkit | morgan2_2048 | parse_inclusive | 8.659085 | [8.634145, 8.683611] | 24 | yes |
| 0 | cosmolkit | morgan2_2048 | parse_inclusive | 7.604886 | [7.580635, 7.624546] | 24 | yes |
| 0 | rdkit | morgan2_chiral | first_use | 1.308470 | [1.295722, 1.313640] | 24 | yes |
| 0 | cosmolkit | morgan2_chiral | first_use | 0.829285 | [0.825749, 0.832607] | 0 | no |
| 0 | rdkit | morgan2_chiral | repeated_use | 5.210014 | [5.177800, 5.237111] | 24 | yes |
| 0 | cosmolkit | morgan2_chiral | repeated_use | 3.168404 | [3.127564, 3.205150] | 24 | yes |
| 0 | rdkit | morgan2_chiral | parse_inclusive | 4.540472 | [4.523944, 4.551168] | 24 | yes |
| 0 | cosmolkit | morgan2_chiral | parse_inclusive | 3.993054 | [3.982285, 4.013794] | 24 | yes |
| 0 | rdkit | morgan3_2048 | first_use | 3.433136 | [3.421191, 3.446713] | 24 | yes |
| 0 | cosmolkit | morgan3_2048 | first_use | 2.240135 | [2.228068, 2.250686] | 24 | yes |
| 0 | rdkit | morgan3_2048 | repeated_use | 4.963280 | [4.937500, 4.990371] | 24 | yes |
| 0 | cosmolkit | morgan3_2048 | repeated_use | 3.120455 | [3.105211, 3.133631] | 24 | yes |
| 0 | rdkit | morgan3_2048 | parse_inclusive | 8.084817 | [8.035051, 8.115111] | 24 | yes |
| 0 | cosmolkit | morgan3_2048 | parse_inclusive | 7.039901 | [7.019973, 7.055891] | 24 | yes |
| 5000 | rdkit | morgan2_2048 | first_use | 1.602271 | [1.594506, 1.614703] | 24 | yes |
| 5000 | cosmolkit | morgan2_2048 | first_use | 1.003052 | [0.995873, 1.007608] | 15 | no |
| 5000 | rdkit | morgan2_2048 | repeated_use | 5.072073 | [5.014813, 5.101004] | 24 | yes |
| 5000 | cosmolkit | morgan2_2048 | repeated_use | 3.008237 | [2.967766, 3.047394] | 24 | yes |
| 5000 | rdkit | morgan2_2048 | parse_inclusive | 4.627948 | [4.595165, 4.650411] | 24 | yes |
| 5000 | cosmolkit | morgan2_2048 | parse_inclusive | 4.455830 | [4.441843, 4.482588] | 24 | yes |
| 5000 | rdkit | morgan2_chiral | first_use | 1.629865 | [1.606947, 1.640573] | 24 | yes |
| 5000 | cosmolkit | morgan2_chiral | first_use | 1.017537 | [1.010288, 1.026289] | 21 | no |
| 5000 | rdkit | morgan2_chiral | repeated_use | 5.036167 | [4.993925, 5.084774] | 24 | yes |
| 5000 | cosmolkit | morgan2_chiral | repeated_use | 2.986397 | [2.963902, 3.004101] | 24 | yes |
| 5000 | rdkit | morgan2_chiral | parse_inclusive | 4.640696 | [4.605481, 4.665608] | 24 | yes |
| 5000 | cosmolkit | morgan2_chiral | parse_inclusive | 4.476007 | [4.457360, 4.496151] | 24 | yes |
| 5000 | rdkit | morgan3_2048 | first_use | 1.815069 | [1.807615, 1.828009] | 24 | yes |
| 5000 | cosmolkit | morgan3_2048 | first_use | 1.176760 | [1.165832, 1.179761] | 24 | yes |
| 5000 | rdkit | morgan3_2048 | repeated_use | 4.832287 | [4.765927, 4.863940] | 24 | yes |
| 5000 | cosmolkit | morgan3_2048 | repeated_use | 2.996025 | [2.977826, 3.012395] | 24 | yes |
| 5000 | rdkit | morgan3_2048 | parse_inclusive | 4.579369 | [4.560603, 4.614404] | 24 | yes |
| 5000 | cosmolkit | morgan3_2048 | parse_inclusive | 4.358562 | [4.353114, 4.377636] | 24 | yes |

Recompute numeric consistency without chemistry SDKs:

```sh
python3 scripts/check_morgan_numeric_packet.py benchmarks/2026-10-17-morgan-inline-adjacency.json.gz
```

The complete local packet is retained independently of the public numeric export. Only the measurement record and reusable benchmark/audit scripts are proposed for main. The failed runtime optimization remains excluded. Source correctness checks are historical evidence for the named candidate, not new performance claims for current main.
