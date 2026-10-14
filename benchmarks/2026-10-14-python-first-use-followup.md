# Python first-use perception follow-up — 2026-10-14

The additional-slice read-plus-ring-count workflow improves 1.076x
[1.067, 1.086] over the previous candidate, winning all 20 blocks. All 27
eligible parse-inclusive comparator lanes beat RDKit and COSMolKit. First-use
TPSA/ASA and initial ring perception still have substantial deficits; complete
speed superiority is not established.

This packet compares a narrow perception candidate with the previous
[scalar/stereo candidate](2026-10-14-python-scalar-stereo.md), RDKit 2026.09.1
and COSMolKit 0.5.0-rc.22. It retains the same seven operations, both fixed
1,000-row slices, and separate first-use, repeated-use and parse-inclusive
boundaries. A general speed-superiority claim requires more than this packet.

## Implementation and rejected development

The production change is limited to three files:

- The RDKit model gate tracks whether its existing atom scan encountered an
  aromatic atom. With none, there is no aromatic component to inspect, so the
  separate Hückel traversal and its allocations are unnecessary. Valence,
  radical, element, cleanup and cyclic-triple checks remain in place.
- The aromaticity identity decision reuses an existing forest proof, or proves
  that a graph is entirely saturated, before allocating cyclic-component
  scratch data. Both shortcuts require no aromatic atom flags and a passing
  model gate. No ring computation or descriptor moves into SMILES parsing.
- Rank-two ring-system traversal uses inline neighbor, branch and path-length
  buffers, retaining a spilling fallback for larger neighborhoods. Ring
  selection and counting rules are unchanged.

Three added regression tests check that closing a conjugated ring or editing
saturated bonds invalidates the earlier identity result, and that unusual
radical, cleanup and coordination inputs retain their model guards.

A broader conversion of scratch vectors to stack buffers was tried first.
Exploratory results were mixed or worse, so it was removed before the formal
comparison. Its source patch and raw observations remain in the
[parity/artifact record](2026-10-14-python-first-use-followup-parity.json).
The narrower candidate was also probed before formal timing; those observations
are retained separately and are not counted as formal blocks.

## Measurement protocol

Each slice has exactly 1,000 nonempty SMILES rows, at offsets 0 and 5000, from
`validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi`, SHA256
`f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`.
Both offsets were inherited unchanged from the earlier packet. These are
exposed development corpora, not sealed holdouts.

Twenty sequential fresh-process blocks rotate four arms evenly. Timed loops
retain native return objects; normalization and comparisons occur afterward.
First use excludes parsing and includes lazy work not done by the parser.
Repeated use performs one untimed call per fresh molecule before timing.
Parse-inclusive timing actually parses and calls the operation in one loop.
Eager/lazy preparation differences remain visible, and all three boundaries
are reported. No input or block is removed to obtain a speed win.

Tables use paired comparator / candidate median ratios and deterministic
10,000-resample paired bootstrap 95% intervals. A strict win requires equal
outputs, a lower confidence bound above 1 and every one of 20 blocks faster.
Output differences suppress ranking. Whole-process peak RSS includes imports,
parsing, preparation and all seven operations; it is not operation-only or
streaming memory evidence.

## Formal results

Bold cells pass the strict all-block win gate. Other cells remain measured results, not strict wins.

### First use

| Operation | RDKit, slice 0 | COSMolKit, slice 0 | RDKit, slice 5000 | COSMolKit, slice 5000 |
|---|---|---|---|---|
| tpsa | 1.460x [1.427, 1.488]; 19/20 | 0.487x [0.481, 0.505]; 0/20 | 0.317x [0.315, 0.335]; 0/20 | 0.085x [0.084, 0.096]; 0/20 |
| labute_asa | 0.705x [0.687, 0.714]; 1/20 | 0.328x [0.325, 0.334]; 0/20 | 0.156x [0.152, 0.161]; 0/20 | 0.067x [0.063, 0.073]; 0/20 |
| num_rings | 0.021x [0.021, 0.022]; 0/20 | Not ranked: outputs differ | 0.096x [0.095, 0.104]; 0/20 | 0.081x [0.074, 0.100]; 0/20 |
| morgan2_2048 | **3.160x [3.124, 3.193]; 20/20** | **2.045x [2.030, 2.056]; 20/20** | **1.402x [1.390, 1.415]; 20/20** | 0.910x [0.901, 0.936]; 1/20 |
| morgan2_chiral | **1.265x [1.255, 1.272]; 20/20** | 0.830x [0.819, 0.840]; 0/20 | **1.437x [1.406, 1.450]; 20/20** | 0.926x [0.918, 0.949]; 1/20 |
| morgan3_2048 | **3.354x [3.325, 3.370]; 20/20** | **2.233x [2.202, 2.271]; 20/20** | **1.623x [1.595, 1.638]; 20/20** | 1.087x [1.067, 1.122]; 19/20 |
| maccs | **1.885x [1.878, 1.892]; 20/20** | **2.169x [2.158, 2.175]; 20/20** | **1.629x [1.615, 1.641]; 20/20** | **1.756x [1.740, 1.764]; 20/20** |

### Repeated use

| Operation | RDKit, slice 0 | COSMolKit, slice 0 | RDKit, slice 5000 | COSMolKit, slice 5000 |
|---|---|---|---|---|
| tpsa | **10.006x [9.864, 10.104]; 20/20** | **1.671x [1.633, 1.765]; 20/20** | **9.871x [9.597, 9.944]; 20/20** | 1.378x [1.356, 1.533]; 19/20 |
| labute_asa | **2.724x [2.683, 2.761]; 20/20** | **1.437x [1.384, 1.531]; 20/20** | **2.672x [2.641, 2.716]; 20/20** | 1.210x [1.156, 1.291]; 17/20 |
| num_rings | **1.857x [1.802, 1.916]; 20/20** | Not ranked: outputs differ | **1.731x [1.686, 1.809]; 20/20** | 1.070x [1.033, 1.105]; 17/20 |
| morgan2_2048 | **5.227x [5.204, 5.275]; 20/20** | **3.226x [3.188, 3.267]; 20/20** | **4.824x [4.769, 4.901]; 20/20** | **2.968x [2.907, 3.039]; 20/20** |
| morgan2_chiral | **5.031x [4.981, 5.075]; 20/20** | **3.093x [3.072, 3.132]; 20/20** | **4.795x [4.695, 4.847]; 20/20** | **2.946x [2.907, 2.989]; 20/20** |
| morgan3_2048 | **4.853x [4.812, 4.889]; 20/20** | **3.130x [3.095, 3.179]; 20/20** | **4.556x [4.468, 4.638]; 20/20** | **2.888x [2.858, 2.938]; 20/20** |
| maccs | **1.950x [1.947, 1.964]; 20/20** | **2.254x [2.248, 2.261]; 20/20** | **1.724x [1.712, 1.742]; 20/20** | **1.860x [1.848, 1.874]; 20/20** |

### Parse inclusive

| Operation | RDKit, slice 0 | COSMolKit, slice 0 | RDKit, slice 5000 | COSMolKit, slice 5000 |
|---|---|---|---|---|
| tpsa | **12.953x [12.898, 13.130]; 20/20** | **12.187x [12.033, 12.249]; 20/20** | **4.083x [4.041, 4.118]; 20/20** | **4.220x [4.197, 4.290]; 20/20** |
| labute_asa | **12.306x [12.241, 12.401]; 20/20** | **11.647x [11.589, 11.698]; 20/20** | **3.902x [3.880, 3.938]; 20/20** | **4.191x [4.155, 4.234]; 20/20** |
| num_rings | **8.821x [8.731, 8.872]; 20/20** | Not ranked: outputs differ | **11.636x [11.504, 11.707]; 20/20** | **12.613x [12.539, 12.778]; 20/20** |
| morgan2_2048 | **8.645x [8.538, 8.712]; 20/20** | **7.625x [7.499, 7.676]; 20/20** | **4.269x [4.243, 4.301]; 20/20** | **4.116x [4.102, 4.186]; 20/20** |
| morgan2_chiral | **4.463x [4.448, 4.526]; 20/20** | **3.949x [3.898, 3.987]; 20/20** | **4.280x [4.244, 4.295]; 20/20** | **4.121x [4.097, 4.185]; 20/20** |
| morgan3_2048 | **7.988x [7.894, 8.050]; 20/20** | **6.993x [6.940, 7.034]; 20/20** | **4.196x [4.130, 4.255]; 20/20** | **4.022x [3.944, 4.068]; 20/20** |
| maccs | **2.107x [2.100, 2.118]; 20/20** | **2.382x [2.371, 2.393]; 20/20** | **1.795x [1.787, 1.803]; 20/20** | **1.929x [1.920, 1.947]; 20/20** |

### Incremental change from the previous candidate

Only additional-slice parse-inclusive ring count passes the incremental strict gate. First-use ring-count medians improve 4.0% and 3.4%, with favorable intervals, but win only 19/20 and 18/20 blocks. Original-slice TPSA/ASA medians are approximately 0.8%/1.2% slower and their intervals reach or cross parity. First-use plain R2 on the additional slice is 0.6% slower [0.984, 0.999]. These mixed results are retained.

| Operation / mode | Previous / candidate, slice 0 | Previous / candidate, slice 5000 |
|---|---|---|
| tpsa/first_use | 0.992x [0.981, 1.001]; 6/20 | 1.011x [1.000, 1.017]; 14/20 |
| tpsa/repeated_use | 0.997x [0.975, 1.014]; 9/20 | 1.008x [0.998, 1.029]; 13/20 |
| tpsa/parse_inclusive | 0.995x [0.989, 1.002]; 6/20 | 1.029x [1.022, 1.037]; 17/20 |
| labute_asa/first_use | 0.988x [0.974, 1.000]; 6/20 | 1.011x [1.005, 1.021]; 16/20 |
| labute_asa/repeated_use | 0.992x [0.986, 1.011]; 8/20 | 1.004x [0.983, 1.028]; 10/20 |
| labute_asa/parse_inclusive | 0.999x [0.989, 1.004]; 8/20 | 1.033x [1.023, 1.045]; 18/20 |
| num_rings/first_use | 1.040x [1.030, 1.048]; 19/20 | 1.034x [1.030, 1.041]; 18/20 |
| num_rings/repeated_use | 1.019x [1.002, 1.036]; 14/20 | 1.002x [0.970, 1.023]; 10/20 |
| num_rings/parse_inclusive | 1.025x [1.018, 1.030]; 18/20 | **1.076x [1.067, 1.086]; 20/20** |

### Memory boundary

| Slice | Candidate median whole-process peak MiB | RDKit MiB | COSMolKit MiB |
|---|---:|---:|---:|
| 0 | 46.60 | 106.56 | 105.33 |
| 5000 | 41.16 | 88.38 | 82.84 |

All 27 eligible parse-inclusive comparator lanes pass the strict gate; original-slice COSMolKit ring count remains unranked because six records differ. This maintains the earlier workflow advantage and does not establish first-use or toolkit-wide superiority. COSMolKit remains faster for first-use TPSA/ASA on both slices, first-use chiral Morgan on both slices, and plain R2 on the additional slice. Initial ring counting remains slower than RDKit on both slices and COSMolKit on the comparable additional slice. Short prepared scalar all-block gates on the additional slice remain open.


## Correctness and provenance

The exact previous source wheel is the baseline, SHA256
`b3a7b1a050f6517420a4f3e627a0f19ab407b310914c2848ad2461af87cd1408`.
The new candidate wheel SHA256 is `8952cc18021c845ba2c834187b4148f07eb0011685cdf0ff1b2cfad59703630f`.
Base commit is `ffe682c32cb4c76de0d309649427202a7d8e956d`. Apply the prior
speed-frontier and scalar/stereo patches, then this candidate patch, in an
isolated checkout. The recorded source hashes distinguish this candidate from
other wheels carrying source version 1.1.3. Unrelated dirty workspace work,
including the distance-matrix, Open Babel, Python API and 3D changes, is excluded.

The release wheel uses Rust 1.97.0, opt-level 3, LTO, one codegen unit and
native-inchi, on macOS arm64 / CPython 3.13.6 / 10 logical CPUs. Changed crate
dependencies were rebuilt before linking; the runner verifies the installed
extension against the exact supplied wheel in every block. A wrapper verifies
frozen build sources, changed workspace sources, benchmark tools, both wheels,
the manifest and comparator native libraries before and after timing. Compiler
load was allowed to finish before formal timing began. This is source-wheel
evidence, not a public release, hosted CI, commit or merge.

Normal 10,000 plus unusual 4,072 inputs produce 112,576 first/repeated
descriptor outcomes and 56,288 Morgan outcomes byte-identical to the previous
candidate. Float bits, errors and nonfinite outcomes are retained; matching
nonfinite outcomes are not successful chemical agreement. Existing lossless
dumps are reused after their hashes were verified. The broader RDKit residuals
reported in the scalar/stereo packet remain unchanged; retaining the baseline
does not erase them.

Isolated release tests for core, perception, chemistry, fingerprints and SMILES:
2,371 passed, zero failed, ten existing ignored cases. The three new edit/guard
tests pass. Release Clippy for both changed crates and all targets passes with
warnings denied. The selected Python molecule, descriptor, RDKit, fingerprint,
pickle and cross-binding suites pass 577 tests. Initial isolated test builds
needed additional pinned fixtures restored from the base commit before the
successful run. No production change was made to resolve those missing files.

Five unchanged benchmark-gate tests pass. Changed Rust files pass rustfmt,
and the workspace diff has no whitespace errors. All raw paired summaries
were independently recalculated and every retained artifact hash was verified.

- [Original-slice timing blocks](2026-10-14-python-first-use-followup-slice0-paired20.json)
- [Additional-slice timing blocks](2026-10-14-python-first-use-followup-slice5000-paired20.json)
- [Source, wheel and comparator manifest](2026-10-14-python-first-use-followup-source.json)
- [Exact parity, rejected development and validation logs](2026-10-14-python-first-use-followup-parity.json)
- [Reconstructable candidate patch](2026-10-14-python-first-use-followup-candidate.patch)
- [Unchanged benchmark runner](../scripts/bench_python_speed_frontier.py)

## Reproduce and remaining scope

Use the exact previous wheel as baseline, reconstruct the new candidate from
the three recorded patches, and build with the recorded settings. Install the
two wheels separately. Run `bench_python_speed_frontier.py` with the same
corpus, offsets, `--limit 1000 --blocks 20`, wheel paths, installed paths and
this source manifest. The parity record retains the wrapper used to add source
and comparator checks around the unchanged runner. Reproduce descriptor and
Morgan dumps using the existing frontier verification scripts.

Toolkit-wide superiority remains unproven. The original first-use, repeated-use
and parse-inclusive results are all retained. Broader common operations,
reactions and 3D workflows, malformed-input behavior, memory-specific work,
sealed corpora, other hosts and exact published artifacts require separate
equivalent-output comparisons.
