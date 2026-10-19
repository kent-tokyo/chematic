# COSMolKit speed audit: strict rotors and Layered fingerprints — 2026-10-19

Layered fingerprints were slower than COSMolKit in both operation-only and
parse-plus-operation screening. Inline path bookkeeping removes this loss on
the fixed 1,000-row lane: parse-inclusive source speed improves
2.149x and is
1.762x as fast as COSMolKit.
Strict rotor first use improves 1.025x,
but remains slower than COSMolKit (0.846x).
MQN first use improves 1.015x against source;
its external comparison is unranked because of existing output differences.

These are source wheels on one macOS arm64 host, not published-package,
cross-platform, 41-operation, or toolkit-wide superiority results.

## Source identity and protocol

- Baseline main: `8b953768882a003fddca22f099f3e17583213718`.
- Frozen candidate: `23f774e1b087e500d0a3c325517b1a4a5bc3a760`.
- CPython 3.13; RDKit `2026.09.1`;
  COSMolKit `0.5.0`.
- Both source wheels used opt-level 3, full LTO, one codegen unit,
  `--release --features native-inchi`, the same compiler and build flags.
  The baseline started from an empty target; the candidate reused Cargo's
  dependency target with incremental compilation disabled. Wheel/native hashes
  and source-file digest are retained in the packet. No binary repair was used.
- Corpus: first 1,000 nonempty rows, offset zero, of the pinned
  `validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi`;
  complete-file SHA-256 `f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`.
- Before timing, three source-primary conditions and all eight operations were
  registered. Eight alternating-source pilot blocks were diagnostic only.
  The formal run used fresh processes and all 24 four-arm permutations once,
  shuffled with seed 20261019. Every observation is retained, including losses.
- First use starts with fresh parsed molecules. Repeated use primes separate
  fresh objects once outside timing. Parse-inclusive calls parse and evaluate
  inside the same timed loop. Native return objects are retained; bit/vector
  normalization and every-row agreement checking occur outside timing.
- Paired ratios are comparator time divided by candidate time. The table shows
  median, deterministic 10,000-resample paired-bootstrap 95% interval, and
  faster-block count. A strict external win requires interval lower bound > 1
  and all 24 blocks faster. Source adoption requires all three registered
  primary interval lower bounds > 1; every secondary lane remains visible.

## Implementation and output preservation

Strict rotor predicates test the existing degree and methyl exclusions before
recursive multi-atom exclusions. The conjunctions, endpoint order, perceived
view, and explicit-H filtering are unchanged. MQN and QED use this same rotor
calculation. A frozen old SMARTS/pair-filter oracle compares every endpoint
pair on all 14,072 normal/stress inputs and their Kekule views.

Layered hashing uses inline small vectors for path-local degrees, atom indices,
bond-neighbor counts, and the six layer buffers instead of per-feature hash
maps, B-trees and heap vectors. The default seven-bond paths fit inline;
configured longer paths spill dynamically. Enumeration, sorted layer values,
hash folding, bit numbering, atom counters, and resource limits are unchanged.
A frozen old bookkeeping/hash implementation checks exact seeds and atom
membership, including paths exceeding all inline capacities.

The COSMolKit benchmark adapter also resolves five final-0.5 fingerprint method
names outside timing, while retaining the rc22 names as fallbacks. Before that
repair, the obsolete calls raised AttributeError on every row and were excluded
from rankings. They were never valid fast fingerprint results.

Full source-wheel audits retain all 14,072 rows: rotor counts, complete 42-entry
MQNs, complete Layered on-bit sets, and QED. All **56,288 values exactly match
the source baseline**, including QED floating-point values.
On the full normal 10k corpus both chematic and COSMolKit match RDKit on every
rotor count, Layered bit set and QED value (QED oracle tolerance 1e-9).
MQNs in both engines differ from RDKit on normal rows 9, 23, 28, 29, 30 and 34;
no rows are removed to enable a speed claim.

Stress residuals are retained: chematic has four MQN and 12 Layered differences;
RDKit refuses QED on one row. COSMolKit cannot parse six stress rows and
has four further MQN differences; QED has 279 source-unavailable states and
one both-unavailable state. These are correctness/availability boundaries,
not timing wins. The packet includes every non-equal row index and dump digest;
raw molecular streams, wheels and process-local paths remain private.

## Secondary regressions and adoption tradeoff

The three pre-registered primary confidence bounds pass. This adoption rule
is not a claim of no regressions. QED first-use and parse-inclusive source
speed ratios are 0.988x and 0.986x; MACCS first/repeated/parse-inclusive ratios
are about 0.989–0.990x. Their intervals are below 1. Chiral Morgan first/repeated
source ratios are about 0.992–0.993x. Repeated QED has a 0.981x median but its
interval includes 1, so its apparent 1.9% loss is inconclusive. These secondary
losses are retained explicitly as a tradeoff for the roughly 2.15x Layered
pipeline gain and smaller rotor/MQN gains. Existing QED and MACCS external
wins remain; chiral Morgan first-use is tied with COSMolKit on this 1k slice.
No claim that every operation improved is made.

## All formal lanes

| Operation/boundary | Source baseline / candidate | RDKit / candidate | COSMolKit / candidate |
|---|---:|---:|---:|
| rotatable_bonds/first_use | 1.025x [1.023, 1.030]; 24/24 | 1.406x [1.395, 1.416]; 24/24 | 0.846x [0.843, 0.848]; 0/24 |
| rotatable_bonds/repeated_use | 1.030x [1.027, 1.032]; 24/24 | 1.705x [1.696, 1.715]; 24/24 | 1.022x [1.021, 1.029]; 24/24 |
| rotatable_bonds/parse_inclusive | 1.016x [1.012, 1.022]; 24/24 | 2.943x [2.936, 2.959]; 24/24 | 2.382x [2.376, 2.389]; 24/24 |
| mqn/first_use | 1.015x [1.009, 1.018]; 24/24 | unranked | unranked |
| mqn/repeated_use | 1.034x [1.024, 1.038]; 24/24 | unranked | unranked |
| mqn/parse_inclusive | 1.015x [1.012, 1.019]; 23/24 | unranked | unranked |
| fp_layered/first_use | 2.162x [2.148, 2.176]; 24/24 | 1.575x [1.568, 1.589]; 24/24 | 1.683x [1.674, 1.691]; 24/24 |
| fp_layered/repeated_use | 2.171x [2.165, 2.191]; 24/24 | 1.591x [1.585, 1.602]; 24/24 | 1.704x [1.688, 1.711]; 24/24 |
| fp_layered/parse_inclusive | 2.149x [2.137, 2.152]; 24/24 | 1.673x [1.657, 1.682]; 24/24 | 1.762x [1.754, 1.781]; 24/24 |
| qed/first_use | 0.988x [0.986, 0.992]; 1/24 | 3.297x [3.287, 3.304]; 24/24 | 2.787x [2.780, 2.793]; 24/24 |
| qed/repeated_use | 0.981x [0.949, 1.020]; 10/24 | 13249.554x [12857.109, 13507.274]; 24/24 | 11252.376x [10860.065, 11452.728]; 24/24 |
| qed/parse_inclusive | 0.986x [0.985, 0.988]; 0/24 | 3.492x [3.486, 3.500]; 24/24 | 2.984x [2.977, 2.992]; 24/24 |
| tpsa/first_use | 1.010x [1.004, 1.021]; 17/24 | 8.033x [7.891, 8.081]; 24/24 | 2.694x [2.660, 2.721]; 24/24 |
| tpsa/repeated_use | 0.999x [0.986, 1.012]; 11/24 | 10.057x [9.886, 10.197]; 24/24 | 1.595x [1.566, 1.649]; 24/24 |
| tpsa/parse_inclusive | 0.994x [0.991, 0.998]; 5/24 | 10.058x [9.996, 10.115]; 24/24 | 9.360x [9.334, 9.377]; 24/24 |
| labute_asa/first_use | 1.009x [1.007, 1.019]; 21/24 | 3.423x [3.349, 3.468]; 24/24 | 1.575x [1.551, 1.619]; 24/24 |
| labute_asa/repeated_use | 0.992x [0.984, 1.004]; 8/24 | 2.894x [2.813, 2.918]; 24/24 | 1.258x [1.230, 1.267]; 24/24 |
| labute_asa/parse_inclusive | 0.997x [0.990, 1.000]; 8/24 | 9.845x [9.740, 9.876]; 24/24 | 9.268x [9.216, 9.347]; 24/24 |
| morgan2_chiral/first_use | 0.992x [0.987, 0.997]; 7/24 | 1.547x [1.532, 1.560]; 24/24 | 0.999x [0.994, 1.005]; 11/24 |
| morgan2_chiral/repeated_use | 0.993x [0.991, 0.997]; 7/24 | 5.324x [5.296, 5.335]; 24/24 | 3.252x [3.235, 3.262]; 24/24 |
| morgan2_chiral/parse_inclusive | 0.997x [0.989, 1.001]; 8/24 | 4.430x [4.410, 4.445]; 24/24 | 3.878x [3.868, 3.902]; 24/24 |
| maccs/first_use | 0.989x [0.986, 0.993]; 1/24 | 1.788x [1.784, 1.796]; 24/24 | 2.063x [2.056, 2.067]; 24/24 |
| maccs/repeated_use | 0.990x [0.988, 0.993]; 0/24 | 1.835x [1.832, 1.840]; 24/24 | 2.117x [2.112, 2.123]; 24/24 |
| maccs/parse_inclusive | 0.990x [0.987, 0.992]; 2/24 | 1.985x [1.982, 1.987]; 24/24 | 2.244x [2.240, 2.248]; 24/24 |

Parse-only medians (ms per 1,000): source baseline 5.686,
candidate 5.677, RDKit 52.932, COSMolKit
53.011. Peak RSS is retained per worker as a whole-process
measurement, including imports and all operations; it is not per-operation RSS.

## Remaining scope and reproduction

The broader 20-operation, three-block screening is retained as diagnostic
numeric data. It identified rotor, MQN and Layered operation-time gaps.
Screening alone does not establish output equivalence; only the formally
validated conditions above receive a speed ranking. Its pipeline values are
parse-plus-operation sums; the formal table uses actual timed pipelines.
MMFF/UFF, tautomer enumeration, file interchange and bio-structure APIs were
not rebenchmarked here. This record does not supersede published-wheel or
full-10k Morgan timing records, and gains do not imply every API is faster.

- [All 72 comparisons (CSV)](2026-10-19-rotor-layered-lanes.csv).
- [Numeric packet, builds, plan, pilots, screening and correctness](2026-10-19-rotor-layered-numeric.json.gz).
- `scripts/bench_python_speed_frontier.py` accepts all eight operations in the
  packet's plan; supply the two frozen source wheels and installed paths,
  `--blocks 24 --balanced-permutations --offset 0 --limit 1000`.
- `scripts/verify_python_scalar_frontier.py --ops rotatable_bonds,mqn,fp_layered,qed`
  dumps every normal/stress output or error with a selected installed wheel.
- `python3 scripts/check_python_rotor_layered_packet.py benchmarks/2026-10-19-rotor-layered-numeric.json.gz --require-adoption`
  rechecks source/native identity, permutations, every raw duration, all
  intervals, correctness accounting and primary adoption.

Local validation: 991 chem and 353 fingerprint tests passed (one existing chem
test ignored), Clippy with all targets and warnings denied passed, 562 focused
Python tests passed, and 22 method/packet corruption tests passed. The older
packet/method suite also passed 46 tests. The sparse local benchmark-index
check verified tracked path consistency; hosted CI validates real full-checkout
files. Subsequent evidence/tooling commits leave the measured runtime source
hashes unchanged.
