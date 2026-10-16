# Python Labute ASA loop locality — 2026-10-16

Local adoption gate: **FAIL**. Toolkit-wide superiority remains unproved. The complete measured candidate and every observation are retained.

The two first-use ASA ratios are 1.004 [0.993, 1.014] and 1.003 [1.000, 1.011], so neither establishes the required improvement. Offset-0 repeated ring count has median 0.982 [0.962, 1.013], below the 0.99 floor; its interval is inconclusive and does not establish a causal ring regression. The candidate remains unpromoted and must not be merged as a demonstrated optimization.

## Source change and correctness

Baseline: `96a1addb742b9b4cb3c8258e123563cee0d7055c`. The candidate fuses the final hydrogen and atom-surface loops in `labute_asa_core`. Radii, bond overlap, RDKit bond traversal, formulas, independent atom-index accumulation orders, model preparation and APIs are unchanged. It also iterates the existing scratch array directly. No new dependency or unsafe code.

The regression oracle freezes the entire original core from the named base and compares every per-atom value, the hydrogen pool and total by exact float bits, with and without per-atom collection, cold and repeated calls on 14,072 pinned inputs. Empty, dummy/nonfinite, explicit-H, dative, disconnected, aromatic and inline/spill-capacity boundaries are checked separately. Public descriptor outcomes, Morgan bytes and MACCS bytes are also compared between the two measured wheels on all 14,072 inputs. These are development inputs, not an independent holdout. Preserved nonfinite outcomes are not accuracy passes.

## Protocol

Two unfiltered 1,000-input slices (offsets 0 and 5,000) of the pinned 10,000-row corpus. Twenty rotating four-arm fresh-process blocks per slice. Seven operations and three preparation modes yield 42 source-relative lanes. Native return objects remain retained during timing; output normalization is outside the timer. First use excludes parsing and includes lazy work; repeated use warms the operation once; parse inclusive runs actual parsing inside the loop. Published comparator parsing may eagerly perform perception.

Builds: Rust 1.97.0, release opt3/full LTO/one codegen unit/one job, native-inchi, explicit strip=none. Source wheels are local builds, not published chematic artifacts. Comparators: RDKit 2026.09.1 and COSMolKit 0.5.0, CPython 3.13.6 on one macOS arm64 host. Package-local native libraries are hashed before and after timing; this excludes system libraries and is not a full runtime closure. Host load averages and whole-process peak RSS are retained; RSS includes imports/preparation/warmup and is not per-operation memory. External host activity is uncontrolled.

Ratios are reference time / candidate time. Brackets are paired-bootstrap 95% intervals (existing 10,000-resample procedure, three-decimal output). A strict comparator win requires equivalent outputs, CI lower >1 and all 20 blocks faster. COSMolKit ring differences are unranked. The local adoption gate was fixed before timing: first-use ASA CI lower >1 on both slices, exact source outputs, and every source-relative median >=0.99. All negative intervals are disclosed; no adaptive timing reruns or pooled observations.

## Complete matrix: offset 0

| Operation / mode | Baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 0.994 [0.987, 1.007]; 8/20 | 1.506 [1.473, 1.527]; 19/20 | 0.551 [0.540, 0.565]; 0/20 |
| tpsa/repeated_use | 1.003 [0.949, 1.104]; 10/20 | 9.515 [9.227, 10.022]; 20/20 **strict win** | 1.823 [1.662, 2.030]; 20/20 **strict win** |
| tpsa/parse_inclusive | 0.996 [0.986, 1.001]; 7/20 | 13.023 [12.879, 13.192]; 20/20 **strict win** | 12.177 [12.087, 12.296]; 20/20 **strict win** |
| labute_asa/first_use | 1.004 [0.993, 1.014]; 11/20 | 0.737 [0.713, 0.752]; 0/20 | 0.370 [0.358, 0.392]; 1/20 |
| labute_asa/repeated_use | 1.018 [0.961, 1.041]; 12/20 | 2.702 [2.619, 2.765]; 19/20 | 1.562 [1.343, 1.767]; 19/20 |
| labute_asa/parse_inclusive | 0.994 [0.987, 1.006]; 8/20 | 12.160 [12.075, 12.411]; 20/20 **strict win** | 11.708 [11.523, 11.764]; 20/20 **strict win** |
| num_rings/first_use | 1.002 [0.990, 1.014]; 12/20 | 0.037 [0.037, 0.039]; 0/20 | unranked: output difference |
| num_rings/repeated_use | 0.982 [0.962, 1.013]; 8/20 | 1.729 [1.692, 1.765]; 19/20 | unranked: output difference |
| num_rings/parse_inclusive | 0.993 [0.986, 1.008]; 9/20 | 11.530 [11.411, 11.673]; 20/20 **strict win** | unranked: output difference |
| morgan2_2048/first_use | 1.002 [0.996, 1.017]; 13/20 | 3.135 [3.064, 3.172]; 20/20 **strict win** | 2.029 [2.002, 2.058]; 20/20 **strict win** |
| morgan2_2048/repeated_use | 1.014 [0.999, 1.036]; 14/20 | 5.222 [5.128, 5.258]; 20/20 **strict win** | 3.228 [3.195, 3.307]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 1.004 [0.996, 1.009]; 11/20 | 8.668 [8.578, 8.757]; 20/20 **strict win** | 7.595 [7.510, 7.695]; 20/20 **strict win** |
| morgan2_chiral/first_use | 0.999 [0.995, 1.020]; 9/20 | 1.256 [1.243, 1.273]; 19/20 | 0.825 [0.814, 0.840]; 1/20 |
| morgan2_chiral/repeated_use | 1.003 [0.990, 1.022]; 11/20 | 4.988 [4.917, 5.046]; 20/20 **strict win** | 3.099 [3.057, 3.131]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 1.000 [0.991, 1.006]; 10/20 | 4.432 [4.383, 4.451]; 20/20 **strict win** | 3.896 [3.844, 3.957]; 20/20 **strict win** |
| morgan3_2048/first_use | 0.998 [0.991, 1.023]; 10/20 | 3.349 [3.298, 3.388]; 20/20 **strict win** | 2.254 [2.234, 2.278]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 0.995 [0.983, 1.013]; 9/20 | 4.895 [4.808, 5.019]; 20/20 **strict win** | 3.148 [3.110, 3.251]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 0.999 [0.989, 1.017]; 10/20 | 7.991 [7.882, 8.214]; 20/20 **strict win** | 6.974 [6.892, 7.087]; 20/20 **strict win** |
| maccs/first_use | 0.998 [0.994, 1.012]; 8/20 | 1.803 [1.789, 1.818]; 20/20 **strict win** | 2.074 [2.058, 2.113]; 20/20 **strict win** |
| maccs/repeated_use | 1.002 [0.983, 1.007]; 11/20 | 1.857 [1.853, 1.863]; 20/20 **strict win** | 2.142 [2.131, 2.198]; 20/20 **strict win** |
| maccs/parse_inclusive | 0.999 [0.991, 1.006]; 8/20 | 2.018 [2.000, 2.043]; 20/20 **strict win** | 2.278 [2.266, 2.316]; 20/20 **strict win** |

## Complete matrix: offset 5000

| Operation / mode | Baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 1.002 [0.995, 1.008]; 11/20 | 0.473 [0.463, 0.482]; 0/20 | 0.138 [0.134, 0.144]; 0/20 |
| tpsa/repeated_use | 1.058 [0.964, 1.111]; 13/20 | 9.276 [8.892, 9.792]; 20/20 **strict win** | 1.529 [1.363, 1.644]; 19/20 |
| tpsa/parse_inclusive | 0.997 [0.988, 1.012]; 7/20 | 5.262 [5.179, 5.277]; 20/20 **strict win** | 5.393 [5.349, 5.450]; 20/20 **strict win** |
| labute_asa/first_use | 1.003 [1.000, 1.011]; 14/20 | 0.228 [0.214, 0.239]; 0/20 | 0.108 [0.102, 0.111]; 0/20 |
| labute_asa/repeated_use | 1.015 [0.972, 1.048]; 11/20 | 2.529 [2.394, 2.623]; 20/20 **strict win** | 1.268 [1.220, 1.349]; 20/20 **strict win** |
| labute_asa/parse_inclusive | 0.995 [0.988, 1.004]; 8/20 | 4.857 [4.831, 4.940]; 20/20 **strict win** | 5.203 [5.177, 5.256]; 20/20 **strict win** |
| num_rings/first_use | 1.001 [0.984, 1.011]; 10/20 | 0.105 [0.101, 0.109]; 0/20 | 0.124 [0.112, 0.132]; 0/20 |
| num_rings/repeated_use | 1.001 [0.984, 1.051]; 10/20 | 1.740 [1.710, 1.826]; 20/20 **strict win** | 1.129 [1.069, 1.151]; 18/20 |
| num_rings/parse_inclusive | 0.995 [0.982, 1.003]; 8/20 | 11.777 [11.616, 11.927]; 20/20 **strict win** | 12.770 [12.652, 12.904]; 20/20 **strict win** |
| morgan2_2048/first_use | 1.003 [0.992, 1.013]; 12/20 | 1.545 [1.519, 1.562]; 20/20 **strict win** | 1.025 [1.010, 1.034]; 16/20 |
| morgan2_2048/repeated_use | 1.008 [1.000, 1.022]; 14/20 | 4.643 [4.600, 4.706]; 20/20 **strict win** | 2.918 [2.866, 2.950]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 1.001 [0.985, 1.014]; 10/20 | 4.632 [4.541, 4.679]; 20/20 **strict win** | 4.453 [4.403, 4.491]; 20/20 **strict win** |
| morgan2_chiral/first_use | 0.999 [0.995, 1.009]; 10/20 | 1.600 [1.586, 1.623]; 20/20 **strict win** | 1.057 [1.045, 1.074]; 19/20 |
| morgan2_chiral/repeated_use | 0.997 [0.990, 1.009]; 9/20 | 4.758 [4.699, 4.924]; 20/20 **strict win** | 3.010 [2.915, 3.141]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 0.998 [0.979, 1.006]; 9/20 | 4.613 [4.520, 4.679]; 20/20 **strict win** | 4.437 [4.326, 4.469]; 20/20 **strict win** |
| morgan3_2048/first_use | 0.997 [0.988, 1.017]; 9/20 | 1.780 [1.771, 1.799]; 20/20 **strict win** | 1.207 [1.181, 1.212]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 1.001 [0.983, 1.003]; 11/20 | 4.492 [4.410, 4.545]; 20/20 **strict win** | 2.975 [2.887, 3.017]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 0.997 [0.960, 1.005]; 9/20 | 4.481 [4.352, 4.519]; 20/20 **strict win** | 4.304 [4.120, 4.383]; 20/20 **strict win** |
| maccs/first_use | 1.002 [0.997, 1.005]; 12/20 | 1.572 [1.565, 1.585]; 20/20 **strict win** | 1.686 [1.681, 1.695]; 20/20 **strict win** |
| maccs/repeated_use | 1.001 [0.997, 1.005]; 11/20 | 1.655 [1.641, 1.673]; 20/20 **strict win** | 1.777 [1.770, 1.790]; 20/20 **strict win** |
| maccs/parse_inclusive | 0.999 [0.991, 1.002]; 9/20 | 1.738 [1.734, 1.745]; 20/20 **strict win** | 1.868 [1.849, 1.888]; 20/20 **strict win** |

## Every negative source interval

```json
{
  "0": {},
  "5000": {}
}
```

## Validation and retention

- `native-tests.log`: test result: ok. 960 passed; 0 failed; 1 ignored; 0 measured; 0 filtered out; finished in 0.78s; test result: ok. 3 passed; 0 failed; 1 ignored; 0 measured; 0 filtered out; finished in 5.24s; test result: ok. 4 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s; test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s; test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s; test result: ok. 4 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s; test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s; test result: ok. 8 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s; test result: ok. 5 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s; test result: ok. 4 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s; test result: ok. 2 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s; test result: ok. 4 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s; test result: ok. 6 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.01s; test result: ok. 11 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s; test result: ok. 16 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.01s; test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.03s; test result: ok. 9 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.01s; test result: ok. 2 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 8.59s; test result: ok. 2 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s; test result: ok. 122 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.01s; test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s; test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s; test result: ok. 12 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s; test result: ok. 8 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s; test result: ok. 12 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s; test result: ok. 13 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.22s; test result: ok. 14 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s; test result: ok. 15 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 0.00s; test result: ok. 7 passed; 0 failed; 3 ignored; 0 measured; 0 filtered out; finished in 0.00s
- `python-tests-final.log`: ======================= 2302 passed, 2 skipped in 7.48s ========================
- `runner-tests.log`: 9 passed in 0.09s

Release Rust tests use LTO=false and 16 codegen units in a separate target; those are correctness checks, not the full-LTO timing executables. All-target chemistry Clippy, formatting, refactor budget and the frozen-oracle extraction check passed.

An initial default-strip wheel failed to import with the macOS LINKEDIT alignment error. Both source arms were rebuilt with explicit strip=none and imports validated before timing. A candidate rebuild also stopped on ENOSPC before timing. Failed logs and unusable-binary hashes remain; the unused regenerable failed binaries and completed owned build caches were cleaned with preservation/process/open-file checks. No failed or partial build contributed timing samples.

[2026-10-16-python-labute-locality-artifacts](2026-10-16-python-labute-locality-artifacts/README.md) contains complete raw matrices, every negative interval, the predeclared plan, source/library/wheel manifests, output hashes, reconstructable candidate patch, test/build logs, and exact launch/validation drivers. Full source wheels, installed extensions and complete output streams remain under `/private/tmp/chematic-labute-locality-evidence-20261016`; hashes are versioned. Source reconstruction does not imply bit-reproducibility of a rebuilt wheel.

The offset-5000 checkpoint initially stopped on ENOSPC after seven saved blocks (the eighth in-memory observation was lost). That incomplete epoch and its logs are retained separately. Before inspecting performance summaries, one administrative restart was declared: a fresh complete 20-block offset-5000 epoch, with no pooling or sample substitution. The complete original offset-0 matrix was retained. The recovery reserves 64MiB and retries only the same checkpoint after ENOSPC, without rerunning a worker or timer.

The scope excludes additional hosts, held-out chemistry, other APIs, 3D/reaction throughput, cold CLI startup and published artifacts. An immediate-baseline result does not resolve historical regressions against older candidates, and a source optimization does not establish universal RDKit/COSMolKit superiority.
