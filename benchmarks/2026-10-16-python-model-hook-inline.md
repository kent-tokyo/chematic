# Inline aromatic model-gate scratch — 2026-10-16

Local adoption gate: **FAIL**. RDKit/COSMolKit-wide superiority remains unproved. Every measured lane is retained.

This candidate must remain unpromoted: the fixed primary or retention gate below did not pass. CI correctness alone cannot override this performance decision.

## Source change and preservation

Baseline: `3e7f04b20593af4b07569d886b56afc1184468fe` (source version 1.1.5). In `rdkit_model_hook.rs`, common aromaticity-model validation uses inline scratch storage for the visited state, electron counts, component traversal, local adjacency and BFS. The FIFO queue uses a contiguous buffer and a read cursor. Atom/edge visitation order, shortest-path parents, electron counting, rejection rules, sanitization, caches and public APIs are unchanged. Larger buffers spill to the heap. No dependency, unsafe code or new chemical shortcut.

The inline buffers trade bounded stack storage for fewer heap allocations. Whole-process RSS is recorded below, rather than presented as per-operation or large-input memory proof.

The test oracle freezes three complete original gate functions from the named base and checks their saved tokens after restoring renamed identifiers. All 14,072 pinned inputs preserve cold/repeated gate decisions. Directed empty, dummy, radical, charged, dative, disconnected, aromatic and macrocycle cases are covered, along with fused graphs around the inline/spill boundaries, alternative electron assignments and non-component spectator atoms.

The two measured source wheels preserve all descriptor float bits/integers, Morgan bytes and MACCS bytes on the same 14,072 inputs: 112,576 descriptor outcomes, 56,288 Morgan outputs and 28,144 MACCS outcomes. Six pre-existing nonfinite descriptor outcomes are preserved, not scored as accuracy passes. These are exposed development inputs, not an independent holdout.

## Fixed protocol and gates

Two unfiltered 1,000-row slices (offsets 0 and 5,000) of the pinned 10,000-input corpus. Twenty rotating four-arm fresh-process blocks per slice, seven operations and three preparation modes: 42 source-relative lanes. No adaptive timing reruns or pooled epochs. Every worker retains native return objects during timing and normalizes values afterwards. First use parses fresh objects outside the timer and includes lazy work; repeated use warms each object once; parse inclusive runs actual parsing and the operation together. Competitor parsing may eagerly perceive chemistry. The timed TPSA lane uses RDKit-default N/O contributions.

Both source wheels were freshly built with Rust 1.97.0, opt3/full LTO/one codegen unit/one job, native-inchi and strip=none. Comparators: RDKit 2026.09.1, COSMolKit 0.5.0 and CPython 3.13.6 on one shared macOS arm64 host. Package-local native files are hashed before and after each profile; system/transitive libraries and native data files outside those package directories are excluded, so this is not a full runtime closure. Startup/import time is excluded. External host activity is uncontrolled; load averages and whole-process peak RSS are retained.

Ratios are reference time / candidate time. Brackets are paired-bootstrap 95% intervals using the unchanged deterministic 10,000-resample procedure, rounded to three decimals. Local adoption requires exact source outputs, positive interval lower bounds (>1) for all four first-use TPSA/ASA primary lanes, and all 42 source-relative medians >=0.99. Strict competitor wins additionally require output equivalence and all 20 blocks faster. Mismatches remain unranked. Negative intervals and budget failures are disclosed below.

## Primary decision

| Offset | Operation | Source ratio |
|---|---|---:|
| 0 | tpsa/first_use | 1.000 [0.998, 1.008]; 10/20 |
| 0 | labute_asa/first_use | 0.999 [0.996, 1.013]; 9/20 |
| 5000 | tpsa/first_use | 1.007 [0.994, 1.011]; 13/20 |
| 5000 | labute_asa/first_use | 1.002 [0.991, 1.005]; 11/20 |

## Complete matrix: offset 0

| Operation / mode | Baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 1.000 [0.998, 1.008]; 10/20 | 1.501 [1.489, 1.512]; 20/20 **strict win** | 0.531 [0.524, 0.540]; 0/20 |
| tpsa/repeated_use | 1.004 [0.999, 1.021]; 12/20 | 10.124 [9.913, 10.274]; 20/20 **strict win** | 1.704 [1.661, 1.753]; 20/20 **strict win** |
| tpsa/parse_inclusive | 0.995 [0.983, 1.005]; 7/20 | 13.035 [12.814, 13.122]; 20/20 **strict win** | 12.131 [12.062, 12.192]; 20/20 **strict win** |
| labute_asa/first_use | 0.999 [0.996, 1.013]; 9/20 | 0.735 [0.715, 0.755]; 0/20 | 0.354 [0.348, 0.360]; 0/20 |
| labute_asa/repeated_use | 1.005 [0.999, 1.034]; 13/20 | 2.733 [2.680, 2.771]; 20/20 **strict win** | 1.421 [1.386, 1.535]; 20/20 **strict win** |
| labute_asa/parse_inclusive | 0.993 [0.988, 0.999]; 6/20 | 12.272 [12.171, 12.359]; 20/20 **strict win** | 11.704 [11.588, 11.743]; 20/20 **strict win** |
| num_rings/first_use | 1.003 [0.999, 1.010]; 13/20 | 0.037 [0.036, 0.039]; 0/20 | unranked: output difference |
| num_rings/repeated_use | 1.002 [0.995, 1.025]; 11/20 | 1.750 [1.737, 1.788]; 20/20 **strict win** | unranked: output difference |
| num_rings/parse_inclusive | 0.995 [0.988, 1.000]; 6/20 | 11.556 [11.450, 11.632]; 20/20 **strict win** | unranked: output difference |
| morgan2_2048/first_use | 0.999 [0.996, 1.006]; 9/20 | 3.139 [3.120, 3.161]; 20/20 **strict win** | 2.031 [2.021, 2.046]; 20/20 **strict win** |
| morgan2_2048/repeated_use | 1.000 [0.997, 1.014]; 10/20 | 5.185 [5.156, 5.247]; 20/20 **strict win** | 3.189 [3.170, 3.211]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 1.002 [0.989, 1.010]; 10/20 | 8.575 [8.514, 8.637]; 20/20 **strict win** | 7.524 [7.483, 7.595]; 20/20 **strict win** |
| morgan2_chiral/first_use | 1.006 [1.001, 1.013]; 14/20 | 1.260 [1.254, 1.267]; 20/20 **strict win** | 0.827 [0.818, 0.835]; 0/20 |
| morgan2_chiral/repeated_use | 1.001 [0.996, 1.007]; 11/20 | 5.011 [4.960, 5.037]; 20/20 **strict win** | 3.104 [3.068, 3.119]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 1.003 [1.000, 1.006]; 14/20 | 4.460 [4.435, 4.489]; 20/20 **strict win** | 3.908 [3.889, 3.918]; 20/20 **strict win** |
| morgan3_2048/first_use | 1.001 [0.994, 1.008]; 12/20 | 3.382 [3.355, 3.406]; 20/20 **strict win** | 2.246 [2.218, 2.257]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 1.000 [0.994, 1.007]; 10/20 | 4.883 [4.851, 4.910]; 20/20 **strict win** | 3.106 [3.092, 3.123]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 1.001 [0.999, 1.007]; 13/20 | 8.002 [7.937, 8.047]; 20/20 **strict win** | 6.968 [6.914, 7.016]; 20/20 **strict win** |
| maccs/first_use | 0.998 [0.996, 1.002]; 7/20 | 1.801 [1.797, 1.805]; 20/20 **strict win** | 2.063 [2.057, 2.066]; 20/20 **strict win** |
| maccs/repeated_use | 1.000 [0.998, 1.004]; 10/20 | 1.856 [1.853, 1.859]; 20/20 **strict win** | 2.133 [2.129, 2.136]; 20/20 **strict win** |
| maccs/parse_inclusive | 1.001 [0.998, 1.003]; 13/20 | 2.025 [2.019, 2.031]; 20/20 **strict win** | 2.274 [2.270, 2.281]; 20/20 **strict win** |

## Complete matrix: offset 5000

| Operation / mode | Baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 1.007 [0.994, 1.011]; 13/20 | 0.444 [0.442, 0.453]; 0/20 | 0.120 [0.118, 0.123]; 0/20 |
| tpsa/repeated_use | 1.014 [1.000, 1.019]; 13/20 | 9.960 [9.757, 10.067]; 20/20 **strict win** | 1.345 [1.318, 1.388]; 20/20 **strict win** |
| tpsa/parse_inclusive | 1.000 [0.987, 1.004]; 10/20 | 5.186 [5.160, 5.223]; 20/20 **strict win** | 5.333 [5.302, 5.371]; 20/20 **strict win** |
| labute_asa/first_use | 1.002 [0.991, 1.005]; 11/20 | 0.203 [0.199, 0.207]; 0/20 | 0.086 [0.084, 0.087]; 0/20 |
| labute_asa/repeated_use | 1.011 [0.989, 1.032]; 11/20 | 2.632 [2.598, 2.684]; 20/20 **strict win** | 1.184 [1.161, 1.219]; 20/20 **strict win** |
| labute_asa/parse_inclusive | 0.996 [0.990, 1.001]; 6/20 | 4.807 [4.784, 4.865]; 20/20 **strict win** | 5.175 [5.136, 5.221]; 20/20 **strict win** |
| num_rings/first_use | 1.000 [0.993, 1.008]; 9/20 | 0.098 [0.096, 0.101]; 0/20 | 0.084 [0.080, 0.089]; 0/20 |
| num_rings/repeated_use | 1.015 [0.993, 1.027]; 13/20 | 1.687 [1.646, 1.729]; 20/20 **strict win** | 1.149 [1.102, 1.177]; 19/20 |
| num_rings/parse_inclusive | 0.990 [0.983, 1.009]; 9/20 | 11.555 [11.479, 11.680]; 20/20 **strict win** | 12.663 [12.592, 12.715]; 20/20 **strict win** |
| morgan2_2048/first_use | 0.999 [0.992, 1.005]; 10/20 | 1.549 [1.541, 1.556]; 20/20 **strict win** | 0.994 [0.987, 1.007]; 8/20 |
| morgan2_2048/repeated_use | 0.999 [0.988, 1.007]; 10/20 | 4.654 [4.616, 4.690]; 20/20 **strict win** | 2.857 [2.788, 2.905]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 1.000 [0.991, 1.010]; 10/20 | 4.592 [4.539, 4.627]; 20/20 **strict win** | 4.446 [4.426, 4.461]; 20/20 **strict win** |
| morgan2_chiral/first_use | 0.999 [0.992, 1.005]; 10/20 | 1.598 [1.580, 1.601]; 20/20 **strict win** | 1.024 [1.007, 1.033]; 15/20 |
| morgan2_chiral/repeated_use | 0.996 [0.994, 1.001]; 8/20 | 4.672 [4.648, 4.703]; 20/20 **strict win** | 2.844 [2.822, 2.866]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 0.996 [0.989, 1.000]; 5/20 | 4.614 [4.577, 4.640]; 20/20 **strict win** | 4.434 [4.415, 4.468]; 20/20 **strict win** |
| morgan3_2048/first_use | 0.997 [0.992, 1.003]; 9/20 | 1.797 [1.776, 1.809]; 20/20 **strict win** | 1.176 [1.169, 1.187]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 0.996 [0.992, 1.002]; 7/20 | 4.537 [4.491, 4.570]; 20/20 **strict win** | 2.872 [2.847, 2.894]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 0.997 [0.991, 1.004]; 8/20 | 4.510 [4.466, 4.550]; 20/20 **strict win** | 4.268 [4.216, 4.293]; 20/20 **strict win** |
| maccs/first_use | 1.002 [0.998, 1.005]; 12/20 | 1.577 [1.568, 1.585]; 20/20 **strict win** | 1.681 [1.679, 1.689]; 20/20 **strict win** |
| maccs/repeated_use | 1.002 [1.000, 1.005]; 14/20 | 1.666 [1.657, 1.671]; 20/20 **strict win** | 1.775 [1.768, 1.778]; 20/20 **strict win** |
| maccs/parse_inclusive | 1.000 [0.998, 1.002]; 10/20 | 1.739 [1.733, 1.745]; 20/20 **strict win** | 1.861 [1.850, 1.869]; 20/20 **strict win** |

## Every negative source interval

```json
{
  "0": {
    "labute_asa/parse_inclusive": {
      "median": 0.993,
      "ci95": [
        0.988,
        0.999
      ],
      "winning_blocks": 6,
      "blocks": 20
    }
  },
  "5000": {}
}
```

## Every median budget failure

```json
{
  "0": {},
  "5000": {}
}
```

## Absolute batch times and whole-process memory

Milliseconds below are medians for a complete 1,000-input batch.

| Offset / operation / mode | Baseline ms | Candidate ms | RDKit ms | COSMolKit ms |
|---|---:|---:|---:|---:|
| 0 / tpsa/first_use | 1.325 | 1.323 | 1.999 | 0.702 |
| 0 / tpsa/repeated_use | 0.035 | 0.035 | 0.357 | 0.060 |
| 0 / tpsa/parse_inclusive | 4.289 | 4.317 | 56.139 | 52.178 |
| 0 / labute_asa/first_use | 1.388 | 1.386 | 1.019 | 0.488 |
| 0 / labute_asa/repeated_use | 0.036 | 0.035 | 0.096 | 0.051 |
| 0 / labute_asa/parse_inclusive | 4.462 | 4.483 | 55.019 | 52.309 |
| 0 / num_rings/first_use | 1.697 | 1.692 | 0.063 | 0.066 |
| 0 / num_rings/repeated_use | 0.032 | 0.031 | 0.055 | 0.037 |
| 0 / num_rings/parse_inclusive | 4.660 | 4.687 | 54.121 | 51.801 |
| 0 / morgan2_2048/first_use | 4.943 | 4.941 | 15.547 | 10.077 |
| 0 / morgan2_2048/repeated_use | 2.920 | 2.902 | 15.093 | 9.259 |
| 0 / morgan2_2048/parse_inclusive | 8.162 | 8.167 | 70.223 | 61.857 |
| 0 / morgan2_chiral/first_use | 12.460 | 12.383 | 15.611 | 10.258 |
| 0 / morgan2_chiral/repeated_use | 3.070 | 3.064 | 15.314 | 9.464 |
| 0 / morgan2_chiral/parse_inclusive | 15.828 | 15.758 | 70.318 | 61.739 |
| 0 / morgan3_2048/first_use | 5.845 | 5.845 | 19.769 | 13.133 |
| 0 / morgan3_2048/repeated_use | 3.996 | 3.991 | 19.578 | 12.437 |
| 0 / morgan3_2048/parse_inclusive | 9.364 | 9.350 | 74.767 | 65.179 |
| 0 / maccs/first_use | 226.635 | 226.697 | 408.954 | 467.485 |
| 0 / maccs/repeated_use | 219.290 | 219.481 | 407.259 | 467.853 |
| 0 / maccs/parse_inclusive | 230.525 | 230.490 | 466.638 | 524.402 |
| 5000 / tpsa/first_use | 3.574 | 3.563 | 1.593 | 0.429 |
| 5000 / tpsa/repeated_use | 0.036 | 0.035 | 0.352 | 0.048 |
| 5000 / tpsa/parse_inclusive | 4.898 | 4.912 | 25.489 | 26.237 |
| 5000 / labute_asa/first_use | 3.394 | 3.402 | 0.696 | 0.295 |
| 5000 / labute_asa/repeated_use | 0.036 | 0.036 | 0.093 | 0.042 |
| 5000 / labute_asa/parse_inclusive | 5.031 | 5.048 | 24.367 | 26.134 |
| 5000 / num_rings/first_use | 0.569 | 0.568 | 0.056 | 0.048 |
| 5000 / num_rings/repeated_use | 0.032 | 0.031 | 0.052 | 0.036 |
| 5000 / num_rings/parse_inclusive | 2.033 | 2.041 | 23.673 | 25.857 |
| 5000 / morgan2_2048/first_use | 5.429 | 5.439 | 8.412 | 5.413 |
| 5000 / morgan2_2048/repeated_use | 1.730 | 1.734 | 8.094 | 4.933 |
| 5000 / morgan2_2048/parse_inclusive | 7.041 | 7.052 | 32.373 | 31.328 |
| 5000 / morgan2_chiral/first_use | 5.297 | 5.300 | 8.456 | 5.411 |
| 5000 / morgan2_chiral/repeated_use | 1.752 | 1.753 | 8.204 | 5.000 |
| 5000 / morgan2_chiral/parse_inclusive | 7.032 | 7.080 | 32.590 | 31.356 |
| 5000 / morgan3_2048/first_use | 5.836 | 5.840 | 10.492 | 6.845 |
| 5000 / morgan3_2048/repeated_use | 2.241 | 2.252 | 10.188 | 6.481 |
| 5000 / morgan3_2048/parse_inclusive | 7.663 | 7.682 | 34.759 | 32.798 |
| 5000 / maccs/first_use | 136.913 | 136.757 | 215.191 | 229.726 |
| 5000 / maccs/repeated_use | 129.592 | 129.340 | 214.908 | 229.359 |
| 5000 / maccs/parse_inclusive | 138.631 | 138.358 | 240.533 | 257.536 |

Peak RSS includes imports, preparation, warmup and every operation in the worker. It is not isolated operation memory.

| Offset | Baseline MiB | Candidate MiB | RDKit MiB | COSMolKit MiB |
|---|---:|---:|---:|---:|
| 0 | 46.6 | 46.6 | 106.6 | 106.1 |
| 5000 | 40.0 | 40.0 | 88.4 | 83.8 |

## Validation, retention and limits

Rust perception/chemistry: 1,529 passed, 0 failed, 7 ignored across 34 suites. Release all-target perception Clippy, formatting and the refactor budget passed. Python: 2,302 passed, 2 skipped. Native correctness uses opt3/LTO=false/16 codegen units in a separate target; these are not the timed full-LTO wheels.

Pre-timing failures are retained: the initial baseline build snapshot lacked a workspace tools directory; the initial SmallVec array sizes were unsupported by the existing dependency features; and the new empty-input test incorrectly asked the parser to accept an empty SMILES. These were corrected before formal timing. No failed build/test supplied observations.

[2026-10-16-python-model-hook-inline-artifacts](2026-10-16-python-model-hook-inline-artifacts/README.md) retains complete raw matrices, all slice observer outputs and mismatch rows, decisions, the frozen plan, source/library/wheel manifests, exact reconstructable candidate patch, launch/validation drivers and build/test logs. Full wheels, installed extensions, original source snapshots and all 14,072-input output streams remain under `/private/tmp/chematic-model-hook-inline-evidence-20261016` with versioned hashes. Reconstructable source does not imply bit-reproducible rebuilt binaries.

Additional hosts, sealed holdout inputs, other operations, 3D/reaction workflows, cold CLI startup, Python/WASM publication and toolkit-wide superiority remain open. Source-relative interval movement in an untouched API does not establish a causal code regression. A source result also does not establish a published-package win.

