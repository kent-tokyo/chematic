# Neutral boron initial placement: exposed development comparison

The input is `../2026-11-27-selenium-crystal-stereo/source-development-input.json` (SHA-256 `e3c2cf4aa5c7c5740fd3d23c670e26535515b97f16a02ee2bf387661513726a1`). It contains 1,345 records from the exposed COD development partition. No reserved holdout records were selected. The probe parses every input and selects neutral boron with three explicit neighbours and no implicit hydrogen: 24 molecules, 36 centers, zero parse failures. It does not cover implicit-hydrogen boron centers.

`before.json` is from commit `40e6feb9`; `after.json` is from that commit plus the change in `crates/chematic-3d/src/dg_connectivity_ordered.rs` (diff SHA-256 `14523f372a7e5e93c7d3493cd5a326f36aa5da4f14e4da8ea0d9296490d318fc`). Both use the same `crates/chematic-3d/examples/boron_development_geometry.rs` (SHA-256 `3600f9f5e035ffa32897c608e523fbf20152d69901890eeb25747f0a668732f9`), macOS arm64, and rustc 1.100.0. Use separate Cargo target directories for the two checkouts: sharing a target directory between different checkouts produced a stale comparison during this investigation.

From each checkout, after copying the probe into the baseline checkout, run:

```sh
CARGO_TARGET_DIR=/path/to/unique-target cargo run --locked --quiet -p chematic-3d --example boron_development_geometry > result.json
```

The two stages are rule-based placement and ETKDG, respectively. For each center the result records its three pair angles in degrees, all molecular bond lengths in deterministic atom-pair order, declared stereo statuses, and the count of nonbonded atom pairs closer than 0.8 Å. These are diagnostic measurements, not a crystal-structure fit or an independent stereo oracle.

Only record 2877 changed in the reported measurements. Its boron angles changed from 120.0°, 120.0°, 97.2° to 120.0°, 120.0°, 120.0° at placement, and from 96.6°, 102.6°, 97.2° to 120.0°, 120.0°, 120.0° after ETKDG. Initial bond lengths stayed unchanged. After ETKDG two aromatic C–C bonds changed from 1.870 and 1.166 Å to 1.400 Å each; other bond lengths changed by less than 0.05 Å. The nonbonded pair counts stayed at 4 for placement and 0 after ETKDG. Its one declared tetrahedral status remains `Violated` in both stages, before and after; this change does not fix that separate issue. The reported measurements for the other 23 molecules are unchanged.

This fixed development sample tests a narrow subset, and its crystal coordinates were not used as an oracle here. It does not establish general geometric accuracy, performance, or parity with RDKit.
