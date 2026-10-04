# RDKit 2026.09.1 independent C++ source rebaseline

This is an **output differential**, not a speed benchmark. It compares RDKit
`Release_2026_03_6` (`0e0d85f4ca34aeae15dfc0f7cf5503bdb0a8e985`) with
`Release_2026_09_1` (`fece8caa860bdf6c9c82bdb9253f22a0cf87150c`). Both
were built from pinned source with the same flags in
[`rdkit_native_rebaseline_build.json`](../validation/rdkit_native_rebaseline_build.json),
without Python wrappers. They are **not** distributed Python/nanobind wheels or
native binary releases. The CI build used Ubuntu 24.04, GCC 13.3.0 and CMake
3.31.6. Fresh builds and the first output packet succeeded in
[CI #37171205566](https://github.com/kent-tokyo/chematic/actions/runs/37171205566);
the final packet with the old-Python cross-check and parser-failure accounting
succeeded in [CI #37173214297](https://github.com/kent-tokyo/chematic/actions/runs/37173214297).

The shared exposed corpus contains 10,000 SMILES (SHA-256
`f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`)
and 31 SMARTS queries (SHA-256
`250a241a8905d013eb58fff0ccfb1e0dbe3e97c08ff56e021828fbfec957b341`).
No sealed accuracy cohort was used. Both C++ lanes emitted every row. The
older C++ oracle also agreed on canonical SMILES, atom/bond CIP labels and
Morgan 2048-bit fingerprints with the pinned published RDKit 2026.03.6 Python
wheel on all 10,000 rows. This checks the C++ runner's operation settings; it
does not compare every SMARTS cell against Python because the old Python packet
did not preserve every oracle target set.

| Old to new C++ output | Changed |
|---|---:|
| SMILES parse status / canonical spelling | 0 / 0 rows |
| CIP atom / bond labels | 0 / 0 rows |
| Morgan radius-2 2048-bit fingerprint | 0 rows |
| SMARTS target-atom sets | 12 of 310,000 cells, across 6 rows |

The changed input indices are **9, 23, 28, 29, 30 and 34**. Each has one
`[R2]` and one `[R3]` change. Their old and new target-atom sets exactly match
the independently published `@rdkit/rdkit` npm/WASM delta after numeric atom-set
normalization. This is an **RDKit oracle change** in these ring-count queries,
not a measured CheMatic regression or a claim of complete SMARTS parity.

The checked-in [summary](../validation/results/rdkit-native-rebaseline-2026-03-6-to-2026-09-1-2026-10-04-summary.json),
[old rows](../validation/results/rdkit-native-rebaseline-2026-03-6-to-2026-09-1-2026-10-04-old.jsonl.gz),
[new rows](../validation/results/rdkit-native-rebaseline-2026-03-6-to-2026-09-1-2026-10-04-new.jsonl.gz)
and [delta rows](../validation/results/rdkit-native-rebaseline-2026-03-6-to-2026-09-1-2026-10-04-delta.jsonl.gz)
retain the raw evidence. The summary records source, compiled-probe and
uncompressed-row SHA-256 digests, build flags and toolchain. The Linux binary
and linked RDKit libraries themselves were not checked in; the source-build
lane is therefore marked partial-provenance in the manifest, separately from
unavailable distributed binary lanes.

Recheck the versioned evidence, including every row and the npm delta:

```bash
python3 scripts/check_rdkit_native_rebaseline_evidence.py
```

The RDKit 2026.09.1 **distributed Python/nanobind** lane remains unmeasured:
its PyPI version endpoint was still HTTP 404 when this packet was assembled.
Do not transfer the C++ or npm result to that binding or infer Python call
overhead from this record.
