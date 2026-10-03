# RDKit 2026.09.1 rebaseline — published npm lane (2026-10-04)

This is a version-pinned **partial** rebaseline. Official `@rdkit/rdkit@2026.9.1`
is published and its browser runtime reports `2026.09.1`. At the collection time,
`rdkit==2026.9.1` returned PyPI 404; the GitHub release had no attached binary.
Therefore Python/nanobind call overhead, Python CIP/SMARTS parity, and an
independent native-C++ lane are **not measured** here. Do not carry over their
2026.03.6 results or call the entire rebaseline complete.

## Fixed inputs and artifacts

- Exposed, non-sealed 10,000-SMILES corpus:
  `validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi`, SHA-256
  `f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`.
- Official RDKit npm tarballs: `2026.3.6` and `2026.9.1`; published CheMatic
  npm tarball: `@kent-tokyo/chematic@1.0.33`. Exact tarball hashes, npm SRI,
  package metadata, WASM hashes, runtime versions, Node/macOS environment, and
  corpus hash are in the two `rdkit-rebaseline-npm-provenance-v1.0.33-*` records.
- Chrome on the same macOS arm64 host, 20 isolated processes per arm, alternating
  CheMatic/RDKit order within each version's run, 20 warm-up rows. `single_batch`
  timing avoids sub-millisecond per-row clock quantization. The two RDKit
  versions were run **sequentially**, not directly counterbalanced against each
  other. No confidence interval or cross-host speed claim is made.

## Chemistry output

| Check | RDKit 2026.03.6 | RDKit 2026.09.1 |
|---|---:|---:|
| CheMatic/RDKit radius-2, 2048-bit Morgan exact, supported rows | 9,999 / 9,999 | 9,999 / 9,999 |
| Typed unsupported coordination row | 1 / 10,000 | 1 / 10,000 |

The unchanged unsupported row is input index 8341, an Fe(II) coordination
structure. It remains visible in the direct browser gate, not removed from the
corpus. This does not prove parity for other fingerprint options or inputs.

The **RDKit-version-to-RDKit-version** direct comparison covered all 10,000
inputs and the 31 fixed SMARTS queries (310,000 cells):

| Oracle-version difference | Count |
|---|---:|
| Parse acceptance | 0 / 10,000 |
| Canonical SMILES spelling | 0 / 10,000 |
| Morgan bit string | 0 / 10,000 |
| Atom/bond JSON graph representation | 0 / 10,000 |
| CIP stereo tags, after graph correspondence | 0 / 10,000 |
| SMARTS atom-set result | 12 / 310,000 |

All 12 SMARTS cells are `[R2]` or `[R3]`, six cells each, on six complex
polycyclic inputs (indices 9, 23, 28, 29, 30, 34). The CIP comparison uses
MinimalLib `get_stereo_tags()` and compares tags only after the ordered atom
and bond JSON graph agrees. This is an **oracle change**,
not yet adjudicated as an RDKit fix or regression, and it is not a new CheMatic
CIP/SMARTS parity measurement. Every row is retained with status, output hashes,
and explicit changed atom sets in the compressed JSONL packet.

## Browser measurements

WASM size uses raw bytes and gzip level 9, only the `.wasm` file; package
tarball and JavaScript sizes are independently recorded in raw evidence.

| Metric | RDKit.js 2026.03.6 | RDKit.js 2026.09.1 | CheMatic 1.0.33 in new-version run |
|---|---:|---:|---:|
| WASM raw / gzip-9 (bytes) | 7,333,095 / 2,379,975 | 7,379,879 / 2,405,173 | 4,379,739 / 1,610,304 |
| Cold initialization p50 (ms) | 74.45 | 84.60 | 24.45 |
| Parse p50 of per-process means (ms/mol) | 0.159495 | 0.188905 | 0.004940 |
| Parse + canonical write | 0.224310 | 0.252370 | 0.075065 |
| Parse + configured Morgan | 0.216037 | 0.244619 | 0.025818 |
| Prepared fingerprint first use | 0.052145 | 0.052275 | 0.006601 |
| Prepared fingerprint reuse | 0.051580 | 0.051485 | 0.006541 |

CheMatic timing is the same published tarball in both runs. The later RDKit.js
run's higher parse-inclusive medians are **observations**, not a causal speed
regression claim: the versions were measured in separate runs. Parsed-molecule
cache/perception boundaries also differ between implementations. The prepared
lanes exclude preparation from the timed operation; RSS was not measured.

## Evidence and reproduction

- Direct bit gate:
  `validation/results/rdkit-rebaseline-npm-morgan-parity-v1.0.33-vs-{2026.03.6,2026.09.1}-2026-10-04.json`.
- Isolated 20-run browser records:
  `validation/results/rdkit-rebaseline-npm-browser-runtime-v1.0.33-vs-{2026.03.6,2026.09.1}-2026-10-04.json`.
- Oracle delta and **all 10,000 row records**:
  `validation/results/rdkit-rebaseline-npm-oracle-delta-2026.03.6-to-2026.09.1-2026-10-04.{json,jsonl.gz}`.
- Artifact availability snapshot:
  `validation/results/rdkit-2026-09-1-artifact-availability-2026-10-04.json`.
- Runners: `scripts/compare_rdkitjs_release_chemistry.cjs`,
  `scripts/check_browser_rdkit_ecfp4_parity.py`,
  `scripts/bench_browser_wasm_vs_rdkit_isolated.py`, and the argv-safe
  `validation/rdkit_rebaseline_execution.json` command catalog.

Install exact npm versions into separate scratch prefixes, then invoke the
runners with those extracted package directories, the fixed corpus/query
paths above, and an explicit Chrome executable. Verify the packet with
`python3 scripts/check_rdkit_2026_09_rebaseline_evidence.py`.

**Remaining to close the full rebaseline:** obtain and hash a released
2026.09.1 Python artifact; identify Boost.Python versus nanobind from that
artifact rather than the release note; rerun the 10k complete-row Python
SMILES/CIP/SMARTS/Morgan and typed-error/overhead lanes against the same
2026.03.6 baseline; classify each new residual; run a separately pinned C++
lane only if a reproducible binary/build packet is available. Recheck the
package endpoints before those runs.
