# Parse + compatible Morgan performance plan

Updated 2026-09-23. The optimization phase is complete for the measured
v1.0.20 browser lanes. This document now defines the maintenance and rebaseline
contract; detailed optimization history remains in the dated benchmark record
and Git history.

## Current result

Against official `@rdkit/rdkit@2026.3.6` on the fixed exposed 10,000-row
browser corpus, registry-installed `@kent-tokyo/chematic@1.0.20` records:

| Lane | Speedup | 95% lower bound | Correctness |
|---|---:|---:|---|
| Parse-inclusive compatible Morgan | 1.398x | 1.363x | 9,999/9,999 supported rows bit-exact |
| Prepared compatible Morgan | 3.511x | 3.407x | 9,999/9,999 supported rows bit-exact |

One Fe(II) coordination input remains a typed refusal and is not counted as a
match. The result applies to the named package versions, browser protocol,
radius-2/2048-bit options, corpus, and failure policy. It is not a claim about
every browser, fingerprint setting, native ECFP4, or similarity-search workload.

Evidence:

- [v1.0.20 public-package fingerprint and 3D record](../benchmarks/2026-09-23-public-package-fingerprint-3d-v1.0.20.md)
- [source optimization and multi-browser record](../benchmarks/2026-09-20-parse-morgan-rdkitjs.md)
- [benchmark methodology](benchmark.md)

## Profiles must remain separate

- **Compatible Morgan:** the measured profile is radius 2, 2,048 folded bits,
  chirality off, bond types on, and redundant environments off. Supported rows
  must be bit-identical to the pinned RDKit profile.
- **Native ECFP4:** uses a different definition and hash space. It may be faster
  or useful for native-native retrieval, but cannot substitute for compatible
  Morgan parity.
- **Detail APIs:** sparse counts, raw identifiers, folded/raw bitInfo, chirality,
  other radii or widths, and search are independent contracts. A bit-only
  optimization does not prove these paths unchanged without their regressions.

## Required gate when code or dependencies change

### Correctness first

1. Record chematic source/package identity and RDKit npm/runtime versions.
2. Account for every input as success, failure, refusal, or skipped.
3. Require exact packed-bit equality on every supported row.
4. Recheck Rust/Python/Node/WASM parity where the shared kernel or binding
   conversion changed.
5. Run sparse-count, bitInfo, atom-order, aromatic/Kekulé, isotope, explicit-H,
   edit-invalidation, and search regressions affected by the change.

### Equivalent timing

The primary timed operation is:

`SMILES parse -> compatible preparation -> Morgan fingerprint -> consume the same 256-byte packed output`

Prepared timing keeps construction outside the interval for both engines.
Parse-inclusive timing keeps it inside for both. Options are built outside the
loop, warm-up is symmetric, and unsupported rows remain visible. Startup,
download, memory, and search are separate measurements.

### Statistics

- use fresh-process paired A/B or B/A runs;
- keep raw pairs and the predetermined stopping rule;
- report paired geometric speedup and a 95% interval;
- pass only when the lower bound is greater than 1.0;
- keep browsers, hosts, sessions, and corpora as separate strata;
- do not rerun until a favorable sample appears or drop outliers after seeing
  the result.

## Rebaseline triggers

Rerun the public-package gate when any of these changes:

- the compatible-Morgan kernel or its preparation path;
- SMILES parsing, aromaticity, ring membership, or default sanitization used by
  that profile;
- Rust/WASM binding conversion or packed-byte output;
- chematic's published npm version;
- the official RDKit.js stable version;
- browser/runtime behavior material to the timing protocol.

A source-candidate result remains source evidence until the published tarball
is installed and measured. An old result remains historical; it is never
silently relabeled as evidence for a new release.

## Resource and regression limits

Performance work must not:

- change native ECFP defaults or supported compatible-Morgan bits;
- weaken typed refusals, computation limits, cancellation, or input accounting;
- introduce unbounded caches or stale prepared state after molecule edits;
- increase WASM size, startup, or measured memory materially without a recorded
  tradeoff;
- replace ring or aromaticity semantics with an approximation that merely fits
  the benchmark corpus.

The previous extra 1.10x SMILES-only stretch target remains abandoned. The
active requirement is to preserve correctness and the scoped published
Parse + compatible Morgan advantage, not to optimize indefinitely.

## Maintained entry points

- browser runner: `scripts/bench_browser_wasm_vs_rdkit_isolated.py`
- speed checker: `scripts/check_parse_morgan_rdkit_speed_gate.py`
- bit parity: `scripts/check_browser_rdkit_ecfp4_parity.py`
- cross-binding parity: `scripts/rdkit_ecfp4_cross_binding_parity.py`
- hosted workflow: `.github/workflows/parse-morgan-rdkitjs-gate.yml`
- dated records: [`benchmarks/README.md`](../benchmarks/README.md)

The hosted performance workflow is an environment-specific gate. A successful
run proves only its recorded browser, host, packages, corpus, and operation.
