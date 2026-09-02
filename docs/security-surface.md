# Security surface inventory

This inventory records the public input boundaries reviewed for the v0.89.0
security baseline. It is a scope document, not a claim that every hostile input
has been exhaustively fuzzed.

## Trust model

| Attacker | Boundary | Primary risk |
| --- | --- | --- |
| Untrusted caller | Rust parser and chemistry APIs | panic, allocation exhaustion, CPU exhaustion, incorrect state |
| Bulk/service caller | SMARTS, MCS, reaction, 3D, and batch APIs | combinatorial CPU or output growth |
| Untrusted document | file-format readers and serializers | malformed records, oversized fields, numeric edge cases |
| Untrusted web/client input | Python, WASM, Node, MCP, and CLI bindings | boundary panic, oversized request/response, path or URL misuse |
| Dependency/CI actor | Cargo dependencies, Actions, build and generated files | supply-chain compromise or artifact mismatch |
| Local privileged user | explicit filesystem/process operations | overwrite, path policy, or host-side privilege mistakes |

## Public boundary inventory

| Area | Representative entrypoints | Current controls | Residual risk / next gate |
| --- | --- | --- | --- |
| SMILES / SMARTS / TDT | `parse_with_limits`, `parse_smi_file_with_limits`, `SmilesRecordReader`, `TdtRecordReader`, `parse_smarts_with_config`, bounded match/MCS APIs | input/atom limits, table/stream line/record/field/tag limits, match budgets and typed outcomes on bounded paths | legacy convenience APIs remain less explicit; S1 should make policy uniform |
| Molecule formats | SDF, MOL, CML, MolJSON, CJSON, mmCIF, CIF, XYZ, PDBQT, ORCA, Gaussian, Cube, OpenDX, LAMMPS, InChI | format-specific byte, line, record, depth, atom/bond, or voxel limits where documented; 3D XYZ and pure-Rust InChI now have explicit limits | continue inventory for formats without a `*ParseLimits` API |
| KET / MRV / XML-like data | KET and MRV readers | KET input/atom/bond limits; MRV scanner depth/attribute/input checks | malformed corpus and parser-wide panic gate |
| Reactions | reaction parsing, transform, enumeration, retro, mapped I/O | RXN input/component limits, typed errors, and selected match/output budgets | unify generation/output limits and cancellation |
| 3D / force fields | conformer embedding, distance geometry, UFF, descriptors, volumetric grids | stage/config limits in several APIs; grid point and coordinate validation | audit every public O(n²)/search path and add allocation/time budgets |
| Rust file/process surface | CLI and explicit file APIs | core library is offline by default; CLI policy is host-controlled | S3 path, symlink, overwrite, and atomic-write tests |
| Python / WASM / Node | PyO3 and WASM format/descriptor bindings | typed binding errors and browser sandbox for WASM | boundary size, panic-unwind, buffer, and response-size tests |
| MCP | molecule/format/search tools | URL encoding and selected timeout/atom guards | strict schemas, request/output/concurrency policy, SSRF/path audit |
| Supply chain | Cargo.lock, workflows, release scripts | local audit, Clippy, fmt, build, cargo-deny workflow definitions | clean-room artifact, SBOM, provenance, and live repository-setting verification |

The optional `native-inchi` feature is the sole C FFI boundary. Its unsafe
calls are confined to `chematic-inchi/src/native`, use fixed `repr(C)` layouts,
validate non-empty input and signed count ranges before crossing the boundary,
and release every library-owned output through `FreeStdINCHI`. The default
pure-Rust and WASM paths do not enable this feature.

## Verification performed for this baseline

- `cargo audit --no-fetch` scanned the lockfile with the locally available
  advisory database; no vulnerability advisory was reported.
- `cargo clippy -p chematic-mol --all-targets -- -D warnings` passed.
- MolJSON bounded parser regression tests passed.
- Unsafe inventory found no unsafe code outside the optional native-InChI FFI
  module; all other reviewed public crates forbid unsafe code.
- `cargo deny check` passed advisories, bans, licenses, and sources against the
  current locked dependency graph.
- GitHub Linux AddressSanitizer, LeakSanitizer, and ThreadSanitizer all passed
  on the core/parser scope; focused Miri tests also passed in run 33588355817.

The audit reported two unmaintained transitive packages (`rustybuzz` and
`ttf-parser`) used by the depiction dependency chain. This is an explicitly
accepted residual dependency risk, not a security vulnerability: the packages
are not directly exposed as a public API, no vulnerable advisory is reported,
and replacing the rendering chain would change artifact behavior. The release
owner must re-evaluate this decision when upstream `usvg`/`resvg`/`svg2pdf`
moves to maintained font components, or when an advisory affects the chain.

## Classification rule

An accuracy or compatibility difference is a correctness issue unless it can
cause confidentiality, integrity, availability, memory-safety, or supply-chain
impact. Resource exhaustion, panic on attacker-controlled input, arbitrary
filesystem/process access, unintended network access, and secret leakage are
security issues.
