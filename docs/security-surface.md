# Security surface inventory

This inventory records the public input boundaries reviewed for the v0.83.0
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
| SMILES / SMARTS | `parse_with_limits`, `parse_smarts_with_config`, bounded match/MCS APIs | input/atom limits, match budgets and typed outcomes on bounded paths | legacy convenience APIs remain less explicit; S1 should make policy uniform |
| Molecule formats | SDF, MOL, CML, MolJSON, CJSON, mmCIF, CIF, XYZ, PDBQT, ORCA, Gaussian, Cube, OpenDX, LAMMPS | format-specific byte, line, record, depth, atom/bond, or voxel limits where documented | continue inventory for formats without a `*ParseLimits` API |
| KET / MRV / XML-like data | KET and MRV readers | bounded records and MRV scanner depth/attribute/input checks | malformed corpus and parser-wide panic gate |
| Reactions | reaction parsing, transform, enumeration, retro, mapped I/O | typed errors and selected match/output budgets | unify generation/output limits and cancellation |
| 3D / force fields | conformer embedding, distance geometry, UFF, descriptors, volumetric grids | stage/config limits in several APIs; grid point and coordinate validation | audit every public O(n²)/search path and add allocation/time budgets |
| Rust file/process surface | CLI and explicit file APIs | core library is offline by default; CLI policy is host-controlled | S3 path, symlink, overwrite, and atomic-write tests |
| Python / WASM / Node | PyO3 and WASM format/descriptor bindings | typed binding errors and browser sandbox for WASM | boundary size, panic-unwind, buffer, and response-size tests |
| MCP | molecule/format/search tools | URL encoding and selected timeout/atom guards | strict schemas, request/output/concurrency policy, SSRF/path audit |
| Supply chain | Cargo.lock, workflows, release scripts | local audit, Clippy, fmt, build, cargo-deny workflow definitions | clean-room artifact, SBOM, provenance, and live repository-setting verification |

## Verification performed for this baseline

- `cargo audit --no-fetch` scanned the lockfile with the locally available
  advisory database; no vulnerability advisory was reported.
- `cargo clippy -p chematic-mol --all-targets -- -D warnings` passed.
- MolJSON bounded parser regression tests passed.
- `cargo deny` could not complete locally because the required crate archive was
  unavailable while crates.io access was restricted.

The audit reported two unmaintained transitive packages (`rustybuzz` and
`ttf-parser`) used by the depiction dependency chain. They are recorded as
known residual dependency risk, not as security vulnerabilities. Re-evaluate
when the upstream `usvg`/`resvg`/`svg2pdf` chain moves to maintained font
components.

## Classification rule

An accuracy or compatibility difference is a correctness issue unless it can
cause confidentiality, integrity, availability, memory-safety, or supply-chain
impact. Resource exhaustion, panic on attacker-controlled input, arbitrary
filesystem/process access, unintended network access, and secret leakage are
security issues.
