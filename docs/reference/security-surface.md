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
pure-Rust crates and the WASM build do not enable this feature; the published
Python wheels do (`[tool.maturin] features` in `crates/chematic-py/pyproject.toml`),
so that `standard_inchi` returns IUPAC InChI instead of an error.

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

The original audit reported two informational unmaintained advisories for
`rustybuzz` and `ttf-parser` in the `svg2pdf` PDF depiction chain. Neither
advisory reported an exploit or patched release. The `svg2pdf` repository is
archived, and its suggested successor `krilla-svg` 0.8.1 requires Rust 1.92,
above chematic's Rust 1.88 minimum.

The PDF depiction feature now uses `harumi` 1.23.8 with its SVG converter and
an embedded OFL-licensed Lato font. This release fixes the text placement,
three-digit hex colors, marker definitions, SVG path commands, and empty
transform panics reported in harumi issues
[#79](https://github.com/kent-tokyo/harumi/issues/79),
[#80](https://github.com/kent-tokyo/harumi/issues/80),
[#82](https://github.com/kent-tokyo/harumi/issues/82),
[#85](https://github.com/kent-tokyo/harumi/issues/85), and
[#86](https://github.com/kent-tokyo/harumi/issues/86). Updating its `lopdf`
dependency also removes the transitive `ttf-parser` reported in
[#81](https://github.com/kent-tokyo/harumi/issues/81). Harumi 1.23.7 removes
the `ttf-parser` dependency reintroduced in 1.23.6
([#115](https://github.com/kent-tokyo/harumi/issues/115)). The locked PDF
depiction graph contains neither `rustybuzz` nor `ttf-parser`. The PDF export
contract tests cover vector geometry, labels, dashed bonds, reaction and grid
SVGs, page dimensions, and SVG text alignment. Version 1.23.2 also fixes
SVG subpath separation and polygon fill/stroke colors reported in
[#92](https://github.com/kent-tokyo/harumi/issues/92) and
[#91](https://github.com/kent-tokyo/harumi/issues/91). Version 1.23.5 also
accepts non-rendering SVG metadata and CSS named/`rgb()` colors reported in
[#107](https://github.com/kent-tokyo/harumi/issues/107) and
[#108](https://github.com/kent-tokyo/harumi/issues/108). Version 1.23.8 also
includes root SVG style inheritance and nested opacity handling for
individual primitives reported in [#119](https://github.com/kent-tokyo/harumi/issues/119)
and [#120](https://github.com/kent-tokyo/harumi/issues/120). The `pdf` feature
remains opt-in. Harumi pins `aes` to 0.9.2 for Rust 1.88 consumers,
resolving [#89](https://github.com/kent-tokyo/harumi/issues/89); chematic no
longer needs its own temporary direct dependency.

Upstream status: [svg2pdf](https://github.com/typst/svg2pdf),
[krilla-svg 0.8.1](https://docs.rs/crate/krilla-svg/0.8.1),
[resvg changelog](https://github.com/linebender/resvg/blob/main/CHANGELOG.md),
[harumi v1.23.8](https://github.com/kent-tokyo/harumi/releases/tag/v1.23.8).

## Classification rule

An accuracy or compatibility difference is a correctness issue unless it can
cause confidentiality, integrity, availability, memory-safety, or supply-chain
impact. Resource exhaustion, panic on attacker-controlled input, arbitrary
filesystem/process access, unintended network access, and secret leakage are
security issues.
