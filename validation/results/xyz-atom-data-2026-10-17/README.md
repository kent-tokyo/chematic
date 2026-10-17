# XYZ atom-data local evidence — 2026-10-17

Baseline `35b932f9` has no public parser/record for headerless XYZ+data.
Python baseline files were compared with Git; runtime probes use the reused
native extension and are not claimed as a fresh baseline native build.
`provenance.json` records source/artifact hashes and exact commands.

- `pytest.log`: 234 passed, zero failures/skips: 72 new / 162 existing.
  Tests cover widths 1/3/9, original text/CRLF, immutable values, unsupported
  schemas, extra/missing data, finite binary64 boundaries and explicit caps.
- `crossread.json`: ASE 3.26.0 and native extxyz each read 150 synthetic
  records (2,670 atom rows / 11,570 payload scalars). Coordinates/payloads
  agree as binary64 bytes, including negative zero, subnormal values and
  maximum finite payloads. A separate ASE check covers 118 element symbols.
  Opaque comments are JSON-tagged to prevent ASE scalar/type inference;
  native metadata retains the marker and is decoded explicitly in the check.
- The Open Babel 3.2.1 diagnostic retains original input/output: seven-column
  XYZ+vec is rewritten as four-column XYZ. No other Open Babel API/format
  is assessed by that observation.
- `typing.json` / logs: mypy 1.8.0 strict helper passes, valid consumer
  passes and four intentional public misuse cases yield four errors.
  Consumers use `--follow-imports=silent --ignore-missing-imports`; full
  legacy-package strict checking is not claimed.
- `source-gates.json`: refactor, source stub, native binding inventory,
  cross-binding manifest and unsafe checks pass. Native inventory is
  unchanged (213 Python functions / 317 WASM exports).
- `mkdocs.log`: strict build passed; existing unlisted pages are informational.
  `product-name.log`: canonical spelling passed.

The parser and serializer are standard-library Python. The native artifact
is reused for package import/old tests/cross-read only. ASE is installed in
a separate temporary directory, not a dependency. No fresh native build,
release, graph/physical-unit/tensor-layout validation, completed hosted CI or
speed result is represented by this evidence. Fixtures are original project
data under MIT OR Apache-2.0. See `docs/xyz-atom-data.md` for the checked profile.
