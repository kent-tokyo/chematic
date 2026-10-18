# MOPAC keyword continuations — 2026-10-18

Native implementation of the standalone `+` portion of Open Babel #1231.
The official OpenMOPAC specification says a `+` announces another keyword
line, followed eventually by two description lines. `++` and `&` have distinct
semantics and remain unsupported outside double-quoted values. Implementation
and synthetic fixtures are original; no Open Babel/MOPAC code or attachments
were copied.

## Evidence

- `before.log`: all 11 initial regression tests fail at the unmodified native
  baseline. The exact byte-limit assertions were strengthened afterward.
- `rust-tests.log`: 840 passing tests (815 library, 7 existing Cartesian,
  7 existing archive, 11 continuation). The old blanket '+' rejection assertion
  is removed; incomplete blocks remain rejected by the new tests.
- `after.log`: 36 freshly compiled native cases, 1–12 keyword lines with '+'
  at the beginning, middle or end of each nonfinal line. Original input and
  serialized output retain keyword text, title/comment and two-atom geometry.
- `oracle-summary.json`: 72 explicit single-line normalization controls are
  cross-read by Open Babel 3.2.1, checking C/O identity and 432 coordinate
  values. Six direct ordinary single-line reads also pass.
- `direct-openbabel-diagnostics.json`: all 66 direct multiline input/output
  reads yield no molecule with an invalid-geometry diagnostic. All still exit
  0, so process success alone would give false conversion-success evidence.
  These failures are not counted as positive format interoperability.
- `clippy.log`, `source-gates.log`, `mkdocs.log`: warnings-denied Clippy,
  formatting, complexity/source checks and strict docs pass.
- `guide-and-binding-check.log`: the guide's Rust example is compiled against
  the new native rlib and executed. Python test/stub syntax parses as 3.9;
  source inspection confirms the existing binding delegates unchanged.

Limits are 64 physical keyword lines and 65,536 combined LF-block bytes, plus
existing per-line/input/geometry limits. Tests include exact byte/count limits,
blank final keyword lines, CRLF, labels/flags, quoted markers, signed keyword
values, writer edits and physical error lines. New Python parametrized tests
are included for a newly rebuilt extension, but are not claimed as locally
executed. No Python/WASM binary was rebuilt. This is geometry/header I/O,
not a quantum-chemistry or full-format conformance claim.

## Reproduce

```sh
export CARGO_TARGET_DIR=/private/tmp/chematic-mopac-continuation-target
export CARGO_BUILD_JOBS=1 CARGO_PROFILE_DEV_DEBUG=0 CARGO_PROFILE_TEST_DEBUG=0 CARGO_INCREMENTAL=0
cargo test -p chematic-mol --offline --lib --test mopac_cartesian_contract --test mopac_archive_contract --test mopac_keyword_continuation
cargo test -p chematic-mol --offline --test mopac_keyword_continuation -- --nocapture > validation/results/mopac-keyword-continuation-2026-10-18/after.log 2>&1
python3.13 validation/results/mopac-keyword-continuation-2026-10-18/oracle.py
cargo clippy -p chematic-mol --offline --lib --test mopac_keyword_continuation --no-deps -- -D warnings
python3.13 -m mkdocs build --strict --site-dir /private/tmp/chematic-mopac-keyword-site
```

Rust dependencies must be available locally for `--offline`. The reference
script uses `/opt/homebrew/bin/obabel`; substitute the matching local binary
if necessary. Source/tool/artifact identities are in `provenance.json`.
Hosted CI is a separate integration gate.
