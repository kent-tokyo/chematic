# Python stub packaging gate

[RDKit #9540](https://github.com/rdkit/rdkit/issues/9540) reports that malformed
shipped `.pyi` files can stop downstream type checking before user code is
analysed. The current chematic source stubs are valid; the missing capability
was a source **and packaged-wheel** gate in CI and the PyPI workflow.

`check_python_stubs.py` now fails if a type stub is malformed, omitted or stale.
It uses only Python's standard library and does not import chematic, NumPy,
compiled native modules or the stubs under test. It never rewrites broken stubs,
drops declarations or replaces signatures with permissive `*args` fallbacks.

```sh
python scripts/check_python_stubs.py
python scripts/check_python_stubs.py --wheel "dist/chematic-*.whl" --report stub-gate.json
```

Quoted glob patterns work on Linux, macOS and Windows. Every matched wheel is
checked; an unmatched pattern is a failure. A nonzero exit stops the calling
CI or publication step. Errors include path, line and column where available;
optional JSON reports are written on both success and failure.

## Checks

- Require the source package directory, `__init__.pyi` and `py.typed` marker.
- Parse every source `.pyi` with Python 3.9 grammar, then compile its AST to a
  code object without executing it. Compilation catches additional errors such
  as duplicate argument names that `ast.parse` alone can accept.
- Inspect every `.pyi` archive member, including any outside `chematic/`.
  Require `chematic/py.typed`, identical marker bytes, and exactly the source
  stub inventory under `chematic/`, with byte-identical contents.
- Reject corrupt or duplicate-member archives. Read stub/marker members without
  extracting a wheel or importing a binary extension. Individual files are
  capped at 4 MiB and total stub bytes at 32 MiB. The wheel SHA-256 is streamed.

The syntax policy uses the package's declared minimum Python version (3.9).
`ast.parse(feature_version=(3,9))` on a newer interpreter is a grammar check,
not an emulation of Python 3.9. Therefore CI runs the source checks and negative
controls on **actual Python 3.9 and 3.13** too. The two diagnostic phases,
`parse` and `compile`, are reported separately.

This checks syntax and package contents. It does not establish type-checker
semantic correctness, stub/runtime signature parity, native ABI compatibility,
wheel authenticity, or chemical behavior. The compiler check is deliberately
stricter than AST parsing; no claim is made that every downstream type checker
rejects every compiler-policy violation. The source distribution archive is
outside the wheel-content check.

## Workflow placement

`ci.yml` has a dedicated Python 3.9/3.13 source job. Existing built-source-wheel
lanes run the wheel gate immediately after maturin, before installing the wheel:
Linux Python 3.12, Linux Python 3.9, macOS and Windows.

`publish-pypi.yml` checks every Linux/macOS/Windows built wheel before uploading
its artifact. The final publication job independently checks all downloaded
wheels against the same selected source ref, immediately before the PyPI upload
action. These checks do not execute foreign-platform extensions. No release,
tag or PyPI publication is performed by this change.

## Local validation

- All **17 standard-library tests** pass on actual Python **3.9.6** and
  **3.13.6**. Tests cover malformed identifiers/default ordering/duplicate
  arguments, invalid encodings, newer grammar, missing or mismatched stubs and
  markers, corrupt and duplicate archives, budgets, glob handling, diagnostics
  and CLI exit codes. A side-effect sentinel confirms that stub code is never
  executed.
- An existing retained **chematic 1.1.5 macOS arm64 wheel** passes on both
  interpreters, with source-identical stub bytes. This is an artifact retained
  from an earlier build, not a freshly built or newly published candidate.
- Eight temporary mutations of that actual wheel are rejected with CLI exit 1:
  illegal annotation target, default ordering, duplicate argument, missing
  stub, missing marker, stale valid stub, changed marker and duplicate member.
  The original wheel is read only; temporary mutations are removed.
- A read-only audit of installed **RDKit 2026.09.1** finds two AST parse failures
  among 284 stubs. The stricter compiler phase finds five more violations:
  one duplicate-argument definition and four misplaced future imports. These
  additional diagnostics are local policy observations, not a claim that all
  five stop every type checker. No RDKit files are changed or copied.

Reproduce the CLI/archive experiment with an existing wheel and, optionally,
a pinned installed RDKit stub directory:

```sh
python -m unittest scripts.tests.test_check_python_stubs -v
python scripts/validate_python_stub_gate.py \
  --wheel /path/to/chematic-1.1.5-platform.whl \
  --rdkit-stubs /path/to/rdkit-stubs --rdkit-version 2026.09.1
```

Reports, source hashes and raw CLI diagnostics are retained in
`validation/results/python-stub-gate-2026-10-16`. Cross-platform hosted CI and
future release execution remain separate from these local observations.

The syntax/AST distinction follows the
[Python 3.9 AST documentation](https://docs.python.org/3.9/library/ast.html#ast.parse);
packaged type information follows the
[Python typing distribution specification](https://typing.python.org/en/latest/spec/distributing.html).
