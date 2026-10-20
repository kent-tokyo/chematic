# Static XBS input — October 19, 2026

Open Babel [#2576](https://github.com/openbabel/openbabel/issues/2576), confirmed
open, requests the XBS dialect, which differs from Open Babel's `bs` format.
Source base `b7004777` already writes XBS scenes but has no reader. This adds
Python `read_xbs` / immutable `XbsScene`, independently from the author's
[README](https://github.com/Montmorency/XBS#readme) and the
[manual](https://manpages.debian.org/bookworm/xbs/xbs.1.en.html).

The bounded static gray-style profile retains exact decoded source and atom,
style and rule order with 1-based physical line numbers. Unsupported commands
are refused, not skipped; drawing distance rules are never chemical bonds.

## Results

- **271 pytest cases pass**: 62 new reader cases, 47 existing writer cases,
  162 existing compatibility cases. Tests include source/CRLF/line retention,
  deep immutability, canonical readback, forward references, duplicate rules,
  invalid syntax/commands, numeric underflow and positive/negative resource caps.
  Exact 4 MiB and 20,000-line inputs, 10,000 atoms, 128 styles and 1,024 rules
  are exercised without checking in large fixtures.
- **144 independently authored generated scenes / 8,016 atom rows** compare
  against hofmann **0.21.1**. Both original source and normalized output are
  independently parsed, plus a separate 10,000-atom scene: **289 cross-reads**,
  **865 style checks**, **576 rule checks**. All coordinates match binary64
  bytes, including negative zero and minimum subnormal. The packet tests both
  atom orders, LF/CRLF, general/scientific notation, comments, interleaving and
  forward style references. Source and original atom-line arrays also match.
- Scoped strict mypy passes for both XBS helpers and a public typed consumer.
  Four invalid calls/assignments yield four expected errors. The package root's
  81 pre-existing generic-parameter annotation errors are excluded from this
  scoped check; a whole-package strict pass is not claimed.
- Black, source-stub gate, product naming and pinned MkDocs strict build pass.

`independent.json` records versions, seed, packet digest and candidate source
hashes. New controls are authored project data, not copied manual coordinates;
`controls.json` records provenance, license, minimization, oracle and hashes.

## Reproduction and boundaries

With the current Python sources importable:

```sh
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.13 -m pytest -p no:cacheprovider -q crates/chematic-py/tests/test_xbs_input.py crates/chematic-py/tests/test_xbs_output.py crates/chematic-py/tests/test_rdkit_compat.py
python3.13 validation/validate_xbs_input.py --output /tmp/xbs-input-independent.json
python3.13 -m mkdocs build --strict --site-dir /tmp/xbs-input-docs
```

For local runtime verification the full current Python package was copied to
`/private/tmp/chematic-xbs-input-python`, excluding pre-existing `.so` and cache
files. A previously built 1.1.9 native extension was linked only for package
import and the unchanged compatibility tests; the new reader calls no native
functions. Exact native hash is in `environment.json`. The oracle package was
the existing `/private/tmp/chematic-xbs-oracle` installation, NumPy 2.5.3;
Python 3.13.6 / pytest 9.1.1. No native build or published wheel is claimed.

The mypy configuration uses `strict = True`, `follow_imports = skip`, and a
`[mypy-chematic] ignore_errors = True` section for the existing package-root
stub. Both `_xbs.py` and `_xbs_read.py` are explicit checked source targets;
the new reader and public consumer are not ignored. The good consumer assigns
all scene fields to their tuple types and serializes to `str`; negative calls
are `read_xbs(12)`, `XbsScene(None)`, integer `preserve_source`, and an assignment
to frozen `labels`. Logs retain the actual checks.

No generic `convert_format` alias, chemical graph, automatic bond detection,
unit conversion, movie/RGB/wildcard support, GUI/rendering, Rust/WASM API,
performance or release result is implied. Hosted CI/main integration are
separate gates. Open Babel's existing `bs` dialect is unchanged.
