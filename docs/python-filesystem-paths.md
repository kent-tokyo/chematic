# Python filesystem paths

Python file readers and writers accept `str`, `bytes`, `os.PathLike[str]` and
`os.PathLike[bytes]`, following the filesystem-path contract proposed in
[RDKit #9601](https://github.com/rdkit/rdkit/issues/9601).

```python
from pathlib import Path
import chematic

records = list(chematic.iter_sdf(Path('molecules.sdf')))
with chematic.SDWriter(Path('selected.sdf')) as writer:
    writer.write(records[0].mol)
```

## Covered APIs

| File operation | Accepted path parameter |
| --- | --- |
| `iter_sdf`, `iter_sdf_checked`, `iter_sdf_batched` | `path` |
| `iter_xyz_batched`, `iter_extxyz_batched` | `path` |
| `SdfDataset`, `SdfDataset.open` | `path` |
| `SDMolSupplier`, `SmilesMolSupplier`, `TDTMolSupplier` | `path` |
| `SDWriter`, `SmilesWriter`, `TDTWriter` | `path` |
| `Report.save` | `path` |
| `report` | `output` (or `None` for no file output) |

These are 15 entry points, with matching stub annotations. Existing string
paths retain their behavior. Chemistry algorithms and native formats are
unchanged; no Python native exports are added.

The shared converter uses Python `os.fsdecode`, including its `os.fspath`
protocol, then PyO3's native `PathBuf` conversion. Bytes use the filesystem
encoding, rather than UTF-8 chemistry-text decoding. POSIX surrogate-escaped
names can be represented losslessly, provided the filesystem permits them.
Paths are not expanded or resolved, and arbitrary objects are not stringified.

Invalid types and invalid `__fspath__` return values raise `TypeError`.
Exceptions from `__fspath__` propagate, and it is called once per path. Missing
parents and NUL paths remain `OSError`/`IOError`. Error strings use readable
path displays, which need not reproduce arbitrary filename bytes losslessly.
`report(output=None)` still means no file output.

## Source validation — 2026-10-16

Base: `2251a9f5d2e706087efaf55cb0f8ab1ff20b8016`. The before-change module was
retained from `b066f0d70fb952995ebe019ce7f47d9054fa66ea`; its four affected
pre-change I/O/report source files match the base exactly.

The old module passed 30 string acceptance cases and failed 75 bytes/Path/
PathLike cases. The new module passed all 105, including Japanese directories
and spaces. The selected run passed 390 tests and skipped 30, with no failures:
301 new path/stub checks and 89 existing I/O, SDF selection and input-error
checks. Controls cover callback exceptions, wrong types, one protocol call,
file contents, missing parents, NUL paths and stub annotations for all 15 APIs.

The 30 raw non-UTF-8 filename controls require a supporting POSIX filesystem.
On this macOS host, Python's own file creation raises `EILSEQ`; those controls
verify an OS error at the chematic boundary and skip actual file round trips.
Linux and Windows behavior has not been locally validated. Strict Clippy,
formatting, binding inventory, shared manifest, refactor budget, unsafe/naming
checks and strict MkDocs passed.

See `crates/chematic-py/tests/test_filesystem_paths.py` and the commands, raw
logs, hashes and limitations in
`validation/results/filesystem-paths-2026-10-16/`. This implements the filesystem
portion of RDKit #9601, not all wrapper requests or a package release.
