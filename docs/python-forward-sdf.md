# Python SDF readers for binary streams

`chematic.ForwardSDMolSupplier` and
`chematic.rdkit_compat.ForwardSDMolSupplier` read SDF records from a binary
object with a callable `read(size)` method. This adds the file-object input
missing from chematic's path-based suppliers, motivated by the unreadable
stream report in [RDKit #9601](https://github.com/rdkit/rdkit/issues/9601) and
its proposed [upstream fix #9615](https://github.com/rdkit/rdkit/pull/9615).

```python
import gzip
from chematic import rdkit_compat as Chem

with gzip.open("compounds.sdf.gz", "rb") as source:
    for mol in Chem.ForwardSDMolSupplier(source, strictParsing=False):
        if mol is not None:
            print(Chem.MolToSmiles(mol), mol.GetProp("_Name"))
```

The supplier validates `read` during construction without reading input.
Missing or non-callable `read` raises `ValueError`. Each read requests at most
8192 bytes; short reads, including one-byte reads, are accepted. `read(size)`
must return immutable `bytes`, with `b""` indicating EOF. Text streams and
other return types raise `TypeError`; oversized results raise `ValueError`.
Original Python lookup/read exceptions propagate with their identity intact.

Iteration reuses the native SDF reader and its default limits: 1 GiB input,
16 MiB record and physical line, and 100,000 nonempty records. Reads are
incremental with up to 8192 bytes of read-ahead. Native limit checks occur
after line/record collection; they are not a hard peak-memory ceiling for
an adversarial overlong line. No background reader or whole-file collection
is introduced. V2000/V3000, UTF-8 titles and SD fields use the diagnostic-complete
native SDF reader, including wedge and E/Z perception. The existing path-based
`SDMolSupplier` retains its separate lightweight mode; its V2000 wedge
perception is not a parity target for this addition.

The shared SD-field parser now accepts optional header information after the
field name, including RDKit's `> <Activity> (1)` record identifiers. Previously
that header was ignored and its value disappeared. The correction also applies
to existing path and in-memory SDF readers; the suffix is not part of the key.

With `strictParsing=True`, malformed graph records raise `ValueError`; the
next call can read the following record. With `False`, each malformed record
yields `None` in its original position. IO, invalid UTF-8 and resource-limit
errors always raise. A Python read failure ends iteration after the exception.

The supplier retains its source and read callable, including when the caller
releases its reference. Python cyclic garbage collection can release a source
that refers back to its supplier. The supplier never closes the caller's
stream. An emitted molecule remains usable after the supplier is destroyed.

This is a forward-only iterator: no restart, indexing, length or full RDKit
supplier-method compatibility is promised. As with chematic's existing
`SDMolSupplier`, `sanitize` and `removeHs` are accepted compatibility no-ops;
the native parser determines hydrogen/aromaticity behavior. Compression is
handled by the supplied file object, rather than filename detection.

## Local validation

The regression suite is `crates/chematic-py/tests/test_forward_sdf.py`.
It covers both APIs, missing/non-callable reads, exception identity,
short reads, Unicode/CRLF, EOF, gzip, source lifetime/cycles, malformed-record
positions, binary return types, V3000 parity and lazy construction.

`validation/rdkit_issues/forward-sdf/probe.py` isolates the missing-read RDKit
reproduction in a child with a two-second timeout. It also compares 80
molecule/property results across ten graphs, V2000/V3000, BytesIO/gzip and
native/compatibility APIs. Graph comparisons normalize both outputs through
RDKit, including isotope, formal charge, tetrahedral and E/Z controls.
These are local source checks, not wheel-publication or full RDKit parity
evidence. The upstream tracking issue includes other work outside this change.

The dated outputs are retained under
`validation/results/forward-sdf-stream-2026-10-17/`.
Local results: 255 Python tests (66 new stream cases and 189 existing I/O/
compatibility cases), 788 Rust library tests, 17 SDF integration tests and
80 successful RDKit-normalized molecule/property comparisons. Strict Clippy,
format, stub, binding-manifest and documentation checks also passed.
