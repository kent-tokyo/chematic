# Position-preserving Standard InChIKey batches

`inchi_keys_checked` returns one indexed record for each input in an accepted
batch, including failures. Open Babel
[#2231](https://github.com/openbabel/openbabel/issues/2231) reports that failed
InChI conversions disappear from output, making input/output rows ambiguous.

```python
import chematic as ch

sources = ["InChI=1S/H2O/h1H2", "not-an-inchi", "InChI=1S/CH4/h1H4"]
rows = ch.inchi_keys_checked(iter(sources))
assert [r.index for r in rows] == [0, 1, 2]
assert [r.source_text for r in rows] == sources
assert [r.key for r in rows] == [
    "XLYOFNOQVPJJNP-UHFFFAOYSA-N", None, "VNWKTOKETHGBQD-UHFFFAOYSA-N"
]
assert rows[1].status == "invalid"
assert rows[1].code == "inchi_prefix"
```

The tuple and its `InchiKeyRecord` objects are immutable. Each record has
`index`, exact `source_text`, `status`, optional `key`, optional `code` and
optional `message`. Order and duplicates are preserved; no normalization or
whitespace stripping is applied. Messages have at most 1024 characters.
For a tabular export, `[r.key for r in rows]` has the original length, and
applications can explicitly render `None` as `NA` while retaining status/code.

| Status | Meaning |
| --- | --- |
| `ok` | Exact Standard InChI readback and native IUPAC key generation succeeded. |
| `invalid` | The identifier fails the accepted text/formula checks or the native parser reports an input error. |
| `unsupported` | Nonstandard identifiers, changed readback or native unsupported/ambiguous input; not a declaration that the source is chemically invalid. |
| `resource_limit` | The individual identifier exceeds byte, numeric or formula limits. |
| `conversion_error` | Native InChI generation/key output failed or the feature is unavailable. |

## Identity boundary

The API reads with `from_inchi`, regenerates `Mol.standard_inchi`, and requires
that string to equal the original input exactly before accessing
`Mol.standard_inchikey`. This matters because the existing reader may fall back
to the pure-Rust parser when the native reader declines an input. Successful
parsing alone therefore does not prove that the original identifier survived.

Only Standard `InChI=1S/` identifiers are supported. A valid source may be refused
because of noncanonical spelling or a limitation in the existing molecular
readback. No key for an altered identifier is returned. The API does not use
`Mol.inchikey`, which is the separate nonstandard Rust approximation. Native
InChI support is required for successful standard-key generation.

## Limits and failure handling

`max_records` is an integer from 1 through 10000, default 1000; bool is refused.
Sources must yield plain strings. All source items and aggregate input size are
checked before native conversion. Wrong types, iterator errors, aggregate input
above 16 MiB or a record-count overflow abort the batch without a partial tuple.
Count overflow consumes at most one extra source item. Iterator side effects
are not rolled back.

Each identifier is limited to 8192 UTF-8 bytes. Numeric fields may not exceed
4096, and the formula may contain at most 4096 atoms, including hydrogens and
component multipliers. These conservative checks run before native parsing.
Unknown formula syntax, whitespace/control/non-ASCII data and invalid Unicode
are explicit failure rows. The formula check is a resource guard, not a new
chemical parser. Native parser input errors and native `RuntimeError` become
row results; `MemoryError` and unexpected exceptions propagate. No files are
written and no molecule objects are retained in the result.

## Evidence

The mixed reference packet preserves all 219 positions. Of 144 valid identifier
rows, 142 produce keys matching RDKit 2026.09.1's direct `InchiToInchiKey` call.
Two rows, duplicates of one source, are refused because readback changed.
All failures remain in the denominator. The separate Open Babel 3.2.1 control
emits two keys for three input lines even with `-e` to continue after the
malformed identifier, reproducing the positional-loss behavior.
191 Python tests pass, including 39 new batch tests. Reproduction and artifact
boundaries are in `validation/results/positioned-inchi-keys-2026-10-19/`.
