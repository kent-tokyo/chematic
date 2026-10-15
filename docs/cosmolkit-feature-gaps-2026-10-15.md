# Selected COSMolKit feature gaps: indexed SDF workflows

Compared the published **COSMolKit 0.5.0 CPython 3.13 macOS universal wheel**
with the current isolated chematic source candidate, preserving earlier work
through GZMAT precision. Public docs and upstream source were inspected at
`d99c96c8abb21c45698826ef3d3d506d82d2f34a`; that inspected source revision is
not asserted to be the wheel's build commit.

This is a selected workflow audit, not an exhaustive feature count or a
claim of overall COSMolKit/RDKit parity. The
[COSMolKit SDF guide](https://kit.cosmol.org/python/io)
documents advanced dataset indexing and position-preserving batches; those
operations were also run using the pinned public wheel.

## Findings and implementation

| Capability | chematic finding | This candidate |
| --- | --- | --- |
| Indexed seekable SDF access, titles and byte ranges | Already present in Rust/Python `SdfDataset`; molecular records are parsed lazily. | Reused the existing file index and locked handle. No duplicate index. |
| Slices, integer lists, boolean-mask dataset selection | Python `__getitem__` accepted only a scalar index. Rust had only scalar and contiguous-range access. | Added checked ordered Rust `read_indices` and Python selection with negative indices, steps and duplicates. |
| Python indexed-file batches retaining invalid positions | Rust streaming `SdfBatchReader` already retains errors. Existing Python streaming batches skip rejected records, and the dataset had no batch method. | Added `SdfDataset.batches()` and `SdfSelectionBatch`, retaining source indices and failed slots. Existing streaming behavior remains compatible. |
| Mixed V2000/V3000 indexed reads | A native V3000 parser existed, but SDF graph/record/file readers called only V2000 parsing. | Dispatch the version on the fourth header line; send only the CTAB to the strict V3000 reader and retain SD fields separately. |
| Pickle, default Layered fingerprints, scalar SDF access | Already implemented. | No replacement or duplicate feature. |
| General Layered parameter/provenance controls | COSMolKit exposes `LayeredFingerprintParams` and `fingerprint_layered_with_params`; chematic's exact Layered function is a fixed default 2048-bit profile. | Remains a separate gap; not silently equated with the native `layered_fp` family. |
| Unified complete `BioStructure` hierarchy and mmCIF serialization | COSMolKit has the unified value. chematic has protein projection plus general mmCIF atom/residue records and explicit adapters, which are narrower interfaces. | Remains separate; no full structural-hierarchy equivalence claimed. |

The implementation is original; COSMolKit is used as a behavioral reference,
not copied source. Existing SD property parsing was moved into `sdf/fields.rs`
to keep the refactoring budget; no new debt allowance was added.

## Python usage

```python
import chematic

dataset = chematic.SdfDataset.open("library.sdf")
record = dataset[-1]                  # existing scalar result: SdfRecord
selected = dataset[[8, 2, 8]]         # order and duplicate index retained
reverse = dataset[::-1]              # slice and step semantics
masked = dataset[[True, False, True]] # mask must match the dataset length

# Explicitly retain rejected inputs; each slot keeps its original source index.
batch = dataset.select([8, 2, 8], errors="keep")
print(batch.source_indices())
print(batch.valid_mask())
print(batch.to_smiles_list())         # None where parsing failed
print(batch.errors())                 # position, source_index, code, message
for record in batch.records():
    if record is not None:
        print(record.name, record.get("supplier_id"))

for batch in dataset.batches(size=1024, errors="keep"):
    indices = batch.source_indices()
    valid = batch.valid_mask()
    # Keep indices/valid beside downstream feature values.
```

The example mask illustrates its shape requirement; it is only valid for a
three-record dataset. Indexing with a container defaults to raising on any
malformed record, matching the measured COSMolKit default. `select()` defaults
to `errors="keep"`. Both `select()` and `batches()` accept `"keep"` or `"raise"`;
`"skip"` is refused because it loses the positional contract.

`SdfSelectionBatch` is a separate result type rather than an alias for
COSMolKit's `MoleculeBatch` or chematic's existing valid-only `MoleculeBatch`.
It retains source indices, `SdfRecord` values, fields and parse diagnostics.
`len(batch)` includes rejected slots. `batch[i]` and `records()` return
`None` for a rejected slot; `errors()` contains one entry for each failed
occurrence, including duplicates. Returned index/mask/error/list containers
are independent snapshots. Parsed duplicate occurrences are separate records;
changing a molecule property in one does not change the other or a future read.

Selectors accept Python integers, slices, lists and tuples. Scalar booleans
and mixed integer/bool selectors are rejected. Boolean masks must have exactly
one value per source record. NumPy arrays must be explicitly converted with
`.tolist()`; no generic iterator or implicit ndarray coercion is implemented.
Python's slice normalization determines clipping, negative steps and empty
ranges. All indices are normalized and validated before file I/O.

A dataset-batch iterator owns a shared reference to the indexed file and
survives deletion of the Python dataset wrapper. It creates indices only for
the next chunk. `len(iterator)` is the remaining planned batch count, including
malformed records. `cancel()` stops subsequent reads. Exhaustion, cancellation,
I/O failure or a raised malformed batch leaves the iterator exhausted.

## Rust and limits

`SdfDataset::read_indices(indices)` returns `SdfIndexedBatch`:

- `indices: Vec<usize>` preserves the exact requested ordering and duplicates.
- `records: Vec<Result<SdfRecord, MolParseError>>` retains a result per request.
- `len()`, `is_empty()` and `valid_mask()` describe the returned positions.

`read_indices_with_limits(indices, &SdfSelectionLimits)` accepts explicit
record-count and cumulative source-byte limits. Defaults are **10,000 selected
positions** and **64 MiB of source-record bytes**, counting duplicates each
time. All bounds and cumulative sizes are checked before any seek/read.
An invalid index or exceeded budget returns `SdfSelectionError`; actual I/O
and UTF-8 read errors abort the operation. Molecular parse errors stay in the
per-position result vector. The existing scalar and range APIs retain their
contracts.

Python selected reads retain the 10,000-position cap. Batch size must be
1–10,000. Both selection and batching expose `max_bytes` (default 64 MiB)
with strict positive-integer validation. A batch exceeding its byte budget
fails; it is not silently split or truncated. Source-byte limits constrain
input work, not exact process RSS or the allocator overhead of parsed objects.

The backing file must remain unchanged while indexed. The file handle uses
the existing mutex, so concurrent reads do not interfere with each other's
seek offsets. This is not a filesystem snapshot or a parallel parser.

## CTAB boundary handling

V3000 MOL parsing is strict about its terminator and accepts no trailing SD
fields. The SDF adapters now identify the CTAB boundary after the three header
lines/counts line, pass that portion to the V3000 parser, and parse SD fields
from the remaining record. Literal `M  END` in the title or comment does not
end the CTAB. V2000's existing permissive terminator recognition is retained.

This dispatch covers `SdfReader`, `SdfRecordReader`, normal/fast
`SdfFileReader` and the dataset paths that use them. It does not assert new
parameter parity for every other MOL/SDF convenience function. Existing
version-specific coordinate and stereo diagnostics are reused.

## Validation and reproduction

Original fixture provenance, redistribution terms, expected values, oracle
settings and normalized payload hash are in
`validation/cosmolkit_comparison/sdf_selections/manifest.json`.
The minimal corpus has four one-atom valid records (including V3000) and one
malformed record, plus SD fields. It does not reuse a competitor example or
third-party molecular dataset.

`python scripts/validate_sdf_selections.py` compares the public COSMolKit
wheel with this extension across LF/CRLF and terminated/unterminated variants.
It exhausts the chosen slice grid, short integer-list products, all boolean
masks and all five-record permutations. Python list/slice semantics and the
explicit fixture expectations provide separate selection controls. The shared
indexing boundary raises on malformed selected input in both engines; the
additional chematic `select(errors="keep")` result is checked against the
explicit slot expectations. Keeping batches are compared directly in both
engines. Engine-specific exception messages are not required to match.

The run passed **12,828** selection cases: **9,028** successful shared
indexing results and **3,800** matching malformed-input rejections. It checked
**22,312** selected positions, **56** keeping batch values and **16** scalar
metadata/field controls across the four physical-file variants. Order, masks,
SMILES, fields and input bytes had no mismatches.

Local native validation passed **817 Rust tests** (753 library, 10 new
selection tests, and 54 coordinate/stereo/semantic-interchange regressions).
Python passed **68 tests** (41 new selection tests and 27 existing I/O tests).
Clippy with warnings denied, formatting, strict MkDocs, stub syntax, static
binding inventory and refactoring budget checks passed. This does not claim
a full-workspace or hosted CI run.

The evidence is saved as
`validation/results/cosmolkit-sdf-selections-2026-10-15.json`, with source,
compiled-extension and COSMolKit package hashes. These are local source
measurements, not release, speed or general chemical-parity evidence.

```sh
cargo test -p chematic-mol --lib --test sdf_selections \
  --test v3000_semantic_interop --test conformer_3d_io \
  --test stereo_reader_integration --offline --quiet
cargo clippy -p chematic-mol -p chematic-py --all-targets --offline --no-deps -- -D warnings
cargo fmt --all -- --check
python -m pytest crates/chematic-py/tests/test_sdf_selections.py crates/chematic-py/tests/test_io.py -q
python scripts/validate_sdf_selections.py
python scripts/check_refactor_budget.py
python scripts/check_binding_surface.py
python scripts/check_product_names.py
mkdocs build --strict
```

The validator requires COSMolKit 0.5.0 and an extension built from this source.
The review patch is
`validation/cosmolkit_comparison/sdf-selections-2026-10-15.patch`.
It requires the inherited candidate context through GZMAT precision. Main
integration, feature commit, push and publication are separate steps.
