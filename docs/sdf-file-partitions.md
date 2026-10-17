# Source-preserving SDF file partitions

`chematic.plan_sdf_split` creates exactly the requested number of nonempty,
contiguous SDF groups, balanced by record count. It implements the workflow in
[Open Babel #2147](https://github.com/openbabel/openbabel/issues/2147) without
rebuilding molecular records or loading the input file in full.

```python
from pathlib import Path
import chematic

plan = chematic.plan_sdf_split("large.sdf", 400, prefix="batch")
print(plan.record_count, plan.source_sha256)
# batch001.sdf through batch400.sdf; earlier groups receive any extra record.

output = Path("parts")
output.mkdir(exist_ok=True)
for index, part in enumerate(plan.parts):
    with (output / part.filename).open("xb") as destination:
        for block in plan.iter_part(index, chunk_bytes=65536):
            destination.write(block)
```

The API creates no output files. The example uses exclusive file creation and
writes one file at a time. Reading, verification or filesystem errors can leave
partial output; it is not a multi-file transaction. Consume each iterator fully
to complete its verification, or stage output in an application-managed temporary
directory before publishing it.

## Partition and preservation contract

For `N` records and `P` requested parts, `N // P` records go into each part and the
first `N % P` parts receive one extra. `1 <= P <= min(N, 10000)` is required;
requesting more parts than records raises instead of creating empty files.
Numbering starts at one, padded to the digit count of `P`. The prefix is a
1–64 character ASCII stem starting with a letter/digit and continuing with
letters, digits, underscores, periods or hyphens. It cannot be a path.

The public frozen result objects contain:

| Type | Fields |
| --- | --- |
| `SdfSplitPlan` | `source_path`, `source_bytes`, `source_sha256`, `record_count`, `unterminated_final_record`, `parts` |
| `SdfSplitPart` | `filename`, `byte_offset`, `byte_len`, `first_record`, `record_count`, `sha256` |

Record indices are zero-based. Byte intervals are contiguous, cover the whole
source and are disjoint. Concatenating all part streams in order reproduces the
original bytes exactly. Source CRLF/LF, titles (including blank titles), UTF-8 or
opaque SD values, duplicate properties, CTAB text, atom order, coordinates and
stereo annotations are never decoded or rewritten.

A record boundary is an exact `$$$$` physical line terminated by LF, CRLF or EOF,
consistent with the indexed native dataset reader. Embedded `$$$$` substrings
and whitespace-padded dollar lines remain data. Empty framed records and an
empty/whitespace-only file are refused. A nonempty final record without a delimiter
is retained and reported; trailing whitespace is attached to the last part.
This is **raw record framing**, not molecular validation. Malformed MOL blocks
and unsupported CTAB content are retained rather than skipped; use the existing
SDF readers to validate chemistry separately. Counts refer to framed records,
not necessarily successfully parsed molecules.

## Memory, source identity and limits

Planning scans the regular file once to find boundaries and compute its SHA-256,
then reads the planned ranges in a second pass to compute each part's SHA-256.
The two whole-source digests and source identity/stat data must agree. Planning
retains one 8-byte offset per record plus the part descriptors; it retains no
molecule graphs or complete file contents. `iter_part(index)` opens on first
iteration, seeks to the part offset and yields at most the requested chunk size.
The absolute source path is retained, so changing the working directory does
not retarget a plan. String, bytes and `os.PathLike` inputs are accepted.

Source replacement/observable mutation is checked before reading, and the part
hash and source stat are checked at iterator exhaustion. A late mismatch can
raise after some chunks have been yielded. Early iterator closure releases the
file but does not verify unread bytes; this API does not provide a file snapshot
or synchronization with external writers.

| Bound | Default and hard maximum |
| --- | --- |
| Input bytes | 1 GiB |
| Framed records | 1,000,000 |
| Physical line bytes, including its ending | 1 MiB |
| Record span bytes, including its delimiter | 64 MiB |
| Requested parts | 10,000 |
| Per-yield chunk bytes | 1 MiB (default); callers can choose a smaller positive size |

`max_input_bytes`, `max_records`, `max_line_bytes` and `max_record_bytes` may lower
the planning bounds. Integer parameters accept index integers and exclude Python
bool. Gzip input is refused explicitly; decompress before constructing a plain
file plan. Direct output writing, CLI splitting and compressed-file splitting
are outside this addition.

## Evidence

The dated record `validation/results/sdf-split-2026-10-17/` contains source/fixture
hashes, commands, tests and comparison results. An independent mixed V2000/V3000
corpus covers 52 input/partition configurations and 1,194 part streams. All source
bytes reassemble exactly, with 9,148 successful record reads each by native
chematic and RDKit 2026.09.1. Open Babel 3.2.1 reads all three candidate groups in
order; the reported `--splitinto 3` command still emits a single output on the
same seven-record source.

A separate 35,418,112-byte, 4,096-record input is partitioned and streamed with a
measured Python allocation peak of about 2.1 MB and maximum yielded blocks of
65,536 bytes. This is a `tracemalloc` observation, not process RSS or a general
performance claim. All 47 new contracts pass, including 400-way partitioning,
source mutation/replacement/hash failures and finite bounds.

This is a standard-library Python implementation. A recorded existing native
extension is used for package import and reader checks; no native rebuild or
Rust/WASM split implementation is claimed. The new helper passes strict mypy in
isolation; valid public usage and four expected invalid-usage diagnostics are
recorded without a package-wide strict typing claim.
