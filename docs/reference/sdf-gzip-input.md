# Per-input gzip selection for SDF

`iter_sdf_checked` accepts an explicit compression policy:

```python
import chematic

for index, record, error in chematic.iter_sdf_checked(
    "library.data", compression="auto", max_input_bytes=64 * 1024 * 1024,
):
    if error is not None:
        raise ValueError(f"record {index}: {error}")
    print(record.smiles, record.name, record.properties())
```

- `plain` is the default and reads bytes unchanged, preserving existing calls.
- `gzip` requires gzip input.
- `auto` recognizes the two-byte gzip magic independently for each input.
  A `.gz` suffix is neither required nor used for detection.

Rust callers wrap any `Read` source with `SdfInput::new(source, policy)` and
pass it directly to `SdfFileReader::new` or `with_limits`. `SdfCompression`
selects `Plain`, `Gzip` or `Auto`; `SdfInput::compression()` returns the resolved
mode. Auto probing consumes at most two bytes and replays them to the reader,
including short or non-seekable sources. Explicit modes do not probe.

Each wrapper owns its compression state. Reading another source constructs
another wrapper; full reads, partial reads and errors cannot carry a prior
gzip detection into it. There is no global flag or shared conversion session.

The decoder supports concatenated gzip members as one byte stream, including
a physical line or UTF-8 character split between members. It validates headers,
payloads and trailers. A corrupt/truncated member or trailing non-gzip bytes
yields a terminal IO-error entry through `iter_sdf_checked`. Header-probe IO
failures raise `IOError` at construction. Unknown modes raise `ValueError`, and
non-string modes raise `TypeError`.

Read through EOF to verify the entire compressed input. A trailer error can
arrive after earlier records have been emitted; these records are not rolled
back. Stopping early does not establish gzip integrity.

The native SDF parser retains source order, names, SD fields and its existing
stereo perception. Malformed records occupy their own error positions; later
valid records can still be read. Input, record, physical-line and record-count
limits apply to decoded data. Default limits remain 1 GiB total input, 16 MiB
per record/line and 100,000 records. The inherited line/record checks happen
after collection and are not a hard peak-memory ceiling. This addition does
not provide ZIP, zlib-wrapped input, compressed output, or compression support
for other file APIs.

## Upstream report and scope

[Open Babel #3030](https://github.com/openbabel/openbabel/issues/3030) reports
an auto-detected gzip flag persisting when an `OBConversion` replaces its input.
An original C++ probe reproduced this on Open Babel **3.2.1**: the first gzip
SMILES stream read successfully, but subsequent plain string/stream inputs
failed; equivalent fresh conversions succeeded.

chematic has no matching shared `OBConversion` session, and the prior checked
SDF API had no gzip policy. This adds the missing compressed-input capability
with independent per-source selection; it does not claim to fix an existing
chematic sticky-session defect or change Open Babel itself. The prior API's
source matched the retained baseline binding and rejected the new compression
keyword in all six initial path/protocol cases.

The probe is `validation/openbabel_issues/sdf-gzip/probe.cpp`; dated results,
source/artifact hashes and test logs are retained under
`validation/results/sdf-gzip-input-2026-10-17/`. These are local source checks,
not package-publication or completed hosted-CI evidence.

Local validation passed 879 Rust tests (including nine new compression tests),
74 Python tests (including 33 new cases), and 80 plain/compressed record
comparisons. One existing optional OGB reference test was skipped because the
`ogb` package was absent. Strict Clippy, formatting, source stubs, refactor
budgets, binding/manifests, naming, unsafe-surface and strict MkDocs checks passed.
