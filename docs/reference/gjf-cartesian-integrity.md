# Gaussian Cartesian geometry integrity

chematic provides a bounded Cartesian GJF/COM profile. It preserves explicit
coordinates and optional freeze codes, refuses unsupported geometry instead of
partially parsing it, and offers exact source retention when no edit is needed.
It is not a complete Gaussian route parser or a guarantee that a job will run.

## API overview

| Need | Rust | Python / CLI / WASM |
|---|---|---|
| Read one Cartesian geometry | `parse_gjf_geometry` | `chematic.parse_gjf_geometry(text)` |
| Read with custom limits | `parse_gjf_geometry_with_limits` | Python byte/atom options |
| Write explicit Cartesian geometry | `write_gjf_cartesian` | `chematic.write_gjf_cartesian(...)` |
| Retain one source exactly | `parse_gjf_source_record` / `write_gjf_source_record` | `roundtrip_gjf_record`, common self-conversion, CLI `convert` |
| Retain and inspect Link1 jobs | `parse_gjf_source_document` / `write_gjf_source_document` | `roundtrip_gjf_document`, `parse_gjf_document_json`, CLI `gjf-document` |

The legacy `parse_gjf` tuple uses the corrected geometry parser but omits
freeze codes. Use `GjfGeometry` when those codes matter.

## Geometry reader contract

Each atom row must have exactly one of these forms:

```text
element x y z
element integer-freeze-code x y z
```

Whitespace and comma delimiters, symbols, atomic numbers, numeric suffix labels
such as `Cl12`, LF/CRLF, and `!` comments are accepted. The result contains:

- atoms in input order, without inferred bonds or per-atom formal charges;
- three finite coordinates per atom in the input numeric units;
- one optional signed freeze code per atom;
- job charge and positive multiplicity.

```python
import chematic

text = "# HF/STO-3G\n\ngeometry\n\n-1 2\nC -1 1.25 2.5 3.75\n\n"
g = chematic.parse_gjf_geometry(text)

assert g["coords"] == [[1.25, 2.5, 3.75]]
assert g["freeze_codes"] == [-1]
assert (g["charge"], g["multiplicity"]) == (-1, 2)
```

The reader rejects short or extra-field atom rows, invalid freeze codes,
non-finite coordinates, non-positive multiplicity, Z-matrices, unsupported
atom annotations, fragment state, and multiple Link1 jobs. Use the document
API for multiple jobs.

Route keywords such as `Units=Bohr` are retained but not interpreted by the
geometry reader. Coordinates are not converted. Job charge is not distributed
over atoms.

## Checked Cartesian writer

`write_gjf_cartesian` accepts explicit elements, coordinates, optional freeze
codes, charge, multiplicity, and declared units:

```python
output = chematic.write_gjf_cartesian(
    g["elements"],
    g["coords"],
    charge=g["charge"],
    multiplicity=g["multiplicity"],
    units="angstrom",
    freeze_codes=g["freeze_codes"],
    route="# HF/STO-3G",
    title="Explicit Cartesian geometry",
)

back = chematic.parse_gjf_geometry(output)
assert back["freeze_codes"] == g["freeze_codes"]
```

Finite coordinates are written with enough precision to reproduce their f64
bits, including signed zero and finite subnormals. Element, coordinate, and
flag counts must match. Missing coordinates are never replaced with zeros.

The `units` argument describes the supplied numbers and emits exactly one
`Units=Ang` or `Units=AU` keyword; it does not rescale coordinates. Routes must
be ASCII, begin with one `#`, and must not add another units declaration.
Headers, arrays, and output are bounded.

The older GJF writer remains available for compatibility, but rounds
coordinates, omits freeze codes, and may fill missing coordinates. Prefer the
checked writer for integrity-sensitive output.

## Exact source retention

`GjfSourceRecord` pairs an immutable source string with its parsed Cartesian
projection. Writing the record returns the original UTF-8 bytes instead of
rebuilding the job.

```python
retained = chematic.convert_format(text, ".GJF", "com")
assert retained.encode("utf-8") == text.encode("utf-8")
```

The same path is available through CLI `convert` and WASM common conversion.
It retains Link0 directives, wrapped routes, titles, labels, comments, inline
sections, line endings, BOM, whitespace, and final-newline state.

Source retention validates the leading geometry and text framing. It refuses
malformed geometry, unsupported atom syntax, bare CR, forbidden controls,
Unicode line separators, and multiple jobs. Opaque route/additional sections
are copied without validating calculation semantics or opening external files.

## Link1 source documents

`GjfSourceDocument` validates every explicit Cartesian job separated by a
standalone `--Link1--` line and then retains the whole source exactly. No
earlier job is returned if a later job fails.

```python
import json

retained = chematic.roundtrip_gjf_document(text)
jobs = json.loads(chematic.parse_gjf_document_json(text))["jobs"]

for job in jobs:
    source = text.encode("utf-8")[job["source_start"]:job["source_end"]]
    print(job["index"], job["charge"], job["multiplicity"], source)
```

The inspection schema is `chematic.gjf-document.v1`. Each job reports ordered
elements, coordinates, freeze codes, charge, multiplicity, and byte range.
Checkpoint-only geometry, empty jobs, malformed delimiters, unsupported
geometry, and invalid framing are rejected. Routes and relationships between
jobs remain uninterpreted.

```bash
chematic gjf-document --input jobs.com --output retained.gjf
chematic gjf-document --inspect --input jobs.com --output jobs.json
```

## Default operation limits

| Surface | Input | Atoms | Jobs |
|---|---:|---:|---:|
| Rust | 64 MiB | 1,000,000 | 128 for documents |
| Python | 16 MiB | 100,000 | 128 |
| CLI | 64 MiB | 10,000 | 128 |
| WASM | 1,000,000 bytes | 10,000 | 128 |

Native limits also bound physical lines, line length, and section count.
Python and CLI expose selected overrides. These are per-operation limits, not
a process memory guarantee.

## Verification and scope

The validation packets cover numeric boundaries, freeze codes, unit labels,
exact source bytes, Link1 jobs, malformed input, Python, CLI, and real
Node/WASM transport. Reports are stored in:

- `validation/results/gjf-write-2026-10-22.json`;
- `validation/results/gjf-source-2026-10-22.json`;
- `validation/results/gjf-document-2026-10-22.json`.

Reproduction commands:

```bash
python3 scripts/validation/validate_gjf_cartesian.py --obabel /path/to/obabel
python3 scripts/validation/validate_gjf_write.py --obabel /path/to/obabel
python3 scripts/validation/validate_gjf_source.py --python-root /path/to/package
python3 scripts/validation/validate_gjf_document.py --python-root /path/to/package
cargo test -p chematic-mol --test gjf_cartesian_write
```

Open Babel 3.2.1 advertises GJF as write-only and its local reader rejects GJF
input. The comparison evidence therefore covers bounded writer controls, not
full-job interoperability or toolkit-wide superiority. Gaussian execution,
external resources, arbitrary route validity, and chemistry holdout quality
are outside this contract.
