# Checked PDB model splitting

`chematic.split_pdb_models` partitions an in-memory PDB coordinate document into
complete model documents and names them using the digit count of the total.
Ten models yield `split01.pdb` through `split10.pdb`; one hundred yield
`split001.pdb` through `split100.pdb`. It implements the naming workflow in
[Open Babel #2476](https://github.com/openbabel/openbabel/issues/2476).

```python
from pathlib import Path
import chematic

text = Path("ensemble.pdb").read_text(encoding="ascii")
result = chematic.split_pdb_models(
    text, prefix="nmr", allow_metadata_loss=True,
)
print([(loss.line, loss.record_name) for loss in result.losses])

# The complete result has already passed validation. File creation is explicit;
# exclusive mode refuses an existing destination. Use a separate output folder.
output = Path("models")
output.mkdir(exist_ok=True)
for item in result.files:
    with (output / item.filename).open("x", encoding="ascii", newline="\n") as file:
        file.write(item.pdb)
```

The API itself performs no file I/O. The example's writes are individual file
operations; an I/O failure does not roll back files written earlier.

## Result and identity

The public frozen dataclasses are:

| Type | Fields |
| --- | --- |
| `PdbModelSplit` | `source_text: str`, `files: tuple[PdbModelFile, ...]`, `losses: tuple[PdbSplitLoss, ...]` |
| `PdbModelFile` | `filename: str`, `model_num: int \| None`, `pdb: str`, `atom_count: int` |
| `PdbSplitLoss` | `line: int` (one-based original input), `record_name: str` |

Outputs follow source order. Their numbering starts at one, independently of
source MODEL numbers. Distinct positive four-column MODEL numbers, including
nonconsecutive numbers, are retained verbatim in each output. An implicit single
model has `model_num=None`. All explicit MODEL records require ENDMDL; an empty
model, duplicate number, nesting, mixed implicit/explicit coordinates, a model
after global CONECT, or nonblank content after END is refused. A final input END
is optional; outputs contain exactly one final END.

ATOM, HETATM, ANISOU and TER line contents, serials, atom names, alternate locations,
residues, insertion codes and coordinates are retained. LF replaces CRLF, blank
lines are omitted, and a missing terminal newline is supplied. No coordinate
rounding, serial renumbering, hydrogen addition or bond perception occurs.
Models need not have matching atom counts or row order; the result does not assert
that they are conformers of the same graph.

Global CONECT rows after all explicit models are copied into every file. Every
reference must exist in every model. Missing references are refused rather than
filtered. Model-local CONECT is outside this profile. Each output passes the
existing native [PDB selection validator](pdb-atom-selection.md), including finite
coordinates, atom serial uniqueness, ANISOU identity and TER consistency.

## Metadata and bounds

Only coordinate records, MODEL/ENDMDL, global CONECT and END are emitted.
Other records—including HEADER, NUMMDL, SEQRES, CRYST1 and MASTER—require
`allow_metadata_loss=True`. Each omitted record is listed once using its original
line number, and the exact input remains in `source_text`. Omitting them avoids
copying stale model totals, record counts or whole-entry metadata into every
single-model file. This is a coordinate-only splitter, with explicit metadata-loss
consent, rather than a general PDB metadata rewrite.

The prefix is a 1–64 character ASCII stem: an initial letter/digit, followed by
letters, digits, underscores, periods or hyphens. Paths and hidden-file stems are
refused. Generated names are distinct within a result.

Input is bounded at 64 MiB, 1,000,000 records and 1,024 printable ASCII bytes per
line. Native validation limits each model to 99,999 atoms. `max_models` defaults
to 9,999 and `max_output_bytes` to 64 MiB; callers may lower these positive integer
limits. Bool limits are refused. Aggregate output size includes repeated CONECT
rows and is checked before constructing output strings. All models must pass
before a result is returned; partial results are never yielded.

This addition is a typed Python orchestration API around the existing Rust
validator. It adds no Rust/WASM splitter or command-line split option.

## Validation

The dated record in `validation/results/pdb-model-split-2026-10-17/` contains the
upstream issue snapshot, native artifact provenance, test/type logs and comparison
results. The independent validator covers 80 synthetic inputs, 7,530 retained atom
lines and eight Open Babel XYZ cross-reads. The issue's RCSB 1NER input yields 30
models of 1,187 atoms, named `split01.pdb` through `split30.pdb`; all 35,610 ATOM
and 30 TER lines match the input exactly. Its 222 metadata records require explicit
loss consent. Open Babel 3.2.1 still produces unpadded names for that same input.
The fixture is fetched for validation and identified by hash; it is not vendored.
