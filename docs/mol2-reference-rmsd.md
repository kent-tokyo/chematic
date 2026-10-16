# All-reference MOL2 RMSD and source-preserving records

[Open Babel #2671](https://github.com/openbabel/openbabel/issues/2671) requests
RMSD for every reference conformation against one test structure. chematic's
single-record MOL2 readers rejected such concatenated input. The new Rust and
Python APIs make that workflow explicit and retain each record's identity.

```python
import json
from pathlib import Path
import chematic

references = Path("references.mol2").read_text()
probe = Path("probe.mol2").read_text()
report = json.loads(chematic.mol2_reference_rmsds_json(references, probe))
for row in report["results"]:
    print(row["reference_index"], row["reference_name"], row["rmsd"])

# Opt into rigid alignment; the default compares the original coordinates.
aligned = json.loads(chematic.mol2_reference_rmsds_json(references, probe, align=True))
records = json.loads(chematic.parse_mol2_records_json(references))["records"]
assert "".join(record["source"] for record in records) == references
```

## Correspondence and metric

Every reference is compared with exactly one probe, in file order. Results
include a zero-based reference index, both names, heavy-atom count and RMSD.
Hydrogen positions are excluded from the metric. Coordinates use the input
MOL2 unit, normally Å; no unit conversion is inferred.

Correspondence is deliberately checked: all atoms, including hydrogen atoms,
must have the same row order, source IDs, names, types, elements and formal
charges. Source bond endpoint IDs and bond types must match; bond row order
and reversed undirected endpoint spelling may differ. Partial charge values,
record names, substructure labels and opaque metadata do not define atom
correspondence. Identical element sequences alone are insufficient. Different
atom ordering, hydrogen counts or typing requires a separate explicit mapping.

Raw RMSD is `sqrt(sum(||reference_i - probe_i||²) / heavy_atom_count)`.
`align=True` uses the existing Kabsch implementation for a proper rotation
and translation over the paired heavy atoms. Reflections are excluded.
This API performs no symmetry/isomorphism search, graph normalization, docking,
energy optimization or selection of a best conformation.

Empty references, hydrogen-only probes, mismatches, malformed records,
resource violations and non-finite calculations are errors. Python raises
`ValueError` and returns no partial JSON. Rust returns `Mol2RmsdError`.

## Record boundaries, preservation and limits

Rust `Mol2SourceReader::new(text)` lazily parses a borrowed in-memory input.
`with_limits(text, Mol2StreamLimits)` bounds total input bytes, record count,
and each record through `Mol2ParseLimits`. This is not a `BufRead` file reader:
the full input already exists in memory and each successful record owns its
source. `mol2_reference_rmsds(reference_text, probe_text, align, &limits)`
materializes the result list, while holding one parsed reference at a time.

Only a complete, trimmed, case-insensitive `@<TRIPOS>MOLECULE` line starts a
new record. Preamble belongs to the first record; intervening comments and
whitespace belong to the preceding record. CRLF, opaque sections, atom types,
charges and source text are retained. Empty/whitespace input produces zero
inspection records; nonempty input without a marker is invalid. Single-record
APIs retain their rejection of multi-record input.

The Rust iterator emits an indexed error for a malformed record, then permits
the caller to continue at the next known boundary. Global input/count errors
end iteration. The Python inspection function and RMSD function stop on every
error rather than silently omitting a record. Inner parse line numbers are
record-local, and the record index is zero-based.

Python accepts keyword-only `max_input_bytes` (64 MiB), `max_records` (10,000),
`max_record_bytes` (16 MiB), `max_atoms` (100,000) and `max_bonds` (200,000).
Atom/bond limits apply per record, including the probe. Other native record
limits retain their defaults. Both Python APIs materialize their JSON output.
These are UTF-8 text APIs, not filesystem-path or binary-input APIs.

## Upstream diagnostic and validation

The original five-record fixture has translations of lengths 0, 1, 2, 3 and
5. Thus every paired heavy atom has a hand-derived raw RMSD of 0, 1, 2, 3 or
5 respectively, and rigid alignment removes each translation.

Open Babel 3.2.1 produced one line for `obrms references.mol2 probe.mol2`.
However, its [own manual](https://github.com/openbabel/openbabel/blob/master/doc/obrms.1)
documents the default as pairing the nth records of both files. Reversing the
inputs and using `obrms -f probe.mol2 references.mol2` returned all five values;
adding `-m` returned five aligned zeros. The one-line diagnostic therefore
does **not** establish an upstream parser defect. This addition implements the
issue's requested workflow through an explicit chematic API; it is not a
claim that Open Babel lacks a working invocation.

Original fixtures and the independent runner are under
`validation/openbabel_issues/mol2-rmsd/`. Local results, baseline/candidate
extension hashes and source hashes are recorded under
`validation/results/mol2-rmsd-2026-10-16/`. Tests cover all records, malformed
middle records, identity/order/charge/bond mismatches, resource boundaries,
hydrogen exclusion, source retention, rigid transforms and finite arithmetic.
No upstream source or attachments were copied. Native source validation is
separate from hosted CI and package publication; no speed or universal
chemical-parity claim follows from these small fixtures.

The local candidate passed 43 Rust tests (20 new and 23 existing) and 78 Python
tests (32 new and 46 existing). Thirty-nine generated rigid/distorted/mirrored
four-heavy-atom controls matched direct raw-distance calculations within
`1e-10` and independent NumPy 2.5.3 SVD alignment within `1e-9` (observed maximum
aligned difference `8.46e-14`). RDKit 2025.09.3 alignment with explicit atom
maps, `reflect=False` and `maxIters=1000` was checked separately at `1e-5`;
its observed maximum difference was `8.71e-7`, including a small iterative
residual on a near-zero rigid control. Twenty-three changed atom-order
permutations were rejected. These are fixed-fixture numerical checks, not
general equivalence to either tool.

```sh
PYTHONPATH=/path/to/candidate-package python \
  validation/openbabel_issues/mol2-rmsd/validate.py \
  --baseline-package /path/to/baseline-package --obabel /path/to/obabel \
  --obrms /path/to/obrms --output /tmp/mol2-rmsd-validation.json
```
