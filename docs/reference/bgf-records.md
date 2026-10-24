# BGF atom and connectivity records

chematic reads and writes a bounded **nonperiodic BIOGRF 200** profile through
editable records. These APIs preserve source fields independently of the common
molecular graph. They do not infer elements, formal charges, missing bonds,
hydrogens, or force-field types.

## Supported fields

The standard 80-column `ATOM` / `HETATM` layout retains the record kind, source
ID, atom name, residue name, chain, residue identifier, Cartesian coordinates
in angstroms, force-field type, maximum covalent-bond capacity, lone-pair count,
and partial charge. Residue IDs are A5 strings, including nonnumeric IDs.
Maximum bond capacity is independent of the actual number of neighbors.

`DESCRP`, repeated `REMARK`, and `FORCEFIELD` values are retained, including
whether optional metadata was absent. Atom and connectivity row order and
neighbor order are retained. Source IDs need not be contiguous.

`CONECT` uses `(a6,12i6)`: one central ID and at most eleven neighbors. An optional
immediately following `ORDER` has one integer code per neighbor. Codes 0 and 1
both mean a single bond, 2 a double bond, and 3 a triple bond. The API preserves
explicit 0, explicit 1, and an absent `ORDER` distinctly. One-sided connectivity
is accepted; supplied reciprocal rows must agree on bond meaning. Unknown
references, self bonds, duplicate IDs/neighbors/centers, and conflicting orders
are refused.

The atom declaration, when supplied, must be:

```text
FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)
```

Omitting the two standard FORMAT declarations uses these defaults. The writer
emits both declarations. It normalizes padding and declaration placement;
this is a field-preserving writer, not a byte-identical source return route.

## APIs

```python
import chematic

record = chematic.parse_bgf(open("input.bgf", "rb").read())
record["atoms"][0]["partial_charge"] = -0.12345
record["atoms"][0]["lone_pairs"] = 2
text = chematic.write_bgf(record)
assert chematic.parse_bgf(text) == record
```

Rust exposes `parse_bgf`, `parse_bgf_with_limits`, `BgfRecord`, `write_bgf`, and
`bgf_record_from_json`. Python accepts text or bytes and returns a dictionary.
WASM/Node exposes `bgf_record_json(text)` and `bgf_from_record_json(json)`.
All JSON/dictionary writers refuse unknown fields at every record level.

```sh
chematic bgf -i input.bgf -o normalized.bgf
chematic bgf -i input.bgf --inspect -o record.json
chematic bgf -i record.json --from-json -o edited.bgf
```

## Refusals and limits

This profile excludes periodic `XTLGRF`, `CRYSTX`, fractional coordinates,
symmetry, cell-displacement records, restraints, arbitrary FORTRAN layouts,
extended atom columns, noninteger/aromatic bond-order dialects, and multiple
structures. Unknown or misplaced records refuse instead of being skipped.
Thus this addition does not cover every BGF dialect supported by other tools.

Atom IDs must be unique in 1..99999; text fields must fit their printable ASCII
columns; chains are uppercase or blank. Coordinates and partial charges must
be finite, fit their F10.5/F8.5 fields, and retain their exact numerical value
after formatting. Overflow and edits needing more than five decimal places
refuse. Implied decimal-point and exponent conventions are not guessed.

Default Rust limits are 16 MiB of input/output, 4096 bytes per input line, and
99999 atoms/connectivity rows. JSON input is bounded at 64 MiB. WASM additionally
limits text input to 1000000 bytes, atoms to 10000, and JSON/output to 16 MiB.
CLI `--max-input-bytes` can lower the input cap; writer field/output limits remain.

## Local comparison

The original packet in `validation/bgf/cases.json` contains 15 valid cases and
26 refusal cases, shared by Rust, Python, WASM, and CLI tests. It covers typed
atoms, source metadata, sparse IDs, source order, absent/explicit/zero ORDER,
one-sided connections, default formats, CRLF, and fixed-column boundaries.
It is a small diagnostic packet, not a representative survey of BGF files.

With Open Babel **3.2.1**, `obabel -ibgf -obgf -l1` produced observable BGF for
12/15 inputs; two inputs produced no molecule and one was rejected. Failures
to produce output are recorded separately from field changes. On the 12 common
observable cases, chematic retained all represented fields. Open Babel retained
the exact coordinates in 12/12; it changed the supplied atom names, lone-pair
counts and partial charges in all 12. For example, the first case's charges
`+0.12345/-0.12345` became `+0.27945/-0.27945` under this command.

The observer permits generated remarks and normalizes reciprocal connection
entries and ORDER 0/1 when comparing the indexed bond graph. It reports source
IDs and ORDER encoding separately. Source-absent metadata is excluded from
preservation denominators. Full command receipts, field-level observations,
binary/input/source hashes, and refusals are in
`validation/bgf/openbabel-3.2.1-comparison.json`. Reproduce with:

```sh
python3 scripts/validation/check_bgf_preservation.py \
  --cli target/debug/chematic --obabel /path/to/obabel \
  --output /tmp/bgf-comparison.json
```

This result establishes the declared profile and this packet only; it does not
establish Open Babel-wide format coverage or superiority.

## Format references

- [MSI Cerius² Modeling Environment, BGF records](https://www.chem.cmu.edu/courses/09-560/docs/msi/modenv/D_Files.html)
  (archival vendor manual; the original URL may now redirect).
- [ReaxFF user manual, examples 2.1–2.3](https://www.engr.psu.edu/ADRI/Upload/reax_um.pdf)
  documents BIOGRF 200 and distinguishes periodic cells and restraints.
- [VMD Biograf plugin notes](https://www.tcbg.illinois.edu/Research/vmd/plugins/molfile/bgfplugin.html)
  describes force-field types and the separate valence/lone-pair fields.
- [Open Babel's BGF format entry](https://openbabel.org/docs/FileFormats/MSI_BGF_format.html).
