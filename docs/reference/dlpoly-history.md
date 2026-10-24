# DL-POLY Classic HISTORY trajectories

The direct HISTORY APIs read, edit and write a bounded **formatted Classic
layout**: an 80-character title, a three-integer header, and one or more
`timestep` frames. All frames are retained, including their source atom order.
This is separate from CONFIG/REVCON, binary HISTORY and newer extended layouts.

## Retained fields

| Record | Fields |
|---|---|
| File header | Title, `keytrj`, `imcon`, atom count |
| Frame | Step number, integration timestep, optional ordered a/b/c cell vectors |
| Atom in each frame | Opaque label, ID, mass, source atomic charge, position, optional velocity and force |

`keytrj=0` carries positions, 1 adds velocities, and 2 adds forces. Presence is
validated exactly; absent data is never zero-filled. Massless shells are
allowed. Labels need not be elements or unique, and a label equal to `timestep`
does not become a frame delimiter inside an atom record. Every frame's IDs
uniquely cover 1..N, but source order can change.

Frame counts and keys must agree with the Classic file header. Repeated or
decreasing step numbers are retained in file order. Integration timesteps may
vary by frame; they are not interpreted as elapsed time, and no elapsed time
is reconstructed by multiplying step and timestep.

`imcon=0` has no cell records. Modes 1 and 2 require positive axis-aligned
cubic/orthorhombic cells. Mode 3 accepts a finite nonsingular parallelepiped,
retaining orientation, handedness and changes between frames. Coordinates are
not wrapped and cell vectors are not reconstructed from lengths/angles.

All scalar/vector values remain in **raw file units**. The layout itself does
not encode a unit scheme. FIELD/CONTROL context is needed to distinguish the
DL_POLY family's molecular and DPD conventions. No element, bond, formal
charge or unit conversion is inferred from labels, masses or atomic charges.

## APIs

```python
import chematic

trajectory = chematic.parse_dlpoly_history(open("HISTORY", "rb").read())
trajectory["frames"][1]["atoms"][0]["charge"] = -0.7654321
trajectory["frames"][1]["atoms"][0]["force"] = [3.0, 4.0, 5.0]
text = chematic.write_dlpoly_history(trajectory)
assert chematic.parse_dlpoly_history(text) == trajectory
```

Rust exposes `DlpolyHistory`, `DlpolyHistoryFrame`, `DlpolyHistoryAtom`,
`DlpolyHistoryLimits`, `parse_dlpoly_history`, `parse_dlpoly_history_with_limits`,
`write_dlpoly_history` and `dlpoly_history_from_json`. Python accepts text or
bytes and returns a dictionary. WASM/Node exposes `dlpoly_history_json` and
`dlpoly_history_from_json`. Writers validate current edited data and reject
unknown JSON fields at file, frame and atom levels.

```sh
chematic dlpoly-history -i HISTORY -o normalized.HISTORY
chematic dlpoly-history -i HISTORY --inspect -o trajectory.json
chematic dlpoly-history -i trajectory.json --from-json -o edited.HISTORY
```

The canonical writer retains Classic field widths: 10 columns per integer,
12 per real, eight per label; frame headers are 60 bytes, atom headers 42,
and vector records 36. Real fields use explicit decimal points and optional
exponents. An exact f64 representation must fit twelve columns. It can use
more digits than the original program's default display precision when those
digits fit. Values needing wider fields are refused, never silently rounded.
The reader also accepts whitespace-separated equivalents and Fortran D
exponents. Implied decimal points are refused because guessing them can change
the numerical value. This is field-preserving normalization, not byte identity.

## Refusals and limits

The profile excludes unformatted/binary and netCDF files; five-field modern
file headers; elapsed-time frame extensions; displacement/RSD atom extensions;
boundary modes 4–7; varying per-frame atom counts or header keys; name-only
atom headers; comments/unknown trailers; malformed/truncated records; and
nonfinite or underflowing values. Mass must be nonnegative and timestep positive.
Names are 1..8 printable ASCII nonspace bytes; the title is at most 80 bytes.
Integer step values must fit ten columns.

Default Rust limits are 64 MiB input/output, 4096 bytes per input line, 100000
atoms per frame, 10000 frames and 1000000 atom records across all frames.
JSON input is capped at 64 MiB. WASM further caps text input/output at 1000000
bytes, atoms at 10000, frames at 1000, total atom records at 100000, and JSON
input/output at 16 MiB. Expanded fixed-width text exceeding the reader cap
refuses before it is returned. The public record API materializes a bounded trajectory; it is not a
streaming interface for arbitrarily large simulation files.

## Evidence and comparison scope

The original packet in `validation/dlpoly_history/cases.json` shares 18 valid
and 47 refusal cases across Rust, Python, WASM and CLI. It covers all three
data levels, changing cells, massless shells, opaque/repeated labels, changed
atom order, repeated steps and variable timestep, adjacent fixed fields,
precision refusals, and truncated/extended data. Tests exercise edits and
unknown JSON fields; native checks cover exact record widths and signed zero.
WASM checks include the 1000/1001-frame and 10000/10001-atom boundaries.

Open Babel 3.2.1 advertises HISTORY as **read-only**. The committed comparison
therefore measures coordinate observation rather than a HISTORY round trip.
Both tools read the same seekable input file. chematic emits inspection JSON;
Open Babel uses `-iHISTORY input.HISTORY -ocjson -O frame.cjson -m`, producing
one CJSON file per frame. Six original nonperiodic C/O trajectories yielded
all-frame raw coordinates equal to the expected values in **6/6 cases for both
tools**. Frame output is collected in numeric order; no first-frame selection
is used. This validates bounded coordinate interoperability, not superiority.

The observer scores only raw coordinates, canonical C/O element identity and
frame count. It does not isolate parser behavior from CJSON serialization or
measure unit interpretation, masses, charges, velocities, forces, step metadata
or periodic cells. Those source fields are checked by chematic's separate
same-format record tests; absence in a chosen CJSON representation is not
claimed as an Open Babel HISTORY reader defect.

Commands, output-file receipts, input/source/binary hashes and case results
are in `validation/dlpoly_history/openbabel-3.2.1-comparison.json`. Reproduce:

```sh
python3 scripts/validation/check_dlpoly_history.py \
  --cli target/debug/chematic --obabel /path/to/obabel \
  --output /tmp/history-comparison.json
```

An independent GNU Fortran 16.1.0 probe also reads 961 real fields from all 18
valid cases plus three numerical boundary variants using each of the manual's
F12.6, E12.4 and G12.4 input descriptors. All 2883 reads retain the expected f64
bits, including signed zero. The probe, runner and receipt are committed as
`validation/dlpoly_history/real_field_probe.f90`,
`scripts/validation/check_history_fortran_fields.py` and
`validation/dlpoly_history/fortran-real-fields.json`. Reproduce with the runner's
`--cli`, `--compiler` and `--output` arguments.

These are toolkit I/O and real-field numeric tests; no DL_POLY simulation or
postprocessor executable was run to validate an end-to-end workflow.

## Primary references

- [DL_POLY Classic manual 1.10, section 4.2.1.1](https://www.nist.gov/system/files/documents/2025/04/29/DL_POLY_Classic_User_Manual.pdf)
- [Current DL_POLY output layouts](https://ccp5.gitlab.io/dl-poly/UserManual/data_files/output_files.html)
- [DL_POLY DPD unit conventions](https://ccp5.gitlab.io/dl-poly/UserManual/DPD/units.html)
- [Open Babel HISTORY format direction](https://openbabel.org/docs/FileFormats/DL-POLY_HISTORY.html)
