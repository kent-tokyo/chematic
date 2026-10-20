# FHI-aims geometry.in

This independently implemented Rust/Python reader and writer fills a format
gap confirmed in the previous source candidate. Open Babel exposes the format
as [`fhiaims`](https://openbabel.org/docs/FileFormats/FHIaims_XYZ_format.html).
The grammar reference is FHI-aims manual **261002**, particularly its
[geometry definitions](https://fhi-aims.org/uploads/manual/Ch3/S2.html),
[relaxation constraints](https://fhi-aims.org/uploads/manual/Ch3/S12.html), and
[initial density](https://fhi-aims.org/uploads/manual/Ch3/S5.html) sections,
retrieved on 2026-10-15. No Open Babel, ASE or FHI-aims source is incorporated
into the implementation, and the core gains no external runtime dependency.

## Python

```python
import chematic

text = chematic.write_fhiaims(
    ["Si", "Cl"],
    [[0.5, 0.25, -0.5], [1, 2, 3]],
    coordinate_kinds=["fractional", "cartesian"],
    lattice_angstrom=[[4, 1, 0], [-1, 5, 2], [1, 0, -6]],
    constraints=[[False, True, False], [False, False, False]],
    initial_charges=[0.25, None],
    initial_moments=[-1.5, None],
)
geometry = chematic.parse_fhiaims(text)
assert geometry["positions_angstrom"] == [[1.25, 1.75, 3.5], [1.0, 2.0, 3.0]]
```

The parser returns a dictionary with:

| Key | Meaning |
|---|---|
| `elements` | Canonical element symbols in atom order |
| `coordinates` | Stored source coordinate triples |
| `coordinate_kinds` | `cartesian` or `fractional` per atom |
| `positions_angstrom` | Derived Cartesian coordinates |
| `lattice_angstrom` | Zero to three ordered Cartesian lattice row vectors |
| `periodicity` | Number of lattice vectors; no padding or cell inference |
| `constraints` | Fixed Cartesian x/y/z components per atom |
| `lattice_constraints` | Fixed x/y/z components per lattice vector |
| `initial_charges` | Optional initial-density charge values in +e |
| `initial_moments` | Optional initial collinear N_up − N_down values |
| `source_text` | Exact original text sidecar |

The writer takes `elements`, `coordinates` and the six corresponding optional
keyword arrays shown in the signature. It consumes source coordinate triples,
not the derived `positions_angstrom` field. Default kinds are Cartesian;
omitted constraints are false, lattice vectors default to an empty list, and
initial density values default to absent. Supplied parallel arrays must have
matching lengths. A value of `0.0` stays explicitly present and distinct from
`None`. The writer does not infer molecular formal charge or spin multiplicity.

Round-trip a parsed dictionary with:

```python
fields = ("coordinate_kinds", "lattice_angstrom", "constraints",
          "initial_charges", "initial_moments", "lattice_constraints")
output = chematic.write_fhiaims(
    geometry["elements"], geometry["coordinates"],
    **{key: geometry[key] for key in fields},
)
```

## Rust and supported boundary

The crate root exports `FhiaimsGeometry`, `FhiaimsAtom`,
`FhiaimsCoordinates::{Cartesian, Fractional}`, `FhiaimsLatticeVector`,
`FhiaimsError`, `FhiaimsParseLimits`, `parse_fhiaims`,
`parse_fhiaims_with_limits`, and `write_fhiaims`.
`positions_angstrom()` derives Cartesian positions; `validate()` checks
mutable public data before serialization.

- `atom`, `atom_frac`, `lattice_vector`, `constrain_relaxation`,
  `initial_charge`, and `initial_moment` are supported. Unknown directives
  fail explicitly rather than disappearing from structured output.
- Cartesian geometry supports zero to three independent lattice vectors.
  Fractional coordinates require exactly three. The basis is preserved,
  including handedness and order. Atom coordinates are neither wrapped nor
  rotated. Mixed coordinate kinds and lattice declarations after atoms work.
- Constraints attach to the most recently declared atom or lattice vector.
  Initial charge/moment must follow an atom, with no intervening lattice
  declaration. Repeated initial-value directives are refused.
- Constraint values are `.true.`, `.false.`, `x`, `y`, and `z`, one per line.
  Multiple active axes accumulate. Conflicting `.false.` after an active
  constraint is deliberately refused. Canonical writing emits `.true.` for
  three fixed axes, individual axis lines otherwise, and omits all-false
  defaults. Their structured meaning is retained.
- Custom species names, ghost/empty sites, isotope aliases, noncollinear
  moments, occupation blocks, velocities, external files, symmetry settings,
  Hessians and relaxation restart controls are outside this API. Canonical
  element names are required; no external `control.in` mapping is guessed.
- This is geometry interchange. Initial charge and collinear moment remain
  starting-density metadata. The code does not validate total system charge,
  electronic occupations or `control.in`, and does not run FHI-aims.
- All numbers and derived Cartesian positions must be finite. Lattice vector
  lengths must be finite and positive. Normalized two-vector cross-product
  length or three-vector determinant magnitude must exceed `1e-6`.
- Defaults are 64 MiB input, 16,384 bytes per line, and 1,000,000 atoms.
  Python exposes smaller read limits through keyword arguments. The atom
  maximum is a hard cap. Canonical output is limited to 64 MiB.
- LF/CRLF, `#` comments and Fortran D/d exponents are accepted. Canonical
  output uses LF and scientific notation with 17 decimal places. Original
  text is retained as a sidecar; writing uses current structured values.
- Rust errors include a one-based source line for rejected fields. Whole
  geometry consistency errors use zero. Python maps these to `ValueError`.
- The delivery exposes Rust and Python APIs; WASM, CLI and the common
  conversion bridge have not been extended.

## Validation

`scripts/validation/validate_fhiaims.py` generates 192 original synthetic documents from
eight small geometry recipes, six dimensionality/frame choices, metadata
on/off and LF/CRLF. Every structured field must match the original input
arrays, and every document must retain its fields through an exact internal
round trip. Fractional Cartesian positions are checked independently with
NumPy matrix multiplication.

| Local check, 2026-10-15 | Result |
|---|---|
| Exact structured round trips | 192 passed |
| ASE 3.27.0 comparisons | 128 passed |
| Open Babel 3.2.1 Cartesian comparisons | 64 passed |
| Fractional or mixed-frame cases | 64 passed |
| Maximum Open Babel CJSON lattice error | 8.882e-16 Å |
| Maximum Open Babel rounded output error | 4.762e-6 Å |
| Unexpected mismatches | 0 |
| Rust library + FHI-aims + EXYZ tests | 753 + 8 + 9 passed |
| Python FHI-aims + EXYZ + TURBOMOLE + graph-query tests | 138 passed |

Common Cartesian coordinates and cells use a `2e-12` Å comparison tolerance.
Open Babel's five-decimal output uses `5.00001e-6` Å. Native structured fields
are checked exactly. Rust separately covers all eight constraint masks,
partial dimensionality, invalid/near-dependent cells, fractional overflow,
source association, explicit zero/absence, limits and mutable-data edits.

The ASE comparison covers homogeneous input or a Cartesian projection of
mixed-frame input; lattice constraint directives are omitted from its input.
It verifies positions, element order, atom constraints, initial values and
three-dimensional cells. Partial periodicity and lattice constraints have
exact internal contracts, without an external compatibility claim.

Measured Open Babel differences are retained separately in the report:

- The tested `atom_frac 0.5 0.25 -0.5 Si` in the skewed example cell is returned
  as `[0.5, 0.25, -0.5]` Cartesian. The official formula, NumPy and ASE give
  `[1.25, 1.75, 3.5]`. Fractional input is outside the Open Babel parity check.
- The tested one-dimensional lattice is absent from Open Babel CJSON output.
- Open Babel's writer omits atom/lattice constraints and initial charge/moment
  lines. Across the common-geometry cases it drops 128 explicit initial-value
  fields. chematic retains those fields and their absence/presence.

All fixture inputs are original project data under `MIT OR Apache-2.0`.
No external molecular corpus or restricted FHI-aims source is redistributed.
The generated report records provenance, minimization, oracle versions and
settings, expectations and normalized LF UTF-8 hashes. The static contract
fixture is in `validation/openbabel_issues/fhiaims/contract-fixtures.json`.

```sh
python scripts/validation/validate_fhiaims.py --obabel /path/to/obabel
cargo test -p chematic-mol --lib --test fhiaims_contract --test exyz_contract
```

The Python gate requires the built candidate extension, ASE 3.27.0 and NumPy
2.5.4. These are local source-candidate results, not package publication,
calculation validation, performance evidence or general Open Babel parity.
The isolated candidate inherits previous EXYZ and RDKit graph-query work;
its portable patch and evidence are delivered without committing or pushing.
