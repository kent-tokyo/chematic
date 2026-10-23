# Checked DMol3 .incoor geometry

The isolated candidate adds standalone DMol3 `.incoor` reading/writing to
Rust and Python. Open Babel exposes this coordinate grammar through its
[`dmol` / `outmol` format](https://openbabel.org/docs/FileFormats/DMol3_coordinates_format.html).
The new API retains raw Cartesian Bohr values, atom order and an optional
oriented three-vector cell. Angstrom coordinates are a separate projection.
It does not parse complete `.outmol` calculation logs, `.car` or `.arc` files.

The implementation is independently written from the public
[ASE format description](https://docs.ase-lib.org/_modules/ase/io/dmol.html)
and observed format examples generated from original synthetic inputs.
No Open Babel or ASE implementation code or upstream corpus was copied.

## API

```python
import chematic

text = chematic.write_dmol3(
    ["Si", "Cl"],
    [[0.125, -1.25, 2.5], [-1.5, 0.125, 0.5]],  # Bohr
    cell_bohr=[[4, 1, 0], [-1, 5, 2], [1, 0, -6]],
)
g = chematic.parse_dmol3(text, max_atoms=1000)
assert g["periodicity"] == 3
assert g["positions_bohr"][0] == [0.125, -1.25, 2.5]
assert g["positions_angstrom"][0][0] == 0.125 * 0.529177210544
```

The dictionary contains `elements`, `positions_bohr`, `cell_bohr`,
`positions_angstrom`, `cell_angstrom`, `bohr_angstrom`, `periodicity` and
`source_text`. `cell_bohr=None` and `cell_angstrom=None` mean nonperiodic;
otherwise they contain exactly three Cartesian row vectors. The writer accepts
Bohr values explicitly. Element and coordinate counts must agree.

Rust exports `Dmol3Atom`, `Dmol3Geometry`, `Dmol3Error`, `Dmol3ParseLimits`,
`DMOL3_BOHR_ANGSTROM`, `parse_dmol3`, `parse_dmol3_with_limits` and
`write_dmol3`. `Dmol3Geometry::positions_angstrom()` and `cell_angstrom()`
validate current fields before projecting them. The core is pure Rust;
there is no external chemistry runtime dependency. No WASM or common-format
converter adapter is added in this increment.

## Grammar and preservation

```text
$cell vectors
4 1 0
-1 5 2
1 0 -6
$coordinates
Si .125 -1.25 2.5
Cl -1.5 .125 .5
$end
```

The optional `$cell vectors` block precedes `$coordinates` and has exactly
three rows. The coordinate block contains at least one canonical element
symbol followed by exactly three finite Bohr values. `$end` is required.
Blank lines and leading/trailing whitespace are accepted; header spelling
and internal spacing are exact. Both LF and CRLF input are supported.

Raw binary64 Bohr values are authoritative. Writing uses current validated
fields with 17 fractional digits in scientific notation, canonical element
symbols, LF and a final newline. It preserves negative zero, subnormal and
largest-finite coordinates through native read/write. The exact original
input remains in `source_text`; canonical writing never replays stale text.

Cells retain all components, order, orientation and handedness. Positions
are not wrapped, rotated or transformed into a standardized cell. Validation
requires positive finite vector norms, normalized absolute determinant above
`1e-6`, and finite positive volume. Dependent or nearly dependent vectors and
overflowing/underflowing cell volumes are refused.

## Units and measured engine differences

The native Angstrom projection multiplies Bohr by **0.529177210544**, the
central value from [CODATA 2022](https://physics.nist.gov/cuu/pdf/wall_2022.pdf).
This conversion does not alter raw storage. Normal floating-point rounding
still applies to the derived Angstrom values.

The pinned engines use different constants:

| Operation | Pinned factor |
|---|---|
| chematic Bohr → Angstrom | `0.529177210544` |
| Open Babel 3.2.1 Bohr → Angstrom | `0.529177249` |
| Open Babel 3.2.1 Angstrom → Bohr | `1.889725989` |
| ASE 3.27.0 Bohr → Angstrom | `0.5291772105638411` |

The Open Babel read/write factors are not exact reciprocals. Its raw-Bohr
roundtrip therefore scales values by their product, approximately
`1.0000000002228242`, before 14-decimal output rounding. This is reported
as a known difference; the oracle checks the explicitly pinned factors
instead of silently changing chematic's units to match that drift.

## Explicit boundaries

The reader refuses output-log preambles, repeated jobs, unknown blocks,
comments, trailing nonblank data and missing terminators. It also refuses
extra atom fields, isotope/atom-label spellings, atomic-number identifiers,
noncanonical symbols, non-finite values, partial periodicity and invalid cells.
It performs no bond perception, charge/spin inference, optimizer execution or
validation of DMol3 element parameterization.

Fortran `D/d` exponents are accepted natively. They are outside the common
Open Babel/ASE lane: in a separate three-dimensional diagnostic, chematic
reads `1.25D-1` as `0.125` Bohr; Open Babel reads it as `1.25` Bohr before its
Angstrom conversion; ASE rejects it with `ValueError`.

ASE's `.incoor` interface is restricted to 3D periodic systems. On the isolated
nonperiodic input, ASE 3.27.0 reading raises `UnboundLocalError`; its writer
explicitly rejects nonperiodic `Atoms`. These are recorded diagnostics and
excluded from the ASE common comparison. Nonperiodic documents are checked
natively and with Open Babel. See also the
[ASE DMol3 interface](https://docs.ase-lib.org/ase/calculators/dmol.html).

Hard input/output bounds are 64 MiB, 1024 bytes per physical input line and
one million atoms. Rust limits and Python `max_input_bytes`, `max_line_len`
and `max_atoms` can tighten these bounds; larger overrides do not raise them.
Errors identify source lines where local, or line zero for whole-geometry
validation. Python format/validation failures raise `ValueError`.

## Local validation, 2026-10-15

Oracles: **Open Babel 3.2.1** (`Jul 11 2026 -- 19:21:27`) and **ASE 3.27.0**.
The corpus contains original synthetic recipes with sixteen seeds,
nonperiodic/3D cells, `E`/`D` exponents and LF/CRLF. Inputs are authored
independently of the production writer. Every manifest record carries
provenance, license, redistribution, minimization, oracle settings, expected
values and a normalized UTF-8/LF SHA-256.

| Check | Result |
|---|---|
| Generated exact native roundtrips | 128 / 128 |
| Original analytic manifest fixture | 1 / 1 |
| Open Babel read / write, E-exponent 0D/3D lane | 64 / 64 each |
| ASE read / write, E-exponent 3D lane | 32 / 32 each |
| Unexpected mismatches within declared normalization/tolerance | 0 |
| Max Open Babel read position / cell residual | 0 / `8.881784197001252e-16` Angstrom |
| Max ASE read position / cell residual | 0 / 0 Angstrom |
| Max Open Babel write residual after its pinned scale | `5.329070518200751e-15` Bohr |
| Max ASE raw write residual | `4.884981308350689e-15` Bohr |
| Max Open Babel raw roundtrip drift | `1.336950106178847e-9` Bohr |
| Max Open Babel/native Angstrom unit difference, positions | `1.2086171419234404e-7` Angstrom |

Read comparison uses `obabel -idmol INPUT -ocjson -ab`; write comparison
uses `obabel -idmol INPUT -odmol -ab`. Bond inference is disabled. ASE uses
`ase.io.read/write(format="dmol-incoor")` with default `bohr=True`.
Declared tolerances are `2e-12` Angstrom for reads against each engine's own
factor and `2e-12` Bohr for writes, after the documented Open Babel scale.
The native structured roundtrip is exact rather than tolerance-based.

Rust validation passed 753 unit tests plus 22 integration tests (775 total).
The new Rust contracts also roundtrip all 118 canonical elements. Python
130 tests passed, including MOPAC, FHI-aims and RDKit path regressions.
Scoped Clippy with warnings denied, formatting, the unchanged refactor
allowance budget, binding inventory, product naming and strict MkDocs are
recorded in the evidence. These local results do not establish complete
DMol3 support, scientific validity, hosted CI or a published release.

## Reproduction and delivery

```sh
cargo test -p chematic-mol --lib --test dmol3_contract \
  --test mopac_cartesian_contract --test fhiaims_contract
python -m pytest crates/chematic-py/tests/test_dmol3.py
python scripts/validation/validate_dmol3.py --obabel /opt/homebrew/bin/obabel
```

Build/install the candidate extension before running Python. The oracle
script checks Open Babel 3.2.1 and ASE 3.27.0 explicitly. Fixtures are in
`validation/openbabel_issues/dmol3/contract-fixtures.json`; evidence is in
`validation/results/openbabel-dmol3-2026-10-15.json`.

The source candidate is `/private/tmp/chematic-dmol3`, based on
`/private/tmp/chematic-paths`. The incremental patch is
`validation/openbabel_issues/dmol3-2026-10-15.patch`. Only this document,
patch and evidence are copied to the primary checkout, preserving unrelated
source work. The feature is uncommitted, unpushed and unpublished.
