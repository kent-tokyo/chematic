# Checked Maestro output

Base `540086f978612df1e21baa731c0884dea7441f46`.
[Open Babel #2774](https://github.com/openbabel/openbabel/issues/2774) remains
OPEN and asks about SDF-to-Maestro availability in a Windows installation.
No chematic Maestro writer existed. This adds Python `write_mae_checked`,
`v2000_to_mae_checked` and frozen `MaeExport`, over the existing checked native
CDJSON graph projection. No Rust source or bindings change.

## Validation

- **239 Python tests pass**: 77 new Maestro cases plus 162 compatibility cases.
  Covers explicit H/ions, numeric and size bounds, invalid/profile chemistry,
  quotes/backslashes, frozen result/source immutability, 3D V2000 extraction,
  opaque SD metadata, multiple-record rejection, raw query/stereo/isotope/radical
  flags, duplicate bonds, charge-reference/overflow/duplicate errors and numeric
  underflow/nonstandard tokens. Strict mypy accepts the helper/public consumer;
  five wrong calls/mutations produce five expected errors.
- **1404 V2000 MOL/SD conversions and 1404 direct graph exports** independently
  read by RDKit **2026.09.1** (its Maestro/maeparser consumer), with 39 authored
  ordinary nonaromatic explicit-H graphs/ions, four coordinate scales, three
  translations and three title styles. No corpus case is skipped/refused.
- Another **five coordinate-edge controls** and **two hypervalent direct graphs**
  pass. **18959 atoms / 16216 bonds** are checked for exact indexed elements,
  formal charges, zero isotope/implicit-H counts, bond endpoints/order, title
  and XYZ values. Observed coordinate difference is **zero**. Does not compare
  stereolabels, source file metadata, physical validity or entire Maestro schema.
- Source stub, binding, cross-binding, unsafe, refactor and product-name gates;
  strict MkDocs, executed guide and three API anchors. Source syntax also parses
  under Python 3.9 grammar; runtime/typechecking use Python 3.13.6 / NumPy 2.5.3.

`water.sdf` / `water.mae` are original generated reference/control samples.
`oracle.py` independently generates its SD inputs via RDKit rather than using
chematic's own SD writer. The retained warning about not removing an isolated
hydrogen comes from RDKit's proton source construction; counts are still checked.

## Output contract

One ASCII m2io 2.0.0 file / f_m_ct structure, indexed atom/bond rows with atomic
numbers, formal charges, original order and caller XYZ in angstroms, using .17g
numbers and quoted/escaped titles. Requires `allow_metadata_loss=True`: source
records/SD properties, residue information, labels and force-field types are
explicitly omitted. It does not add H, infer bonds/coordinates or write files.
Returned source text is retained for the V2000 adapter. Output and report are
independent strings/frozen metadata. Title defaults to `mol.GetName()`.

Profile: fully explicit H, no isotopes, aromatic/query/dummy/mapped atoms,
declared stereo or enhanced groups; plain single/double/triple bonds. A small
explicit ordinary-valence table covers H/C/N/O/halogens/P/S and selected ions;
all other charge/valence/element combinations refuse. The API only validates
this valence profile, not electronic states or physical geometry. Readers can
infer stereo or 2D/3D status from coordinates; label/dimensional-flag parity is
not promised. No automatic kekulization or hydrogen completion.

Direct bounds: 1..4096 atoms, <=8192 bonds, <=2 MiB output, title <=256 printable
ASCII characters, finite coordinate magnitude <=1e8 and nonzero >=1e-100.
V2000 adapter: single MOL/SD record, <=999 atoms/bonds, <=16 MiB ASCII input,
<=20000 lines, <=4096 bytes/line, optional checked unique M CHG references with
charges -3..3 and at most 8 entries per row. Atom flags are zero except ordinary
charge codes / valence 15 (no implicit H); unsupported CTAB records refuse.
Raw coordinate decimal grammar and nonzero underflow are checked. SD fields
are opaque omitted data, not a complete SDF schema validation. V3000, batch
streams, isotope/stereo export, compression, reader, Rust/WASM/CLI surfaces,
Maestro GUI/force-field execution and speed claims are outside this addition.

## Open Babel control

The local **Open Babel 3.2.1 -- Jul 11 2026 -- 19:21:27** format registry lacks
`mae`; the `-isdf water.sdf -omae` control returns 1, reports `cannot write output
format`, and prints usage to stdout. `openbabel-stdout.log` is that diagnostic,
not an exported Maestro file. Registry/help, version/command/status and stderr
are retained. This is an installation-specific missing-format observation,
not proof of absent Maestro support in all upstream builds or an installer fix.

## Native fixture and reproduction

The existing **ca7e14e8** native extension is used without modification as an
import/projection/parser fixture. Its SHA-256 matches that task's recorded
provenance. Selected CDJSON module/binding paths are unchanged at this base;
other core/parser dependencies have subsequent changes. This is **not** a fresh
current-main native build, published-wheel or cross-platform validation.
Source/native hashes and pinned upstream source revision: `provenance.json`.
Only the independent reader uses RDKit; production code requires no RDKit.

Set `CANDIDATE_PY` to this source package plus the verified existing extension,
and `RDKIT_PY` to RDKit 2026.09.1. Python 3.13.6 / pytest 9.1.1.

```sh
PYTHONPATH="$CANDIDATE_PY" python3.13 -m pytest \
  crates/chematic-py/tests/test_maestro_output.py \
  crates/chematic-py/tests/test_rdkit_compat.py -q -p no:cacheprovider
PYTHONPATH="$CANDIDATE_PY:$RDKIT_PY" python3.13 \
  validation/results/maestro-output-2026-10-18/oracle.py
MYPYPATH="$CANDIDATE_PY" python3.13 -m mypy --no-incremental \
  --config-file validation/results/maestro-output-2026-10-18/mypy.ini \
  crates/chematic-py/python/chematic/_maestro.py \
  validation/results/maestro-output-2026-10-18/typing_valid.py
# typing_invalid.py: expect five errors, exit 1.
python3.13 -m mkdocs build --strict --site-dir /tmp/chematic-mae-site
```

Schema sources: [maeparser README](https://github.com/schrodinger/maeparser/blob/b8106cc17c4857fb2b0072159bf42d9df2d8ef34/README.md),
[block/property constants](https://github.com/schrodinger/maeparser/blob/b8106cc17c4857fb2b0072159bf42d9df2d8ef34/MaeConstants.hpp).
Original working-tree document changes are preserved.
