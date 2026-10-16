# Explicit single/double bond encoding in MOL2

[Open Babel #2668](https://github.com/openbabel/openbabel/issues/2668) requests
Kekulé bond output when converting SDF to MOL2. The
[documented MOL2 write options](https://open-babel.readthedocs.io/en/latest/FileFormats/Sybyl_Mol2_format.html)
do not include that choice. chematic provides an explicit, bounded MOL2 record
operation in Rust, Python and WASM:

```python
import chematic

with open("aromatic.mol2", encoding="utf-8", newline="") as file:
    source = file.read()
output = chematic.kekulize_mol2_record(source)
with open("explicit-bonds.mol2", "w", encoding="utf-8", newline="") as file:
    file.write(output)
```

The JavaScript/Node API is `wasm.kekulize_mol2_record(text)`. Rust exposes
`chematic_mol::kekulize_mol2_record(&str) -> Result<String, Mol2KekuleError>`.
For SDF input, first produce a MOL2 record using an existing conversion
workflow, then apply this operation. This API does not add a new SDF converter.

## Exact source contract

Only the fourth token of an aromatic row in `@<TRIPOS>BOND` changes:
`ar` becomes `1` or `2`. A pyridine ring obtains three single bonds and three
double bonds. Original IDs, row order, endpoint order, atom names, coordinates,
types, residue fields, charge-method labels, partial-charge columns, formal
charge attributes, status bits, comments, unknown sections, whitespace and
line endings remain byte-for-byte identical. There is no coordinate rounding.

Atom types such as `C.ar` and `N.ar` deliberately remain. This is a bond
encoding choice; it does not retype atoms or remove aromatic chemical identity.
Unknown metadata is retained as uninterpreted source, not chemically validated.
Non-aromatic bond tokens, including `am`, are unchanged. Applying the operation
again returns exactly the same text. Records without `ar` bonds are bounded,
validated exact no-ops. Existing writers and exact source round-trips retain
their original defaults.

The core matching algorithm chooses an assignment using the supplied atom
order. Different row orders may select different valid resonance forms.
No canonical Kekulé form or numbering-independent byte output is promised.

## Conservative accepted profile

One record is limited to **8 MiB, 512 atoms and 1,024 bonds**. Conversion accepts
neutral aromatic **C/N/O/S** atoms with degree two or three, ordinary covalent
graph bonds, and aromatic bonds belonging to cycles of aromatic bonds.
Each converted atom must retain its implicit hydrogen count and stay within
the profile's neutral valence (C: 4, N: 3, O/S: 2). Every converted bond must
be recognized as aromatic again by chematic's Hückel model after the explicit
assignment. A complete matching alone is insufficient: antiaromatic and
saturated donor cycles are refused.

Pyrrole-like hydrogens must be explicit atom nodes; an omitted N-H is not
guessed from a type string. N-substituted pyrroles with explicit substituents
are supported. Charged aromatic atoms, other aromatic elements, ambiguous
query/zero bonds, malformed records, unsuccessful matching and aromaticity or
hydrogen-count changes return errors. Python raises `ValueError`; WASM throws.
There is no partial output. Charged atoms outside the aromatic system may be
retained unchanged, as demonstrated by the benzene/Na/Cl fixture.

This bounded profile is not a general aromaticity validator, aromaticity-model
parity claim, charge predictor, atom-typing engine or coordinate generator.
Source coordinates and metadata remain authoritative for their own meaning.

## Validation and reproduction

The dated packet records **204 matching Rust/Python/WASM outputs** and exact
token-preservation proofs, **204 RDKit graph/explicit-valence checks**, and
**17 Open Babel before/after identity checks**, with zero mismatches. These
counts apply to the listed synthetic controls and reference settings.

Original fixtures cover pyridine, benzene, fused/biphenyl rings, explicit-H and
N-substituted pyrroles, furan, thiophene, substituents, an amide, an ionic spectator,
implicit carbon hydrogens, an unchanged non-aromatic control, and configured
tetrahedral/geometric sidechains. Fixture provenance is recorded in
`validation/openbabel_issues/mol2-kekule/cases.json`; no upstream attachment or
implementation was copied.

The independent runner permutes atom and bond rows with three fixed seeds and
retains an identity order. LF, CRLF and absent final newline variants exercise
source preservation. It checks actual Rust/Python/Node-WASM output equality,
exact token changes, idempotence, and independent RDKit molecular identity and
explicit per-atom valence. RDKit uses an independently normalized graph
projection with sequential IDs and explicit UNITY formal charges; its reader
does not accept the sparse-ID, comment/status-rich source profile directly.
This is a graph check, not raw-record RDKit interoperability. Open Babel
independently cross-reads the base fixtures with LF line endings; its stdin
reader rejects CRLF on this host. Sparse IDs and comments inside bond rows are
covered by the separate native/binding source-preservation control, not that
external-reader lane. The Open Babel lane reparses its valid output SMILES
before comparing identity, because its MOL2 reader retains explicit bond
encoding without necessarily flagging it aromatic. Its default versus `-xk`
pyridine output is recorded separately
as an upstream diagnostic, not used as the transformation oracle.

```sh
cargo build -p chematic-mol --example mol2_kekule_dump
python scripts/validate_mol2_kekule.py \
  --native-dump target/debug/examples/mol2_kekule_dump \
  --wasm-module /path/to/node-package/chematic_wasm.js
```

The environment needs the matching candidate Python extension, RDKit 2026.09.1,
Node and Open Babel. The runner writes compressed raw comparisons and source,
fixture, packet and runtime hashes under
`validation/results/mol2-kekule-2026-10-16/`. Local validation and hosted CI
status are recorded separately; neither constitutes a package release.

Reference interpretation is bounded: Open Babel reads the fixture's sodium
spectator as neutral even before conversion, while RDKit treats hydrogen-omitted
MOL2 carbons as radicals. Their respective before/after interpretations remain
the same. These controls establish preservation, not agreement on absolute
formula or charge across engines. Exact UNITY bytes and native implicit-H
counts are checked independently; use explicit hydrogens for portable records.
