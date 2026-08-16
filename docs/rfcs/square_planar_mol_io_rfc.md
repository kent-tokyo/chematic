# Square-planar stereo in MOL/SDF (V2000 + V3000)

## 1. The question

`Chirality::SquarePlanar(SquarePlanarPermutation)` (PR #313/#326) already round-trips
through SMILES's `@SP1`/`@SP2`/`@SP3` tags. MOL/SDF (V2000 and V3000) has no
such symbolic tag at all. Before writing any code, this RFC answers: is there
*any* way to preserve a declared square-planar configuration through a MOL/SDF
round trip, and if so, what is it?

## 2. What the MDL/CTfile format actually has

Neither V2000 nor V3000 has a field that stores a permutation value
(1/2/3, or any other discrete SP-class marker) for a stereocenter:

- **V2000's atom-line "stereo parity" field** exists in the spec, but this
  codebase never reads it for *any* geometry (see `mol2000.rs`'s atom-block
  parser — the field is skipped entirely), matching RDKit's own primary
  `MolFromMolBlock` path. It has no permutation-value semantics regardless
  (0/1/2/3 = none/odd/even/either, a 2-state-plus-unknown tetrahedral parity
  flag, not a 3-way class selector).
- **V3000's atom-line `CFG`** is likewise unread today (explicitly out of
  scope, `mol3000.rs`'s own comment). Its defined values (0/1/2/3) are again
  a tetrahedral-shaped parity/either encoding, not a permutation-class field.
- **V3000's `COLLECTION`/`STEABS`/`STEOR`/`STEAND` block** groups atoms into
  absolute/OR/AND *certainty* classes for already-defined tetrahedral
  centers. It carries no per-atom geometry-class or permutation information
  and cannot be repurposed for this without inventing semantics MDL never
  defined.

**Conclusion: there is no lossless, standard, symbolic MOL/CTfile field for
`@SPn`.** This is Tier 1 in the design space below, and it is a dead end —
confirmed independently against RDKit 2026.03.4 (see §4) and consistent with
public RDKit documentation/source.

## 3. What DOES exist: 3D-coordinate-derived reperception

Square-planar geometry, unlike tetrahedral parity, is a real physical shape a
3D conformer can encode directly: 4 ligands and a center, coplanar, with two
ligand pairs at ~180° (trans) and the rest at ~90° (cis). A MOL/SDF atom
block already carries real `(x, y, z)` coordinates. If those coordinates
geometrically encode a square-planar arrangement, the SP1/SP2/SP3 class can
be **reperceived from geometry**, the same way this crate's existing 2D
wedge/hash pathway reperceives tetrahedral parity from a drawing rather than
storing it symbolically.

## 4. RDKit oracle findings (2026.03.4, empirical, this session)

RDKit was available locally (`.venv`, `rdkit==2026.03.4`) and used as a live
oracle, matching this repo's established methodology (PR #313's own square
planar RFC used the same approach).

**Finding 1 — no symbolic round trip.** A molecule built with
`ChiralType.CHI_SQUAREPLANAR` + `_chiralPermutation`, written via
`Chem.MolToMolBlock` with **no existing 3D conformer**, gets *auto-generated
2D coordinates* that geometrically encode the correct trans-pairing (verified
by hand: the emitted `(x, y)` layout puts SP1's declared trans-pairs at
antipodal positions around the center). But reading that exact block back
with plain `Chem.MolFromMolBlock` (default arguments) returns
`CHI_UNSPECIFIED` — the tag is **not** recovered automatically, from either
V2000 or V3000, from a 2D-declared or a flat-z record.

**Finding 2 — 3D reperception is real.** Calling
`Chem.AssignStereochemistryFrom3D(mol)` explicitly, on a parsed molecule
whose real 3D coordinates encode a square-planar shape, *does* recover
`CHI_SQUAREPLANAR` with the correct `_chiralPermutation`, for **both** V2000
(`forceV3000=False`, confirmed with a hand-built RWMol using plain covalent
bonds, no dative bonds involved) and V3000. Round-tripping that block again
through `MolToMolBlock` reproduces the identical coordinates. This holds even
when every z is exactly `0.0` (a genuinely planar 3D structure is
legitimately flat) — RDKit's own `MolFromMolBlock`/`AssignStereochemistryFrom3D`
split treats "is this really 3D data" as a question of *which function the
caller invokes*, not a geometry heuristic on z-values.

**Finding 3 — element gating, not pure geometry.** The exact same coplanar
`(F, Cl, Br, I)` square arrangement around a center is perceived as
`CHI_SQUAREPLANAR` for Pt/Pd/Ni/Cu/Au/Fe/Xe/K/Ca/Zn/Ga/Ge/Sn/Pb/Bi, but
**not** for C/Si/Na/Mg/Al — RDKit does not assign non-tetrahedral chirality
to an element whose periodic table entry carries a small, fixed default
valence. Pure coordinate geometry is not the only signal; element identity
gates eligibility too (§6).

**Finding 4 — angle tolerance is real but bounded.** Distorting a square
arrangement's trans-pair angle away from 180° (equivalently, the cis angle
away from 90°) by up to ~30° still perceives correctly; by 40° it fails
(`CHI_UNSPECIFIED`). RDKit's own cutoff sits strictly between those two
values.

These four findings settle the design question: **Tier 2 (3D-coordinate-derived
reperception) is real, RDKit-precedented, and implementable** — not a
hypothetical.

## 5. Design decision: Tier 2, read + write, both formats

**Read**: for every atom whose element has no fixed default valence (§6) and
exactly 4 heavy-atom neighbors, when the record has genuine (non-flat, i.e.
not all-`z≈0`) 3D coordinates, classify the local 5-point geometry (center +
4 neighbors). If it resolves unambiguously to one of SP1/SP2/SP3 (coplanar,
within tolerance, with exactly one of the three trans-pairings both
≥135°), set `Chirality::SquarePlanar` and `stereo_neighbor_order`. If the
center is coplanar-*ish* but the classification is ambiguous or degenerate,
surface a typed diagnostic (`SquarePlanarPerceptionDiagnostic`) rather than
silently leaving `Chirality::None` — the same "explain why, don't guess"
discipline `StereoDiagnostic`/`Stereo3DDiagnostic` already use elsewhere in
this crate. If the center simply isn't coplanar (the overwhelmingly common
case — most 4-coordinate transition-metal centers are tetrahedral or
octahedral, not square-planar), nothing is reported: this mirrors
`StereoDiagnostic`'s existing "no wedge present -> nothing to say" precedent.

**Write**: only permitted when a 3D conformer already exists, and the
existing coordinates are validated against the declared `@SPn` tag before
writing — never fabricated from nothing, and never silently trusted to
match. See §8 for why fabricating coordinates was considered and rejected.

**Both formats, verified not assumed**: V2000 and V3000 atom lines both
carry independent `(x, y, z)` fields at the same fixed columns, and both
formats share this crate's existing `classify_geometry_rank`/`parse_dimension_code`
(mol3000.rs already calls into mol2000.rs's implementations directly — see
that file's existing `wedge_vs_3d_conflicts` reuse). The new perception and
validation logic is added once, in `mol2000.rs`, and reused by `mol3000.rs`
the same way. This was confirmed, not assumed: RDKit finding 2 above was
independently reproduced against *both* a V2000 block (hand-built,
`forceV3000=False`, plain covalent bonds) and a V3000 block.

## 6. Element eligibility gate: reusing `Element::normal_valences()`

RDKit's oracle behavior (§4, finding 3) shows non-tetrahedral perception is
gated by element, not geometry alone: a real sp3 carbon is never square
planar even if handed adversarial/synthetic coordinates that happen to look
that way. This project already has a fitting primitive:
`Element::normal_valences() -> &'static [u8]`, documented as "Empty =
undefined (transition metals, etc.)". The new geometric classifier only
considers atoms whose element has `normal_valences().is_empty()`.

This was RDKit-tested directly (same coplanar-square coordinate fixture, per
element, `AssignStereochemistryFrom3D`) for two groups: elements this
crate's `normal_valences()` gives a defined list for and therefore *excludes*
(**C, Si** — both correctly rejected by RDKit too, `CHI_UNSPECIFIED`), and
elements it returns empty (undefined) for and therefore *includes*
(**Pt, Pd, Ni, Cu, Au, Fe, Xe, K, Ca, Zn, Ga, Ge, Sn, Pb, Bi** — all correctly
perceived by RDKit as `CHI_SQUAREPLANAR`). The remaining entries in this
crate's 15-element `normal_valences()` table (H, B, N, O, F, P, S, Cl, As,
Se, Br, Te, I) were **not** individually RDKit-tested; they are assumed, not
verified, to follow the same defined-valence-excludes-them rule by the same
logic as C/Si (all are ordinary main-group organic-subset elements with
small, fixed valences in RDKit's own periodic table too), consistent with
this project's policy against unverified "fixed"/"matches" claims.

**One confirmed, honest divergence**: RDKit's periodic table also assigns
Na/Mg/Al a small fixed valence and excludes them (RDKit-tested directly,
§4 finding 3), but this crate's `normal_valences()` table has no entry for
Na/Mg/Al and returns empty for them too, so this crate's gate is *more
permissive* than RDKit's for those 3 specific, RDKit-confirmed-divergent
elements. In practice this divergence is very unlikely to matter: it can
only be reached by real (or adversarially-constructed) coordinates in which
an Na/Mg/Al center is genuinely coplanar with 4 ligands at square angles —
chemically implausible for these elements' actual bonding, so the
coplanarity gate (§7) would almost always exclude them anyway even without
the element check. Replicating RDKit's full internal default-valence table
for every element was judged out of scope for this PR; this is recorded
here as a known,
low-risk, documented gap rather than silently overclaiming exact parity.

## 7. Geometry classification: coplanarity + trans-pair angle

Given a candidate center (`Point3`) and its 4 neighbor positions
(`[Point3; 4]`, ordered per `Molecule::neighbors()`/`stereo_neighbor_order`):

1. **Degenerate-bond-vector check first**: reject (typed reason
   `DegenerateBondVector`) if any neighbor's distance from the center is
   below `1e-3` Å — many orders of magnitude below any real bond length,
   purely a corrupt-data guard. Checked before coplanarity so a coincident
   point doesn't trivially "pass" a plane fit.
2. **Coplanarity**: reuse `mol2000.rs::classify_geometry_rank` directly
   (not reimplemented) on the 5-point set `{center, n0, n1, n2, n3}`. Only
   `Coplanar`/`FlatZero` pass; `ThreeD` is rejected (`NotCoplanar`) — this is
   the check that makes a real tetrahedral center's ~109.5° angles
   self-reject: 4 tetrahedral substituents plus their center are never
   coplanar within `COPLANAR_EPS` (`1e-2` Å), so no separate "is this
   tetrahedral-shaped" check is needed. Reuses `COPLANAR_EPS`, explicitly
   **not** `VOLUME_EPS` — `VOLUME_EPS` is a signed-tetrahedral-volume
   epsilon from a different check with different failure semantics; using it
   here was the RFC's own stated anti-pattern to avoid (matching PR #326's
   deferred-scope note).
3. **Trans-pair angle**: for each of the 3 candidate tags (SP1/SP2/SP3),
   read its `trans_pairs()` (existing, oracle-verified — not
   reimplemented) and test whether *both* pairs have
   `cos(angle) <= -1/sqrt(2)` (i.e. angle ≥ 135°) as measured from the
   center. 135° is the geometric midpoint between an ideal cis angle (90°)
   and an ideal trans angle (180°) — a self-justifying threshold requiring
   over 45° of distortion from either ideal to misclassify, in the same
   family as RDKit's own empirically-observed tolerance (finding 4: real
   cutoff between 30° and 40°) without copying RDKit's number. Exactly one
   tag matching resolves the class; zero or more than one is
   `AmbiguousTransPairing` (typed, surfaced, never guessed).

This function (`classify_square_planar_geometry`) is the single shared
implementation used by both the reader (candidate scan) and the writer
(pre-write validation) — not two independent copies.

## 8. Why the writer never fabricates coordinates

An earlier draft of this design considered having the writer synthesize a
fresh coordinate layout (mirroring RDKit's `MolToMolBlock` auto-2D-layout
behavior, §4 finding 1) whenever a `Chirality::SquarePlanar` atom had no
existing conformer, so that e.g. a SMILES-parsed (coordinate-less) cisplatin
molecule could still be written to MOL. This was rejected: a writer that
invents geometry to satisfy its own serialization is exactly the kind of
silent, unverifiable behavior this project's fail-closed policy exists to
prevent — nothing downstream could distinguish "these coordinates are real
structural data" from "these coordinates are a fabrication that happens to
decode back to the right tag". If the caller has no conformer, that is a
**Tier 1** situation (no lossless standard encoding exists) and must return
a typed error, not trigger geometry fabrication. Round-trip test fixtures in
this PR therefore construct their own explicit `Coords3D` (the same way
existing tetrahedral wedge/hash round-trip tests in
`stereo_reader_integration.rs` construct explicit 2D `coords` rather than
having the writer invent a layout) — see `docs/rfcs/square_planar_stereo_rfc.md`'s
sibling precedent.

## 9. Write-side validation, not silent trust

Because MOL/SDF's only mechanism for this is geometry, a caller could in
principle hand the writer a `Chirality::SquarePlanar(SP1)` atom next to a
conformer whose real coordinates encode SP2 — writing that out verbatim
would silently manufacture a file that reperceives as the *wrong* tag on the
next read. The writer path added by this PR
(`write_mol_with_conformer_checked` / `write_mol_v3000_with_conformer_checked`,
backed by the shared, `pub` `validate_square_planar_for_write`) runs the
exact same `classify_square_planar_geometry` classifier used on read against
every `Chirality::SquarePlanar` atom's `stereo_neighbor_order` positions
before emitting anything, and fails closed with a typed
`MolStereoWriteError` (carrying the offending atom, the declared
`chematic_core::StereoGeometry`, the target `MolFormat`, and one of:
no conformer supplied; missing `stereo_neighbor_order`; an implicit-H slot in
that order (§11); a missing conformer position; geometry that classifies to
a *different* tag than declared; or geometry that doesn't classify at all)
on any mismatch, missing data, or ambiguity — never a silent overwrite.

It also rejects a whole-molecule-flat (`GeometryRank::FlatZero`/
`Indeterminate`) conformer outright, for the same reason the reader's
`conformer` field is `None` in that case (§10) — a conformer this function
would happily validate as geometrically SP1-shaped at z≈0 exactly would be
silently unrecoverable by this crate's own reader, which is precisely the
"writer says encoded, reader says nothing" failure this design exists to
prevent.

## 10. The z=0 ambiguity, and why fixtures use a nonzero z

`GeometryRank::FlatZero` (every atom's z within `1e-4` of exactly 0) is this
crate's existing, pre-PR#326-vintage signal for "this all-zero-z record is
indistinguishable from an ordinary 2D wedge-only depiction" — the existing
`conformer: Option<Coords3D>` field is already `None` in that case (see
`mol2000.rs`'s `read_mol_with_diagnostics`), and this PR's square-planar
perception is gated on `conformer.is_some()`, i.e. it **never runs** against
a flat-z record. This is deliberate, not an oversight: it is the direct
implementation of §2's conclusion that a pure-2D, wedge-only MOL block has
no way to encode square-planar stereo at all — flat z=0 coordinates are
exactly what a 2D depiction looks like, indistinguishable from "no 3D data
was ever provided". A genuinely-planar real 3D square-planar structure whose
plane happens to coincide with z=0 exactly is an unavoidable, accepted
false-negative of this same ambiguity (documented limitation, not fixed by
this PR — matches this file's own `GeometryRank::Coplanar` doc comment,
which already states a flat real molecule "legitimately lands" in the
zero-evidence bucket).

Because of this, every round-trip fixture in this PR's test suite places its
whole molecule off the z=0 plane (e.g. z=1.5 for every atom) — this is not a
workaround for a test-only quirk, it is what *any* real embedded 3D
conformer looks like in practice (an embedder essentially never lands a
structure exactly on z=0), and it is exactly the geometry the write-side
validator (§9) now also requires.

## 11. Explicitly out of scope

- **Pure-2D, wedge-only square-planar encoding.** No such mechanism exists
  in MDL/CTfile (§2); not invented here.
- **3-heavy-neighbor + implicit-H square-planar centers.** SMILES syntax
  permits `@SPn` with `STEREO_H_SENTINEL` in `stereo_neighbor_order` (no
  geometry-specific restriction in the parser), but there is no real spatial
  position for an implicit hydrogen to validate or write against, and this
  shape is chemically unusual for the d8 transition-metal chemistry
  square-planar stereo describes. The writer returns a typed
  `MolStereoWriteError` (reason: implicit-H slot present) rather than
  attempting anything; the reader's candidate scan requires exactly 4
  explicit heavy neighbors and never considers implicit-H centers.
- **Full RDKit-parity element eligibility table** (§6) — Na/Mg/Al divergence
  accepted and documented, not closed.
- **A chematic-specific lossless SDF property-block extension** (Tier 3 —
  e.g. a `<CHEMATIC_NONTETRAHEDRAL_STEREO>` SD field). Not built: Tier 2
  already gives genuine, real-3D-data round-tripping for the realistic case
  (a molecule that actually has a conformer), and this PR's scope is
  intentionally kept to the light-gate policy in effect for this work (no
  new SD-property-block format to design, document, and test). Left for
  future work if a pure-2D or coordinate-free use case is ever prioritized.
- **Retrofitting every pre-existing 2D-only MOL/SDF writer
  (`write_mol`/`write_mol_with_coords`/`write_mol_v3000`/the `write_sdf*`
  family) to a fallible signature.** These functions have never consulted
  `Atom.chirality` for *any* geometry (tetrahedral round-trips via the
  wedge/hash `BondOrder`, which they do write faithfully; only
  square-planar has no such always-written channel). Changing their return
  type ripples through `chematic-py`, `chematic-wasm`, and every existing
  example/test in this crate for a class of molecule vanishingly rare in
  those call sites today. Instead, this PR adds new, additive, fully
  fail-closed entry points (`*_checked`) for the one path that can actually
  represent this stereo class (the conformer-aware writers), and leaves the
  pre-existing infallible 2D writers' documented behavior as "does not
  encode non-tetrahedral stereo, use the checked conformer writer instead" —
  an explicit scope decision, flagged here for the maintainer to overrule if
  a broader breaking change is actually wanted.
- **`write_sdf_record_with_conformer` has no `_checked` sibling.** Unlike
  the 2D-only writers above, this one *does* carry a real `Coords3D`
  conformer, so the gap here is narrower: a caller can still get fail-closed
  behavior for it today via `validate_square_planar_for_write(mol,
  Some(conformer), MolFormat::V2000)` before calling it (the function is
  `pub` precisely so this composes), just not through a single wrapped call.
  Not added because no fixture or test in this PR needed the SDF-record
  (`$$$$`-terminated) framing specifically — every round-trip test operates
  at the plain MOL-block level, which the `*_checked` writers already
  cover.

## 12. `wedge_vs_3d_conflicts` fix

`wedge_vs_3d_conflicts` (`mol2000.rs`) gated on `atom.chirality ==
Chirality::None`, i.e. it ran its tetrahedral-signed-volume comparison
against *any* non-`None` chirality. Before this PR, `Chirality::SquarePlanar`
was never produced by a MOL/SDF reader, so the bug was latent — this PR makes
it reachable. Fixed by switching the gate to `!atom.chirality.is_tetrahedral()`,
using the existing allowlist-shaped `Chirality::is_tetrahedral()` (`atom.rs`,
PR #326) rather than enumerating known non-tetrahedral variants by name. This
is the same "positive match on the known-good case, not negative match on the
known-bad case" fix shape as the two prior instances of this exact bug class
in this project's history (`chematic-3d/src/stereo_constraints.rs` and
`chematic-chem/src/cip.rs`, per project history) — the fix is
exhaustive-match-safe by construction: any future non-tetrahedral geometry
(trigonal-bipyramidal, octahedral, sketched but unimplemented in
`stereo_geometry.rs`) is automatically excluded by this same gate without
requiring a new match arm here, because `is_tetrahedral()` only ever returns
`true` for the two known-tetrahedral variants.

## 13. Sibling-gate audit (not fixed, and why)

A grep for `Chirality::None ==`-shaped gates elsewhere in `chematic-mol`/
`chematic-perception` found two more production sites
(`chematic-perception/src/stereo_validation.rs`, `validate_stereo` and
`stereo_centers`). Both were checked for the same false-positive shape
`wedge_vs_3d_conflicts` had (a downstream computation that assumes
tetrahedral geometry and produces a *wrong* answer for a non-`None`,
non-tetrahedral chirality). Neither does: both are generic
neighbor-count/Morgan-rank distinctness checks that do not depend on which
concrete tetrahedral-vs-square-planar shape is present, so a
`Chirality::SquarePlanar` atom passing through them does not produce an
incorrect `StereoError`/`stereo_centers` entry. `stereo_centers`'s doc
comment labels its output "tetrahedral candidates", which is a
pre-existing, orthogonal naming/scope-accuracy gap (a square-planar center
with 4 distinct neighbors would be counted there too) rather than a
correctness bug, and duplicating this PR's new element-eligibility gate into
that module to fix the label would be scope creep beyond the one bug PR #326
named as deferred. Recorded here rather than silently left unmentioned.

## 14. Verification

- RDKit 2026.03.4 oracle (`.venv`), used interactively this session — not a
  pre-existing script, since this is a narrow point-confirmation, not a
  corpus measurement (matches this PR's light-gate scope).
- Cisplatin/transplatin MOL/SDF round-trip differential test against the
  existing, independently oracle-verified SMILES fixtures in
  `chematic-smiles/tests/square_planar_stereo.rs` (source of truth for which
  tag means cis vs trans for this specific structure) — not re-derived from
  scratch.
- Atom-renumbering invariance, coordinate rounding, and near-degenerate
  (distorted-but-still-valid) geometry — see `chematic-mol/tests/`.
- `cargo test -p chematic-mol` run immediately after wiring the reader (before
  writing new tests) to check for regressions in existing 4-coordinate
  fixtures (`platinum_benchmark.rs`, `conformer_3d_io.rs`) that might newly
  acquire a tag.
