# P1-A0: 2D-Stereo Reader Integration Boundary — Diagnosis

**Status:** diagnosis only. No production code changed. Not merged.

**Scope:** this document does not implement a wedge→SMILES stereo converter. It
answers a narrower, prior question: *where should chematic call stereo
perception at all*, given that today none of its coordinate-bearing readers
do. Implementation (P1-S1a/S1b/S2) is deferred to follow-up PRs once this
question is settled.

## 1. Why this exists

IO-3 (MRV, PR #128) RDKit-oracle validation found that 69/206 fixtures lose
tetrahedral/E-Z stereo when converted from a 2D format to SMILES
(`validation/mrv_io_parity_summary.json`,
`known_divergence_breakdown.phase1.tetrahedral_or_ez_stereo_lost_converting_wedge_bonds_to_smiles`).
The original hypothesis (see project memory
`chematic_smiles_bracket_h_and_wedge_chirality_gaps.md`) was that
`chematic_perception::stereo2d::apply_stereo_from_2d` writes the wrong field
(`Atom.cip_code` instead of `Atom.chirality`). Source-level research for this
diagnosis found that framing is incomplete: **no reader in the codebase calls
`assign_stereo_from_2d`/`apply_stereo_from_2d` at all.** Fixing the field
those functions write to would not fix anything, because nothing calls them
from a reader in the first place. The real open question is an integration
one: *when and under what conditions should a coordinate-bearing reader
invoke stereo perception?*

## 2. Reader inventory

Audited: MOL V2000, SDF, MOL V3000, CDXML, CML (all on `main`), plus MRV
(pending, unmerged PR #128, read via `git show feat/io-mrv:...` without
checking out the branch). A quick pass over MOL2 (Tripos), Ketcher (KET),
MolJSON, CJSON, and CIF determined whether they carry stereo-relevant data at
all.

| Format | 2D coords | 3D coords | Wedge/dash → | E/Z | Writes `chirality`/`stereo_neighbor_order`/`cip_code`? | Calls `chematic_perception` stereo2d? |
|---|---|---|---|---|---|---|
| MOL V2000 (`mol2000.rs`) | yes, returned as separate `Vec<(f64,f64)>` | no (Z bytes never parsed) | `BondOrder::Up`/`Down` | not represented | no (zero hits) | no |
| SDF (`sdf.rs`, wraps mol2000) | same as V2000, per-record | same as V2000 | same as V2000 | not represented | no | no |
| MOL V3000 (`mol3000.rs`) | yes, separate `Vec<(f64,f64)>` | no (z token parsed position exists but never read) | **broken asymmetrically**: writer emits bond `CFG=1/6`, but the reader never parses bond-line `CFG`/`KEY=VALUE` tokens at all — a V3000 file with a real wedge/dash bond silently round-trips to plain `Single`, and the writer's own `CFG` output can't be read back by this same reader | not represented | no | no |
| CDXML (`cdxml.rs`) | yes (ChemDraw Y-down convention) | not supported by this reader | `WedgeBegin`/`WedgedHashBegin` → Up, `Hash`/`Dash`/`WedgeEnd` → Down | **yes — actively perceived**: `parse_cdxml_all` calls `chematic_perception::assign_ez_from_2d` at parse time (`cdxml.rs:159`), writing `Atom.cip_code` | `cip_code` only, via the call above | **yes** — the one existing precedent, E/Z only |
| CML (`cml.rs`) | yes (`x2`/`y2`, Y-up) | doc comment claims x3/y3/z3 parsed but not returned — **stale**, zero occurrences of x3/y3/z3 in the actual code | **not implemented at all** — module doc states outright "stereochemistry attributes are ignored"; CML's own `<bondStereo>` element (same mechanism MRV uses) is never scanned for | none | no | no |
| MRV (`mrv.rs`, **pending PR #128, unmerged**) | yes (`x2`/`y2`) | yes (`x3`/`y3`/`z3`) | `<bondStereo>` `"W"`/`"H"` (+ MDL `convention="1"/"6"`) → Up/Down | **explicitly discarded** — cis/trans values matched and dropped with comment "matching RDKit" | no | no (only calls `kekulize_inplace`, unrelated) |
| Tripos MOL2 | no | yes | no — format has no wedge notation | n/a | n/a | n/a |
| Ketcher (KET) | yes (projects 3D→2D) | yes | **yes** — bond `"stereo"` field `1`/`6` → Up/Down, same convention family as MOL V2000/MRV. Contrary to the "probably out of scope" assumption, KET *is* in scope for a future stereo-perception design | not represented | no | no |
| MolJSON | no | no | no (pure connectivity format, no coordinate fields in the schema at all) | n/a | n/a | n/a |
| CJSON | no | yes | no | n/a | n/a | n/a |
| CIF | no | yes (fractional→orthogonal) | n/a (crystallographic, no wedge concept) | n/a | n/a | n/a |

Also relevant: `chematic_core::BondOrder::Up`/`Down` (`crates/chematic-core/src/bond.rs`)
is a single representation overloaded for two unrelated stereo concepts —
MDL/CDXML/MRV/KET wedge/dash depiction, and OpenSMILES `/`/`\` E/Z
directional-bond markers. `BondOrder::smiles_token()` maps both unconditionally
to `"/"`/`"\\"`. This is load-bearing for the pathway trace in §3.

Python/WASM bindings for MOL V2000/V3000/CDXML/CML mostly **discard the
returned coordinates** at the binding layer (e.g. `chematic-py`'s `from_cdxml`,
`from_cml` map away the `coords` tuple element); only MOL V2000/V3000 and SDF
expose coordinates through to Python/WASM callers today.

## 3. Current stereo pathway — verified end-to-end (empirically run, not just read)

Traced for MOL V2000, the richest-on-`main` format, using a real fixture run
through the actual reader and writer (not a hypothetical):

```
file: wedge bond (MDL stereo flag 1) + 2D coordinates
        |
        v
mol2000::parse_mol_with_coords
        |  sets bond.order = BondOrder::Up on the wedge bond   [verified: mol2000.rs:244]
        |  atom.chirality stays Chirality::None (default)      [zero writes anywhere in mol2000.rs]
        |  atom.cip_code stays None
        |  Molecule::stereo_neighbor_order has no entry
        v
  <-- nothing calls chematic_perception::stereo2d from here -->
        |
        v
chematic_smiles::write(&mol)   (naive call, exactly what a caller gets today)
        |
        +-- emit_atom reads atom.chirality (None) -> no @/@@ printed
        |
        +-- write_chain's bond-elision match treats BondOrder::Up/Down as
        |   "never implicit" (writer.rs: `_ => false` catch-all), so the
        |   wedge bond is written explicitly via smiles_token() -> "/"
        v
OUTPUT: e.g. "C(F)(Cl)/Br"  -- NOT the same as dropping the stereo silently.
```

**This is worse than a silent drop, and was confirmed by actually running it**
(not just reading the writer's source): the wedge annotation is re-encoded as
an OpenSMILES *directional bond marker*, which is only semantically meaningful
adjacent to a double bond for E/Z. Here there is no double bond at all. Fed
back through chematic's own parser, the token round-trips as `BondOrder::Up`
(self-consistent), but fed to RDKit (`Chem.MolFromSmiles`), it parses without
error and is silently dropped from the canonical output (`FC(Cl)Br`, no `/`,
verified against the venv's `rdkit==2026.03.3`) — i.e. **any spec-compliant
consumer discards it identically to the case where nothing was written**, so
the extra token is not a mitigant, just a second way to reach zero information
transfer while looking superficially like something survived.

The "obvious fix" — wire `apply_stereo_from_2d` into readers — would **not**
close this gap by itself: that function computes the correct R/S answer
(verified: it returns `S` for a real 4-neighbor test fixture) but writes
only `Atom.cip_code`. `chematic_smiles::write`/`canonical_smiles` read only
`Atom.chirality` (plus `stereo_neighbor_order` for the canonical writer's
DFS-order correction) — `cip_code` has zero readers in `writer.rs`/`canonical.rs`.
Making a stereocenter round-trip through SMILES requires setting **both**
`chirality` and a consistent `stereo_neighbor_order` together (the only place
in the repo that currently does this correctly is stereoisomer enumeration,
`chematic-chem/src/stereo.rs:150-199`) — `cip_code` (an absolute R/S label) and
`chirality` (a neighbor-order-relative local parity) are answers to two
different questions, and today's 2D-perception code only computes the first.

CDXML is the one existing precedent for wiring 2D perception into a reader at
parse time (`assign_ez_from_2d` called from `cdxml.rs:159`) — but it's E/Z
only, and even that write is a dead end for SMILES purposes since `cip_code`
still has no writer-side reader.

## 4. Frozen diagnostic fixtures

14 hand-built MOL V2000 fixtures, one per required mechanism, generated
programmatically (not hand-typed fixed-width text) by
`crates/chematic-mol/examples/stereo2d_fixture_dump.rs`. This example calls
only existing public APIs (`mol2000::parse_mol_with_coords`,
`chematic_perception::{assign_stereo_from_2d, apply_stereo_from_2d,
assign_ez_from_2d}`, `chematic_smiles::{write, canonical_smiles}`) exactly as
an external caller would — it does not modify any reader, writer, or the
stereo2d module.

Reproduce:
```bash
cargo run -p chematic-mol --example stereo2d_fixture_dump \
    > validation/results/stereo2d_fixture_dump.jsonl
.venv/bin/python scripts/stereo2d_diagnosis.py
```

The Python script re-parses each fixture's raw MOL block with RDKit (pinned
commit `8afba32ec539dcb2369bc84549d802aca3f7eb39`, matching the installed
`rdkit==2026.03.3`), extracts RDKit's own chiral-tag/CIP/E-Z verdict via
`rdCIPLabeler`, and classifies the (chematic, RDKit) pair into a failure
bucket. Output: `validation/results/stereo2d_diagnosis_summary.json`.

The script is fail-closed by construction, not just by convention: it checks
the fixture-ID set against an explicit expected set (extra or missing IDs
abort the run), rejects duplicate IDs, uses a bucket *whitelist* rather than
a `startswith("unexpected")` string check (so a new, unrecognized bucket name
can't silently slip through as "not unexpected"), verifies the specific
evidence each bucket claims (e.g. the contradictory-wedge bucket asserts two
direction tokens actually appear, the coord-mismatch bucket asserts a
non-empty result actually came back, the degenerate-coordinate bucket checks
RDKit's own side agrees), and runs two self-tests of its own fail-closed
logic before touching any real data. `sys.exit(1)` on any of these. All of
this was verified empirically during this review round, not just written and
trusted: duplicate-ID, missing-ID, unknown-mechanism, and weakened-evidence
inputs were each injected by hand and confirmed to produce exit code 1
before being reverted.

**Result: 14/14 fixtures classified, 0 unexplained, exit code 0.**

| # | Fixture | Mechanism | chematic result | Failure bucket |
|---|---|---|---|---|
| 1 | `tetrahedral_3heavy_implicit_h` | 3 heavy + implicit H | `assign_stereo_from_2d` → `[]` (RDKit assigns a tag) | `rs_not_computed_3heavy_implicit_h_gap` |
| 2 | `tetrahedral_4neighbors_explicit_h` | 3 heavy + 1 *explicit* H = 4 neighbors (NOT 4 heavy atoms — see fixture #3) | `[{"cip_code":"S"}]`, but `chirality` never set; naive SMILES emits a stray `/` on the wedge bond | `rs_computed_but_writer_emits_meaningless_bond_direction_token` |
| 3 | `tetrahedral_4heavy_no_h` | genuinely 4 heavy atoms, C(F)(Cl)(Br)(I), zero H anywhere (implicit or explicit) — added in review to stop conflating this with #2 | `[{"cip_code":"R"}]`, same shape as #2 (chirality never set, stray `/` token) | `rs_computed_but_writer_emits_meaningless_bond_direction_token` |
| 4 | `solid_wedge_only` | solid wedge, 3 heavy + implicit H | `[]` | `rs_not_computed_despite_rdkit_success` |
| 5 | `dashed_wedge_only` | hash wedge, 3 heavy + implicit H | `[]` | `rs_not_computed_despite_rdkit_success` |
| 6 | `wedge_atom_order_reversed` | bond line atom1=substituent, atom2=center (non-standard) | `[{"cip_code":"S"}]`, matches RDKit's own reading of the *same* non-standard file | `wedge_atom_order_reversed_agrees_with_rdkit_on_same_file` |
| 7 | `multiple_stereocenters` | 2,3-dibromobutane, 2 independent centers | `[]` for both (same 3-heavy gap as #1) | `rs_not_computed_3heavy_implicit_h_gap` |
| 8 | `no_wedge_negative_control` | same skeleton as #7, no wedges | `[]`, RDKit also assigns nothing | `correctly_no_stereo_both_agree` |
| 9 | `cip_priority_tie` | C(CH3)(CH3)(F)(Cl), 2 branches tie | `[]` (tie correctly detected pre-geometry) | `correctly_no_stereo_both_agree` |
| 10 | `degenerate_2d_coordinates` | same graph as #1, all atoms at (0,0) | `[]` (CIP ranks fine, geometry degenerate → `v≈0`); RDKit's own side independently confirmed to also assign nothing | `degenerate_coords_correctly_yields_no_stereo` |
| 11 | `ez_geometry_2butene` | defined cis 2D layout, no wedge | `assign_ez_from_2d` → `Z` correctly, but naive SMILES has zero `/`/`\` (no wedge bonds exist to tokenize) | `ez_computed_but_no_bond_direction_for_writer` |
| 12 | `terminal_alkene_propene` | terminal `=CH2` | `[]`, RDKit also assigns nothing | `correctly_no_stereo_both_agree` |
| 13 | `contradictory_wedge_annotations` | 2 wedges (both "up") from the same center | both silently tokenized independently as `/` (verified: exactly 2 direction tokens present, not just "the mechanism ran") — no consistency check anywhere in the pipeline | `no_consistency_check_both_wedges_silently_tokenized` |
| 14 | `coord_atom_count_mismatch` | 3 coords passed for a 5-atom molecule (API misuse, not a file case — a truncated MOL *file* just fails to parse before perception is reached) | **no panic, but not a safe no-op either**: out-of-range neighbors silently fall back to the *center's own* (x,y) via `unwrap_or(*center_pos)` in `assign_rs`, and a CIP code (`S`) is still returned from this corrupted geometry instead of `None`/an error (verified: the result is checked to be genuinely non-empty, not assumed). That it matches the true answer here is coincidental to this fixture's geometry, not a property of the fallback | `silent_result_from_corrupted_fallback_positions_not_error` |

Fixture #2 was originally named `tetrahedral_4heavy_explicit_h`, which was
inaccurate — its center has 3 heavy neighbors plus 1 explicit H atom, not 4
heavy atoms — and its result had been (incorrectly) cited elsewhere in an
earlier draft of this document as evidence about a "genuinely 4-heavy-atom"
center. Renamed to `tetrahedral_4neighbors_explicit_h`, and fixture #3
(`tetrahedral_4heavy_no_h`, C(F)(Cl)(Br)(I)) added so the two shapes — "H
present as a real atom that RDKit's `removeHs` would later strip" vs "no H
anywhere, `removeHs` is a no-op" — are tested and reported separately rather
than conflated under one claim.

Two additional findings surfaced only by *running* the fixtures, not visible
from reading source alone:

- **Coordinate-symmetry trap in `assign_rs`.** An early, more "natural"-looking
  coordinate choice for fixture #2 (methyl/wedge substituents placed at mirror
  positions) produced an *exactly* zero signed volume once fixture #6's
  atom-order reversal flipped one wedge's z-sign — not a chematic bug, but a
  reminder that `assign_rs`'s coplanarity check is exact-algebraic, not just a
  numerical-noise guard, and trivially-symmetric fixture geometry can
  accidentally land on it. Fixed by using an asymmetric layout (see the code
  comment on fixture #2 in `stereo2d_fixture_dump.rs`).
- **RDKit rejects a layout the original chematic-perception unit test uses.**
  The existing test `stereo2d.rs::test_r_s_bromochlorofluoromethane` places two
  wedged substituents exactly 180° apart through the center. Feeding the
  equivalent MOL file to RDKit produces `Warning: ambiguous stereochemistry -
  opposing bonds have opposite wedging - at atom 0 ignored` and RDKit assigns
  **no** chiral tag at all — even though chematic's simpler geometric model
  computes one confidently. This is a genuine algorithm-behavior divergence
  (RDKit source audit §5, below, confirms this is a deliberate, documented
  guard in `atomChiralTypeFromBondDirPseudo3D`), not a bug in either engine,
  and is exactly the kind of case design question (j) below is about.

## 5. RDKit source audit (pinned commit `8afba32ec539dcb2369bc84549d802aca3f7eb39`)

Two independent research passes, all claims cite exact file/function/line at
the pinned commit (fetched via `raw.githubusercontent.com`).

### 5a. MolFromMolBlock flow, sanitize ordering, 3-heavy+implicit-H handling

**The single most decision-relevant fact in this whole audit**: RDKit
perceives wedge-based chirality **before** sanitization — deliberately, with
an explicit comment explaining why — and that perception is **not** gated by
the `sanitize` flag at all. Only the later CIP-labeling step is.

**Ordered call sequence**, `FileParserUtils::finishMolProcessing`
(`MolFileParser.cpp:3457-3533`), the single code path every public entry
point (`MolFromMolBlock`, `MolFromMolFile`, the legacy `MolDataStreamToMol`
shim) converges on:

1. `atom->calcExplicitValence(false)` for every atom (need degree/valence known).
2. `ProcessMolProps` (misc mol-file flag postprocessing).
3. **Chirality *tag* assignment — unconditional**, gated only by
   `chiralityPossible || conf.is3D()` (a structural check: are there any
   wedge/dash bonds, or is this a 3D conformer at all — **not** `params.sanitize`):
   `MolOps::assignChiralTypesFromBondDirs(*res, conf.getId(), true)` for the
   2D case (the function chematic's future Stage 1 must mirror), or
   `MolOps::assignChiralTypesFrom3D` for a real 3D conformer.
4. `Atropisomers::detectAtropisomerChirality`.
5. `MolOps::clearSingleBondDirFlags` — the wedge annotation is dropped from
   the bond *only after* the tag is safely captured on the atom.
6. **Only if `params.sanitize`** (default `true`, but a real, user-facing
   toggle): `sanitizeMol`, `detectBondStereochemistry` (E/Z), optional
   `removeHs`, then `MolOps::assignStereochemistry(*res, true, true, true)`
   (CIP `_CIPCode` labeling + bond-stereo refinement).

**Why this order, in RDKit's own words** (comment directly above step 3,
`MolFileParser.cpp:3480-3484`):
> "we detect the stereochemistry before sanitizing/removing hydrogens because
> the removal of H atoms may actually remove the wedged bond from the
> molecule. This wipes out the only sign that chirality ever existed..."

This is a concrete, load-bearing design constraint for chematic, not just a
RDKit implementation detail: **if a future chematic Stage 1 (wedge→parity)
step is ever deferred until after some other mutation (H removal,
kekulization, a lazy/on-demand call triggered by the caller at an arbitrary
later point) that mutation can silently destroy the very wedge bond the
perception step needs.** Whatever chematic builds must run at, or immediately
after, parse time — before anything else gets a chance to touch the H atoms
or bonds carrying the wedge annotation.

**`sanitize=False` behavior**: chiral *tags* are still set (step 3 sits
outside the `if (params.sanitize)` block entirely); what's skipped is CIP
`_CIPCode` labeling, bond-stereo refinement/cleanup, and (for double bonds)
`detectBondStereochemistry` moves from the `sanitize` branch to a bare `else`
call — so E/Z direction detection *always* runs too, just via a different
branch. Net effect: **with `sanitize=False`, a parsed RDKit mol has chiral
tags and (single-bond) `BondDir` values, but no R/S/E-Z human-readable
labels.** The `v1` legacy API's own doc comment describing `sanitize` as
gating "sanitization and stereochemistry perception" is measurably
imprecise — it's really CIP-labeling that's gated, not tag/direction
assignment.

**3-heavy-neighbor + implicit-H handling — first-class, not a fallback.**
`atomChiralTypeFromBondDirPseudo3D`'s entry gate is `nNbrs < 3 || nNbrs > 4`
(`Chirality.cpp:537-545`, comment: *"we can implicitly add a single H to
3 coordinate atoms"*), with a dedicated ordering/parity code path for
`nNbrs == 3` (lines 644-654) and explicit, IUPAC-guideline-cited
(`ST-1.2.10`/`ST-1.2.12`) rejection of ambiguous 3-coordinate geometries
(T-shaped/collinear substituents, conflicting wedge directions,
lines 718-751) — i.e. RDKit's own version of design question (j) for this
specific case is already answered in its source: reject with the same
warn-and-skip pattern as the 4-neighbor case (§5b), not a special error.
After the tag is set, RDKit **materializes the implicit H as an explicit H
count** so CIP ranking downstream has a concrete neighbor to work with
(`assignChiralTypesFromBondDirs`, `Chirality.cpp:3801-3810`:
`atom->setNumExplicitHs(1)` when `getDegree()==3 && getNumImplicitHs()==1`)
— and the *same* 3-neighbor-plus-H case is independently special-cased a
second time in CIP-code assignment itself (`assignAtomChiralCodes`,
lines 1797-1800: `if (nbrIndices.size() == 3 && atom->getTotalNumHs() == 1)
++nSwaps;`). Two independent code paths treating this shape as fully
supported, not an edge case, is strong confirmation for design question (h):
chematic's `assign_rs` needs an equivalent "treat the implicit H as a
synthetic 4th neighbor" step, not a bail-out.

**Naming note**: the function name guessed going into this audit,
`DetectAtomStereoChemistry`, does exist but is a deprecated one-line
pass-through to `MolOps::assignChiralTypesFromBondDirs`
(`MolFileStereochem.cpp:201`, header-annotated as deprecated) — the current
canonical name is `assignChiralTypesFromBondDirs`, used throughout this
section.

**One more confirmed, previously-suspected divergence** (not one of the
original 5 questions, found incidentally while reading `ParseMolFileBondLine`,
`MolFileParser.cpp:1778-1803`): MDL bond-stereo code `4` ("either"/unknown
single bond) maps to `Bond::BondDir::UNKNOWN` in RDKit — a **third**,
distinct state from a definite wedge (`1`) or hash (`6`). chematic's
`mol2000.rs` (`stereo_raw` match, `1 | 4 => BondOrder::Up`) currently
collapses codes `1` and `4` into the same `Up` value, losing the
"direction not actually specified" signal code `4` is meant to carry. Noted
as a known gap for whoever picks up P1-S1a/S1b — not fixed here (this PR
changes no reader code).

**Not fully verified** (flagged rather than guessed, per the audit's own
open-questions list): `Chirality::assignChiralTypesFromMolParity` reads the
V2000 atom-line "parity" field but has zero callers found in the four files
searched — it may be reachable from some other entry point not traced here.
`MolOps::assignStereochemistry`'s own three `bool` parameters
(`cleanIt`/`force`/`flagPossibleStereoCenters`) were not traced in detail,
since they weren't required to answer the five audit questions.

### 5b. 2D E/Z detection, terminal alkenes, chiral-tag-vs-CIP-code, Python/C++ parity, degenerate/contradictory-wedge handling

**RDKit's model is a genuine two-stage pipeline — geometry→parity is a
different function, run at a different time, than priority→label:**

1. **Stage 1 (always runs): geometry → direction/parity, no CIP involved.**
   - Double bonds: `MolOps::detectBondStereochemistry` → `setDoubleBondNeighborDirections`
     → `updateDoubleBondNeighbors` (`Chirality.cpp:164`) reads 2D coords,
     computes the dihedral, and sets `Bond::BondDir` (`ENDUPRIGHT`/`ENDDOWNRIGHT`)
     on the *adjacent single bonds* — never sets a final E/Z label on the
     double bond itself (only `STEREOANY` for ambiguous cases).
   - Atoms: `assignChiralTypesFromBondDirs` → `atomChiralTypeFromBondDirPseudo3D`
     (`Chirality.cpp:428-846`) computes `Atom::ChiralType` (`CHI_TETRAHEDRAL_CW`/`CCW`)
     from wedge z-offset + 2D xy + raw neighbor-bond order — **never reads a
     CIP rank**. This is exactly chematic's missing "produce something a
     SMILES writer can consume" step.
   - Both are called **unconditionally** in `MolFileParser.cpp`'s
     `finishMolProcessing` (lines 3505-3527), in every sanitize branch
     including `sanitize=False`.
   - The SMILES writer consumes `BondDir` directly:
     `SmilesWrite.cpp::GetBondSmiles` maps `ENDDOWNRIGHT → "\\"`, `ENDUPRIGHT → "/"`.

2. **Stage 2 (only runs when sanitizing): priority → label.**
   - `assignBondStereoCodes` (E/Z) and `assignAtomChiralCodes` (R/S,
     `Chirality.cpp:1742-1822`) run only from `MolOps::assignStereochemistry`,
     require CIP ranks, and are what sets `Bond::STEREOE`/`STEREOZ` and the
     `_CIPCode` property.
   - **CIP-priority ties only suppress Stage 2's label, never Stage 1's
     tag/direction**: `isAtomPotentialChiralCenter` sets `hasDupes=true` on a
     rank collision (`Chirality.cpp:1729`); `assignAtomChiralCodes` skips
     setting `_CIPCode` when `hasDupes` (line 1778) but the chiral tag from
     Stage 1 is untouched.

3. **Terminal alkenes**: two independent degree-based guards —
   `isBondCandidateForStereo` (`Chirality.cpp:382-390`, requires
   `getDegree() > 1` on both double-bond atoms) and `assignBondStereoCodes`'s
   own degree-2-or-3 check (`Chirality.cpp:1850-1851`) — both silently skip,
   no exception.

4. **Python vs C++ defaults: no divergence.** `Code/GraphMol/Wrap/rdmolfiles.cpp`'s
   `MolFromMolBlock` binding calls the identical `MolDataStreamToMol` used by
   the C++ v1 API, with byte-identical defaults
   (`sanitize=true, removeHs=true, strictParsing=true`,
   `FileParsers.h:47-51` vs `rdmolfiles.cpp:1008-1013`). The only Python-side
   difference is catching `FileParseException` into a logged warning + `None`
   instead of propagating a C++ exception — not a stereo-behavior change.

5. **Degenerate/contradictory handling — extensive, explicit, warn-and-skip
   (never an exception):** all inside `atomChiralTypeFromBondDirPseudo3D`:
   zero-length bond (`Chirality.cpp:518-523`), overlapping neighbor directions
   (552-557), **opposing bonds with opposite wedging** (709-712 — exactly the
   case chematic's own unit-test fixture triggers, §4 above), an explicit
   3-coordinate wedge-contradiction check citing IUPAC ST-1.2.10 (745-749),
   and zero chiral volume at both the 4-coordinate (806-809) and 3-coordinate
   (838-841) stages. The double-bond analogue, `isLinearArrangement`
   (`Chirality.cpp:100-113`), treats a near-zero-length substituent vector as
   linear/ambiguous (covers all-identical/all-zero coordinates) and sets
   `STEREOANY` rather than guessing.

**Design-relevant takeaway**: RDKit's Stage 2 (labels) is defined *in terms
of* Stage 1's output (`BondDir`/`ChiralType`), never the reverse. chematic
currently only has something like a Stage-2-flavored output (`cip_code`, an
absolute R/S/E/Z label) and has no Stage-1 equivalent a SMILES writer could
consume. Building the Stage-1 equivalent (parity/direction, in a
neighbor-order-relative frame matching what `@`/`@@` and `/`/`\` need) is the
correct next building block — not a fix to what `apply_stereo_from_2d`
writes its answer *into*, since its answer (R/S) is the wrong *kind* of value
regardless of field name.

## 6. Design questions — answers

| # | Question | Answer | Rationale |
|---|---|---|---|
| a | Should low-level parse functions permanently retain raw wedge info? | **Yes** | Already true today (`BondOrder::Up`/`Down` on bonds survive parse); no reader currently discards it. The gap is downstream perception, not raw retention. |
| b | Should default parse auto-perceive stereo? | **Yes, for a Stage-1-equivalent (parity/direction) step — once it exists** | Revised after §5a: RDKit's real precedent is that Stage 1 (`assignChiralTypesFromBondDirs`) is gated only by *structural* possibility (`chiralityPossible \|\| conf.is3D()`), never by a sanitize-like flag — it is unconditional whenever there's a wedge bond to read. It also must run at/immediately-after parse time specifically because a later mutation (H removal, in RDKit's case) can destroy the wedge annotation before a deferred/lazy call would ever see it (§5a's `finishMolProcessing` comment). This overrides an earlier draft answer of "no" — Stage 1 being cheap, CIP-free, and structurally-gated (not caller-opt-in) is exactly why RDKit runs it unconditionally. |
| c | Is an opt-in option needed? | **Only for CIP-priority labeling (a Stage-2-equivalent), not for the Stage-1 parity/direction step** | Matches RDKit's real split precisely: `sanitize` (default `true`, a genuine opt-out) gates only `assignStereochemistry`/CIP labeling; the chirality-tag step it gates nothing. chematic should expose the Stage-2-equivalent (CIP R/S/E-Z labeling) as toggleable — mirroring `assign_stereo_from_2d`'s current existence as an explicit call — while Stage 1 becomes reader-default behavior once built. |
| d | Should only an RDKit-compatible Python wrapper auto-run it? | **No** — the Rust core should own the Stage-1 step so every binding (Python, WASM, MCP) gets it for free; a Python-only wrapper would leave Rust callers and WASM with the gap | Consistent with the project's "pure-Rust RDKit alternative" positioning (per `chematic-mol`'s own `Cargo.toml` description) — behavior shouldn't depend on which binding a caller uses. |
| e | Should the Rust API use an explicit conversion function? | **Yes for the Stage-2/CIP-labeling call (mirrors today's `assign_stereo_from_2d`); Stage 1 (parity/direction) should be reader-internal, called from within each format's parse function itself, not a separate opt-in call the caller must remember to make** | Directly follows from (b)/(c): if Stage 1 must run before any mutation can destroy the wedge annotation (§5a), it cannot be a function a caller invokes later at their discretion — CDXML's existing `assign_ez_from_2d` call from inside `cdxml.rs:159` (§2/§3) is already the right shape for this, just needs a Stage-1 (parity, not just E/Z-label) equivalent and needs porting to the other wedge-bearing readers (V2000, SDF, MRV once merged, KET). |
| f | Should the original wedge survive after conversion? | **Yes** | `BondOrder::Up`/`Down` costs nothing to retain and is the only way to detect/debug a conversion disagreement later (as this diagnosis itself needed to). |
| g | CIP code or SMILES-native chirality as source of truth? | **`chirality` (+ `stereo_neighbor_order`) for round-tripping; `cip_code` remains the human-facing label, derived, never the writer's input** | Directly matches RDKit's own architecture (§5): the writer-consumable value is the parity/direction (Stage 1), not the CIP label (Stage 2). Verified in chematic: `cip_code` has zero readers in `writer.rs`/`canonical.rs` today. |
| h | What must be extended to include implicit-H centers in P1-S1a? | The parity computation (today's `assign_rs`'s neighbor collection, `stereo2d.rs:223-226`, currently `if nbs.len() != 4 { return None }`) needs a `STEREO_H_SENTINEL`-based substitution path when `nbs.len() == 3`, mirroring the existing pattern in `chematic-cip/src/assign.rs`. Critically, this belongs in the **parity** step (S1a), not the labeling step (S1b) — RDKit's own 3-neighbor handling (§5a) sits entirely inside `atomChiralTypeFromBondDirPseudo3D`/`assignChiralTypesFromBondDirs`, before any CIP ranking is computed. | `STEREO_H_SENTINEL` (`chematic-core/src/molecule.rs:43`) already exists and is used exactly this way elsewhere; this is a known, proven pattern to port, not new design. Confirmed by fixtures #1/#4/#5/#7 in §4: RDKit resolves these, chematic currently returns `None` for all of them. |
| i | Should tetrahedral and E/Z be separate PRs? | **Yes** | RDKit itself treats them as fully separate functions/data (`detectBondStereochemistry` vs `atomChiralTypeFromBondDirPseudo3D`, different Bond/Atom fields). chematic's own gap is structurally separate too: fixture #11 shows E/Z perception (`assign_ez_from_2d`) already computes the right label but has *no* bond-direction output step at all (a different missing piece than the tetrahedral case, which computes an R/S answer but writes it to the wrong-shaped field). |
| j | No-coords / degenerate-coords / contradictory-annotation: error, warning, or unspecified? | **Warning (log) + unspecified (`None`/no assignment), matching RDKit's own "ambiguous/conflicting stereochemistry ... ignored" pattern — never a hard error** | RDKit (§5b point 5) never throws for any of these; it warns and proceeds with unassigned stereo. chematic's current behavior already matches this for degenerate coordinates and CIP ties (fixtures #9, #10) but has two real gaps to close: (1) the coord/atom-count-mismatch fallback (fixture #14) produces a *silent, possibly-wrong* answer instead of `None`+warning — this should change to return `None` for any atom whose own or neighbors' coordinates are out of range, not degrade to a corrupted position; (2) contradictory wedges (fixture #13) are currently accepted with zero detection anywhere in the pipeline — a future S1a should port something like RDKit's "opposing bonds have opposite wedging" / "bond wedging contradiction" checks (`Chirality.cpp:709-712`, `745-749`). |

## 7. Recommended final split

**This section originally had S1a/S1b in the wrong order** (R/S-first, then
derive `chirality` from it) — backwards relative to this same document's own
§5a/§5b RDKit audit, which shows RDKit computes the chiral *tag* (parity)
first, with **zero** CIP involvement, and only *afterward*, as a separate
optional stage, ranks CIP priority to produce an R/S *label*. Caught in
review before any implementation started. The corrected order below computes
parity before labeling, matching RDKit and matching this document's own §3
pathway trace and §5 audit. Getting this backwards would have been a real
implementation defect, not just a naming issue: an R/S-first design cannot
produce `Atom.chirality` at all for a CIP-priority tie, a CIP-unresolvable
case, a pseudoasymmetric center (a label computed in a later CIP rule than
chematic implements today), or simply a caller who only wants the drawn
local parity preserved losslessly without paying for CIP ranking — R/S would
be undefined or unavailable in every one of those cases, and a
"derive-chirality-from-R/S" step would have nothing to derive from. Parity
never needed CIP ranking in the first place (§5b point 1: `atomChiralTypeFromBondDirPseudo3D`
"never reads a CIP rank"), so building it first avoids the dependency
entirely.

1. **P1-S1a** — wedge + 2D coordinates + the neighbor order already available
   at reader-parse time → `Atom.chirality` + `Molecule::stereo_neighbor_order`.
   CIP-independent (mirrors `assignChiralTypesFromBondDirs`/
   `atomChiralTypeFromBondDirPseudo3D`, §5a/§5b). Must include: the
   3-neighbor+implicit-H case (design question h), a degenerate-coordinate
   guard, and a contradictory-wedge guard (design question j) — all of these
   are parity-stage concerns in RDKit's own model, not labeling-stage ones.
   Runs automatically from inside each wedge-bearing reader's own parse
   function (design questions b/c/e), immediately after the wedge bond is
   read, before anything else can mutate it away.
2. **P1-S1b** — local parity (from S1a) + CIP ranking → `Atom.cip_code`
   (R/S label). Mirrors `assignAtomChiralCodes` (§5b point 2): a labeling
   stage built **on top of** S1a's already-computed parity, consuming it,
   never replacing or gating it. On a CIP-priority tie, `chirality` (from
   S1a) is preserved exactly as-is and `cip_code` stays `None` — this
   directly matches RDKit's `hasDupes` behavior (§5b point 2) and is only
   possible because parity never depended on CIP ranking to begin with. If
   chematic wants to keep `assign_stereo_from_2d`'s current public shape (an
   R/S-label-producing function), it should become S1b's implementation,
   internally consuming S1a's parity result rather than recomputing geometry
   independently from scratch, as it does today.
3. **P1-S2** — 2D E/Z → SMILES `/`/`\` bond directions. Independent of
   S1a/S1b (different data: `assign_ez_from_2d` already computes the right
   *label*; what's missing is a direction-setting step analogous to RDKit's
   `setDoubleBondNeighborDirections`, not a labeling fix).

Per design questions (b)-(e): S1a alone is enough to make wedge-derived
tetrahedral centers round-trip through SMILES (it produces exactly what
`chematic_smiles::write`/`canonical_smiles` already read). S1b adds
human-readable R/S labels on top and is not required for round-tripping.
Once S1a exists, the readers in §2 that carry wedge bonds (V2000, SDF, V3000
once its bond-`CFG` reading gap is fixed, KET, and MRV once merged) should
call it from *inside* their own parse function, unconditionally whenever a
wedge bond is present — mirroring both RDKit's `finishMolProcessing`
ordering and chematic's own existing CDXML precedent (`cdxml.rs:159`) — not
expose it as a separate function callers must remember to invoke after the
fact. S1b (CIP-label perception, if chematic wants an
`assignStereochemistry`-shaped API at all) can remain an explicit,
caller-invoked, opt-in function, matching `assign_stereo_from_2d`'s current
shape.

None of P1-S1a/S1b/S2 should also fix the `BondOrder::Up`/`Down` →
meaningless-`/`-token writer bug found in §3/§4 as a side effect — that's a
distinct, narrower writer-side bug (the writer should only emit a directional
token when the bond is actually adjacent to a resolved double bond) worth its
own tiny follow-up, flagged here so it isn't lost, not folded into the
stereo-perception work.

## 8. What this PR does not do

- Does not call stereo perception from any reader.
- Does not modify `Atom.chirality`, `Molecule::stereo_neighbor_order`, or any
  SMILES writer.
- Does not modify `assign_stereo_from_2d`/`apply_stereo_from_2d`/`assign_ez_from_2d`.
- Does not rebase/restack `#126`/`#127`/`#128`.
- Is not merged.
