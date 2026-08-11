# Platinum coordination-chemistry compatibility: feasibility, baseline, and fix

**Scope.** This is a cheminformatics representation/interoperability
benchmark, not an anticancer-activity project. It does not implement or
claim anything about IC50, cisplatin resistance, toxicity, pharmacokinetics,
or clinical efficacy. The question is narrower and purely structural:

> Can chematic correctly represent, parse, round-trip, characterize, and
> compare anticancer platinum coordination complexes, without silently
> corrupting their coordination chemistry?

**Killer claim measured by this benchmark:** cisplatin-family structures
round-trip with preserved formula, charge, and coordination topology.
Cis/trans identity preservation is **measured and reported as a real
gap**, not silently claimed as solved (see "P0 finding not fixed" below).

**Do not confuse this corpus with `validation/manifests/dataset_provenance.json`'s
existing `platinum_diverse_dataset` entry** — that is an unrelated,
previously-referenced third-party dataset (Genentech's "Platinum Diverse
Dataset" of small-molecule drug conformers, used for a completely different
3D-conformer benchmark). This document and `pt_corpus.jsonl` are about
platinum *coordination complexes* (the metal, not a dataset name), and share
no data or code path with that entry.

## 1. Corpus

**Target:** the task requested 100–500 structures if reachable. This
benchmark shipped with **18** hand-verified entries and explicitly does not
reach that target. Reason, per the task's own explicit permission (§2, "don't
fabricate weak-ground-truth structures for the sake of count"): public
structured databases essentially never provide platinum coordination bonds
as real bonds (see §2 below) — every usable source for this benchmark had to
be built compound-by-compound, cross-checking formula/charge/mass against
PubChem's numeric fields (reliable) while independently authoring the actual
bonded connectivity from established inorganic chemistry (crystal structure
literature / textbook coordination geometry), because PubChem's own
connectivity fields are not usable. Scaling this to 100+ compounds at the
same verification bar was not achievable in the time available without
either (a) resorting to unverified/scraped structures, which the task
explicitly forbids, or (b) a dedicated CSD (Cambridge Structural Database)
license and query, which is out of scope here (see "Not attempted" below).

**File:** `validation/platinum/pt_corpus.jsonl` — one JSON object per line.
Fields: `id`, `name`, `category` (list, matching §13's requested buckets:
`pt2_classical`, `pt2_trans`, `pt2_chelate`, `pt4_prodrug`, `charged`,
`neutral`, `sulfur_ligand`, `oxygen_ligand`, `nitrogen_ligand`,
`carbon_ligand`, plus a `generalization_gate` bucket for the two non-Pt
sanity-check entries), `identity_group` (which structures must/must-not
canonicalize to the same identity), `cis_trans`, `pt_oxidation_state`,
`pt_coordination_number`, `geometry`, `ligand_donor_atoms`,
`smiles_dative` (the actual chematic-parseable input), `formula_expected`,
`charge_expected`, `exact_mass_expected` (`null` where not independently
sourced — see below), and a `source` object recording exactly where each
number came from and what was corrected and why.

**Composition:**

| category | entries |
|---|---|
| Classical Pt(II) drugs (cis) | cisplatin, carboplatin, oxaliplatin, nedaplatin, lobaplatin, picoplatin, dicycloplatin |
| Pt(II) trans | transplatin |
| Pt(IV) prodrugs | satraplatin, iproplatin, tetraplatin, oxoplatin |
| Charged species | tetrachloroplatinate(II) anion, cisplatin's diaqua activation dication |
| Non-drug reference compounds (S/C-ligand diversity, §13) | cis-[PtCl2(DMSO)2] (sulfur donor), a minimal Pt–C sigma-bond reference compound (carbon donor) |
| Generalization-gate (non-Pt, §16) | [Co(NH3)6]3+, [Fe(H2O)6]2+ |

Every Pt(II)/Pt(IV) drug requested by the task's §1 "classical Pt(II)" list
is present **except heptaplatin**, and neither `phenanthriplatin` nor
`pyriplatin` (task §1's suggested research-compound examples) could be
resolved: PubChem's compound-name search and its autocomplete API returned
no match for any of the three under those names (checked directly, not
inferred) in the time available. Rather than fabricate structures for them
from memory, they are omitted and reported as gaps here. `heptaplatin` in
particular is one of the 7 explicitly requested compounds — if this
benchmark is extended, resolving its correct name/CID (or sourcing it from
its original literature/patent) is the most direct way to close that gap.

**"≥1 Pt(IV) oxaliplatin-derived prodrug" (task §1) was not resolved
either** — the corpus's 4 Pt(IV) entries are all cisplatin-lineage
(oxoplatin, satraplatin) or isopropylamine/DACH-based (iproplatin,
tetraplatin), not oxaliplatin-derived specifically. A literature search for
a well-documented, PubChem-resolvable oxaliplatin→Pt(IV) prodrug was not
completed in the time available; this is a corpus gap, not a chematic
result.

**Not attempted / explicitly out of scope for this pass:**
- ChemRxiv/journal Supporting Information mining for research-stage Pt(II)/Pt(IV)
  complexes (the task's §1 "research compounds" tier). Time-boxed effort went
  to PubChem-verifiable structures instead, which gave a higher
  confidence-per-hour yield. This is the most likely place to find the
  remaining structures needed to approach the 100–500 target.
- Cambridge Structural Database (CSD) — the actual primary source for
  crystallographically-confirmed Pt coordination geometry — requires a
  license not available in this environment. All coordination topology in
  this corpus is instead sourced from well-established, multiply-cited
  textbook/review inorganic chemistry (see each entry's `source.note`), which
  is appropriate confidence for the classical drugs used here (their
  structures have been confirmed by crystallography for decades and are not
  scientifically controversial), but would not be an adequate bar for
  research-stage / structurally-uncertain compounds — which is exactly why
  none of those were included.
- Redistributable dataset gathering — everything in this corpus is either
  small-molecule structural data derived from public PubChem records (US
  government/NIH data, redistributable) or hand-authored by this benchmark
  (original, redistributable). No CSD-derived or other non-redistributable
  data is vendored, per the task's explicit instruction.

## 2. A major, independently-confirmed finding about existing "canonical" sources

Before any chematic code was touched, cross-checking PubChem's own records
for these compounds surfaced a finding that reshaped the whole corpus
methodology (§3's "don't treat RDKit/PubChem parsing success as ground
truth" warning, borne out concretely):

**PubChem's default SMILES/InChI, and even its 2D SDF record, do not encode
Pt–ligand coordination bonds as bonds for these compounds.** For cisplatin
(PubChem CID 5702198), `PUBCHEM_SMILES` is `N.N.Cl[Pt]Cl` — three
disconnected fragments (`PUBCHEM_COMPONENT_COUNT: 3`) — and the bundled 2D
SDF record's own bond table has **zero** N–Pt bonds (only the two Pt–Cl
bonds and the ammonia N–H bonds are present). This is not a chematic
parsing failure to react to — it is the *starting representation* most
downstream consumers of "the PubChem structure of cisplatin" would inherit,
and it already fails the killer benchmark's condition (B) (coordination
topology preserved) before any parser is even involved, because the topology
was never there to begin with. This is why this benchmark's corpus is not
"PubChem's SMILES for named compound X" — it is "chematic-authored
connectivity, cross-checked at the formula/charge/mass level against
PubChem's numeric property fields (which are computed correctly, unlike the
connectivity fields)".

**A second, independent finding: several PubChem property records for these
exact compounds are internally self-contradictory**, discovered by
cross-checking, not assumed:

- `transplatin` name search resolves to CID 441203, formula
  `"Cl2H6N2Pt+2"`, `Charge: 2`. A neutral trans-[PtCl2(NH3)2] with two
  covalently-drawn (formal charge 0) chlorides cannot have net charge +2 —
  this is an internally inconsistent record, not real transplatin.
- `iproplatin` (CID 155491322): formula `"C6H20Cl2N2O2Pt-4"`, `Charge: -4`.
  Real iproplatin (CHIP, CAS 62802-36-4) is neutral; the base atom
  composition (`C6H20Cl2N2O2Pt`) is correct, only the charge tag is wrong.
- `tetraplatin` (CID 6434704): `Charge: 2`, **and** its own stated
  `ExactMass` (450.952953) does not match its own stated
  `MolecularFormula`/`SMILES` — independently recomputing the exact mass of
  PubChem's own SMILES via RDKit's `Descriptors.ExactMolWt` gives 448.95480,
  not 450.952953 (verified directly, see `pt_corpus.jsonl`'s `tetraplatin`
  entry). Two independent inconsistencies in one record.

None of these three records were used as ground truth; each is corrected in
`pt_corpus.jsonl` with the specific numeric discrepancy documented in its
`source.note` field. This is exactly the scenario the task's §3/§12
warned about ("RDKit/PubChem being parseable is not the same as being
correct") — and it would have silently propagated into this benchmark's
"ground truth" if the corpus had been built by copy-pasting PubChem SMILES,
which is precisely why that shortcut was not taken.

## 3. Baseline (current `main`, before any fix)

Harness: `cargo run --release -p chematic-mol --example platinum_benchmark`
(reads `pt_corpus.jsonl`, writes one JSONL row per entry — see the file's own
doc comment for exactly what each row measures) +
`python scripts/platinum_rdkit_oracle.py` (independent RDKit-only process,
same corpus, never fed by or feeding the chematic run).

Results: `validation/results/platinum_baseline_chematic.jsonl` /
`_summary.json` (chematic, unmodified `main`), `platinum_rdkit_oracle.jsonl`
/ `_summary.json` (RDKit 2026.03.3 oracle). Denominator is 18 (the full
corpus) throughout.

| metric | before |
|---|---|
| SMILES parse success | 18/18 |
| formula matches expected | **1/18** |
| net charge matches expected | 18/18 |
| Pt coordination number matches expected | 16/18 (see note below) |
| MOL V3000 round-trip completes without error | 18/18 |
| MOL V3000 round-trip preserves formula | 18/18 |
| **MOL V3000 round-trip preserves dative bond count** | **1/18** |
| ECFP4: doesn't panic, deterministic across 2 runs | 18/18, 18/18 |
| simple Pt–Cl SMARTS (`[#78]~[#17]`): doesn't panic | 18/18 |

*Coordination-number note:* the harness's coordination-number check looks
specifically for a `Pt` atom; the 2 mismatches are the two non-Pt
generalization-gate entries (Co, Fe), which have no Pt atom at all and
report `None` — a harness scoping limitation (this check is Pt-specific by
design, since this is a Pt benchmark), not a chematic defect on those rows.
All 16 real Pt-containing entries matched.

Exact mass, before (illustrative subset; full data in the results file):

| compound | expected (Da) | chematic (before) |
|---|---|---|
| cisplatin | 298.9556 | 179.9753 |
| carboplatin | 371.0445 | 252.0641 |
| satraplatin | 499.0605 | 380.0802 |
| dicycloplatin | 515.0868 | 396.1063 |

Every single Pt complex's mass is wrong by roughly 119 Da (~30–40%
relative error, worse for smaller molecules), and the error is uniform
across completely unrelated ligand sets — a strong signal of one systemic
mechanism, not compound-specific noise (confirmed below to be two
overlapping systemic bugs).

## 4. Failure taxonomy

**FORMULA_MISMATCH / MASS_MISMATCH (P0, general, not Pt-specific).** Two
independent, compounding root causes, both found via this benchmark and
both fixed (see §6):

1. **`MASS_MISMATCH` root cause — missing periodic-table mass data.**
   `chematic-chem`'s `avg_mass`/`mono_mass` tables
   (`crates/chematic-chem/src/descriptors.rs`) covered only ~24 light
   main-group elements (H..Ca plus As/Se/Br/I) and silently fell back to
   `atomic_number as f64` for every other element — meaning **every
   transition metal, lanthanide, actinide, and heavy post-transition
   element** (platinum included: atomic number 78 was returned as "78.0 Da"
   instead of the real ~195 Da) got a wildly wrong mass with no error, no
   panic, and no diagnostic. This is not connected to coordination bonds at
   all — it would misprice the mass of *any* molecule containing iron,
   copper, zinc, ruthenium, palladium, gold, uranium, or 90+ other elements,
   including ordinary organometallics and simple ionic salts with no
   dative bonds whatsoever. Found only because this benchmark's `exact_mass`
   check compared against an independently-sourced expected value instead of
   trusting chematic's own output.

2. **`FORMULA_MISMATCH` root cause — dative bonds wrongly consume the
   donor's normal valence.** `chematic_core::valence::valence_inferred_hcount`
   summed every bond's `order_int()` (which is `1` for `BondOrder::Dative`,
   same as a plain covalent single bond) to compute implicit hydrogen count.
   For an un-bracketed dative donor atom (e.g. bare `N` in `N->[Pt]Cl`), this
   subtracted the dative bond from the donor's normal valence exactly like a
   real covalent bond would — giving `NH2` (2 implicit H) instead of the
   chemically correct `NH3` (3 implicit H), because a dative bond's whole
   point is that the donor's lone pair is *shared*, not that one of the
   donor's own covalent slots is *consumed*. Verified directly against
   RDKit's identical treatment of the same `->` SMILES syntax (RDKit:
   `[NH3]->[Pt]`, correct; chematic before fix: same string round-trips
   internally as `NH2`). In one case (`cis-[PtCl2(DMSO)2]`) the bug crossed
   a valence-table tier boundary and fabricated a chemically nonsensical
   **extra hydrogen on sulfur** (a spurious S–H on a sulfoxide), not just an
   undercounted one — a sharper illustration that this isn't just "off by
   one," it silently manufactures a different molecule.

**`ROUNDTRIP_TOPOLOGY_CHANGE` / silent corruption (P0, general, externally
validated).** `crates/chematic-mol/src/mol2000.rs` and `mol3000.rs`'s bond
readers mapped MDL bond-type code `9` (dative/coordinate — the exact
convention RDKit itself uses to write `Bond::BondType.DATIVE` in V3000
molfiles, verified by generating one with RDKit directly and feeding it to
chematic) to `_ => BondOrder::Single`, the same catch-all used for genuinely
unknown codes. Reading an **RDKit-generated** V3000 molfile containing a
single dative Pt–N bond, chematic silently returned `BondOrder::Single` —
no error, no warning, connectivity intact but bond semantics quietly
changed. This is precisely the task's named worst case (§6/§7): a
structurally-valid file is accepted and a wrong molecule comes out the other
end with no signal that anything was lost. This is the strongest single
piece of evidence in this benchmark, because it did not require any
self-authored SMILES to trigger — it reproduces on any V3000 file any other
tool (RDKit, and by informal convention several MDL-derived tools) emits for
a dative bond.

**`CIS_TRANS_LOSS` / `CANONICAL_IDENTITY_COLLISION` — P0 finding, measured,
NOT fixed this pass (§10 explicitly scopes this out).** See §5.

**`PARSER_UNSUPPORTED_STEREO_DESCRIPTOR` — explicit failure, not corruption
(no fix needed, matches project policy).** `chematic-smiles`'s parser
rejects `[Pt@SP1]`/`[Pt@SP2]`/`[Pt@SP3]` (RDKit's own extended square-planar
stereo descriptors) outright: `InvalidBracketAtom { detail: "missing ']'" }`.
This is the *good* outcome per the task's explicit policy ("silent wrong
answer is worse than parse failure") — chematic correctly refuses input it
cannot represent, rather than silently dropping the `@SP1` tag and returning
a plausible-looking but stereo-blind molecule. No fix recommended; noted
here only so the killer-benchmark's cis/trans gap (§5) is understood
precisely: it is not that chematic mishandles stereo-tagged input, it is
that chematic has **no representation for square-planar stereo at all**, at
any layer (parser, `Chirality` enum, canonicalizer).

**Not observed in this corpus:** `PARSER_UNSUPPORTED_METAL` (Pt parsed
fine everywhere), `OXIDATION_STATE_AMBIGUOUS` as a chematic defect
(chematic has no oxidation-state field at all — see §7 — so there is
nothing to be ambiguous; the corpus instead carries oxidation state purely
as external metadata, never derived from formal charge), `CANONICAL_FALSE_SPLIT`
(no case where the same structure written two ways got two different
canonical identities — every same-`identity_group` pair with identical
input converged, as expected since they're literally the same string),
`FINGERPRINT_UNSUPPORTED`/`SMARTS_UNSUPPORTED` (both ran to completion,
deterministically, on every entry — no panics observed; **bit-exact RDKit
parity for metal-containing fingerprints was not evaluated and is not
claimed**, per the task's own explicit "not required" scope note), `3D_GEOMETRY_LOSS`
(3D generation was not attempted at all this pass — see §7).

## 5. Killer benchmark: cisplatin-family (task §6)

| condition | before | after |
|---|---|---|
| (A) all structures parse | yes (18/18) | yes (18/18) |
| (B) round-trip preserves Pt coordination topology | yes (connectivity), **no (dative-bond semantics silently lost)** | yes, dative semantics also preserved |
| (C) cisplatin ≠ transplatin | **NO — same canonical SMILES** | **NO — still the same canonical SMILES (not fixed this pass)** |
| (D) formula/charge/MW correct | charge yes, formula/MW **NO** | **yes, all three** |
| (E) same structure ⇒ same identity | yes | yes |
| (F) unsupported features fail explicitly, not silently | partially — MOL bond-type 9 was **silent** (now fixed); `@SP1` etc. was already explicit | MOL bond-type 9 now explicit/correct; `@SP1` etc. still explicit reject (unchanged, correct) |

**Condition (C) is the one P0 finding this benchmark surfaces and does
*not* fix, by design (task §10 explicitly forbids building a full
square-planar/octahedral stereo engine in one pass).** Root cause, confirmed
by reading the type system directly:
`crates/chematic-core/src/atom.rs`'s `Chirality` enum has exactly three
variants (`None`, `CounterClockwise`, `Clockwise`) — tetrahedral only. There
is no `SquarePlanar`/`Octahedral` variant anywhere in the codebase. Cisplatin
and transplatin, given the only SMILES syntax chematic can parse for a
4-coordinate square-planar Pt(II) complex (no stereo bond, no dative-arrow
stereo marker — none exists), are graph-isomorphic and therefore
correctly-but-unhelpfully collapse to one canonical identity:
`[Pt](<-N)(Cl)(Cl)<-N` for both.

**This is not chematic-specific.** RDKit, given the *same plain SMILES*
(no `@SP1`/`@SP2`/`@SP3` tag), also cannot distinguish them — confirmed
directly: RDKit's canonical SMILES for both is
`[NH3]->[Pt](<-[NH3])([Cl])[Cl]`, identical. The difference is that RDKit
*could* distinguish them if the input carried an explicit `@SP1`/`@SP2`/
`@SP3` tag — confirmed directly (`scripts/platinum_rdkit_oracle.py`'s
secondary check): feeding RDKit `N[Pt@SP1](N)(Cl)Cl` /
`N[Pt@SP2](N)(Cl)Cl` / `N[Pt@SP3](N)(Cl)Cl` produces 3 distinct canonical
SMILES. Chematic has no equivalent opt-in path at all — the parser rejects
that syntax outright (§4's `PARSER_UNSUPPORTED_STEREO_DESCRIPTOR`, a safe,
explicit failure, not a silent one).

**Recommended next scope** (design investigation only, not implementation —
see task §10, §19): whether chematic needs its own
`@SP1`/`@SP2`/`@SP3`-equivalent, and whether it can reuse the existing
`Chirality` enum (adding variants) or needs a parallel mechanism, is a real
open design question that touches the SMILES parser, the `Atom`/`Chirality`
type, canonical ranking, and the writer. It was not investigated further
this pass, in line with the task's explicit "cisplatin/transplatin identity
survival is enough scope for one pass" guidance (§10).

## 6. What was fixed (P0/P1, general — see §16 generalization gate)

Two independent root causes, found by this benchmark, both general (not
Pt-specific — see the generalization-gate check below). **They are shipped
as two separate PRs, not one** (see §12 below) — a periodic-table data-gap
fix and a coordination-bond-semantics fix have different blast radii and
different review needs, even though both happened to be found by the same
corpus.

**Fix A — periodic-table mass data (its own PR).**
`crates/chematic-chem/src/descriptors.rs`: `avg_mass`/`mono_mass` rewritten
from a ~24-element partial `match` (silent `atomic_number as f64` fallback
for everything else) to complete 118-element static tables. The 94
previously-*uncovered* elements are sourced from RDKit's own
`PeriodicTable::getAtomicWeight`/`GetMostCommonIsotopeMass` (2026.03.3 — the
same oracle already used throughout this repo's validation suite). **The
~24 previously-covered elements keep their pre-existing values, not
RDKit's** — checked element-by-element, not assumed: two of the 24 differ
between the old hardcoded value and RDKit's value (B: 10.811 vs 10.812,
S: 32.065 vs 32.067 — both sub-0.002 rounding, immaterial either way) and
one differs meaningfully (**Se: pre-existing 78.971 vs RDKit's 78.96** —
RDKit 2026.03.3 ships the pre-2013 IUPAC standard atomic weight for
selenium specifically; the current IUPAC value, already present in
chematic before this benchmark, is kept rather than silently downgraded).
Every other of the 94 gap-filled elements is unaffected by this
distinction (this benchmark did not check whether RDKit's table matches
current IUPAC values element-by-element for all 94; only the pre-existing
24 were cross-checked, since only those risked a silent regression).

**Fix B — dative-bond coordination-chemistry semantics (one combined PR;
these three files share the same corpus evidence and the same test
fallout, see §7's regression note).**

1. **`crates/chematic-core/src/valence.rs`**: `valence_inferred_hcount` now
   treats the donor side of a `BondOrder::Dative` bond (`bond.atom1 == idx`,
   matching `BondOrder::Dative`'s own documented donor→acceptor convention)
   as contributing 0 to the valence sum, not `order_int()`'s 1. Both
   `implicit_hcount` (formula/mass) and `chematic-smiles`'s canonical-writer
   atom invariant (`initial_invariant`, which calls `implicit_hcount`)
   inherit the fix automatically — no second fix site needed there. The
   acceptor side, and any atom that is itself in the organic subset while
   accepting a dative bond (e.g. a hypothetical boron acceptor), is
   deliberately left unchanged — out of scope, since no acceptor in this
   corpus (Pt, Fe, Co) is organic-subset, so nothing here depends on that
   path; noted as a residual limitation, not silently absorbed into the "P0
   fixed" claim. **A separate, independent, `bond_order_sum`-driven function
   (`validate_valence`) still does NOT get this fix** and can produce a
   false-positive `ValenceError` on a *bracketed* dative donor whose only
   listed normal valence is exactly met by its explicit H count (concretely:
   `[OH2]->[Pt]` is flagged invalid, even though it is valid water-donor
   chemistry and `implicit_hcount` agrees) — deliberately not widened,
   since `bond_order_sum` is public and used by `chematic-cip`/tautomer/
   `chematic-ff`, a much larger blast radius than this benchmark measured;
   see §10.
2. **`crates/chematic-mol/src/mol2000.rs` + `mol3000.rs`** (read): MDL bond
   type `9` now maps to `BondOrder::Dative` instead of falling through to
   the "unknown code" `Single` default, in both V2000 and V3000 readers.
   Atom order (`atom1`/`atom2`) is preserved exactly as read, which the
   RDKit-generated test file confirms is already donor→acceptor order.
   Both readers now carry a **committed regression test** built from the
   exact RDKit-generated molblock text that surfaced this finding
   (`mol3000.rs::test_dative_bond_type_9_round_trips`,
   `mol2000.rs::test_parse_bond_type_9_is_dative`) — this is not left as a
   claim in this document alone; `cargo test -p chematic-mol` fails if bond
   type 9 ever silently regresses back to `Single`.
3. **`crates/chematic-mol/src/mol3000.rs`** (write, both writer functions):
   `BondOrder::Dative` now writes MDL bond type `9` instead of collapsing to
   `1` (plain single). **V2000's writer is intentionally left unchanged**
   (still collapses `Dative` to `1`) — RDKit itself cannot write a dative
   bond in V2000 either (confirmed: it auto-upgrades to V3000 the moment a
   dative bond is present), so full dative round-tripping through
   chematic's own MOL writer requires V3000, matching RDKit's own
   constraint rather than a chematic-specific gap.

Fix B's `valence.rs` change also carries 5 new unit tests in
`crates/chematic-core/src/valence.rs` (donor-side N/O keep full valence;
acceptor-side is deliberately unaffected; the fix generalizes across Fe/Co/
Pd/Ru acceptors, not just Pt; and the known `validate_valence` divergence
above, pinned rather than left as prose).

**Explicitly not touched:** `BondOrder`/`Chirality` type design (§8/§10 —
`Dative`/`Zero`/`Quadruple` already existed before this benchmark and needed
no new variants for anything fixed here), CML/MOL2/CJSON/MolJSON writers
(all still collapse `Dative` to a plain bond on write — a pre-existing,
documented, consistent simplification across those formats, not something
this pass changed or regressed), fingerprints, SMARTS semantics, 3D geometry,
oxidation-state inference, square-planar/octahedral stereo.

## 7. Before/after (task §14 — every number, not just the good ones)

Full data: `validation/results/platinum_baseline_chematic.jsonl` (before) vs
`platinum_after_fix_chematic.jsonl` (after), same 18-row corpus, same
harness, only the 4 production files above changed in between.

| metric | before | after |
|---|---|---|
| SMILES parse success | 18/18 | 18/18 (unchanged) |
| formula matches expected | 1/18 | **18/18** |
| net charge matches expected | 18/18 | 18/18 (unchanged) |
| Pt coordination number matches expected | 16/18 | 16/18 (unchanged — the 2 non-Pt rows, see §3 note) |
| exact mass matches expected (within 0.01 Da, 12 entries with an independently-sourced expected value) | 0/12 | **12/12** |
| MOL V3000 round-trip: formula preserved | 18/18 | 18/18 (unchanged) |
| MOL V3000 round-trip: dative bond count preserved | 1/18 | **18/18** |
| ECFP4 panic/determinism smoke check | 18/18 pass | 18/18 pass (unchanged) |
| SMARTS Pt–Cl smoke check | 18/18 pass | 18/18 pass (unchanged) |
| cisplatin ≠ transplatin (killer condition C) | fail | **fail (unchanged, by design — see §5)** |

**Per-fix attribution (measured, not just reasoned about):** the "before"
baseline above stashed all of Fix A + Fix B at once, which conflates their
individual effects. Two additional isolated runs — Fix B applied alone
(Fix A/`descriptors.rs` reverted) and Fix A applied alone (Fix B/`valence.rs`
reverted) — pin down exactly which fix does what
(`validation/results/platinum_attribution_valence_only.jsonl` /
`platinum_attribution_mass_only.jsonl`):

| state | formula matches expected | cisplatin exact mass (expected 298.9556) |
|---|---|---|
| both reverted (baseline) | 1/18 | 179.9753 |
| **Fix B only** (valence.rs; Fix A still reverted) | **18/18** | 181.9910 |
| **Fix A only** (descriptors.rs; Fix B still reverted) | 1/18 | 296.9401 |
| both applied (after) | 18/18 | 298.9556 |

This shows the split is real but **mass is not cleanly attributable to
Fix A alone** — Fix A alone gets cisplatin to 296.9401, not the full
298.9556, because `molecular_weight`/`exact_mass` both call
`implicit_hcount` for their H contribution too, so mass inherits Fix B's
correction as well (the residual 2.0155 Da gap with Fix B still reverted is
almost exactly 2×1.008 Da, i.e. the 2 still-undercounted ammine hydrogens).
**Formula, by contrast, is cleanly attributable to Fix B alone** — Fix A
(the mass table) has zero effect on formula counts, as expected, since
formula only counts atoms/implicit-H, never mass. This is reported plainly
rather than smoothing the two fixes into one "before/after" number that
would overstate Fix A's isolated contribution to mass correctness.

Exact mass, before → after (Da; full 18-row table in the results files):

| compound | expected | before | after |
|---|---|---|---|
| cisplatin | 298.9556 | 179.9753 | 298.9558 |
| carboplatin | 371.0445 | 252.0641 | 371.0446 |
| oxaliplatin | 397.0602 | 278.0798 | 397.0602 |
| satraplatin | 499.0605 | 380.0802 | 499.0607 |
| iproplatin | 417.0550 | 298.0747 | 417.0552 |
| tetraplatin | 448.9548 | 329.9758 | 448.9562 |
| dicycloplatin | 515.0868 | 396.1063 | 515.0868 |

(All "after" values are within 0.0002 Da of expected -- floating-point
summation-order noise, not a remaining defect; every one is well inside
the ≤0.01 Da tolerance the summary table above uses.)

**Regression check:** `cargo test --workspace` run twice — once immediately
after the 4 production fixes (before touching any test), which surfaced one
real, expected fallout (below), and once more after fixing it, both times
green. `cargo clippy -p chematic-core -p chematic-chem -p chematic-mol -p
chematic-smiles -- -D warnings` clean. `cargo fmt --check` clean.

**Fallout found and fixed (not silently absorbed):**
`crates/chematic-smiles/tests/dative_bond_direction.rs`'s
`canonical_writer_flips_the_arrow_when_the_acceptor_is_written_first` test
pinned an exact canonical-SMILES string
(`canonical_smiles(&n_to_fe()) == "[Fe]<-N"`) that depended on N's *old,
wrong* implicit-H invariant (2) to make Fe rank ahead of N in the canonical
DFS start. After the fix, N's invariant correctly becomes 3, which flips the
ranking and the canonical output to `"N->[Fe]"` — a different string, but an
equally valid, still round-trip-correct representation of the same
donor→acceptor pair (verified: re-parsing recovers the identical
donor/acceptor identity either way). The stale assertion was updated to the
new correct value with an explanation tying it to this exact fix, and a
**new** test case (`O donor → Fe acceptor`, which still ranks
acceptor-first post-fix) was added alongside it so the arrow-flip-on-
acceptor-first code path — the actual thing the original test name
describes — stays covered rather than silently losing coverage. This is
the exact risk flagged in advance for this kind of fix (an implicit-H
change moving canonical output elsewhere in the corpus); it is called out
here rather than mentioned only in a commit message.

## 8. RDKit / Open Babel comparison (task §12 — never "RDKit disagrees ⇒
chematic is wrong")

**RDKit 2026.03.3**, used as the primary oracle throughout (already the
project's standard, per `validation/README.md`):
- Supports the same `->`/`<-` dative-bond SMILES syntax as chematic; used
  directly as a second, independent implementation to confirm the
  implicit-H bug (§4) rather than relying on chematic's own before/after
  diff alone.
- Also cannot distinguish cisplatin/transplatin on plain (non-`@SP`-tagged)
  SMILES (§5) — this is not a chematic-vs-RDKit gap on the corpus as
  written, it is a shared limitation of the *input encoding chosen*.
- *Can* distinguish square-planar isomers given explicit `@SP1`/`@SP2`/
  `@SP3` tags, which chematic cannot parse at all (explicit reject, not
  silent corruption).
- Confirmed the two internally-inconsistent PubChem records independently
  (§2) — RDKit parsing PubChem's own `tetraplatin` SMILES gives a *different*
  (and self-consistent) exact mass than PubChem's own stated `ExactMass`
  field for the same record.
- Writes MDL bond type 9 for dative bonds, and *only* in V3000 (auto-
  upgrades from V2000) — chematic's fix (§6) matches this exactly rather
  than inventing its own convention.

**Row-level chematic-vs-RDKit comparison**
(`scripts/platinum_compare_chematic_rdkit.py`, comparing
`platinum_after_fix_chematic.jsonl` and `platinum_rdkit_oracle.jsonl`
row-by-row, not just summary counts): chematic evaluates
`formula_matches_expected` using its own
`chematic_chem::formula::parse_formula`; the RDKit oracle script evaluates
the identical check using a completely independent regex-based formula
parser. **18/18 rows agree** — every row where chematic says the corpus's
`formula_expected` is right, RDKit's independent parser says so too, with
zero disagreements. This means the corpus's expected formulas are
doubly-confirmed, not merely self-consistent with chematic's own parsing of
its own expectation string. For the 6 corpus rows with no
independently-sourced `exact_mass_expected` (the DMSO/methyl-reference/
generalization-gate rows — the rows where the H-count fix is least
externally anchored), chematic's computed mass vs RDKit's computed mass on
the identical input SMILES agrees to **within 0.0019 Da on all 6** — this
is the actual comparison the corpus's `source.note` on those 6 rows
describes, now implemented rather than left as an unfulfilled promise.

**Open Babel 3.1.1** (`openbabel-wheel` PyPI package; no system-level
`obabel` binary was present in this environment, so the Python bindings
were installed and used instead — one bounded `pip install` attempt, as
planned, succeeded here unlike the earlier `obabel` CLI check):
- **Does not support `->`/`<-` dative-bond SMILES syntax at all** — errors
  immediately and explicitly: `*** Open Babel Error in ParseSimple: SMILES
  string contains a character '<' which is invalid`. This is a third data
  point (after chematic's pre-fix silent MOL corruption and RDKit's correct
  handling) on how differently three real tools treat the exact same
  coordination-bond representation choice — underscoring why §8's "design
  investigation before implementation" was the right call, not overcaution.
- On the plain-covalent form (`N[Pt](N)(Cl)Cl`), computes `Cl2H4N2Pt`
  (NH2-style, consistent with RDKit's identical covalent-bond
  interpretation of the same string — not a discrepancy, both tools treat
  a *plain* single bond as consuming the donor's valence, which is the
  documented, correct behavior for that bond type).
- Produces a stable, order-independent canonical SMILES for the
  plain-covalent form (`N[Pt](Cl)(Cl)N` regardless of input atom order) —
  no round-trip/canonicalization defect observed on this limited check.
- MOL/SDF round-trip and cis/trans handling were not evaluated for Open
  Babel beyond the above (time-boxed, per the task's own guidance not to
  chase this indefinitely).

**No case in this benchmark saw RDKit and Open Babel disagree with each
other while chematic's post-fix output matched neither** — the one
place all three genuinely differ (dative-bond SMILES syntax support:
chematic yes, RDKit yes, Open Babel no) is a syntax-support gap, not a
chemistry disagreement, and is reported as such rather than picked as a
"which tool is right" contest.

## 9. Regression corpus (task §13)

`validation/platinum/pt_corpus.jsonl` **is** the permanent regression
corpus — same file used for baseline and after-fix measurement, not a
separate curated subset, so there is no risk of the "regression corpus"
silently drifting from what was actually measured. All 18 entries are
either derived from PubChem's public, redistributable numeric property
fields (with connectivity independently authored, not copied — see §2) or
fully original (the 2 non-Pt generalization-check entries, the DMSO/methyl
reference compounds). Nothing CSD-derived or otherwise non-redistributable
is present, so the whole file is committed directly — no download/build
script indirection was needed for this pass (§13's "don't commit
non-redistributable data" constraint does not currently apply to anything in
this corpus).

## 10. Remaining scientific/engineering limitations (reported honestly, not
buried)

- **Corpus size**: 18, not 100–500 (§1 above — the primary constraint was
  verification time against a real-source bar, not code).
- **Cis/trans (square-planar stereo) is unrepresentable**, confirmed P0,
  not fixed this pass, by design (§5, §10).
- **Octahedral Pt(IV) stereochemistry** (axial/equatorial relationships,
  which axial ligand is "up" vs "down") is likewise entirely
  unrepresented — not even investigated this pass; strictly a superset of
  the square-planar gap above.
- **3D generation was not attempted** for any Pt complex in this corpus.
  Per the task's explicit §11/§15 guidance, "3D generation unsupported" is
  an acceptable, honest conclusion rather than running chematic's organic-
  molecule 3D pipeline on a metal complex and calling whatever coordinates
  come out "supported." (The isolation requirement for this benchmark also
  explicitly forbids touching `crates/chematic-3d/` at all, which was
  respected throughout — confirmed via `git diff --stat` before every
  commit, see PR descriptions.)
- **Oxidation state has no representation in chematic at all** — not a
  bug, a deliberate non-finding: `pt_oxidation_state` in the corpus is
  pure external metadata (drawn from established chemistry, e.g. "Pt(II)"),
  never computed or asserted by chematic, and no chematic code path claims
  to infer it. This matches the task's §9 requirement exactly
  ("`oxidation_state: Unknown` is safer than a confident wrong guess") —
  the safe answer here is simply that the feature doesn't exist, not that
  it exists and returns `Unknown` for these compounds specifically.
- **Fingerprint/SMARTS coverage was only smoke-tested** (panic +
  determinism), not compared bit-for-bit against RDKit, per the task's own
  explicit scope note that bit-exact metal-containing fingerprint parity is
  not required this pass.
- **The acceptor side of a dative bond, when the acceptor is itself an
  organic-subset element** (e.g. a hypothetical boron/nitrogen adduct), is
  untouched by the valence.rs fix — no compound in this corpus exercises
  that path (Pt/Fe/Co are never organic-subset), so it is explicitly
  unverified, not silently assumed fixed.
- **CML/MOL2/CJSON/MolJSON writers** still collapse `Dative` bonds to a
  plain bond order on write (pre-existing behavior, unchanged this pass) —
  only SMILES and MOL V3000 round-trip dative bonds losslessly through
  chematic today.
- **`validate_valence` can still false-positive on a bracketed dative
  donor.** Confirmed directly (pinned as
  `test_known_divergence_bracketed_dative_donor_can_still_false_positive_in_validate_valence`
  in `valence.rs`): `[OH2]->[Pt]` is valid water-donor chemistry —
  `implicit_hcount` correctly agrees (2H) — but `validate_valence` still
  reports a `ValenceError`, because it is built on `bond_order_sum`, a
  *separate* function this fix deliberately did not touch (it is public and
  shared with `chematic-cip`, tautomer enumeration, and `chematic-ff`'s
  bonded-term assembly — widening it was a materially larger blast radius
  than this benchmark measured or scoped). Nitrogen happens to escape this
  in practice only because N's valence list `[3, 5]` has a second tier that
  absorbs the extra count — an element-specific coincidence, not a general
  exemption. This corpus's own `valence_errors` fields are all empty only
  because every donor atom in `pt_corpus.jsonl` is written un-bracketed
  (`N`, `O`, `S`, never `[NH3]`/`[OH2]`), which sidesteps this path
  entirely — worth knowing before relying on `validate_valence` as a
  sanitizer for bracketed-atom coordination SMILES specifically.

## 11. GO / NO-GO and recommended implementation order

**GO**, with the P0 defects found in this pass (mass table, dative-bond
implicit-H, MOL bond-type-9 silent corruption) fixed, and cis/trans
square-planar stereo explicitly flagged as the next real gap rather than
claimed solved.

Recommended order for any follow-on work (highest measured value first):
1. (Done, this pass) Fix the mass-table and implicit-H P0/P1 defects.
2. Resolve the corpus gaps named in §1 (heptaplatin, an oxaliplatin-Pt(IV)
   prodrug, phenanthriplatin/pyriplatin) via literature/ChemRxiv search
   specifically, now that the PubChem-only approach has been exhausted.
3. Design investigation (not implementation) for square-planar/octahedral
   stereo representation — the single highest-value remaining gap measured
   by this benchmark (task §10 explicitly scopes this as a separate,
   later piece of work).
4. Only after (3) has a design: implement the minimal cis/trans-preserving
   mechanism, re-run this exact benchmark, and report before/after again —
   do not fold it into this PR.

## 12. Delivery / PR split

Per the task's explicit "no giant single PR" instruction, this work ships
as three PRs, not one:

- **PR 1 (diagnostic only, no production code changed):** `pt_corpus.jsonl`,
  `crates/chematic-mol/examples/platinum_benchmark.rs`,
  `scripts/platinum_rdkit_oracle.py`,
  `scripts/platinum_compare_chematic_rdkit.py`, this document, and every
  `validation/results/platinum_*` file (baseline, RDKit oracle, and both
  attribution runs — all generated *before* PR 2/3's fixes are applied, so
  this PR's own numbers are reproducible against unmodified `main`).
  `crates/chematic-mol/Cargo.toml`'s two new dev-dependencies
  (`chematic-fp`, `chematic-smarts`, needed only by the new example) are
  the only non-`validation`/`scripts` diff in this PR.
- **PR 2 (Fix B — dative-bond coordination-chemistry semantics):**
  `crates/chematic-core/src/valence.rs` (+ 5 new unit tests),
  `crates/chematic-mol/src/mol2000.rs` (+ 1 new unit test),
  `crates/chematic-mol/src/mol3000.rs` (+ 1 new unit test),
  `crates/chematic-smiles/tests/dative_bond_direction.rs` (1 stale
  assertion updated, 1 new test added — direct, disclosed fallout of the
  `valence.rs` change, see §7). `validation/results/
  platinum_after_fix_chematic*.jsonl` and the two attribution-run result
  files are included here as this PR's own before/after evidence.
- **PR 3 (Fix A — periodic-table mass data):**
  `crates/chematic-chem/src/descriptors.rs` only. Independent of PR 2 —
  reviewable, and revertable, on its own.

None of these three PRs touches anything under `crates/chematic-3d/`
(checked via `git diff --stat` against `main` before each PR was opened,
per this benchmark's isolation requirement) or overlaps the concurrently-
open stereo-aware 3D distance-geometry work. All three are opened Ready,
none merged/tagged/published without separate explicit human approval.
