# Diag/aromaticity-rdkit-parity — chematic vs RDKit aromaticity/kekulization diagnosis

**Status:** diagnosis only. No production code changed. Not merged.

**Branch:** `diag/aromaticity-rdkit-parity`, forked from `main` at
`659baca221f71f135ce0e1780e71245d8770f132`.

**Files touched:**
- `crates/chematic-perception/examples/aromaticity_rdkit_parity_dump.rs` (new) — dumps
  the frozen fixture corpus's chematic-side results as JSONL.
- `scripts/aromaticity_rdkit_parity_diagnosis.py` (new) — RDKit-oracle cross-check,
  fail-closed.
- `validation/results/aromaticity_rdkit_parity_fixture_dump.jsonl` (new) — frozen
  chematic-side dump, reproducible via the example above.
- `validation/results/aromaticity_rdkit_parity_diagnosis_summary.json` (new) — the
  diagnosis script's machine-readable output.
- `docs/aromaticity_rdkit_parity_rfc.md` (this file, new).

**Explicitly out of scope / not touched:** anything under `crates/*/src/**` (no
production code modified — every finding below was confirmed via scratch `examples/`
binaries that were written, run, and then deleted before this PR; only the diagnostic
example above is kept); `feat/io-mrv`, `feat/io-tdt`, `feat/io-smiles-supplier-writer`,
`fix/smiles-bracket-implicit-h`, `diag/stereo-reader-integration-boundary`,
`feat/stereo2d-local-parity` (other agents' active branches — not read, not rebased
onto, not merged from).

**Deliverables:** a deliberately-constructed (not randomly sampled) 40-molecule SMILES
fixture corpus; a Rust dump example; a fail-closed Python diagnosis script cross-
checking against a pinned RDKit; a frozen JSON summary; this RFC.

**Done condition:** all 40 fixtures classify into a named, evidence-checked bucket (no
silent drops), the diagnosis script exits 0 on the frozen baseline, its fail-closed
machinery is verified by hand-injecting failure modes (self-tests plus live
duplicate-ID/missing-ID/baseline-drift injection, see §5), and findings are written up
here with root causes and a recommended (not implemented) fix plan.

## 0. Headline, up front

Read literally, the numbers look almost too good: **40/40 fixtures agree with RDKit on
every atom's aromatic flag; 38/40 also agree bond-by-bond.** That is the honest
top-line number, and it is *not* the main finding — six of those atom-level agreements
are not verified computations. They are chematic's own `apply_aromaticity`/
`apply_aromaticity_ex` silently preserving a **pre-existing, never-independently-
confirmed** flag from the SMILES parser, on molecules where chematic's own Kekulé-
structure algorithm (`chematic_core::kekulize`) **hard-fails** to produce any Kekulé
form at all. The real, actionable finding of this diagnosis is that failure — not a
flag-parity gap — plus the masking mechanism that hides it. Section 1 covers both;
sections 2–4 cover three smaller, independent findings (one already-known false
negative independently reproduced, one stale doc claim, one confirmed-still-fixed
regression guard).

## 1. `kekulize()` hard-fails for 6 of 40 fixtures — and the failure is invisible at the flag level

### 1a. The failure, by fixture

| Fixture | SMILES | `kekulize()` result |
|---|---|---|
| `tropylium_cation` | `c1ccc[cH+]cc1` | `Err`: "atom 6 (C) cannot be assigned a double bond" |
| `imidazolium` | `c1c[nH+]c[nH]1` | `Err`: "atom 3 (C) cannot be assigned a double bond" |
| `pyridinium` | `c1cc[nH+]cc1` | `Err`: "atom 4 (C) cannot be assigned a double bond" |
| `pyrylium` | `c1cc[o+]cc1` | `Err`: "atom 4 (C) cannot be assigned a double bond" |
| `tellurophene` | `c1cc[te]c1` | `Err`: "atom 4 (C) cannot be assigned a double bond" |
| `phosphole` | `c1cc[pH]c1` | `Err`: "atom 4 (C) cannot be assigned a double bond" |

All six parse successfully (chematic's SMILES parser accepts the aromatic notation with
no complaint) and all six are accepted as valid aromatic input by RDKit's sanitizer
(`Chem.MolFromSmiles` + default `sanitize=True`, `rdkit==2026.03.3`). The failure is
specific to `chematic_core::kekulize`'s matching step, and traces to exactly two root
causes in `atom_must_be_matched` (`crates/chematic-core/src/kekulization.rs:571-611`),
independently confirmed for each fixture:

**Root cause A — charge-blind lone-pair-donor rules (`pyridinium`, `pyrylium`,
`imidazolium`).** The O/S/Se rule (line 575, `8 | 16 | 34 => false`) and the
N-with-H rule (line 580, `7 if matches!(atom.hydrogen_count, Some(h) if h > 0) => false`)
both fire unconditionally, regardless of formal charge. A *neutral* pyrrole-type N-H
correctly donates its lone pair and needs no double bond (`false`, no match required).
But a *protonated* ring heteroatom — pyridinium's `[nH+]`, pyrylium's `[o+]` — has that
same H-count/element shape while chemically needing a double bond exactly like neutral
pyridine's bare N does (the added proton occupies what would be the donated lone pair;
the ring still needs its normal alternating pattern). The current rule exempts these
charged atoms from matching anyway, removing one atom from the required-match set and
leaving an **odd** number of atoms that must pair up two at a time — combinatorially
impossible, hence the `Err`. Imidazolium has two N (`[nH+]` and `[nH]`); only the
charged one is mis-classified, but that alone breaks the parity.

**Root cause B — no acceptor-pattern rule for a cationic empty-p-orbital atom
(`tropylium_cation`).** `atom_must_be_matched` has an explicit rule for an *anionic*
ring atom donating its lone pair (line 598, `_ if atom.charge < 0 => false` — this is
what makes `cyclopentadienyl_anion` kekulize correctly, confirmed in the fixture dump).
There is no symmetric rule for a *cationic*, empty-p-orbital atom that should be
excluded from matching the same way (an electron acceptor, contributing 0π, needing no
partner — the same shape aromatic boron already gets, explicitly, at line 578).
Tropylium's `[cH+]` falls through to the generic `_ => true` catch-all, forcing all 7
ring atoms into the must-match set — an odd count, again combinatorially impossible.

**Root cause C — Te(52) missing from the donor list, inconsistent with its own sibling
element (`tellurophene`).** Se(34) is explicitly listed as a lone-pair donor
(line 575); Te(52) is not, despite `AromaticityAlgorithm::RdkitLike`'s own doc comment
(`crates/chematic-perception/src/aromaticity.rs:37-39`) claiming both are supported.
This is a plain cross-module inconsistency between `chematic-core`'s kekulization and
`chematic-perception`'s Hückel model — confirmed by the fact that selenophene's
`kekulize()` *succeeds* while tellurophene's, on an otherwise-identical ring shape,
does not.

**Root cause D — P entirely unhandled, already documented (`phosphole`).** P(15) has no
case in `atom_must_be_matched` (falls to `_ => true`) and no case in
`ring_pi_electrons` either (`crates/chematic-perception/src/aromaticity.rs:976`, falls
to `_ => return None`, "Unsupported element"). Unlike A–C, this one is *not* new: the
`RdkitLike` doc comment already states "P-containing aromatic rings are NOT supported
in this mode (separate sprint)." What *is* new here is the concrete confirmation that
the gap also breaks `kekulize()` outright (not just the Hückel flag), and that RDKit's
own default model — verified live — does treat phosphole as aromatic
(`c1cc[pH]c1` round-trips through `Chem.MolFromSmiles` with `GetIsAromatic()==True` on
every atom).

### 1b. The masking mechanism: why none of this shows up as a flag mismatch

`build_molecule_from_model` (`crates/chematic-perception/src/aromaticity.rs:351-360`)
rebuilds atoms and bonds from a freshly-computed `AromaticityModel`:

```rust
for (idx, atom) in mol.atoms() {
    let mut a = atom.clone();
    if model.is_atom_aromatic(idx) {
        a.aromatic = true;
    }
    builder.add_atom(a);
}
for (bidx, bond) in mol.bonds() {
    let order = if model.is_bond_aromatic(bidx) {
        BondOrder::Aromatic
    } else {
        bond.order
    };
    ...
}
```

Both loops can only **promote** a flag to aromatic when the model confirms it — neither
loop ever **demotes** an already-`true` atom flag (or an already-`Aromatic` bond order)
when the model disagrees. This is technically consistent with `apply_aromaticity`'s own
docstring ("non-aromatic atoms and bonds are unchanged"), but no existing test locks in
this *specific* interaction — every existing `test_apply_aromaticity_*` test starts from
a molecule whose atoms are `aromatic: false` to begin with (Kekulé input), never from an
already-aromatic-flagged molecule the model then disagrees with. This diagnosis is, as
far as could be determined by searching `aromaticity.rs`'s test module, the first case
that actually exercises the "model says no, but the input already said yes" path. Framed
precisely: this is a genuine, previously-unexercised **gap with a surprising
consequence**, not a flatly-documented design guarantee — and it directly explains why
all six kekulization failures above are invisible in a plain "does `atom.aromatic` match
RDKit" check:

- **`tellurophene`, `phosphole`, `tropylium_cation`, `imidazolium`, `pyridinium`,
  `pyrylium`** (bucket `kekulize_fails_atom_bond_flags_survive_coincidentally`, 6
  fixtures): `kekulize()` fails, so the molecule handed to `apply_aromaticity` still
  carries its original `BondOrder::Aromatic` bonds untouched. The Hückel model
  independently computes **zero** aromatic atoms for all six (confirmed via
  `assign_aromaticity_ex` called directly, bypassing the rebuild step — see
  `huckel_model_aromatic_atom_count` in the fixture dump, which reads `0` for every one
  of these six while the final, rebuilt molecule's atom AND bond flags both read
  `true`/`Aromatic`, matching RDKit only by coincidence). Nothing was actually verified;
  the pre-existing flag simply survived.
- **`selenophene`, `azulene`** (bucket
  `kekulize_succeeds_model_disagrees_atom_bond_flags_inconsistent`, 2 fixtures):
  `kekulize()` *succeeds* this time, so the bond loop's fallback (`bond.order`, now
  genuinely Single/Double) correctly reflects the real Kekulé structure — but the atom
  loop's fallback still preserves the stale `true` flag. The result is an **internally
  self-contradictory molecule**: every ring atom flagged `aromatic: true`, every one of
  its ring bonds flagged plain `Single`/`Double` (not `Aromatic`). Bond flags then
  correctly (for the wrong reason: real Kekulé bonds, not real confirmation) disagree
  with RDKit, which is why these two are the only 2 of 40 fixtures with any bond-level
  mismatch at all. Selenophene's root cause is the same "unsupported under default
  Huckel" gap as tellurophene/phosphole (Se(34) *is* handled, but only under
  `AromaticityAlgorithm::RdkitLike`, confirmed: switching to `RdkitLike` gives a fully
  self-consistent, RDKit-matching result for selenophene, atoms AND bonds both). Azulene's
  root cause is unrelated to any unsupported element — see §2.

### 1c. Confirmed downstream impact (traced, not assumed)

Two real production consumers already anticipate `kekulize()` can fail and degrade
gracefully rather than panic or error out — which is precisely why this gap produces no
visible signal to a caller:

- `crates/chematic-depict/src/svg.rs:238-240`, the 2D-depiction renderer's own comment:
  *"kekulize failure silently falls back to aromatic."* `RenderOptions.kekulize: bool`
  (default `false`) is meant to render alternating single/double bonds instead of the
  default aromatic-circle style. **Empirically confirmed**: rendering `pyridine`
  (`c1ccncc1`, unaffected) with `kekulize: true` vs `false` produces two **different**
  SVGs, as expected. Rendering `pyridinium` (`c1cc[nH+]cc1`, one of the six above) with
  `kekulize: true` vs `false` produces **byte-identical** SVGs — the caller's explicit
  request is silently ignored, with no error, no warning, and no way to detect it short
  of diffing the output.
- `crates/chematic-cip/src/assign.rs:117`: `assign_cip_accurate_experimental` computes a
  Kekulé clone once per molecule for its MANCUDE-fractional CIP path; on failure it
  falls back to "the plain, pre-Milestone-3B-1 digraph path... rather than failing the
  whole assignment." Not independently re-run this round (would require a chirally
  interesting molecule among the six above, which none of the 40 fixtures are), but the
  same "degrade silently, no error surfaced" shape as the depict case, confirmed by
  reading the function's own doc comment rather than by execution.

No panics were found anywhere in this trace. The impact is specifically: **a caller who
explicitly asks for Kekulé-style output on one of these 6 molecule shapes gets the
aromatic-style output instead, with no signal that their request was not honored.**

## 2. Azulene: independently reproduced, already-documented false negative — with a purine caveat

Azulene (`c1ccc2cccc-2cc1`, all-carbon, no unsupported elements) lands in the same
"inconsistent flags" bucket as selenophene (§1b), but for a **completely different**
root cause: chematic's per-ring Hückel counting (`ring_pi_electrons`, Pass 1) requires
each SSSR ring to *independently* satisfy 4n+2. Azulene's real system is a 10π
perimeter split as 5+7 (odd/odd) across its two fused rings — neither ring's isolated
count is even, so neither can pass Pass 1 alone, and Pass 2's propagation never seeds
because seeding requires an adjacent ring to have already passed Pass 1 (confirmed via
`assign_aromaticity_ex` directly: `aromatic_atom_count() == 0` under **both**
`Huckel` and `RdkitLike`). This is not a new discovery — `docs/aromaticity_a1_rfc.md`
already lists azulene by name as one of the "known false-negative families." This
diagnosis reproduces it **independently, via a different methodology** (direct
atom/bond SMILES comparison against a live RDKit oracle, rather than the SMARTS corpus
diff that originally found it) — the same kind of independent triangulation
`aromaticity_a1_rfc.md` itself describes for its own bridgehead-N finding.

**Caveat, stated explicitly rather than glossed over**: that same RFC lists azulene
*and* purine together as a "known false-negative family." This diagnosis's bare
`purine` fixture (`c1ncc2[nH]cnc2n1`) **passed cleanly** (`matches_rdkit_exact_kekule`,
no drift, atom/bond flags fully agree). This is not a contradiction this diagnosis
resolves — it suggests the purine gap is specific to some substituent pattern (the
aromaticity.rs source comments mention "N-glycosylated purine" elsewhere) rather than
present in the bare parent structure, but that hypothesis is **not verified** here; no
substituted-purine fixture was built this round.

## 3. CLAUDE.md's indolizine "9-ring" claim is stale

CLAUDE.md states: *"The SSSR algorithm can return a large fundamental cycle (e.g. a
9-ring for indolizine) instead of the two smaller component rings."* This diagnosis's
`indolizine` fixture (`c1ccn2ccccc12`) shows `find_sssr` already returns the correct
`[5, 6]` decomposition **directly**, with zero augmentation needed
(`raw_sssr_ring_count == augmented_ring_count == 2`). This is not a fluke of this one
input: `crates/chematic-perception/src/sssr.rs` already has a dedicated regression test,
`test_indolizine_sssr_minimal`, asserting exactly this `[5, 6]` result on the identical
SMILES — and it passes on current `main`. Cross-referencing project memory, this matches
the "SSSR Horton fix" (2026-07-11), which predates this diagnosis. **CLAUDE.md's prose
was not updated to reflect that fix** — a small, low-stakes documentation staleness, not
a behavioral regression.

This does **not** mean `augmented_ring_set` is now dead code: its own source comments
(`crates/chematic-perception/src/aromaticity.rs`, `count_aromatic_rings`'s
`strip_envelope_rings` machinery) describe 3-ring and 4-ring XOR corrections for
"compact PAHs like pyrene" and "coronene-class PAHs" respectively. **Not verified this
round** — no pyrene/coronene fixture was built; this diagnosis can only confirm the
*specific* indolizine example CLAUDE.md cites no longer reproduces, not that the general
correction mechanism is unused elsewhere.

## 4. Confirmed still fixed: the pyrrole/pyridine Kekulé-erasure regression guard

Project memory records a prior fix in `crates/chematic-perception/src/aromaticity.rs`
for a bug where `apply_aromaticity_ex` erased the pyrrole-vs-pyridine Kekulé distinction
(collapsing "needs 1 implicit H" and "needs 0" to the same post-normalization shape).
**Re-verified on current `main`, not assumed**: the two dedicated regression tests,
`test_apply_aromaticity_preserves_pyrrole_nh_implicit_hydrogen` and
`test_apply_aromaticity_does_not_add_h_to_pyridine_type_n`, both still pass
(`cargo test -p chematic-perception --lib`). Method note, stated plainly per this
diagnosis's own fixture corpus limits: this re-verification used the **existing unit
tests** directly, not the new SMILES fixture corpus built for this round — the fixture
corpus compares aromatic *flags*, not implicit-H counts, so it does not itself exercise
this regression. `n_methylpyrrole` was included in the corpus as a related but distinct
check (substituted pyrrole-type N correctly excluded from the "needs H" question
entirely) and passed (`matches_rdkit_exact_kekule`).

## 5. Fixture corpus and methodology

40 fixtures, one SMILES each, hand-built to deliberately cover (not randomly sample)
every mechanism named in this diagnosis's scope: baseline sanity (1), simple negative
controls incl. an antiaromatic-electron-count case (3), classic Hückel heteroaromatics
with mixed pyrrole-/pyridine-type N in the same ring (13), fused/polycyclic aromatics
incl. the indolizine bridgehead-N SSSR case and the azulene non-alternant case (7),
charged aromatics (5), Se/Te/P/B rings (5), and exocyclic-multiple-bond cases incl. a
sulfoxide (5). Full list and rationale: `crates/chematic-perception/examples/
aromaticity_rdkit_parity_dump.rs`'s `FIXTURES` constant.

**Pipeline** (per fixture, both directions tested from the same in-memory molecule —
no second SMILES string needed): `chematic_smiles::parse` → `chematic_core::kekulize`
+ `apply_kekule` (exercises the kekulization/matching algorithm, §1) →
`chematic_perception::{apply_aromaticity, apply_aromaticity_ex(RdkitLike),
apply_aromaticity_rdkit_parity_experimental}` (exercises the Hückel re-perception
algorithm) → per-atom/per-bond JSON dump, plus `find_sssr`/`augmented_ring_set`/
`count_aromatic_rings` and the **raw** `assign_aromaticity_ex(..).aromatic_atom_count()`
(the model's own verdict, independent of the rebuild step — load-bearing for §1b).

**RDKit oracle**: `scripts/aromaticity_rdkit_parity_diagnosis.py`, pinned to
`rdkit==2026.03.3` / commit `8afba32ec539dcb2369bc84549d802aca3f7eb39` (this repo's
`.venv`, same pin used by the stereo2d diagnosis PR). For each fixture's SMILES,
independently re-parses with `Chem.MolFromSmiles` (default `sanitize=True`), reads
per-atom/per-bond `GetIsAromatic()`, ring info, and a `Chem.Kekulize`-derived Kekulé
bond assignment. Verifies index alignment (chematic atom *i* vs RDKit atom *i*)
element-by-element before trusting any other index-based comparison — added
specifically because every other comparison in this script silently assumes that
alignment holds. Kekulé-structure mismatches are not treated as bugs by default: since
Kekulé structures are not unique (benzene and naphthalene both have >1 valid one), a
raw bond-by-bond diff is independently re-validated by rebuilding RDKit's own molecule
with chematic's bond assignment, clearing aromatic flags, and asking RDKit's own
`SanitizeMol` to re-derive + re-canonicalize — only a structure RDKit itself accepts as
chemically valid *and* identical to the original molecule counts as a legitimate
"different but valid" alternate.

**Fail-closed design**, mirroring `scripts/stereo2d_diagnosis.py`'s style: an explicit
bucket whitelist (`EXPECTED_BUCKETS`, not a `startswith("unexpected")` convention), a
frozen `EXPECTED_BUCKET_BY_ID` per-fixture baseline stricter than the whitelist alone
(several buckets describe mutually exclusive outcomes of what would otherwise look like
the same category), fixture-ID-set verification (missing or extra IDs abort), duplicate-
ID detection, and self-tests of the fail-closed logic itself, run before any RDKit call.
**All of this was verified empirically this round, not just written and trusted**:
baseline drift, a duplicate ID, and a missing ID were each injected by hand into a
temporary copy and confirmed to produce a non-zero exit / `FATAL` message before being
reverted; the checked-in script was re-run clean (exit 0) afterward.

**Result: 40/40 fixtures classified into one of 7 named buckets, 0 unexplained, 0
baseline drift.**

| Bucket | Count | Meaning |
|---|---|---|
| `matches_rdkit_exact_kekule` | 20 | Full agreement, atom+bond flags, single valid Kekulé structure |
| `both_correctly_nonaromatic` | 8 | Both agree non-aromatic (negative controls, benzoquinone, cyclopentadienone, sulfoxide, borole, borazine) |
| `kekulize_fails_atom_bond_flags_survive_coincidentally` | 6 | §1: `kekulize()` hard-fails; flags match RDKit only by coincidence |
| `matches_rdkit_exact_kekule_exocyclic_bond_excluded_correctly` | 2 | 2-pyridone, tropone: aromatic despite exocyclic C=O, correctly |
| `kekulize_succeeds_model_disagrees_atom_bond_flags_inconsistent` | 2 | §1b/§2: selenophene, azulene — atom/bond flags mutually inconsistent |
| `matches_rdkit_kekule_valid_alternate` | 1 | benzene: different (both valid) Kekulé structure than RDKit's choice |
| `sssr_bridge_artifact_not_reproduced_docs_stale` | 1 | §3: indolizine |

Every fixture's full per-atom/per-bond evidence (not just the bucket name) is recorded
in `validation/results/aromaticity_rdkit_parity_diagnosis_summary.json`.

**Ring-membership coverage caveat**: `classify()` compares ring *count* (raw SSSR count,
augmented count, and `count_aromatic_rings`) against RDKit's, plus — for indolizine
specifically — each ring's *size*. It does not compare the actual per-ring atom-*sets*
against RDKit's own SSSR decomposition. This is a deliberate scope limit, not an
oversight: SSSR decomposition is non-unique for symmetric/fused systems (the same
"multiple valid answers" problem the Kekulé-structure check in §5 handles explicitly),
so a raw set-equality check would need the same "independently verify as valid, not just
different" guard already built for Kekulé bonds — not done this round.



## 6. Recommended future fix plan (NOT implemented in this PR)

**Mandatory framing for any future implementation**: the fixes below restore intended,
already-specified behavior — they are bug fixes, not new algorithms, and should ship
as direct corrections, not behind a flag. Two things in this diagnosis are genuinely
algorithm-shaped choices, and for those, the **existing algorithm must remain available
and any new behavior must ship as an explicit opt-in, never a silent replacement** —
exactly the precedent `AromaticityAlgorithm::RdkitLike` and
`apply_aromaticity_rdkit_parity_experimental` already set in this codebase. The two
categories should not be conflated:

**Bug fixes (restore intended behavior, ship directly, no flag needed):**
1. `atom_must_be_matched` (`chematic-core/src/kekulization.rs:571-611`): make the O/S/Se
   and N-with-H lone-pair-donor rules charge-aware (a charge that displaces/consumes the
   lone pair, e.g. protonation, should route the atom back to "must match" — matching
   neutral pyridine's bare-N treatment, not neutral pyrrole's N-H treatment); add a
   symmetric acceptor-pattern rule for a cationic, empty-p-orbital ring atom (mirroring
   the existing `charge < 0 => false` anion rule and the existing boron rule); add Te(52)
   to the existing chalcogen-donor list to match `RdkitLike`'s own documented claim.
2. `build_molecule_from_model` (`chematic-perception/src/aromaticity.rs:351-360`): make
   both the atom and bond rebuild loops fully assign the model's verdict
   (`a.aromatic = model.is_atom_aromatic(idx)`, and the equivalent for bond order)
   instead of only ever promoting. This closes the masking mechanism in §1b so that a
   real kekulization/perception failure becomes *visible* (the flag correctly reads
   `false`) instead of silently echoing a stale, unverified input flag. This is the
   single highest-leverage fix in this diagnosis: it doesn't fix root causes A–D, but it
   turns 6 currently-invisible failures into 6 visible ones, which is a precondition for
   anyone noticing them without re-running a diagnosis like this one. **This is not a
   pure strictly-additive improvement — say so plainly, not just "surfaces hidden
   mismatches":** fixing #2 alone will *demote* any aromatic-notation-parsed molecule
   containing Se/Te/P/B (or any other structurally-unsupported case) from
   `atom.aromatic = true` to `false` under the **default** `Huckel` engine, changing
   the output of any current consumer (descriptors, MW, fingerprints, canonical SMILES)
   that today benefits from the accidental, unverified match. Fix #1 does not rescue
   this at the default-engine level for Se/Te — `AromaticityAlgorithm::RdkitLike`
   remains the migration path for those elements, exactly as it is today; P and the
   azulene-shaped case have no migration path at all until #3/#4 below exist. Anyone
   shipping #2 should expect and communicate this as a deliberate default-engine
   behavior change, not a side effect.

**Algorithm-shaped choices (existing behavior stays the default; new behavior is
opt-in):**
3. P support (root cause D, §1) is already correctly scoped in the existing docs as a
   separate, opt-in-shaped sprint (`RdkitLike`'s own comment already reserves this).
   Any implementation should extend `AromaticityAlgorithm` (a new variant, or an
   additional opt-in flag alongside `RdkitLike`) rather than changing `Huckel`'s
   default behavior.
4. A non-alternant fused-ring-system model for azulene-shaped cases (§2) is a genuinely
   different algorithm from the current per-ring-then-propagate Hückel pass, not a bug
   fix to it — RDKit's own real behavior here was not deeply audited this round (no
   RDKit C++ source citations for this specific mechanism; the live-oracle comparison
   alone was used). Any future implementation of this should be a new, explicitly
   opt-in `AromaticityAlgorithm` variant or a separate function, with `Huckel`'s current
   per-ring behavior remaining the unconditional default — matching this project's
   existing precedent of shipping `apply_aromaticity_rdkit_parity_experimental`
   alongside (never instead of) `apply_aromaticity`.

**Order of operations, recommended but not binding**: fix #2 (the masking bug) first,
specifically *because* fixing it will very likely surface additional currently-hidden
mismatches beyond this diagnosis's 40-fixture corpus (any molecule anywhere with a
pre-existing aromatic-notation flag on a ring chematic's model doesn't independently
confirm) — better to find those with #2 fixed and #1's root causes A–D also fixed in the
same pass, than to fix #2 alone and get a wave of newly-visible-but-still-broken
molecules with no accompanying fix.

## 7. What this diagnosis does not do

- Does not modify any file under `crates/*/src/**`.
- Does not fix `atom_must_be_matched`, `build_molecule_from_model`, or any other
  production code identified above.
- Does not implement a non-alternant aromaticity model or extend `AromaticityAlgorithm`.
- Does not update CLAUDE.md's stale indolizine claim (§3) or add a Te(52)/P doc fix —
  flagged here for whoever picks this up next, not fixed in this PR.
- Does not audit RDKit's C++ source for the azulene/non-alternant mechanism (§2) — the
  finding there rests on live oracle comparison only, not source-level confirmation.
- Does not chase the purine substituent-specific caveat in §2 to a conclusion.
- Does not rebase/restack any other agent's branch or PR.
- Is not merged.
