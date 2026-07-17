# RFC: Aromaticity A1 — component-level Hückel solver (design + diagnosis only, no fix in this round)

Status: **A1-0 landed (diagnosis + observability). No production behavior change.**
A1-1 (the actual component-level solver) is designed below but not implemented.

## Background

`docs/rdkit_compat.md`'s "SMARTS-A0" section root-caused the entire remaining SMARTS
match-set residual (56/80,000 comparisons post-SMARTS-R1) to a single mechanism:
`apply_aromaticity("rdkit_like")` over-extends aromaticity from a genuinely aromatic
ring, across a bridgehead-N ring fusion, into an adjacent ring that RDKit does not
consider aromatic. That diagnosis also found this is **not a new defect** — it is the
same `aromatic_context` propagation mechanism already pinned, with 32 known-wrong
molecules, in `crates/chematic-perception/src/aromaticity.rs`'s
`test_known_regressions_from_bridgehead_n_fix` (a pre-existing "known regressions from
fix #2 (bridgehead-N guard removal)" test, whose own comment already recommends
"removing the bypass in favor of proper ring-system candidate enumeration" — essentially
the same direction as A1-1 below). SMARTS-A0 rediscovered and independently confirmed
this mechanism via a completely different methodology (SMARTS corpus diff vs. a hand-
pinned regression test), and additionally localized it to the exact
`ring_pi_electrons` rule interaction (see "Mechanism" below).

**The defect is not SMARTS-specific.** `apply_aromaticity` is the shared perception
engine behind descriptors, canonical SMILES, fingerprints, and InChI, not just
`chematic-smarts`. Fixing it properly means fixing it once, upstream, not hiding it
behind a SMARTS-local view.

## Decision: fix the shared engine, not a SMARTS-side view

Three options were on the table:

1. **Fix `apply_aromaticity`'s propagation** — broad blast radius, fixes every
   consumer at once.
2. **SMARTS-only aromatic view** — scoped to `chematic-smarts`, but makes a molecule's
   SMARTS aromaticity diverge from its own `atom.aromatic` flags, splitting "truth"
   across `c`/aromatic-bond queries, recursive SMARTS, fingerprint atom typing, and
   canonical SMILES aromatic-bond representation.
3. **Accept 99.93% and defer** — the standing decision already in place for the
   pre-existing 32-molecule false-positive corpus and the azulene/purine false negatives.

**Decision: Option 1.** Rejected Option 2 because it would relocate the residual to a
different surface rather than closing it — the defect is upstream of `chematic-smarts`,
so a SMARTS-local patch creates exactly the "two truths" problem that motivated fixing
`[rN]`/`[kN]` for real in SMARTS-R1 rather than papering over it. Rejected Option 3
because this mechanism has *both* a known false-positive family (this one) and known
false-negative families (azulene, purine, `test_known_order_dependent_regressions`) —
fixing the shared root cause has a chance to close both simultaneously, which a
SMARTS-only view could never do.

**Explicitly not doing a local `if` patch.** A direct patch of the shape "if bridgehead
N and exocyclic C=C, suppress propagation" was explicitly rejected before any code was
written. That would fit the ~33 known false-positive molecules and risk making the
already-known false-negative families (azulene, purine) *worse*, since they need *more*
propagation, not less, and share the exact same `aromatic_context` code path. RDKit's
own model does not propagate a boolean flag atom-by-atom; it evaluates ring/fused-ring-
system *candidates* and sums pi-electron contributions per atom, checking 4n+2 for the
whole candidate. The fix direction is separating **candidate discovery** (which atoms
even form a plausible aromatic system) from **electron-count decision** (does that
candidate satisfy Hückel), not strengthening the existing boolean flood-fill.

## Mechanism, traced to the exact rule interaction

**Corrected from this doc's first version** (which described this as a Pass-2
`aromatic_context`-propagation bug): direct measurement via
`trace_ring_pi_electrons`'s intrinsic (empty-context, Pass-1-only) trace, run against
the actual pinned data in `validation/results/aromaticity_a1_0_diagnosis.jsonl`, shows
this is a **Pass-1 over-count** — the offending ring passes 4n+2 using only its own six
atoms' local contributions, with zero cross-ring context borrowed. The first version's
claim that atoms `2,7` needed `AlreadyAromaticContext` was an artifact of reading the
*full-context* trace column instead of the *intrinsic* one; both atoms independently
qualify via `CarbonEndocyclicDouble` with no context at all, since their shared edge
(`2-7`) is itself one of this ring's own bonds, aromatic-order, with no `Double`-only
gate to trip. This matters concretely: a fix scoped as "correct Pass 2's propagation"
would not touch this ring at all.

On the SMARTS-A0 minimal reproducer `C1=Cc2ccccc2C2=NCCCN12` (benzo ring fused via a
bridgehead N to a third, non-aromatic ring, with an exocyclic-to-benzo `C=C` on the
middle ring), the adjacent ring (atoms `0,1,2,7,8,13`), evaluated **with an empty
context**:

- atoms `0,1,2,7`: `CarbonEndocyclicDouble` → 1π each — all four independently, no
  borrowing (verified: intrinsic and full-context traces agree exactly on these four)
- atom `8`: `CarbonExocyclicHeteroatomDouble` → 0π (`C=N` to atom `9`, exocyclic to
  *this specific ring* — confirmed correct against real RDKit, see "atom-8 rule" below)
- atom `13`: `NitrogenBridgeheadOrSubstitutedLonePair` → 2π (bridgehead N rule: all
  three σ-bonds fill its valence, ring-degree < total-degree *relative to this one ring*)

Total: 1+1+1+1+0+2 = **6π → passes 4n+2 (n=1) → the whole ring gets marked aromatic**,
even though ring C (the third ring this same bridgehead N also belongs to) has three
sp3 CH₂ carbons and cannot itself be part of any real delocalized system. The bridgehead-
N rule's 2π credit is correct for indolizine (both rings genuinely aromatic) but wrong
here, because the rule does not check whether the atom's *other* ring is actually a
valid aromatic partner before granting the lone-pair credit to *this* ring's count.

**The atom-8 rule (`CarbonExocyclicHeteroatomDouble`) is confirmed correct, not the fix
target** — blocking check run against real RDKit before any A1-1 code was written (see
"A1-1a" below): tropone, 2-pyridone, and 4-pyranone all have an in-ring carbon whose
only double bond is exocyclic to a heteroatom (a carbonyl carbon), and RDKit marks the
*whole ring including that carbon* aromatic in all three — matching chematic's existing,
already-passing behavior. The atom-8 rule stays as-is; the fix route is the atom-13
(bridgehead N) side, or candidate generation, not this rule.

## A1-0: fused-component characterization (this round, landed)

Per the 5-step spec, diagnosis/observability only:

**1. Fused candidate components**, built over the *augmented* ring list (the same list
`assign_aromaticity_ex`'s Pass 1/Pass 2 actually iterate — not raw SSSR, which misses
the small XOR sub-rings the augmented list adds, e.g. indolizine's 5-ring). Reused
`chematic-perception`'s existing `find_ring_families` union-find (already classifies
Simple/Fused/Spiro/Bridged) rather than writing a second one — extracted as
`find_ring_families_over(mol, rings: &[Vec<AtomIdx>])`, with `find_ring_families`
becoming a thin wrapper (`find_ring_families_over(mol, sssr.rings())`, byte-identical
behavior, zero-diff for existing callers).

**2. Per-atom pi-electron contribution trace.** `trace_ring_pi_electrons` in
`crates/chematic-perception/src/aromaticity.rs` — a deliberately *separate*
implementation from the production `ring_pi_electrons` (not a wrapper, not a refactor
of it), mirroring its per-atom rules condition-for-condition but returning a
`ContributionReason` per atom instead of early-exiting `None` at the first ineligible
atom. Kept separate specifically so it never touches production behavior; the drift
risk that creates is closed by `trace_matches_ring_pi_electrons_on_corpus`, a new test
asserting `trace_ring_pi_electrons(...).total == ring_pi_electrons(...)` for every ring
in the full false-positive + false-negative + negative-control corpus, in both an empty
context (Pass-1-equivalent) and the model's converged final context (an *observational*
Pass-2 upper bound — not a literal iteration-order replay, since the model doesn't
expose per-iteration context; documented as such in the trace function's doc comment).

**3. Exocyclic-bond contribution changes** are recorded via
`ContributionReason::CarbonExocyclicHeteroatomDouble` (the carbonyl/imine 0π rule) —
already one of the traced reasons, not a separate mechanism.

**4. Per-component electron total and cycle rank**, output by
`crates/chematic-perception/examples/aromaticity_a1_0_report.rs`: for every molecule in
`validation/aromaticity_a1_0_corpus.jsonl`, walks each fused component's rings, reports
`ring_electron_total_{intrinsic,context}`, `ring_aromatic_{intrinsic,context}`, and
`cycle_rank` (= the component's SSSR ring count) alongside the real engine's per-atom
verdict (`assign_aromaticity_ex(mol).is_atom_aromatic`).

**5. Comparison against the RDKit oracle**: `scripts/aromaticity_a1_0_diagnosis.py`
joins the Rust trace output against real RDKit's per-atom aromatic flags (RDKit
bindings are Python-only in this project) and checks that every corpus bucket is
polarized the way its label claims.

### Corpus

`validation/aromaticity_a1_0_corpus.jsonl` (50 molecules, `scripts/gen_aromaticity_a1_0_corpus.py`):

| Bucket | Count | Source |
|---|---|---|
| `false_positive` | 33 | The pre-existing pinned 32-molecule `KNOWN_BRIDGEHEAD_N_FALSE_POSITIVES` regression corpus + PR #86's 14-atom minimal reproducer |
| `false_negative` | 5 | azulene + purine (the two `#[ignore]`d tests) + the pre-existing pinned 3-molecule `KNOWN_ORDER_DEPENDENT_FALSE_NEGATIVES` corpus |
| `negative_control` | 12 | 4-way isolation of the false-positive mechanism's two necessary conditions (bridgehead N alone, exocyclic C=C alone, both-removed, both-present), naphthalene, tetralin, indole, quinoline, 3 carbonyl variants (ruled out as the mechanism in SMARTS-A0), 2 non-aromatic conjugated rings |

The two "known-wrong" buckets are **not new data** — both are the project's own
pre-existing pinned regression tests, hoisted from inline `let cases` locals to named
module-level `const`s (`KNOWN_BRIDGEHEAD_N_FALSE_POSITIVES`,
`KNOWN_ORDER_DEPENDENT_FALSE_NEGATIVES`) so this corpus and those tests share one
source of truth instead of risking a second, silently-drifting copy.

### Polarization check (the "stopped over-propagation, also stopped correct
propagation" guard)

`scripts/aromaticity_a1_0_diagnosis.py`, run after the Rust report:

```
false_positive bucket polarized correctly:   33/33
false_negative bucket polarized correctly:   5/5
negative_control bucket fully agrees:         12/12
All buckets polarized as labeled. 0 problems.
```

Every `false_positive` molecule has at least one atom where chematic says aromatic and
RDKit says not; every `false_negative` molecule has at least one atom the reverse way;
every `negative_control` molecule has 100% per-atom agreement with RDKit. (The
corpus-wide raw per-atom-row agreement is 917/1136 = 80.7% — **not** a general accuracy
figure; this corpus is deliberately loaded with known-wrong cases for the polarization
check, not sampled to be representative.)

**A methodology finding surfaced by building this corpus**: purine's false negative
only reproduces on Kekulized input (matching the pre-existing `test_purine_aromatic`'s
exact input form via `mol_kekulized`). Feeding purine's corpus SMILES through the real
production path (raw aromatic-lowercase parse → `apply_aromaticity`, no explicit
kekulize step) does *not* reproduce it — that path never routes through the
`CarbonExocyclicHeteroatomDouble` rule at all, since that rule requires an explicit
`BondOrder::Double`, which direct aromatic-form parsing never produces (all ring bonds
are `BondOrder::Aromatic` instead). `aromaticity_a1_0_report.rs` therefore kekulizes
every corpus molecule uniformly before tracing, matching the pinned corpus's own
methodology and isolating Pass 1/Pass 2's Hückel-counting logic from the SMARTS-A0-
diagnosed, separately-tracked parser aromatic-flag behavior. This representation-
dependence is itself a real finding for A1-1 to account for, not a corpus bug.

### Files

- `crates/chematic-perception/src/aromaticity.rs`: `ContributionReason`,
  `AtomElectronTrace`, `RingElectronTrace`, `trace_ring_pi_electrons` (new, pure
  addition); `KNOWN_BRIDGEHEAD_N_FALSE_POSITIVES` / `KNOWN_ORDER_DEPENDENT_FALSE_NEGATIVES`
  (hoisted from test-local to module-level `const`); 3 new tests
  (`trace_matches_ring_pi_electrons_on_corpus`, `false_positive_corpus_over_counts_vs_rdkit`,
  `false_negative_corpus_under_counts_vs_rdkit`).
- `crates/chematic-perception/src/ring_family.rs`: `find_ring_families_over` (new,
  `find_ring_families` becomes a thin wrapper — zero behavior change).
- `crates/chematic-perception/examples/aromaticity_a1_0_report.rs` (new).
- `scripts/gen_aromaticity_a1_0_corpus.py`, `validation/aromaticity_a1_0_corpus.jsonl` (new).
- `scripts/aromaticity_a1_0_diagnosis.py`, `validation/results/aromaticity_a1_0_trace.jsonl`,
  `validation/results/aromaticity_a1_0_diagnosis.jsonl` (new).

**No change to `ring_pi_electrons`, `assign_aromaticity_ex`, `apply_aromaticity_ex`, or
any other production decision path.** Confirmed by the full workspace test suite and
`bash scripts/check.sh` passing unchanged, and by `trace_matches_ring_pi_electrons_on_corpus`
proving the new trace never disagrees with the untouched production function it mirrors.

## A1-1a: component model + exact oracle (this round, landed — no production wiring)

Implements the shared per-atom contribution function and component types, plus a
test/diagnostic-only exhaustive-candidate reference oracle. **`ring_pi_electrons`,
`assign_aromaticity_ex`, and `apply_aromaticity_ex` are untouched** — nothing described
here is wired into production. Types match the user's spec (renamed to match this
crate's existing `ContributionReason`/A1-0 naming instead of introducing a second,
parallel vocabulary):

```rust
pub enum PiEligibility { OneElectron, LonePairDonor, ZeroElectron, Ineligible }

pub struct ConjugatedComponent {
    pub atoms: Vec<AtomIdx>,
    pub bonds: Vec<BondIdx>,
    pub source_rings: Vec<usize>,
}

pub struct ContributionDecision {
    pub eligibility: PiEligibility,
    pub reason: ContributionReason,
}

pub fn evaluate_atom_pi_contribution(
    mol: &Molecule,
    atom_idx: AtomIdx,
    component: &ConjugatedComponent,
    algo: AromaticityAlgorithm,
) -> ContributionDecision;
```

**Single source of truth, as specified**: `trace_ring_pi_electrons` (A1-0) now
delegates to `evaluate_atom_pi_contribution` instead of duplicating its own copy of the
per-atom rules — the anti-drift test (`trace_matches_ring_pi_electrons_on_corpus`)
still passes unchanged, confirming the refactor is behavior-preserving.
`evaluate_atom_pi_contribution` is *not yet* called from `ring_pi_electrons` itself
(that wiring, behind a new opt-in `AromaticityAlgorithm` variant, is A1-1b).

**Component construction** (`build_conjugated_components`): connected components of a
"conjugation graph" over each ring family's atoms — nodes are eligible atoms (not
`Ineligible`), edges are *any* bond (single, double, or aromatic) between two eligible
atoms. Two real bugs were found and fixed while building this against the pinned
azulene/indolizine test cases, not guessed:

1. **Connectivity was too narrow.** The first version only bridged single bonds via a
   `LonePairDonor` endpoint (matching the user's original heteroatom-lone-pair framing).
   This left azulene's all-carbon alternating-bond perimeter as 5 disconnected 2-atom
   pairs — the oracle returned an *empty* aromatic set instead of all 10 atoms. Fixed:
   any bond between two eligible atoms connects, matching ordinary conjugation theory
   (butadiene's `C=C-C=C` middle single bond, styrene's vinyl-to-phenyl bond are both
   textbook-conjugated). Two sp2 atoms conjugate across a single bond regardless of
   whether either one is specifically a lone-pair donor.
2. **Degree-sensitive rules broke under a flattened multi-ring context.** The N
   bridgehead/substituted-azole rule tests "does this atom have a bond pointing outside
   *this ring*" (`ring_degree < total_degree`) — evaluated against a flattened
   multi-ring family, a *genuine* bridgehead's every bond is "in-family" by construction,
   so the test never fires and a true bridgehead N (indolizine's) came out `Ineligible`.
   This was a bug in the new A1-1a code, not a pre-existing chematic bug — caught because
   indolizine (a passing, correct chematic test) regressed under the naive component
   model. Fixed via `evaluate_atom_via_home_ring`: evaluate each atom's degree-sensitive
   rule against its own constituent ring, one at a time, not the flattened envelope.

**Test/diagnostic-only exhaustive oracle** (`exhaustive_aromaticity_oracle`): candidates
= every SSSR/augmented ring, individually, **plus** every multi-ring fused-envelope
candidate from `build_conjugated_components`; an atom/bond is oracle-aromatic if *any*
candidate containing it independently satisfies 4n+2, using
`evaluate_atom_pi_contribution` with **no cross-candidate context bootstrapping at all**
(unlike production's Pass 2). Pinned in `exhaustive_oracle_pinned_cases`
(RDKit-atom-index-verified, not guessed).

**Results, run against the corpus expanded to 55 molecules** (tropone, 2-pyridone,
4-pyranone, indolizine, and anthracene added as new `negative_control` entries —
`scripts/gen_aromaticity_a1_0_corpus.py`; joined against real RDKit via
`scripts/aromaticity_a1_0_diagnosis.py`):

| Case | Oracle matches RDKit? | Production (`current_engine`) matches RDKit? |
|---|---|---|
| tropone, 2-pyridone, 4-pyranone | ✅ | ✅ (already correct; confirms the atom-8 rule) |
| indolizine, anthracene | ✅ | ✅ (already correct; confirms A1-1a's own bug fix didn't regress them) |
| **azulene** | ✅ **(fixed)** | ❌ (known false negative) |
| **purine** | ❌ (still wrong) | ❌ (known false negative) |
| **false-positive reproducer (ring B)** | ❌ (still wrong, on purpose) | ❌ (known false positive) |
| All 17 negative-control molecules | ✅ 17/17 | ✅ 17/17 (0 regressions from A1-1a's own code) |
| All 33 false-positive corpus molecules | still over-count (unfixed) | still over-count |
| All 5 false-negative corpus molecules | azulene fixed, purine + 3 order-dependent still wrong | all still wrong |

Overall (informational only — this corpus is deliberately loaded with known-wrong
cases, not representative): oracle vs RDKit 987/1184 (83.36%) vs production vs RDKit
965/1184 (81.50%), oracle vs production 1162/1184 (98.14%) agreement.

**Two honest findings the oracle could NOT resolve, deliberately not chased further
this round** (per explicit scope discipline — solving either is A1-1b, not A1-1a):

- **The false-positive family is not fixable by candidate generation alone.** Ring B's
  *single-ring* candidate (evaluated via `ConjugatedComponent::from_ring`, independent
  of any multi-ring component) still sums to 6π and still passes on its own — adding
  more/better candidates doesn't suppress a bad candidate that's already generated. The
  real fix has to constrain *when* a lone-pair donor's electrons are creditable to a
  ring at all, which needs to know about the atom's other ring's validity — genuinely
  more than a local per-atom, per-candidate rule can decide.
- **Purine's fusion carbons need cross-ring information no single home ring or
  flattened family currently supplies correctly.** Before the home-ring fix, the
  (buggy) flattened evaluation happened to give purine all 9 atoms aromatic *by
  accident* — the same flattening bug that broke indolizine happened to produce the
  right answer here. After the fix, purine regresses to 6/9 on the oracle (still wrong,
  matching production's own pre-existing gap, not a new regression — production's
  Kekulized-input gap was already known from A1-0). This is pinned as an open finding
  in `exhaustive_oracle_pinned_cases`, not silently accepted.

**A1-1b's actual design question, sharpened by this round's evidence**: not "which
existing rule is wrong" (the atom-8 carbonyl rule is confirmed correct; candidate
generation is confirmed necessary-but-insufficient for the false-positive family) but
*how a lone-pair donor's contribution to one ring should depend on whether its other
ring membership is itself valid* — a genuinely different, harder question than local
per-atom or per-candidate evaluation can answer. Studying RDKit's own published
aromaticity algorithm is likely higher-yield for A1-1b's design than further reverse-
engineering rules from cases.

## A1-1b-0: RDKit-parity reference engine (this round, landed — no production wiring)

A1-1a's own conclusion was that local per-atom/per-candidate rule-guessing had hit its
ceiling — the false-positive family and purine both need cross-ring information no
local rule supplies. Per instruction, this round stops guessing new local rules from
cases and instead independently reproduces RDKit's actual default aromaticity algorithm
as a reference engine, source-verified against `Code/GraphMol/Aromaticity.cpp` (fetched
and read in full, not reverse-engineered from black-box test cases).

**The exact divergence point**, confirmed by reading RDKit's real source: RDKit computes
each atom's `ElectronDonorType` **once, globally per molecule**, before any candidate
ring is evaluated. Whether a multiple bond "counts" toward an atom's donor type depends
on whether that bond belongs to *any* SSSR ring in the whole molecule
(`RingInfo::numBondRings(bond) > 0`), not whether it's inside the specific ring/component
currently being scored. chematic's existing `ring_pi_electrons` and A1-1a's own
`evaluate_atom_pi_contribution` are architecturally different: both evaluate per-
candidate-ring/component. This is the precise, source-verified root cause of the
false-positive family — a lone-pair donor's electrons get credited to one ring
independent of whether that same atom's donor type would be disqualified by full-
molecule context.

New module `crates/chematic-perception/src/rdkit_parity.rs`, a faithful line-level port
of RDKit's real functions:

```rust
pub enum ElectronDonorType { Vacant, OneElectron, TwoElectron, OneOrTwo, Any, None }

pub fn get_atom_electron_donor_type(mol: &Molecule, atom_idx: AtomIdx, ring_bonds: &FxHashSet<BondIdx>) -> ElectronDonorType; // ported getAtomDonorTypeArom
pub fn is_atom_candidate_for_aromaticity(mol: &Molecule, atom_idx: AtomIdx, donor_type: ElectronDonorType) -> bool; // ported isAtomCandForArom
pub fn apply_huckel(mol: &Molecule, atoms: &[AtomIdx], donor: &FxHashMap<AtomIdx, ElectronDonorType>) -> bool; // ported applyHuckel

pub fn rdkit_parity_aromaticity(mol: &Molecule) -> (FxHashSet<AtomIdx>, FxHashSet<BondIdx>);
pub fn rdkit_parity_aromaticity_ex(mol: &Molecule, max_num_fused_rings: usize) -> (FxHashSet<AtomIdx>, FxHashSet<BondIdx>);
```

Driver flow, matching RDKit's `aromaticityHelper`: donor type per atom (global,
full-molecule ring-membership context) → candidate ring filter (`isAtomCandForArom`) →
fused-ring adjacency graph (rings sharing ≥1 bond, via union-find, matching
`makeRingNeighborMap`) → for each fused group, enumerate connected ring subsets of size
1..=6 (RDKit's own default `maxNumFusedRings` cap, via a hand-rolled `combinations(n, k)`
matching RDKit's `nextCombination` order — this project has no `itertools` dependency) →
union each subset's atoms, counting an atom only if it appears in exactly 1 or 2 of the
subset's rings (RDKit's own "#2895" fix, avoids double-counting a 3-ring-shared central
atom) → `applyHuckel` on the union → mark the subset's atoms/bonds aromatic on a pass,
matching `applyHuckelToFused`.

**Precondition, matching RDKit's own pipeline**: requires pre-kekulized input (no
`BondOrder::Aromatic`) — RDKit's own sanitization always runs `Kekulize` before
`setAromaticity` as genuinely separate steps, which is also why purine's two-path
representation-dependence needed checking explicitly (see below).

**Calibration battery** (`rdkit_parity.rs`'s own test module, all RDKit-atom-index-
verified before being pinned, not guessed): 15 cases spanning
benzene/pyrrole/furan/thiophene/tropone/2-pyridone/4-pyranone/indolizine/azulene/
naphthalene/anthracene/indole/quinoline/purine/the false-positive reproducer —
`calibration_battery_matches_rdkit`, 15/15 pass. `purine_representation_stable`: both
Kekulization paths (raw aromatic parse's own implicit kekulization vs chematic's
explicit `kekulize()`) agree with each other and with RDKit (all 9 atoms aromatic) —
confirms this round resolves the representation-dependence A1-0 first surfaced.

**Gate result on the 55-molecule diagnosis corpus** (`scripts/aromaticity_a1_1b_0_gate.py`,
joins `rdkit_parity_aromaticity`'s per-atom output against real RDKit per
`(smiles, atom_idx)`, same join methodology as every prior corpus gate in this doc):

```
false_positive corpus fixed:  33/33
false_negative corpus fixed: 5/5
negative_control maintained: 17/17
RDKit atom-flag agreement: 100.00% (1214/1214)
Unexplained differences: 0
=== GATE: PASS ===
```

Full gate met — all conditions the user specified for this round (33 FP fixed / 5 FN
fixed / 17 NC maintained / 100% atom-flag agreement / 0 unexplained diffs) are satisfied.

**Full 5,000-molecule benchmark corpus** (`~/Downloads/SMILES.csv`, the same corpus used
throughout this project's benchmarking — gitignored, user-supplied, not committed): run
via `crates/chematic-perception/examples/rdkit_parity_full_corpus.rs`, joined via
`scripts/aromaticity_a1_1b_0_full_corpus_gate.py` at **set level** (per-atom and
per-bond, not a count comparison — an earlier pass compared only aromatic atom/bond
*counts* per molecule, which is blind to same-count/different-atoms mismatches; this was
caught before being reported and redone as a real per-atom/per-bond join):

**100.0000% atom/bond agreement on all 4,999 comparable molecules** (1/5,000 excluded by
a pre-existing chematic `kekulize()` gap — RDKit parses the excluded molecule fine;
chematic's own `kekulize()` rejects a bridgehead N in a fused purine-like system,
`Cc1cn2c(=O)c3ncn(COCCO)c3nc2n1C` — a separate, narrower, pre-existing limitation not
touched by this round):

```
Comparable molecules:                4,999 / 5,000
Atoms:                              138,635 / 138,635
Bonds:                               150,004 / 150,004
Unexplained aromaticity differences: 0
```

This is a dramatic improvement over the current production baseline of 99.44% atom /
98.82% bond agreement (`docs/rdkit_compat.md`), and — because it is a source-verified
reproduction of RDKit's real algorithm rather than an approximation built from local
rules — appears to resolve both of A1-1a's open findings (the false-positive family and
purine) **by construction**, not by another round of case-specific patching.

**Still not wired into production.** `ring_pi_electrons`, `assign_aromaticity_ex`, and
`apply_aromaticity_ex` are untouched by this round — `rdkit_parity_aromaticity` is a
standalone reference engine, reachable only via its own public API and examples/tests.
Wiring it in behind an opt-in `AromaticityAlgorithm` variant (A1-1b proper) and later
promoting it to `apply_aromaticity("rdkit_like")`'s default (the PR-sequence step after
that) both require a fresh explicit go-ahead, per this initiative's standing diagnosis-
before-wiring discipline — the gate being met (and exceeded) is not itself that
go-ahead.

## A1-1b-1: fallible opt-in production API (this round, landed)

A1-1b-0's gate passed with room to spare (100.0000% set agreement on 4,999/5,000
comparable molecules). This round wires that reference engine into a real, callable
production API — but as a **separate fallible surface**, not a new
`AromaticityAlgorithm` variant, because the existing production API
(`assign_aromaticity_ex`/`apply_aromaticity_ex`) is infallible by contract (`Huckel` and
`RdkitLike` never fail), and this engine's precondition — pre-kekulized input — is not
always satisfiable (the one known `kekulize()` gap). Bending the existing infallible
signature to accommodate one fallible variant would have been the wrong shape; adding a
parallel `Result`-returning API keeps both contracts honest.

```rust
pub enum AromaticityError {
    KekulizationFailed { reason: String },
    InternalInvariantViolation { reason: String },
}

pub fn assign_aromaticity_rdkit_parity_experimental(
    mol: &Molecule,
) -> Result<AromaticityModel, AromaticityError>;

pub fn apply_aromaticity_rdkit_parity_experimental(
    mol: &Molecule,
) -> Result<Molecule, AromaticityError>;
```

Both are always available (no feature flag) — only the low-level `rdkit_parity` engine
internals remain gated. `mod rdkit_parity;` is now unconditionally compiled (it backs
production code), but every item inside it except the two functions above is
`pub(crate)`; the only way to reach the raw engine from outside the crate is through
those two functions, or, for diagnostics/benchmarking, a new
`chematic_perception::diagnostics` module (`#[doc(hidden)]`, gated behind the
`diagnostics` feature) that re-exports just `rdkit_parity_aromaticity`.

**Execution path**, matching the specified shape exactly:

```
input Molecule (&, never mutated)
  -> private clone with every atom's `aromatic` flag reset to false
      (this engine only ever reads BondOrder::Aromatic, never atom.aromatic,
      so this has no effect on the *computation* -- it exists purely so the
      *output* molecule's flags come entirely from this engine's own
      verdict, never leftover from whatever the caller passed in)
  -> kekulize() -- Err(KekulizationFailed { reason }) on failure, `mol` untouched
  -> rdkit_parity_aromaticity(): donor type computed once per molecule,
     candidate rings, fused-ring subsets, Hückel verdict (unchanged from A1-1b-0)
  -> validate_aromaticity_invariants(): every aromatic bond's two endpoints
     must themselves be aromatic atoms -- Err(InternalInvariantViolation) if not
     (0/4,999 fired on the full corpus; exists as a defensive no-panic
     guarantee, not because a violation was ever observed)
  -> AromaticityModel::from_atom_bond_sets() (ring_classifications() and
     antiaromatic_rings() are empty -- this engine, like RDKit's own,
     determines only the aromatic atom/bond sets)
  -> apply(): build_molecule_from_model() -- the SAME shared function
     apply_aromaticity_ex() uses (extracted, not duplicated) for implicit-H
     preservation, bond-direction stashing, and stereo-metadata copying
```

Every step takes `&Molecule` and returns a new value; nothing is mutated in place, so
there is no path that partially rewrites the input before failing (enforced by the type
system, and pinned by `production_api_does_not_mutate_input_on_failure`).

**Verification, not re-derivation.** Since A1-1b-0 already proved the underlying engine
against RDKit, this round's job was to prove the *wiring* is mechanical — that routing
through the new fallible entry points doesn't perturb the already-verified verdicts:

- 55-molecule diagnosis corpus, through `assign_aromaticity_rdkit_parity_experimental`
  (raw aromatic-form input, no manual pre-kekulization) — **byte-identical** JSONL output
  to the A1-1b-0 trace (1,214/1,214 rows).
- Full 5,000-molecule corpus, through the same production function — **byte-identical**
  atom/bond rows to the A1-1b-0 full-corpus trace (138,635 atom rows, 150,004 bond rows),
  transitively inheriting the 100.0000% RDKit set-agreement result without re-joining
  RDKit. Exactly 1 `KekulizationFailed` (the known gap), 0 `InternalInvariantViolation`,
  0 disagreements between `assign_...`/`apply_...` on the same input.

**Downstream checks — additive reading, not bit-identity-with-default.** The engine
changes aromaticity by design on the molecules where current production already
disagrees with RDKit (~0.56% atoms / ~1.18% bonds); demanding bit-identity between
experimental-applied output and default-applied output would fail by construction on
exactly the cases this engine exists to fix. Interpreted additively instead:

1. **Existing default-path suites stay green.** `assign_aromaticity_ex`/
   `apply_aromaticity_ex` are byte-unchanged (only the internals were reorganized into a
   shared `build_molecule_from_model` helper both paths call); `cargo test --workspace
   --lib` is 0 failures across every crate, including the existing 80,000-pair SMARTS
   corpus, descriptor/fingerprint/CIP-MANCUDE suites — none of which touch the new
   experimental path.
2. **Experimental-applied molecules don't crash or corrupt.** New Rust-level check
   (`crates/chematic-smarts/examples/aromaticity_a1_1b_1_downstream_check.rs` — lives in
   `chematic-smarts`, not `chematic-perception`, since the reverse dependency direction
   would be circular), full 5,000-molecule corpus, molecules produced by
   `apply_aromaticity_rdkit_parity_experimental`:
   - SMARTS matching (the same 16-pattern set `scripts/rdkit_compat_diff.py` uses for the
     existing 80,000-pair corpus) completes without panicking on every evaluation:
     **79,984/79,984** (4,999 molecules × 16 patterns). This does not re-run the
     Python-bound 80k corpus itself against the new engine — Python/WASM exposure is
     explicitly out of scope for this PR — it's a proportionate Rust-level substitute
     using the same query set at the same corpus scale.
   - Canonical SMILES round-trip idempotency: **4,912/4,999 (98.26%)** stable. This is
     *not* a new defect: the pre-existing idempotency baseline (documented in
     `docs/rdkit_compat.md`, "large fused-ring-system aromaticity round-trip", a
     canonicalization-tie-break issue upstream of this PR) is already ~98.4% on the same
     corpus under the *default* pipeline (measured same-run baseline: 69/5,000 = 98.62%
     unstable via `apply_aromaticity_ex(.., RdkitLike)`). The two unstable sets overlap
     62/87 (71%); tracing a sample of the 25 molecules unstable *only* under the
     experimental path found their final aromatic atom/bond verdict is **identical** to
     the default path's — the divergence traces instead to `build_molecule_from_model`'s
     implicit-H "Kekule-then-perceive" freeze logic, which behaves slightly differently
     depending on whether its `mol` argument was already Kekule-form (experimental always
     passes its own internally-kekulized clone) or still aromatic-form (default's `mol`
     argument, when the caller passes genuinely aromatic-written SMILES directly, as most
     of this corpus is). Both are internally consistent with the shared function's
     documented intent; this shifts *which* molecules hit the pre-existing
     idempotency gap rather than introducing a new one. Not chased further — root-causing
     the exact shift is out of scope for a mechanical-wiring PR; flagged here rather than
     rounded off.
3. **Representation parity: 4,999/4,999 (100%).** Aromatic-form input and explicitly
   pre-kekulized input produce an identical aromatic atom set through the production API
   — generalizes `purine_representation_stable` (one pinned molecule) to the full corpus.

**Performance — recorded, not optimized** (this PR does not tune either engine for
speed; `crates/chematic-perception/examples/aromaticity_a1_1b_1_perf.rs`, full
5,000-molecule corpus, warm-up pass excluded from timing):

| Engine | mean | p50 | p95 | max |
|---|---|---|---|---|
| `RdkitLike` (current production default) | 408.3µs | 260.2µs | 1021.5µs | 27.7ms |
| `RdkitParityExperimental` | 424.7µs | 275.1µs | 1062.2µs | 27.3ms |

~4% slower on mean/p50/p95 — essentially comparable, not the order-of-magnitude gap seen
elsewhere in this project's opt-in-accuracy tradeoffs (e.g. CIP Accurate's ~10x). Peak
RSS is a coarse whole-process figure (`/usr/bin/time -l` around the full benchmark run,
both engines plus warm-up in one process): ~19MB: not cleanly separable per-engine in a
single process, and not expected to differ meaningfully — both algorithms allocate
`FxHashSet`/`FxHashMap` sized to one molecule's ring/atom count per call, freed
immediately after, so peak RSS is dominated by corpus loading, not algorithm choice.

**Non-goals held**: `RdkitLike` is not replaced, the default is not changed, the one
known `kekulize()` gap is not fixed (it surfaces as `KekulizationFailed`, not a silent
fallback), the reference engine's algorithm is not touched, no SMARTS/canonical-specific
exceptions were added, and this round does not touch Python/WASM bindings or README
default-accuracy numbers.

## A1's gate (for the eventual production-enabling PR, not A1-0)

Per the user's spec — not evaluated in this round, recorded here so A1-1/A1-2's PRs are
measured against a fixed target instead of a moving one:

- SMARTS: 80,000/80,000 = 100.00%
- False-positive corpus (this doc's 33, expandable): 33/33 fixed
- False-negative corpus (azulene, purine, the 3 order-dependent cases): all fixed, or
  at minimum 0 regressions
- Aromatic atom/bond agreement: must not regress below the current 99.44%/98.82%;
  target 99.8%+
- Canonical SMILES: semantic round-trip ≥99.62%, idempotency ≥98.42%, 0 new structural
  breaks
- Descriptors: all existing MW/HBA/HBD/TPSA/LogP gates maintained
- Fingerprints: 0 regressions in existing similarity/ranking
- Performance: aromaticity-batch degradation within 5%
- Representation: identical result across aromatic/Kekulé respelling (already a
  precondition SMARTS-A0 verified for the false-positive family; A1-1 must not
  reintroduce representation-dependence, including the purine gap A1-0 just surfaced)
- README numbers updated only after a full 5,000-molecule corpus regeneration

## PR sequence

1. PR #86 — SMARTS-A0 diagnosis (merged).
2. PR #87 — Aromaticity-A1-0: component trace + corpus, no behavior change (merged).
3. PR #88 — Aromaticity-A1-1a: component model (`PiEligibility`, `ConjugatedComponent`,
   `evaluate_atom_pi_contribution`) + exhaustive reference oracle, no production wiring
   (merged). Fixed 2 bugs found in this round's own new code (connectivity, home-ring
   degree evaluation); confirmed the atom-8 carbonyl rule is correct (not the fix
   target); confirmed candidate generation alone fixes azulene but not the
   false-positive family or purine — both left as open, pinned findings for A1-1b.
4. PR #89 — Aromaticity-A1-1b-0: RDKit-parity reference engine
   (`rdkit_parity_aromaticity`), a source-verified port of RDKit's actual default
   aromaticity algorithm (`ElectronDonorType`, candidate rings, fused-ring adjacency,
   `applyHuckelToFused`-style subset search), no production wiring (merged). Passed the
   full 55-molecule gate (33 FP / 5 FN / 17 NC / 100% atom-flag agreement / 0
   unexplained diffs) and, beyond the stated gate, reached 100.0000% set-level atom/bond
   agreement with real RDKit on all 4,999 comparable molecules of the full 5,000-molecule
   benchmark corpus (1 excluded by a pre-existing, unrelated `kekulize()` gap — vs the
   99.44%/98.82% production baseline) — resolved both of A1-1a's open findings
   (false-positive family, purine) by construction. Ports specific RDKit functions under
   RDKit's BSD 3-Clause license; attribution and full license text recorded in
   `THIRD_PARTY_NOTICES.md`.
5. **This PR — Aromaticity-A1-1b-1**: fallible opt-in production API
   (`assign_aromaticity_rdkit_parity_experimental`/
   `apply_aromaticity_rdkit_parity_experimental`, `AromaticityError`), always available
   (no feature flag); the low-level `rdkit_parity` engine internals are `pub(crate)`,
   reachable externally only through those two functions or, for diagnostics, the new
   `diagnostics` module (gated behind the `diagnostics` feature). Mechanical wiring only
   — no algorithm change, no `RdkitLike` replacement, no default change. Verified
   byte-identical to the A1-1b-0 trace on both the 55-molecule and full 5,000-molecule
   corpora (proving the wiring layer doesn't perturb the already-verified engine);
   downstream checks (SMARTS-doesn't-panic, canonical round-trip, representation parity)
   read additively rather than as bit-identity-with-default, since this engine changes
   aromaticity by design on the molecules where production already disagrees with RDKit.
6. A1-1b-2 — Python/WASM opt-in exposure of the new production API (not started).
7. K0 — diagnose the one remaining `kekulize()` gap
   (`Cc1cn2c(=O)c3ncn(COCCO)c3nc2n1C`) that currently surfaces as
   `AromaticityError::KekulizationFailed` (not started).
8. A1-1c — downstream regression/perf evaluation at a scope beyond this PR's proportionate
   checks, if the additive-reading interpretation above needs revisiting (not started).
9. A1-1d — decide whether to promote the RDKit-parity engine to
   `apply_aromaticity("rdkit_like")`'s default, informed by A1-1b-2/K0/A1-1c (not started).
