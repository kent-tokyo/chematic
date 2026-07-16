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

On the SMARTS-A0 minimal reproducer `C1=Cc2ccccc2C2=NCCCN12` (benzo ring fused via a
bridgehead N to a third, non-aromatic ring, with an exocyclic-to-benzo `C=C` on the
middle ring): Pass 1 correctly marks the benzo ring (atoms `2..7`) aromatic. Pass 2
re-evaluates the adjacent ring (atoms `0,1,2,7,8,13`) using that context:

- atoms `2,7`: `AlreadyAromaticContext` → 1π each (borrowed from the benzo ring, per
  `ring_pi_electrons`'s context rule)
- atoms `0,1`: `CarbonEndocyclicDouble` → 1π each (`0=1`, in-ring for this ring)
- atom `8`: `CarbonExocyclicHeteroatomDouble` → 0π (`C=N` to atom `9`, exocyclic to
  *this* ring)
- atom `13`: `NitrogenBridgeheadOrSubstitutedLonePair` → 2π (bridgehead N rule: all
  three σ-bonds fill its valence, ring-degree < total-degree)

Total: 1+1+1+1+0+2 = **6π → passes 4n+2 (n=1) → the whole ring gets marked aromatic**,
even though ring C (the third ring this same bridgehead N also belongs to) has three
sp3 CH₂ carbons and cannot itself be part of any real delocalized system. The bridgehead-
N rule's 2π credit is correct for indolizine (both rings genuinely aromatic) but wrong
here, because the rule does not check whether the atom's *other* ring is actually a
valid aromatic partner before granting the lone-pair credit to *this* ring's count.

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

## A1-1: component-level solver (designed, not implemented)

Target shape, replacing the current Pass 1/Pass 2 boolean-flood-fill with candidate
discovery separated from electron-count decision:

```rust
struct AromaticCandidateComponent {
    atoms: Vec<AtomIdx>,
    bonds: Vec<BondIdx>,
    cycles: Vec<CycleId>,
}

struct ElectronContribution {
    electrons: u8,
    reason: ContributionReason,
}

fn evaluate_component(
    mol: &Molecule,
    component: &AromaticCandidateComponent,
) -> AromaticityDecision;
```

Flow: ring-candidate extraction → fused-component construction (A1-0's
`find_ring_families_over`, already reusable as-is) → per-atom conjugation eligibility →
reflect exocyclic multiple bonds/charge/valence (A1-0's `trace_ring_pi_electrons` rules
are the starting point, not a rewrite from scratch) → π-electron aggregation *at the
component level*, not just per individual SSSR ring → confirm aromatic atom/bond sets.
Current Pass 2 is demoted to a candidate-discovery input to this solver where possible,
not deleted outright (Pass 1's per-ring evaluation remains valid for simple, unfused
rings — the solver only needs to change how *fused* components are judged).

Not started this round. A1-0 provides the characterization surface (component
construction + trace + 3-bucket regression corpus with a working polarization guard) A1-1
needs to build against and measure regressions with — in particular, the same
`false_positive`/`false_negative`/`negative_control` buckets, expanded, become A1-1's
primary regression gate.

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
2. **This PR — Aromaticity-A1-0**: component trace + corpus, no behavior change.
3. A1-1 — component solver, test-only/opt-in (not started).
4. Enable the solver into `apply_aromaticity("rdkit_like")` + full regression against
   the gate above (not started).
5. Docs/benchmark updates reflecting the final measured numbers (not started).
