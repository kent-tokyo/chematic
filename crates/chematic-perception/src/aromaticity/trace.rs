// ---------------------------------------------------------------------------
// Diagnostic trace (Aromaticity-A1-0) — observational only, no production
// behavior change. `ring_pi_electrons` above is untouched and remains the
// single source of truth for `assign_aromaticity_ex`'s actual decisions;
// this is a parallel, read-only explanation layer for `component/atom/reason`
// tracing, used by `aromaticity_a1_0_report` and the corpus diagnostics in
// `validation/aromaticity_a1_0_corpus.jsonl`. See `docs/rfcs/aromaticity_a1_rfc.md`.
// ---------------------------------------------------------------------------

/// Reason a ring atom contributes (or fails to contribute) pi electrons,
/// mirroring `ring_pi_electrons`'s branches one-to-one. Purely diagnostic.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContributionReason {
    /// Already aromatic from a previous Pass 1/Pass 2 ring: contributes 1π unconditionally.
    AlreadyAromaticContext,
    /// Carbon with an endocyclic double/aromatic bond: 1π.
    CarbonEndocyclicDouble,
    /// Carbon whose only double bond is exocyclic to O/N/S: 0π (e.g. a carbonyl carbon).
    CarbonExocyclicHeteroatomDouble,
    /// Carbanion with no double bond: 2π (lone pair).
    CarbonCarbanionLonePair,
    /// Cationic ring carbon (e.g. tropylium's `[cH+]`): empty p-orbital
    /// electron acceptor, 0π, regardless of representation (Kekule or
    /// aromatic-bond) -- mirrors `CarbonCarbanionLonePair`'s anion rule at
    /// the opposite electron-count extreme.
    CarbonCationVacant,
    /// sp3 carbon (no double bond, not a carbanion): ineligible.
    CarbonSp3Ineligible,
    /// Pyrrole-type N with an H, neutral or anionic: 2π.
    NitrogenPyrroleTypeH,
    /// Pyridine-type N with an explicit double bond (bare, or protonated
    /// N-H+ once it has a Kekule double bond): 1π.
    NitrogenPyridineTypeExplicitDouble,
    /// Bridgehead N (or N-substituted azole N), neutral or anionic:
    /// all-sigma valence, lone pair in p orbital: 2π.
    NitrogenBridgeheadOrSubstitutedLonePair,
    /// N with an in-ring aromatic bond, not a bridgehead (pyridine-type
    /// notation, or a charged N-H+ in aromatic-bond representation): 1π.
    NitrogenAromaticInRing,
    /// N matching none of the above rules: ineligible.
    NitrogenIneligible,
    /// O/S/Se/Te lone-pair donor, neutral or anionic, ring-degree 2: 2π.
    ChalcogenLonePair,
    /// P lone-pair donor in the opt-in RDKit-compatible model: 2π.
    PnictogenOrChalcogenLonePair,
    /// Charged O/S (e.g. pyrylium's `[o+]`): the positive charge consumes
    /// the lone pair, so this atom needs pyridine-type treatment (1π via
    /// its own ring double/aromatic bond) instead of donating 2π.
    ChalcogenCationPyridineType,
    /// O/S/Se/Te with the wrong ring degree, an exocyclic X=O, (Se/Te)
    /// non-RdkitLike mode, or a charged O/S with no ring double/aromatic
    /// bond to fall back on: ineligible.
    ChalcogenIneligible,
    /// Element not supported by the model: ineligible.
    UnsupportedElement,
}

impl ContributionReason {
    /// Whether this reason is an eligible contribution (matches
    /// `ring_pi_electrons` returning `Some`) rather than one that disqualifies
    /// the whole ring (matches it returning `None`).
    pub fn is_eligible(self) -> bool {
        !matches!(
            self,
            ContributionReason::CarbonSp3Ineligible
                | ContributionReason::NitrogenIneligible
                | ContributionReason::ChalcogenIneligible
                | ContributionReason::UnsupportedElement
        )
    }

    /// Coarse `PiEligibility` bucket for this fine-grained reason
    /// (Aromaticity-A1-1a). `AlreadyAromaticContext` has no single fixed
    /// bucket -- it always carries exactly 1π, so it maps to `OneElectron`.
    pub fn eligibility(self) -> PiEligibility {
        use ContributionReason::*;
        match self {
            AlreadyAromaticContext
            | CarbonEndocyclicDouble
            | NitrogenPyridineTypeExplicitDouble
            | NitrogenAromaticInRing
            | ChalcogenCationPyridineType => PiEligibility::OneElectron,
            CarbonCarbanionLonePair
            | NitrogenPyrroleTypeH
            | NitrogenBridgeheadOrSubstitutedLonePair
            | ChalcogenLonePair
            | PnictogenOrChalcogenLonePair => PiEligibility::LonePairDonor,
            CarbonExocyclicHeteroatomDouble | CarbonCationVacant => PiEligibility::ZeroElectron,
            CarbonSp3Ineligible | NitrogenIneligible | ChalcogenIneligible | UnsupportedElement => {
                PiEligibility::Ineligible
            }
        }
    }
}

/// Coarse per-atom pi-eligibility bucket (Aromaticity-A1-1a). A summary view
/// over [`ContributionReason`]'s finer-grained rules -- `electrons()` gives
/// the electron count implied by each bucket.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PiEligibility {
    /// Contributes exactly 1π (e.g. an endocyclic double/aromatic bond).
    OneElectron,
    /// Contributes 2π (a lone pair: pyrrole-type N, chalcogen, bridgehead N, carbanion).
    LonePairDonor,
    /// Contributes 0π but is still sp2 (p-orbital spent on an exocyclic multiple bond).
    ZeroElectron,
    /// Not eligible to be part of any conjugated system (e.g. sp3).
    Ineligible,
}

impl PiEligibility {
    /// Electron count implied by this bucket, or `None` for `Ineligible`.
    pub fn electrons(self) -> Option<u8> {
        match self {
            PiEligibility::OneElectron => Some(1),
            PiEligibility::LonePairDonor => Some(2),
            PiEligibility::ZeroElectron => Some(0),
            PiEligibility::Ineligible => None,
        }
    }
}

/// A candidate conjugated system: some atoms/bonds evaluated together as one
/// pi-electron-counting problem (Aromaticity-A1-1a). Two distinct uses:
/// - a single SSSR/augmented ring, reinterpreted as a trivial one-ring
///   candidate (what `trace_ring_pi_electrons` builds today);
/// - a genuine multi-ring fused envelope, built by
///   [`build_conjugated_components`] as a connected component of the
///   "conjugation graph" (double/aromatic-bonded atoms, plus lone-pair-donor
///   atoms bridging across single bonds) -- the azulene-class candidate
///   `augmented_ring_set`'s own docstring already named as future work
///   ("candidate rings = SSSR ∪ fused envelopes").
#[derive(Debug, Clone)]
pub struct ConjugatedComponent {
    pub atoms: Vec<AtomIdx>,
    pub bonds: Vec<BondIdx>,
    /// Ring indices (into whatever ring list the caller built this from) this
    /// candidate derives from -- one entry for a plain single-ring candidate,
    /// 2+ for a fused envelope spanning multiple rings.
    pub source_rings: Vec<usize>,
}

impl ConjugatedComponent {
    /// Build a trivial single-ring candidate from one ring's atom list (no
    /// bond list needed by [`evaluate_atom_pi_contribution`], which only
    /// consults `atoms` membership).
    fn from_ring(ring: &[AtomIdx], ring_idx: usize) -> Self {
        ConjugatedComponent {
            atoms: ring.to_vec(),
            bonds: Vec::new(),
            source_rings: vec![ring_idx],
        }
    }
}

/// The full per-atom decision from [`evaluate_atom_pi_contribution`]: the
/// coarse eligibility bucket plus the specific rule that produced it.
#[derive(Debug, Clone, Copy)]
pub struct ContributionDecision {
    pub eligibility: PiEligibility,
    pub reason: ContributionReason,
}

impl ContributionDecision {
    pub fn electrons(&self) -> Option<u8> {
        self.eligibility.electrons()
    }
}

/// Per-atom trace entry from [`trace_ring_pi_electrons`].
#[derive(Debug, Clone, Copy)]
pub struct AtomElectronTrace {
    pub atom_idx: AtomIdx,
    /// `None` iff `reason.is_eligible()` is false.
    pub contribution: Option<u8>,
    pub reason: ContributionReason,
}

/// Full per-atom pi-electron trace for one ring — the diagnostic twin of
/// [`ring_pi_electrons`]. Unlike `ring_pi_electrons` (which returns `None` at
/// the first ineligible atom), this always scans every atom so a caller can
/// see exactly which atom(s) disqualify a ring, not just that one did.
#[derive(Debug, Clone)]
pub struct RingElectronTrace {
    pub atoms: Vec<AtomElectronTrace>,
    /// `Some(sum)` iff every atom was eligible — must equal
    /// `ring_pi_electrons(mol, ring, aromatic_context, algo, all_ring_bonds)`
    /// for the same inputs (checked by
    /// `trace_matches_ring_pi_electrons_on_corpus` below).
    pub total: Option<u32>,
}

/// Diagnostic twin of [`ring_pi_electrons`]: identical per-atom rules
/// (delegating to [`evaluate_atom_pi_contribution`], the single source of
/// truth for both this trace and any future experimental production path —
/// see `docs/rfcs/aromaticity_a1_rfc.md`'s A1-1a section), but returns a full
/// trace instead of a single early-exiting `Option<u32>`. Does not call,
/// wrap, or change `ring_pi_electrons` itself — zero effect on
/// `assign_aromaticity_ex`'s behavior. `trace_matches_ring_pi_electrons_on_corpus`
/// is the anti-drift guard that keeps this and `ring_pi_electrons` in sync.
pub fn trace_ring_pi_electrons(
    mol: &Molecule,
    ring: &[AtomIdx],
    aromatic_context: &FxHashSet<AtomIdx>,
    algo: AromaticityAlgorithm,
    all_ring_bonds: &FxHashSet<BondIdx>,
) -> RingElectronTrace {
    let component = ConjugatedComponent::from_ring(ring, 0);
    let mut atoms = Vec::with_capacity(ring.len());
    let mut total: Option<u32> = Some(0);

    for &atom_idx in ring {
        let (contribution, reason) = if aromatic_context.contains(&atom_idx) {
            (Some(1u8), ContributionReason::AlreadyAromaticContext)
        } else {
            let decision =
                evaluate_atom_pi_contribution(mol, atom_idx, &component, algo, all_ring_bonds);
            (decision.electrons(), decision.reason)
        };

        total = match (total, contribution) {
            (Some(t), Some(c)) => Some(t + c as u32),
            _ => None,
        };

        atoms.push(AtomElectronTrace {
            atom_idx,
            contribution,
            reason,
        });
    }

    RingElectronTrace { atoms, total }
}

/// Single source of truth for per-atom pi-electron contribution
/// (Aromaticity-A1-1a): identical rules to `ring_pi_electrons`'s match arms,
/// condition-for-condition, parameterized by an arbitrary candidate
/// [`ConjugatedComponent`] instead of one fixed SSSR ring — the same
/// function evaluates a plain single-ring candidate (via
/// `ConjugatedComponent::from_ring`) or a genuine multi-ring fused envelope
/// (via `build_conjugated_components`) identically. Currently called by
/// `trace_ring_pi_electrons` only — NOT wired into `ring_pi_electrons` or
/// `assign_aromaticity_ex` (that wiring, behind a new opt-in
/// `AromaticityAlgorithm` variant, is Aromaticity-A1-1b, not this round).
pub fn evaluate_atom_pi_contribution(
    mol: &Molecule,
    atom_idx: AtomIdx,
    component: &ConjugatedComponent,
    algo: AromaticityAlgorithm,
    all_ring_bonds: &FxHashSet<BondIdx>,
) -> ContributionDecision {
    let component_atoms: FxHashSet<AtomIdx> = component.atoms.iter().copied().collect();
    let (_electrons, reason) =
        evaluate_atom_pi_contribution_inner(mol, atom_idx, &component_atoms, algo, all_ring_bonds);
    // `reason.eligibility().electrons()` is asserted equal to `_electrons`
    // for every branch by `contribution_decision_electrons_match_inner_on_corpus`.
    ContributionDecision {
        eligibility: reason.eligibility(),
        reason,
    }
}

/// Per-atom contribution logic, mirroring `ring_pi_electrons`'s match arms
/// condition-for-condition, but returning a reason alongside the
/// contribution instead of returning early on `None`.
fn evaluate_atom_pi_contribution_inner(
    mol: &Molecule,
    atom_idx: AtomIdx,
    ring_atom_set: &FxHashSet<AtomIdx>,
    algo: AromaticityAlgorithm,
    all_ring_bonds: &FxHashSet<BondIdx>,
) -> (Option<u8>, ContributionReason) {
    let atom = mol.atom(atom_idx);
    let an = atom.element.atomic_number();
    let context = pi_atom_context(mol, atom_idx, ring_atom_set);
    match an {
        6 => carbon_pi_contribution(mol, atom_idx, ring_atom_set, all_ring_bonds, context),
        7 => nitrogen_pi_contribution(mol, atom_idx, context),
        8 | 16 => chalcogen_pi_contribution(mol, atom_idx, ring_atom_set, context),
        15 | 34 | 52 => extended_lone_pair_contribution(
            mol,
            atom_idx,
            ring_atom_set,
            algo,
            context,
        ),
        _ => (None, ContributionReason::UnsupportedElement),
    }
}

#[derive(Clone, Copy)]
struct PiAtomContext {
    ring_degree: usize,
    total_degree: usize,
    total_coordination: usize,
    has_explicit_double: bool,
    has_double_any: bool,
    has_aromatic_in_ring: bool,
}

fn pi_atom_context(
    mol: &Molecule,
    atom_idx: AtomIdx,
    ring_atom_set: &FxHashSet<AtomIdx>,
) -> PiAtomContext {
    let ring_degree = mol
        .neighbors(atom_idx)
        .filter(|(neighbor, _)| ring_atom_set.contains(neighbor))
        .count();
    let total_degree = mol.degree(atom_idx);
    let has_explicit_double = mol
        .neighbors(atom_idx)
        .any(|(_, bond)| mol.bond(bond).order == BondOrder::Double);
    PiAtomContext {
        ring_degree,
        total_degree,
        total_coordination: total_degree + implicit_hcount(mol, atom_idx) as usize,
        has_explicit_double,
        has_double_any: has_explicit_double
            || mol
                .neighbors(atom_idx)
                .any(|(_, bond)| mol.bond(bond).order == BondOrder::Aromatic),
        has_aromatic_in_ring: mol
            .neighbors(atom_idx)
            .filter(|(neighbor, _)| ring_atom_set.contains(neighbor))
            .any(|(_, bond)| mol.bond(bond).order == BondOrder::Aromatic),
    }
}

fn has_exocyclic_double(
    mol: &Molecule,
    atom_idx: AtomIdx,
    ring_atom_set: &FxHashSet<AtomIdx>,
) -> bool {
    mol.neighbors(atom_idx).any(|(neighbor, bond)| {
        !ring_atom_set.contains(&neighbor) && mol.bond(bond).order == BondOrder::Double
    })
}

fn carbon_pi_contribution(
    mol: &Molecule,
    atom_idx: AtomIdx,
    ring_atom_set: &FxHashSet<AtomIdx>,
    all_ring_bonds: &FxHashSet<BondIdx>,
    context: PiAtomContext,
) -> (Option<u8>, ContributionReason) {
    let charge = mol.atom(atom_idx).charge;
    if charge > 0 {
        return (Some(0), ContributionReason::CarbonCationVacant);
    }
    if !context.has_double_any {
        return if charge == -1 {
            (Some(2), ContributionReason::CarbonCarbanionLonePair)
        } else {
            (None, ContributionReason::CarbonSp3Ineligible)
        };
    }
    let endocyclic_double = mol.neighbors(atom_idx).any(|(neighbor, bond)| {
        ring_atom_set.contains(&neighbor) && mol.bond(bond).order == BondOrder::Double
    });
    let heteroatom_double = mol.neighbors(atom_idx).any(|(neighbor, bond)| {
        !all_ring_bonds.contains(&bond)
            && mol.bond(bond).order == BondOrder::Double
            && matches!(mol.atom(neighbor).element.atomic_number(), 7 | 8 | 16)
    });
    if context.has_explicit_double
        && !context.has_aromatic_in_ring
        && !endocyclic_double
        && heteroatom_double
    {
        (Some(0), ContributionReason::CarbonExocyclicHeteroatomDouble)
    } else {
        (Some(1), ContributionReason::CarbonEndocyclicDouble)
    }
}

fn nitrogen_pi_contribution(
    mol: &Molecule,
    atom_idx: AtomIdx,
    context: PiAtomContext,
) -> (Option<u8>, ContributionReason) {
    let atom = mol.atom(atom_idx);
    if implicit_hcount(mol, atom_idx) > 0 && atom.charge <= 0 {
        (Some(2), ContributionReason::NitrogenPyrroleTypeH)
    } else if context.has_explicit_double {
        (
            Some(1),
            ContributionReason::NitrogenPyridineTypeExplicitDouble,
        )
    } else if context.total_degree == 3
        && context.ring_degree < context.total_degree
        && atom.charge <= 0
    {
        (
            Some(2),
            ContributionReason::NitrogenBridgeheadOrSubstitutedLonePair,
        )
    } else if context.has_aromatic_in_ring {
        (Some(1), ContributionReason::NitrogenAromaticInRing)
    } else {
        (None, ContributionReason::NitrogenIneligible)
    }
}

fn chalcogen_pi_contribution(
    mol: &Molecule,
    atom_idx: AtomIdx,
    ring_atom_set: &FxHashSet<AtomIdx>,
    context: PiAtomContext,
) -> (Option<u8>, ContributionReason) {
    let atom = mol.atom(atom_idx);
    if atom.charge > 0 {
        return if context.has_explicit_double || context.has_aromatic_in_ring {
            (Some(1), ContributionReason::ChalcogenCationPyridineType)
        } else {
            (None, ContributionReason::ChalcogenIneligible)
        };
    }
    let ineligible = context.ring_degree != 2
        || (atom.charge == 0 && context.total_coordination > 2)
        || (atom.element.atomic_number() == 16
            && has_exocyclic_double(mol, atom_idx, ring_atom_set));
    if ineligible {
        (None, ContributionReason::ChalcogenIneligible)
    } else {
        (Some(2), ContributionReason::ChalcogenLonePair)
    }
}

fn extended_lone_pair_contribution(
    mol: &Molecule,
    atom_idx: AtomIdx,
    ring_atom_set: &FxHashSet<AtomIdx>,
    algo: AromaticityAlgorithm,
    context: PiAtomContext,
) -> (Option<u8>, ContributionReason) {
    let atom = mol.atom(atom_idx);
    let ineligible = algo != AromaticityAlgorithm::RdkitLike
        || context.ring_degree != 2
        || (matches!(atom.element.atomic_number(), 34 | 52)
            && atom.charge == 0
            && context.total_coordination > 2)
        || has_exocyclic_double(mol, atom_idx, ring_atom_set);
    if ineligible {
        (None, ContributionReason::ChalcogenIneligible)
    } else {
        (Some(2), ContributionReason::PnictogenOrChalcogenLonePair)
    }
}

/// Evaluate an atom's pi contribution using its "home ring" within a
/// (possibly multi-ring) candidate, instead of the candidate's flattened
/// atom set directly: tries each of `candidate.source_rings` that actually
/// contains the atom, evaluating against *that one ring's own* atom set, and
/// returns the first eligible result found. Falls back to evaluating
/// directly against the flattened `candidate` if `source_rings` is empty or
/// none of them contain the atom (shouldn't happen for well-formed
/// candidates, but keeps this total rather than panicking).
///
/// Needed because degree-sensitive rules (the N bridgehead/substituted-azole
/// rule, `total_degree == 3 && ring_degree < total_degree`) test "does this
/// atom have a bond that points outside THIS ring" -- a genuine multi-ring
/// bridgehead's every bond is "in-family" once the evaluation context is the
/// flattened whole envelope (every neighbor is, by construction, some other
/// family member), which silently defeats that test and makes a real
/// bridgehead N (e.g. indolizine's) look `Ineligible`. Evaluating against
/// one constituent ring at a time preserves the rule's original, correct,
/// per-ring meaning even when the candidate spans multiple rings. This does
/// **not** attempt to resolve whether a bridgehead's lone-pair credit is
/// *legitimately shared* between two rings that are both otherwise valid vs.
/// wrongly borrowed by one ring from another that's actually broken (e.g.
/// by an sp3 atom) -- that is a distinct, harder, open question, deliberately
/// left to Aromaticity-A1-1b (see `docs/rfcs/aromaticity_a1_rfc.md`).
fn evaluate_atom_via_home_ring(
    mol: &Molecule,
    atom_idx: AtomIdx,
    candidate: &ConjugatedComponent,
    rings: &[Vec<AtomIdx>],
    algo: AromaticityAlgorithm,
    all_ring_bonds: &FxHashSet<BondIdx>,
) -> ContributionDecision {
    let mut last = None;
    for &ri in &candidate.source_rings {
        if !rings[ri].contains(&atom_idx) {
            continue;
        }
        let home = ConjugatedComponent::from_ring(&rings[ri], ri);
        let decision = evaluate_atom_pi_contribution(mol, atom_idx, &home, algo, all_ring_bonds);
        if decision.electrons().is_some() {
            return decision;
        }
        last = Some(decision);
    }
    last.unwrap_or_else(|| {
        evaluate_atom_pi_contribution(mol, atom_idx, candidate, algo, all_ring_bonds)
    })
}

/// Build genuine multi-ring conjugated-system candidates (Aromaticity-A1-1a):
/// connected components of the "conjugation graph" over each ring family's
/// atoms -- nodes are atoms whose eligibility (evaluated per-atom against its
/// own home ring, via `evaluate_atom_via_home_ring` -- not the flattened
/// family) is not `Ineligible`; edges are any bond (single, double, or
/// aromatic) between two independently-eligible family atoms: ordinary
/// carbon-carbon single-bond conjugation (butadiene's C=C-C=C middle bond,
/// styrene's vinyl-to-phenyl bond) connects just as directly as a
/// lone-pair-donor heteroatom bridging a sigma bond.
///
/// A pure candidate *generator* -- callers (currently only
/// `exhaustive_aromaticity_oracle`) still run full 4n+2 electron counting on
/// each result. Only components spanning 2+ of a family's rings are
/// returned: a single unfused ring is already covered by
/// `ConjugatedComponent::from_ring`, so this only adds the fused-envelope
/// candidates `augmented_ring_set`'s docstring named as future work
/// ("candidate rings = SSSR ∪ fused envelopes").
pub fn build_conjugated_components(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    ring_families: &[RingFamily],
    algo: AromaticityAlgorithm,
    all_ring_bonds: &FxHashSet<BondIdx>,
) -> Vec<ConjugatedComponent> {
    let mut out = Vec::new();
    for family in ring_families {
        if family.ring_indices.len() < 2 {
            continue;
        }
        out.extend(build_family_components(
            mol,
            rings,
            family,
            algo,
            all_ring_bonds,
        ));
    }
    out
}

fn build_family_components(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    family: &RingFamily,
    algo: AromaticityAlgorithm,
    all_ring_bonds: &FxHashSet<BondIdx>,
) -> Vec<ConjugatedComponent> {
    let family_component = ConjugatedComponent {
        atoms: family.atoms.clone(),
        bonds: Vec::new(),
        source_rings: family.ring_indices.clone(),
    };
    let eligible = family_eligibility(
        mol,
        rings,
        family,
        &family_component,
        algo,
        all_ring_bonds,
    );
    let mut graph = FamilyConjugation::new(family, eligible);
    graph.connect(mol);
    graph
        .groups()
        .into_values()
        .filter_map(|atoms| family_group_component(mol, rings, family, &graph.bonds, atoms))
        .collect()
}

fn family_eligibility(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    family: &RingFamily,
    component: &ConjugatedComponent,
    algo: AromaticityAlgorithm,
    all_ring_bonds: &FxHashSet<BondIdx>,
) -> FxHashMap<AtomIdx, bool> {
    family
        .atoms
        .iter()
        .map(|&atom| {
            let decision = evaluate_atom_via_home_ring(
                mol,
                atom,
                component,
                rings,
                algo,
                all_ring_bonds,
            );
            (atom, decision.electrons().is_some())
        })
        .collect()
}

struct FamilyConjugation {
    atoms: Vec<AtomIdx>,
    atom_set: FxHashSet<AtomIdx>,
    index_of: FxHashMap<AtomIdx, usize>,
    eligible: FxHashMap<AtomIdx, bool>,
    parent: Vec<usize>,
    bonds: Vec<BondIdx>,
}

impl FamilyConjugation {
    fn new(family: &RingFamily, eligible: FxHashMap<AtomIdx, bool>) -> Self {
        let atoms = family.atoms.clone();
        Self {
            atom_set: atoms.iter().copied().collect(),
            index_of: atoms.iter().enumerate().map(|(i, &a)| (a, i)).collect(),
            parent: (0..atoms.len()).collect(),
            atoms,
            eligible,
            bonds: Vec::new(),
        }
    }

    fn root(&mut self, index: usize) -> usize {
        if self.parent[index] != index {
            self.parent[index] = self.root(self.parent[index]);
        }
        self.parent[index]
    }

    fn union(&mut self, left: usize, right: usize) {
        let (left, right) = (self.root(left), self.root(right));
        if left != right {
            self.parent[left] = right;
        }
    }

    fn connect(&mut self, mol: &Molecule) {
        for atom in self.atoms.clone() {
            if !self.eligible[&atom] {
                continue;
            }
            for (neighbor, bond) in mol.neighbors(atom) {
                if !self.atom_set.contains(&neighbor)
                    || !self.eligible.get(&neighbor).copied().unwrap_or(false)
                {
                    continue;
                }
                self.union(self.index_of[&atom], self.index_of[&neighbor]);
                self.bonds.push(bond);
            }
        }
    }

    fn groups(&mut self) -> FxHashMap<usize, Vec<AtomIdx>> {
        let mut groups = FxHashMap::default();
        for atom in self.atoms.clone() {
            if self.eligible[&atom] {
                let root = self.root(self.index_of[&atom]);
                groups.entry(root).or_insert_with(Vec::new).push(atom);
            }
        }
        groups
    }
}

fn family_group_component(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    family: &RingFamily,
    bonds: &[BondIdx],
    atoms: Vec<AtomIdx>,
) -> Option<ConjugatedComponent> {
    let atom_set: FxHashSet<AtomIdx> = atoms.iter().copied().collect();
    let source_rings: Vec<usize> = family
        .ring_indices
        .iter()
        .copied()
        .filter(|&ring| rings[ring].iter().all(|atom| atom_set.contains(atom)))
        .collect();
    (source_rings.len() >= 2).then(|| ConjugatedComponent {
        atoms,
        bonds: bonds
            .iter()
            .copied()
            .filter(|&bond| {
                let bond = mol.bond(bond);
                atom_set.contains(&bond.atom1) && atom_set.contains(&bond.atom2)
            })
            .collect(),
        source_rings,
    })
}

/// Test/diagnostic-only exhaustive-candidate reference oracle
/// (Aromaticity-A1-1a) — **not** used by production or by
/// `trace_ring_pi_electrons`. Evaluates every SSSR/augmented ring AND every
/// multi-ring fused-envelope candidate from `build_conjugated_components`,
/// marking an atom/bond aromatic if ANY candidate containing it
/// independently satisfies 4n+2 via `evaluate_atom_pi_contribution`'s
/// per-atom rules — every candidate is evaluated from a clean slate, with NO
/// `aromatic_context` bootstrapping at all (unlike `assign_aromaticity_ex`'s
/// production Pass 1/Pass 2). Exists to cross-check hypotheses about which
/// per-atom rule needs to change, per the MANCUDE-style bounded-enumeration
/// precedent — see `docs/rfcs/aromaticity_a1_rfc.md`'s A1-1a section.
/// Deliberately simple/slow: O(rings + fused envelopes) candidates, no
/// attempt at Pass-2-style iteration, memoization, or performance tuning.
pub fn exhaustive_aromaticity_oracle(
    mol: &Molecule,
    algo: AromaticityAlgorithm,
) -> (FxHashSet<AtomIdx>, FxHashSet<BondIdx>) {
    let sssr = find_sssr(mol);
    let rings = augmented_ring_set(mol, sssr.rings());
    let families = crate::ring_family::find_ring_families_over(mol, &rings);
    let all_ring_bonds: FxHashSet<BondIdx> =
        rings.iter().flat_map(|r| ring_bond_set(mol, r)).collect();

    let mut candidates: Vec<ConjugatedComponent> = rings
        .iter()
        .enumerate()
        .map(|(i, r)| ConjugatedComponent::from_ring(r, i))
        .collect();
    candidates.extend(build_conjugated_components(
        mol,
        &rings,
        &families,
        algo,
        &all_ring_bonds,
    ));

    let mut aromatic_atoms: FxHashSet<AtomIdx> = FxHashSet::default();
    let mut aromatic_bonds: FxHashSet<BondIdx> = FxHashSet::default();

    for candidate in &candidates {
        let mut total: Option<u32> = Some(0);
        for &atom_idx in &candidate.atoms {
            // Multi-ring candidates evaluate each atom against its home ring
            // (see `evaluate_atom_via_home_ring`'s doc comment); single-ring
            // candidates fall through to the same code path with exactly one
            // source ring, unchanged from evaluating against `candidate` directly.
            let decision = evaluate_atom_via_home_ring(
                mol,
                atom_idx,
                candidate,
                &rings,
                algo,
                &all_ring_bonds,
            );
            total = match (total, decision.electrons()) {
                (Some(t), Some(e)) => Some(t + e as u32),
                _ => None,
            };
        }
        let Some(pi) = total else { continue };
        let (cls, _) = classify_ring_aromaticity(pi);
        if !matches!(cls, RingAromaticity::Aromatic) {
            continue;
        }
        for &a in &candidate.atoms {
            aromatic_atoms.insert(a);
        }
        for &a in &candidate.atoms {
            for (nb, bidx) in mol.neighbors(a) {
                if candidate.atoms.contains(&nb)
                    && matches!(
                        mol.bond(bidx).order,
                        BondOrder::Double | BondOrder::Aromatic
                    )
                {
                    aromatic_bonds.insert(bidx);
                }
            }
        }
    }

    (aromatic_atoms, aromatic_bonds)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------
