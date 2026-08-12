//! MANCUDE (maximum non-cumulated double bonds) ring treatment: production
//! representation plus a reference oracle.
//!
//! # Production: [`MancudeContext`]
//!
//! `assign_cip_accurate_experimental`'s digraph construction needs a Kekulé-invariant
//! fractional atomic number for ring-duplicate nodes in a resonant heteroaromatic system
//! (e.g. pyridine's ring carbons adjacent to N are "6½", not fixed at 6 or 7 depending on
//! which Kekulé form a kekulizer happened to pick). [`MancudeContext::compute`] is a
//! direct, source-verified port of RDKit's actual algorithm
//! (`Code/GraphMol/CIPLabeler/Mancude.cpp`'s `SeedTypes` / `RelaxTypes` / `VisitParts` /
//! `calcFracAtomNums`): classify each ring atom by element + formal charge + bond-order
//! pattern into a small set of seed types, demote any seed with fewer than 2 typed
//! neighbors (resonance needs connectivity), flood-fill the survivors into connected
//! resonance "parts" **along ring bonds only**, then compute each typed atom's fractional
//! atomic number as the mean atomic number of its own same-part ring neighbors (one hop).
//! Linear in atom/bond count -- no combinatorial search, no budget type needed.
//!
//! **Ring-bond membership matters in two places, not one.** A first draft of this port
//! (design/verification work, not shipped) gated seed-type candidacy on "has a ring bond"
//! but let `VisitParts`' flood-fill follow *any* bond to a typed neighbor. That's wrong:
//! two separate aromatic rings joined by one exocyclic single bond (e.g. biphenyl-shaped
//! molecules) would incorrectly merge into a single resonance part through that connecting
//! bond. RDKit's own `VisitPart` explicitly skips non-ring bonds
//! (`if (!mol.isInRing(bond)) continue;`) -- this module does the same, via a local
//! bridge-edge (cut-edge) detection in `ring_bond_set` (an edge lies on some ring iff it is
//! not a bridge -- a standard graph-theory fact, independent of which particular cycle
//! basis a ring-perception algorithm would choose). Originally implemented via
//! `chematic_perception::find_sssr`, which answers a strictly harder question (an actual
//! minimum cycle basis, with real ring membership) than this call site needs (a boolean
//! per bond) -- profiling `prepare_kekule_form` against the full corpus found SSSR
//! dominating its cost entirely (tens of milliseconds on large multi-ring molecules like
//! oligopeptides, vs microseconds for `kekulize` itself), 100% new cost relative to the
//! pre-Milestone-3B-1b engine. Replaced with an O(V+E) DFS bridge search; see the
//! Milestone 3B closeout entry in `docs/rfcs/cip_accurate_rfc.md` for the measurement.
//!
//! **The owner, not the represented atom, carries the fraction.** A
//! [`crate::node::CipNodeKind::MultipleBondDuplicate`] node sitting in atom `i`'s own
//! substituent list represents "whichever atom `i` is double-bonded to" -- and
//! `fractional_atomic_number(i)` (the mean of `i`'s own same-part neighbors) already *is*
//! the resonance-averaged answer to that question. So a duplicate's fractional value must
//! be looked up by its `source_atom` (the owner), never `duplicated_atom` (the specific
//! partner one particular Kekulé form happened to pick) -- see `docs/rfcs/cip_accurate_rfc.md`'s
//! Milestone 3B-1a entry for the full argument (cross-form Kekulé-invariance alone does
//! *not* discriminate this from the wrong reading on every fixture tested; the guard is a
//! concrete asserted value, e.g. quinoline's N-adjacent ring carbon must be exactly 13/2,
//! not 6/1).
//!
//! Real `Atom` nodes always keep their plain integer atomic number -- MANCUDE only ever
//! touches duplicate nodes, never real atoms (confirmed against RDKit: applying a
//! fractional value to a real nitrogen node would make it compare as carbon-valued,
//! silently breaking ordinary, non-resonance-related ranking).
//!
//! Charged seed types (`Nv4D3Plus`, `Nv2D2Minus`, `Cv3D3Minus`, `Ov3D2Plus` in RDKit's
//! naming, plus a secondary charge-relocation fraction pass) are **not implemented** --
//! checked directly against the frozen corpus (`validation/cip_label_corpus.jsonl`): 0 of
//! the 98 MANCUDE-scope cases contain a charged aromatic ring atom. A charged atom that
//! would otherwise seed simply never types, falling back to today's existing (already
//! correct for non-resonant atoms) integer atomic number -- a visible, safe fallback, not
//! a silent wrong average.
//!
//! # Reference model (test-only): [`enumerate_kekule_matchings`] / [`effective_atomic_number`]
//!
//! `chematic_core::kekulization::kekulize` returns *one* canonical Kekulé structure per
//! molecule (a single maximum matching). [`enumerate_kekule_matchings`] instead enumerates
//! *every* valid perfect matching by exhaustive backtracking, reusing the exact same
//! must-match/lone-pair-donor classification `kekulize()` itself uses (via
//! [`chematic_core::kekulization::atom_must_be_matched`], widened to `pub` for this
//! purpose) so the two can never silently disagree about what counts as a valid
//! placement. [`effective_atomic_number`] then computes, per atom, the mean atomic number
//! of whichever neighbor it's double-bonded to *across every one of those matchings* --
//! a literal reading of IUPAC P-92.1.4.4's "mean atomic number over valid Kekulé
//! placements" wording, at the whole-ring-system level rather than RDKit's one-hop
//! same-part-neighbor mean.
//!
//! **These two formulas are not the same function, and the difference is real, not a
//! bug.** Verified this session (Milestone 3B-1a) on a battery of fused/monocyclic/
//! multi-component heteroaromatics: they agree everywhere monocyclic, hydrocarbon-only,
//! or where a heteroatom never seeds at all (furan/thiophene/pyrrole's O/S/NH), and on
//! molecules with multiple separate MANCUDE components joined by a non-ring bond
//! (phenylpyridine, bipyridine). They **diverge** specifically on *fused* systems whose
//! connected resonance part spans the ring fusion and includes a heteroatom seed --
//! quinoline's N-adjacent ring carbon is 13/2 (6.5) under RDKit's one-hop mean but 19/3
//! (≈6.333) under this module's global-Kekulé-enumeration mean; the same pattern repeats
//! in isoquinoline, quinoxaline, and quinazoline. Since the frozen corpus's `modern`
//! column is computed by RDKit (`rdCIPLabeler`) itself, [`MancudeContext`] targets RDKit's
//! formula for production correctness. This does **not** make the enumeration-based
//! formula here "wrong" -- both are source-verified, valid readings of the same IUPAC
//! rule that happen to coincide on simpler systems; this module stays a genuinely useful
//! independent cross-check (exact on monocyclic/hydrocarbon cases), a fused-system
//! divergence detector, and a documented data point for anyone revisiting IUPAC
//! conformance later, not a discredited approach.
//!
//! **Deliberately not a production algorithm.** Full enumeration of every perfect
//! matching is combinatorially explosive on large fused ring systems; [`MancudeBudget`]
//! bounds it and *errors* rather than silently truncating when exceeded (same discipline
//! as [`crate::budget::CipBudget`]).

use std::collections::HashMap;
use std::collections::HashSet;

use chematic_core::kekulization::{KekuleResult, atom_must_be_matched, build_kekule_result};
use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule, implicit_hcount};

use crate::rational::RationalAtomicNumber;

/// A connected resonance ("part") within a molecule's MANCUDE ring system(s). Two atoms
/// share a `MancudeComponentId` iff they're reachable from each other through typed ring
/// atoms connected by ring bonds only -- two separate rings joined by a non-ring bond
/// (e.g. biphenyl's connecting bond) always get distinct ids.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct MancudeComponentId(u32);

/// Per-molecule MANCUDE state: which atoms belong to a resonance component, and each
/// typed atom's RDKit-compatible fractional atomic number. Compute this **once per
/// molecule** and share it across every stereocenter's digraph construction -- it is a
/// whole-molecule quantity, never a per-subtree one.
#[derive(Debug, Clone)]
pub struct MancudeContext {
    fractional_atomic_numbers: Vec<Option<RationalAtomicNumber>>,
    component_ids: Vec<Option<MancudeComponentId>>,
}

impl MancudeContext {
    /// `kekule_mol` must already be in explicit Kekulé form (`Single`/`Double` bond
    /// orders -- e.g. `chematic_core::kekulization::apply_kekule`'s output). Ring-atom
    /// candidacy is seeded from `kekule_mol`'s own bond orders (which requires an explicit
    /// Kekulé form to read at all), **not** gated on `kekule_mol.atom(..).aromatic` -- see
    /// this impl's docs above for why: RDKit's own `SeedTypes` has no aromaticity check
    /// either, and a non-aromatic ring double bond matching the same valence pattern (a
    /// lactone/enone carbon in a fused non-aromatic ring, common in this project's corpus)
    /// is typed identically. An earlier draft asserted "every typed atom must be aromatic"
    /// as a self-check; a full 98-case corpus sweep (`tests/mancude_context.rs`) falsified
    /// it on real molecules, so the check was removed, not the typing. What that sweep
    /// establishes is narrower than "this is correct": the typing is **invariant** (both
    /// Kekulé-form and renumbering sweeps stayed 98/98 clean) and **inert this round** (no
    /// comparator reads this field yet) -- it does not establish that typing these
    /// non-aromatic atoms is *label*-correct, only that doing so structurally mirrors
    /// RDKit's real algorithm shape. That remaining question is Milestone 3B-1b's
    /// corpus-agreement gate to answer, not this round's -- watch specifically for a
    /// non-aromatic typed atom with a *heteroatom* same-part neighbor, since a pure-carbon
    /// one is a no-op today (`Rational(6/1)` compares equal to `Integral(6)`).
    pub fn compute(kekule_mol: &Molecule) -> Self {
        let ring_bonds = ring_bond_set(kekule_mol);
        let mut types = seed_types(kekule_mol, &ring_bonds);
        relax_types(kekule_mol, &mut types);

        let parts = visit_parts(kekule_mol, &types, &ring_bonds);
        let n = kekule_mol.atom_count();
        let mut fractional_atomic_numbers = vec![None; n];
        let mut component_ids = vec![None; n];
        for (idx, _) in kekule_mol.atoms() {
            let i = idx.0 as usize;
            if parts[i] == 0 {
                continue;
            }
            let same_part_neighbors: Vec<u32> = kekule_mol
                .neighbors(idx)
                .filter(|&(nb, _)| parts[nb.0 as usize] == parts[i])
                .map(|(nb, _)| kekule_mol.atom(nb).element.atomic_number() as u32)
                .collect();
            if same_part_neighbors.is_empty() {
                // Typed (survived RelaxTypes' "≥2 typed neighbors, any bond" criterion)
                // but has zero *same-part* (ring-bond-connected) neighbors -- possible
                // only in unusual topologies where the counted neighbors were reached via
                // non-ring bonds. Not a real resonance system; fall back to no fraction
                // rather than risk `RationalAtomicNumber::mean`'s empty-slice panic.
                continue;
            }
            component_ids[i] = Some(MancudeComponentId(parts[i]));
            fractional_atomic_numbers[i] = Some(RationalAtomicNumber::mean(&same_part_neighbors));
        }

        Self {
            fractional_atomic_numbers,
            component_ids,
        }
    }

    /// The RDKit-compatible fractional atomic number for `atom`, or `None` if `atom` is
    /// outside any MANCUDE resonance component.
    pub fn fractional_atomic_number(&self, atom: AtomIdx) -> Option<RationalAtomicNumber> {
        self.fractional_atomic_numbers[atom.0 as usize]
    }

    /// Which resonance component `atom` belongs to, or `None` if it's outside any.
    pub fn component_id(&self, atom: AtomIdx) -> Option<MancudeComponentId> {
        self.component_ids[atom.0 as usize]
    }
}

/// Builds `mol`'s Kekulé-form clone (`AtomIdx`-preserving, per
/// `chematic_core::kekulization::apply_kekule`) and computes its [`MancudeContext`] in one
/// step -- the two are always needed together (`MancudeContext::compute` requires an
/// explicit-Kekulé input) and should be computed **once per molecule**, then shared across
/// every stereocenter's [`crate::digraph::CipDigraph::new_with_mancude`] call, never
/// recomputed per atom or per subtree expansion.
///
/// **Not yet called by [`crate::assign_cip_accurate_experimental`]** (see this module's
/// docs on why switching that entry point's digraph input from aromatic to Kekulé form is
/// deferred to Milestone 3B-1b, alongside comparator wiring: doing so here, before the
/// comparator knows how to weigh the fractional values, would silently change which
/// `MultipleBondDuplicate` nodes exist for *every* aromatic-adjacent stereocenter -- far
/// more than the 98 MANCUDE-labeled corpus cases -- an uncontrolled behavior change this
/// milestone's "byte-identical corpus report" gate forbids).
pub fn prepare_kekule_form(
    mol: &Molecule,
) -> Result<(Molecule, MancudeContext), chematic_core::kekulization::KekuleError> {
    let kekule = chematic_core::kekulization::kekulize(mol)?;
    let kekule_mol = chematic_core::kekulization::apply_kekule(mol, &kekule);
    let ctx = MancudeContext::compute(&kekule_mol);
    Ok((kekule_mol, ctx))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SeedType {
    Cv4D3,
    Nv3D2,
    Other,
}

/// Bond orders SSSR/ring-perception treats as capable of forming a ring at all --
/// mirrors `chematic_perception::sssr`'s private `is_ring_eligible` (kept in sync by
/// inspection; duplicated rather than exposed as new public API for one boolean check).
/// Zero-order/Dative/Query bonds don't participate in ring topology (RDKit PR #9118).
fn is_ring_eligible(order: BondOrder) -> bool {
    matches!(
        order,
        BondOrder::Single
            | BondOrder::Double
            | BondOrder::Triple
            | BondOrder::Quadruple
            | BondOrder::Aromatic
            | BondOrder::Up
            | BondOrder::Down
    )
}

/// A bond is a ring bond iff it lies on some cycle -- equivalently, iff it is *not* a
/// bridge (cut edge) of the molecular graph, a standard graph-theory fact independent of
/// which particular cycle basis a full ring-perception algorithm would pick. Computed via
/// a single O(V+E) DFS (bridge-finding / low-link, Tarjan's algorithm), not
/// `chematic_perception::find_sssr` (a full minimum cycle basis -- see this module's docs
/// for why that was a real, measured perf cost this call site never needed). Topology-only
/// (independent of bond order beyond ring-eligibility), so it's safe to compute on either
/// the aromatic or Kekulé-respelled form of the same molecule.
fn ring_bond_set(mol: &Molecule) -> HashSet<(u32, u32)> {
    let n = mol.atom_count();
    let mut disc: Vec<Option<u32>> = vec![None; n];
    let mut low: Vec<u32> = vec![0; n];
    let mut timer: u32 = 0;
    let mut ring_bonds: HashSet<(u32, u32)> = HashSet::new();

    // Iterative DFS (avoids recursion-depth concerns on long unbranched chains). Each
    // frame tracks the atom being visited, the bond it was reached through (skipped when
    // re-scanning neighbors, so a genuine 2-atom cycle via a *different* bond between the
    // same pair is still detected), and where to resume its neighbor list.
    struct Frame {
        atom: AtomIdx,
        parent_atom: Option<AtomIdx>,
        arrival_bond: Option<BondIdx>,
        neighbors: Vec<(AtomIdx, BondIdx)>,
        next: usize,
    }

    for start in 0..n {
        if disc[start].is_some() {
            continue;
        }
        disc[start] = Some(timer);
        low[start] = timer;
        timer += 1;
        let start_atom = AtomIdx(start as u32);
        let mut stack = vec![Frame {
            atom: start_atom,
            parent_atom: None,
            arrival_bond: None,
            neighbors: mol
                .neighbors(start_atom)
                .filter(|&(_, bidx)| is_ring_eligible(mol.bond(bidx).order))
                .collect(),
            next: 0,
        }];

        while let Some(frame) = stack.last_mut() {
            if frame.next >= frame.neighbors.len() {
                let finished = stack.pop().unwrap();
                if let Some(parent) = finished.parent_atom {
                    let pi = parent.0 as usize;
                    let ui = finished.atom.0 as usize;
                    low[pi] = low[pi].min(low[ui]);
                    // Tree edge (parent, finished.atom) is a ring bond iff it is NOT a
                    // bridge: low[child] <= disc[parent] means the child's subtree can
                    // reach back to `parent` (or higher), so this edge lies on a cycle.
                    if low[ui] <= disc[pi].unwrap() {
                        let a = parent.0.min(finished.atom.0);
                        let b = parent.0.max(finished.atom.0);
                        ring_bonds.insert((a, b));
                    }
                }
                continue;
            }
            let (v, bidx) = frame.neighbors[frame.next];
            frame.next += 1;
            if Some(bidx) == frame.arrival_bond {
                continue;
            }
            let u = frame.atom;
            let vi = v.0 as usize;
            if let Some(vdisc) = disc[vi] {
                // Back edge: undirected-graph DFS has no cross edges, so any already-
                // discovered neighbor (other than the one we arrived through) is
                // necessarily an ancestor on the current path -- always a real cycle.
                let ui = u.0 as usize;
                low[ui] = low[ui].min(vdisc);
                let a = u.0.min(v.0);
                let b = u.0.max(v.0);
                ring_bonds.insert((a, b));
            } else {
                disc[vi] = Some(timer);
                low[vi] = timer;
                timer += 1;
                stack.push(Frame {
                    atom: v,
                    parent_atom: Some(u),
                    arrival_bond: Some(bidx),
                    neighbors: mol
                        .neighbors(v)
                        .filter(|&(_, b)| is_ring_eligible(mol.bond(b).order))
                        .collect(),
                    next: 0,
                });
            }
        }
    }
    ring_bonds
}

fn is_ring_bond(ring_bonds: &HashSet<(u32, u32)>, a: AtomIdx, b: AtomIdx) -> bool {
    ring_bonds.contains(&(a.0.min(b.0), a.0.max(b.0)))
}

/// Direct port of RDKit's `Mancude.cpp` `SeedTypes` (neutral types only -- see module
/// docs' charged-types non-goal). Classifies each ring atom by element, formal charge,
/// and its bond-order pattern (single/double/other counts, seeded from
/// `implicit_hcount` exactly like RDKit's `getTotalNumHs()`-seeded `btypes`, then
/// accumulated over *every* bond, ring and non-ring alike -- an ipso carbon's exocyclic
/// single bond counts toward its own pattern). `ring` requires at least one of the
/// atom's bonds to be a ring bond; an atom with any "other"-multiplicity bond (e.g.
/// triple) never seeds.
fn seed_types(mol: &Molecule, ring_bonds: &HashSet<(u32, u32)>) -> Vec<SeedType> {
    let n = mol.atom_count();
    let mut types = vec![SeedType::Other; n];
    for (idx, atom) in mol.atoms() {
        let mut singles = implicit_hcount(mol, idx) as u32;
        let mut doubles = 0u32;
        let mut other = 0u32;
        let mut in_ring = false;
        for (nb, bidx) in mol.neighbors(idx) {
            match mol.bond(bidx).order {
                BondOrder::Single => singles += 1,
                BondOrder::Double => doubles += 1,
                _ => other += 1,
            }
            if is_ring_bond(ring_bonds, idx, nb) {
                in_ring = true;
            }
        }
        if !in_ring || other != 0 {
            continue;
        }
        let z = atom.element.atomic_number();
        let q = atom.charge;
        match z {
            6 | 14 | 32 if q == 0 && doubles == 1 && singles == 2 => {
                types[idx.0 as usize] = SeedType::Cv4D3;
            }
            7 | 15 | 33 if q == 0 && doubles == 1 && singles == 1 => {
                types[idx.0 as usize] = SeedType::Nv3D2;
            }
            _ => {}
        }
    }
    types
}

/// Direct port of RDKit's `RelaxTypes`: iterative demotion of any typed atom with fewer
/// than 2 typed neighbors (resonance needs connectivity), cascading -- demoting one atom
/// can push a now-under-connected neighbor below the threshold too. Deliberately
/// unfiltered by ring-bond membership here, matching RDKit's own `RelaxTypes` exactly
/// (only `VisitPart`'s flood-fill filters by ring bond -- see module docs).
fn relax_types(mol: &Molecule, types: &mut [SeedType]) {
    let n = types.len();
    let mut counts = vec![0i32; n];
    let mut queue: Vec<AtomIdx> = Vec::new();
    for (idx, _) in mol.atoms() {
        let i = idx.0 as usize;
        for (nb, _) in mol.neighbors(idx) {
            if types[nb.0 as usize] != SeedType::Other {
                counts[i] += 1;
            }
        }
        if counts[i] == 1 {
            queue.push(idx);
        }
    }
    let mut qi = 0;
    while qi < queue.len() {
        let idx = queue[qi];
        qi += 1;
        let i = idx.0 as usize;
        if types[i] != SeedType::Other {
            types[i] = SeedType::Other;
            for (nb, _) in mol.neighbors(idx) {
                let j = nb.0 as usize;
                counts[j] -= 1;
                if counts[j] == 1 {
                    queue.push(nb);
                }
            }
        }
    }
}

/// Direct port of RDKit's `VisitParts`/`VisitPart`: flood-fills connected typed atoms
/// into numbered resonance parts, following **ring bonds only** -- the fix (see module
/// docs) that keeps two separate rings joined by a non-ring bond from merging into one
/// part. Returns one part number (1-based; 0 = untyped) per atom.
fn visit_parts(mol: &Molecule, types: &[SeedType], ring_bonds: &HashSet<(u32, u32)>) -> Vec<u32> {
    let n = types.len();
    let mut parts = vec![0u32; n];
    let mut numparts = 0u32;
    for (idx, _) in mol.atoms() {
        let i = idx.0 as usize;
        if parts[i] == 0 && types[i] != SeedType::Other {
            numparts += 1;
            let mut stack = vec![idx];
            parts[i] = numparts;
            while let Some(cur) = stack.pop() {
                for (nb, _) in mol.neighbors(cur) {
                    let j = nb.0 as usize;
                    if parts[j] == 0
                        && types[j] != SeedType::Other
                        && is_ring_bond(ring_bonds, cur, nb)
                    {
                        parts[j] = numparts;
                        stack.push(nb);
                    }
                }
            }
        }
    }
    parts
}

/// Bounds on [`enumerate_kekule_matchings`]'s search, to keep it a small-ring oracle
/// rather than an accidental production path. Both bounds *error* the whole call rather
/// than returning a truncated/partial result -- an incomplete enumeration would silently
/// corrupt any mean computed from it.
#[derive(Debug, Clone, Copy)]
pub struct MancudeBudget {
    /// Maximum number of must-match atoms in one connected component.
    pub max_atoms: usize,
    /// Maximum number of complete matchings collected.
    pub max_matchings: usize,
    /// Maximum number of backtracking search steps (guards pathological search trees
    /// that explore heavily before finding few, or zero, complete matchings).
    pub max_search_steps: usize,
}

impl Default for MancudeBudget {
    fn default() -> Self {
        Self {
            max_atoms: 24,
            max_matchings: 64,
            max_search_steps: 100_000,
        }
    }
}

/// Errors from [`enumerate_kekule_matchings`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MancudeError {
    /// The must-match component is larger than [`MancudeBudget::max_atoms`].
    TooManyAtoms { count: usize, max: usize },
    /// More valid matchings exist than [`MancudeBudget::max_matchings`] allows.
    TooManyMatchings { max: usize },
    /// The backtracking search exceeded [`MancudeBudget::max_search_steps`].
    SearchBudgetExceeded { max: usize },
}

impl core::fmt::Display for MancudeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            MancudeError::TooManyAtoms { count, max } => {
                write!(
                    f,
                    "mancude component has {count} must-match atoms, over budget {max}"
                )
            }
            MancudeError::TooManyMatchings { max } => {
                write!(
                    f,
                    "mancude component has more than {max} valid Kekulé matchings"
                )
            }
            MancudeError::SearchBudgetExceeded { max } => {
                write!(f, "mancude matching search exceeded {max} steps")
            }
        }
    }
}

impl std::error::Error for MancudeError {}

/// Enumerate *every* valid perfect matching of `mol`'s aromatic must-match subgraph,
/// bounded by `budget`. Returns one [`KekuleResult`] per valid matching, in the exact
/// same format `chematic_core::kekulization::kekulize` produces (every aromatic bond
/// mapped to `Single` or `Double`) -- so `kekulize(mol)`'s own single result is always
/// expected to appear as a member of this function's output (checked in this module's
/// own tests).
///
/// A molecule with no aromatic bonds has exactly one (empty) valid placement, mirroring
/// `kekulize()`'s own no-op convention.
pub fn enumerate_kekule_matchings(
    mol: &Molecule,
    budget: MancudeBudget,
) -> Result<Vec<KekuleResult>, MancudeError> {
    let mut aromatic_bonds: Vec<BondIdx> = Vec::new();
    let mut aromatic_atoms: Vec<AtomIdx> = Vec::new();
    for (bidx, bond) in mol.bonds() {
        if bond.order == BondOrder::Aromatic {
            aromatic_bonds.push(bidx);
            if !aromatic_atoms.contains(&bond.atom1) {
                aromatic_atoms.push(bond.atom1);
            }
            if !aromatic_atoms.contains(&bond.atom2) {
                aromatic_atoms.push(bond.atom2);
            }
        }
    }
    if aromatic_bonds.is_empty() {
        return Ok(vec![KekuleResult::new()]);
    }

    let mut must_match: Vec<AtomIdx> = aromatic_atoms
        .into_iter()
        .filter(|&idx| atom_must_be_matched(mol, idx))
        .collect();
    must_match.sort();

    if must_match.len() > budget.max_atoms {
        return Err(MancudeError::TooManyAtoms {
            count: must_match.len(),
            max: budget.max_atoms,
        });
    }

    let mut adj: HashMap<AtomIdx, Vec<AtomIdx>> = HashMap::new();
    for &bidx in &aromatic_bonds {
        let bond = mol.bond(bidx);
        if must_match.contains(&bond.atom1) && must_match.contains(&bond.atom2) {
            adj.entry(bond.atom1).or_default().push(bond.atom2);
            adj.entry(bond.atom2).or_default().push(bond.atom1);
        }
    }

    let mut results: Vec<HashMap<AtomIdx, AtomIdx>> = Vec::new();
    let mut current: HashMap<AtomIdx, AtomIdx> = HashMap::new();
    let mut steps = 0usize;
    backtrack(
        &must_match,
        &adj,
        &mut current,
        &mut results,
        &budget,
        &mut steps,
    )?;

    Ok(results
        .into_iter()
        .map(|matching| build_kekule_result(&aromatic_bonds, mol, &matching))
        .collect())
}

/// Enumerate every perfect matching of `remaining` under `adj` by always extending the
/// smallest still-unmatched atom -- this visits each complete matching exactly once (the
/// atom's partner in any given valid matching is unique, so branching over its candidate
/// partners partitions the search space without overlap or omission).
fn backtrack(
    remaining: &[AtomIdx],
    adj: &HashMap<AtomIdx, Vec<AtomIdx>>,
    current: &mut HashMap<AtomIdx, AtomIdx>,
    results: &mut Vec<HashMap<AtomIdx, AtomIdx>>,
    budget: &MancudeBudget,
    steps: &mut usize,
) -> Result<(), MancudeError> {
    *steps += 1;
    if *steps > budget.max_search_steps {
        return Err(MancudeError::SearchBudgetExceeded {
            max: budget.max_search_steps,
        });
    }

    let still_free: Vec<AtomIdx> = remaining
        .iter()
        .copied()
        .filter(|a| !current.contains_key(a))
        .collect();
    let Some(&v) = still_free.first() else {
        if results.len() >= budget.max_matchings {
            return Err(MancudeError::TooManyMatchings {
                max: budget.max_matchings,
            });
        }
        results.push(current.clone());
        return Ok(());
    };

    let candidates: Vec<AtomIdx> = adj
        .get(&v)
        .into_iter()
        .flatten()
        .copied()
        .filter(|u| !current.contains_key(u))
        .collect();
    for u in candidates {
        current.insert(v, u);
        current.insert(u, v);
        backtrack(remaining, adj, current, results, budget, steps)?;
        current.remove(&v);
        current.remove(&u);
    }
    Ok(())
}

/// For one atom, the mean atomic number of whichever neighbor it is double-bonded to,
/// across every matching in `matchings` -- the quantity Milestone 3B-1's mancude
/// duplicate node will store. `None` if `atom` is never double-bonded in any of the
/// matchings (a lone-pair donor, e.g. furan's O or pyrrole's `[nH]`, or an atom outside
/// any mancude system) -- such atoms contribute no extra duplicate at all, matching
/// `kekulize()`'s existing single/double-only bond model.
pub fn effective_atomic_number(
    mol: &Molecule,
    atom: AtomIdx,
    matchings: &[KekuleResult],
) -> Option<RationalAtomicNumber> {
    let mut partners = Vec::new();
    for matching in matchings {
        for (nb, bidx) in mol.neighbors(atom) {
            if matching.get(&bidx) == Some(&BondOrder::Double) {
                partners.push(mol.atom(nb).element.atomic_number() as u32);
                break;
            }
        }
    }
    if partners.is_empty() {
        None
    } else {
        Some(RationalAtomicNumber::mean(&partners))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_core::{Atom, Element, MoleculeBuilder};

    fn benzene() -> Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..6)
            .map(|_| b.add_atom(Atom::aromatic(Element::C)))
            .collect();
        for i in 0..6 {
            b.add_bond(atoms[i], atoms[(i + 1) % 6], BondOrder::Aromatic)
                .unwrap();
        }
        b.build()
    }

    /// N1=C2-C3=C4-C5=C6(-N1) around the ring; C2 and C6 are the carbons adjacent to N.
    fn pyridine() -> Molecule {
        let mut b = MoleculeBuilder::new();
        let n = b.add_atom(Atom::aromatic(Element::N));
        let cs: Vec<_> = (0..5)
            .map(|_| b.add_atom(Atom::aromatic(Element::C)))
            .collect();
        let atoms = [n, cs[0], cs[1], cs[2], cs[3], cs[4]];
        for i in 0..6 {
            b.add_bond(atoms[i], atoms[(i + 1) % 6], BondOrder::Aromatic)
                .unwrap();
        }
        b.build()
    }

    fn pyridine_n_and_adjacent_carbon() -> (Molecule, AtomIdx, AtomIdx) {
        let mol = pyridine();
        // atom 0 = N, atom 1 = the carbon adjacent to it (per construction above).
        (mol, AtomIdx(0), AtomIdx(1))
    }

    /// O-C1=C2-C3=C4(-O) around the ring; O is a lone-pair donor, excluded from must_match.
    fn furan() -> Molecule {
        let mut b = MoleculeBuilder::new();
        let o = b.add_atom(Atom::aromatic(Element::O));
        let cs: Vec<_> = (0..4)
            .map(|_| b.add_atom(Atom::aromatic(Element::C)))
            .collect();
        let atoms = [o, cs[0], cs[1], cs[2], cs[3]];
        for i in 0..5 {
            b.add_bond(atoms[i], atoms[(i + 1) % 5], BondOrder::Aromatic)
                .unwrap();
        }
        b.build()
    }

    #[test]
    fn no_aromatic_bonds_is_one_empty_matching() {
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        let mol = b.build();
        let all = enumerate_kekule_matchings(&mol, MancudeBudget::default()).unwrap();
        assert_eq!(all, vec![KekuleResult::new()]);
    }

    /// Consistency check: kekulize()'s own single result must be a member of the full
    /// enumeration, for every fixture -- if it weren't, the two would be using different
    /// notions of "valid matching" and the oracle couldn't be trusted to design against.
    #[test]
    fn kekulize_result_is_a_member_of_the_full_enumeration() {
        for mol in [benzene(), pyridine(), furan()] {
            let single = chematic_core::kekulization::kekulize(&mol).unwrap();
            let all = enumerate_kekule_matchings(&mol, MancudeBudget::default()).unwrap();
            assert!(
                all.contains(&single),
                "kekulize()'s own result must appear in the full enumeration"
            );
        }
    }

    #[test]
    fn benzene_has_exactly_two_kekule_forms() {
        let all = enumerate_kekule_matchings(&benzene(), MancudeBudget::default()).unwrap();
        assert_eq!(all.len(), 2);
    }

    #[test]
    fn furan_has_exactly_one_kekule_form() {
        // O is a lone-pair donor (excluded from must-match); the remaining 4 ring
        // carbons form a PATH (not a cycle, since neither C-O bond is a matching
        // candidate), and a 4-atom path has exactly one perfect matching: {C1-C2, C3-C4}.
        let all = enumerate_kekule_matchings(&furan(), MancudeBudget::default()).unwrap();
        assert_eq!(all.len(), 1);
    }

    /// Hand-derived IUPAC example, matching the design conversation's own "6½" value:
    /// pyridine's ring carbon adjacent to N is double-bonded to N in one Kekulé form and
    /// to its other (carbon) ring neighbor in the other -- effective atomic number
    /// (7 + 6) / 2 = 13/2 = 6½.
    #[test]
    fn pyridine_adjacent_carbon_is_six_and_a_half() {
        let (mol, _n, adjacent_c) = pyridine_n_and_adjacent_carbon();
        let all = enumerate_kekule_matchings(&mol, MancudeBudget::default()).unwrap();
        assert_eq!(all.len(), 2);
        let signature = effective_atomic_number(&mol, adjacent_c, &all).unwrap();
        assert_eq!(signature.numerator(), 13);
        assert_eq!(signature.denominator(), 2);
    }

    /// Furan's O is never double-bonded (lone-pair donor) -- no duplicate contribution.
    #[test]
    fn furan_oxygen_has_no_effective_atomic_number() {
        let mol = furan();
        let all = enumerate_kekule_matchings(&mol, MancudeBudget::default()).unwrap();
        assert_eq!(effective_atomic_number(&mol, AtomIdx(0), &all), None);
    }

    /// The common-signature property M3B-0 exists to demonstrate: two individually-valid,
    /// genuinely different Kekulé forms of the same molecule disagree on a given atom's
    /// immediate double-bond partner, yet both are members of the one enumeration whose
    /// mean is the single MANCUDE signature -- checked on a hydrocarbon fixture (trivial,
    /// can't fail: both forms agree since every ring atom is carbon) and a hetero fixture
    /// (genuinely fractional: the two forms disagree, and the signature averages them).
    #[test]
    fn kekule_form_a_and_b_share_one_common_signature_hydrocarbon() {
        let mol = benzene();
        let all = enumerate_kekule_matchings(&mol, MancudeBudget::default()).unwrap();
        let (form_a, form_b) = (&all[0], &all[1]);
        assert_ne!(
            form_a, form_b,
            "must be two genuinely different resonance structures"
        );

        let atom = AtomIdx(0);
        let partner_in = |form: &KekuleResult| {
            mol.neighbors(atom)
                .find(|&(_, bidx)| form.get(&bidx) == Some(&BondOrder::Double))
                .map(|(nb, _)| mol.atom(nb).element.atomic_number() as u32)
                .unwrap()
        };
        let (pa, pb) = (partner_in(form_a), partner_in(form_b));
        // Hydrocarbon: both partners are carbon, so the two forms happen to agree here --
        // that's expected, not a bug; the fraction only becomes visible with a heteroatom
        // (see the hetero variant of this test below).
        assert_eq!(pa, 6);
        assert_eq!(pb, 6);

        let signature = effective_atomic_number(&mol, atom, &all).unwrap();
        assert_eq!(signature, RationalAtomicNumber::mean(&[pa, pb]));
        assert_eq!(signature, RationalAtomicNumber::integer(6));
    }

    #[test]
    fn kekule_form_a_and_b_share_one_common_signature_hetero() {
        let (mol, _n, adjacent_c) = pyridine_n_and_adjacent_carbon();
        let all = enumerate_kekule_matchings(&mol, MancudeBudget::default()).unwrap();
        let (form_a, form_b) = (&all[0], &all[1]);
        assert_ne!(
            form_a, form_b,
            "must be two genuinely different resonance structures"
        );

        let partner_in = |form: &KekuleResult| {
            mol.neighbors(adjacent_c)
                .find(|&(_, bidx)| form.get(&bidx) == Some(&BondOrder::Double))
                .map(|(nb, _)| mol.atom(nb).element.atomic_number() as u32)
                .unwrap()
        };
        let (pa, pb) = (partner_in(form_a), partner_in(form_b));
        // The two forms genuinely disagree (one has this carbon double-bonded to N, the
        // other to its carbon neighbor) -- this divergence is exactly what a single-form
        // representation (today's digraph) can't average away, and what the common
        // signature below reconciles into one value.
        assert_ne!(
            pa, pb,
            "hetero fixture must show the two forms actually disagreeing"
        );

        let signature = effective_atomic_number(&mol, adjacent_c, &all).unwrap();
        assert_eq!(signature, RationalAtomicNumber::mean(&[pa, pb]));
        assert_eq!(signature.numerator(), 13);
        assert_eq!(signature.denominator(), 2);
    }

    // ---- MancudeContext (production path) --------------------------------------------

    fn kekule_clone(smiles: &str) -> Molecule {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let kekule = chematic_core::kekulization::kekulize(&mol).unwrap();
        chematic_core::kekulization::apply_kekule(&mol, &kekule)
    }

    #[test]
    fn context_pyridine_n_adjacent_carbon_is_six_and_a_half() {
        let kmol = kekule_clone("c1ccncc1");
        let ctx = MancudeContext::compute(&kmol);
        // atom3 = N (per c1ccncc1's atom order); atom2 and atom4 are its ring neighbors.
        let n = AtomIdx(3);
        assert_eq!(kmol.atom(n).element.atomic_number(), 7);
        for adjacent in [AtomIdx(2), AtomIdx(4)] {
            let f = ctx.fractional_atomic_number(adjacent).unwrap();
            assert_eq!(
                (f.numerator(), f.denominator()),
                (13, 2),
                "atom {adjacent:?}"
            );
        }
        // N itself is typed (Nv3D2) but its own *real* atom identity is untouched --
        // MancudeContext still records a fraction for it (mean of its 2 same-part ring
        // neighbors, both carbon = 6/1), but digraph wiring (deliverable 3) must never
        // apply this to N's own `Atom` node, only to duplicates whose owner is N.
        let n_fraction = ctx.fractional_atomic_number(n).unwrap();
        assert_eq!((n_fraction.numerator(), n_fraction.denominator()), (6, 1));
    }

    #[test]
    fn context_benzene_all_atoms_integer_six() {
        let kmol = kekule_clone("c1ccccc1");
        let ctx = MancudeContext::compute(&kmol);
        for i in 0..kmol.atom_count() {
            let f = ctx.fractional_atomic_number(AtomIdx(i as u32)).unwrap();
            assert_eq!((f.numerator(), f.denominator()), (6, 1), "atom {i}");
        }
    }

    #[test]
    fn context_furan_oxygen_never_types() {
        let kmol = kekule_clone("c1ccoc1");
        let ctx = MancudeContext::compute(&kmol);
        let o_idx = (0..kmol.atom_count())
            .map(|i| AtomIdx(i as u32))
            .find(|&idx| kmol.atom(idx).element == chematic_core::Element::O)
            .unwrap();
        assert_eq!(ctx.fractional_atomic_number(o_idx), None);
        assert_eq!(ctx.component_id(o_idx), None);
    }

    /// Hand-derived divergence-table fixture: quinoline's ring carbons directly bonded to
    /// N are 13/2 under RDKit's one-hop same-part mean -- NOT 19/3 or 20/3, which is what
    /// the global-Kekulé-enumeration oracle in this same module computes for the same
    /// atoms (see module docs). This is also the concrete-value assertion that guards the
    /// owner-vs-represented-atom design decision once digraph wiring lands (deliverable 3).
    #[test]
    fn context_quinoline_n_adjacent_carbons_are_six_and_a_half_not_the_oracle_value() {
        let kmol = kekule_clone("n1ccc2ccccc2c1");
        let ctx = MancudeContext::compute(&kmol);
        // atom0 = N; atom1 and atom9 are its two ring neighbors (both plain, non-fusion
        // ring carbons -- see the divergence table in module docs).
        for adjacent in [AtomIdx(1), AtomIdx(9)] {
            let f = ctx.fractional_atomic_number(adjacent).unwrap();
            assert_eq!(
                (f.numerator(), f.denominator()),
                (13, 2),
                "atom {adjacent:?} must be RDKit's 6.5, not the oracle's 6.333/6.667"
            );
        }
        // The two ring-fusion carbons (atom3, atom8) are NOT N-adjacent: all 3 of their
        // same-part ring neighbors are plain carbon, so they reduce to a plain integer
        // (6/1), same as the oracle here -- the divergence table's "isoquinoline atom3 =
        // 19/3" fusion-carbon value comes from isoquinoline's N being adjacent to ITS
        // fusion carbon, which is not the case in quinoline.
        for fusion in [AtomIdx(3), AtomIdx(8)] {
            let f = ctx.fractional_atomic_number(fusion).unwrap();
            assert_eq!((f.numerator(), f.denominator()), (6, 1), "atom {fusion:?}");
        }
        // Every other ring carbon is plain hydrocarbon-typed: 6/1.
        for plain in [AtomIdx(2), AtomIdx(4), AtomIdx(5), AtomIdx(6), AtomIdx(7)] {
            let f = ctx.fractional_atomic_number(plain).unwrap();
            assert_eq!((f.numerator(), f.denominator()), (6, 1), "atom {plain:?}");
        }
    }

    /// Regression test for the ring-bond-only flood-fill fix (see module docs): two
    /// separate pyridine/benzene-type rings joined by one exocyclic single bond must stay
    /// two distinct resonance components, not merge into one through that connecting bond.
    #[test]
    fn context_phenylpyridine_keeps_two_separate_components() {
        let kmol = kekule_clone("c1ccc(-c2ccccn2)cc1");
        let ctx = MancudeContext::compute(&kmol);
        // atom3 = benzene ring's ipso carbon; atom4 = pyridine ring's ipso carbon (per
        // this SMILES's atom order) -- connected by the exocyclic single bond.
        let (ipso_benzene, ipso_pyridine) = (AtomIdx(3), AtomIdx(4));
        let id_benzene = ctx.component_id(ipso_benzene).unwrap();
        let id_pyridine = ctx.component_id(ipso_pyridine).unwrap();
        assert_ne!(
            id_benzene, id_pyridine,
            "the two rings must NOT merge into one component through the connecting bond"
        );
        // The pyridine ipso carbon's 2 ring neighbors are 1 carbon + the ring nitrogen --
        // mean(6,7) = 13/2. If the exocyclic bond were incorrectly treated as same-part
        // (the bug this test guards against), a 3rd same-part neighbor (the benzene ipso
        // carbon, atomic number 6) would pull this to mean(6,6,7) = 19/3 instead.
        let f_pyridine = ctx.fractional_atomic_number(ipso_pyridine).unwrap();
        assert_eq!(
            (f_pyridine.numerator(), f_pyridine.denominator()),
            (13, 2),
            "must be 13/2 (2 same-part ring neighbors only), not 19/3 (3, if the bug regressed)"
        );
    }

    /// Negative control: an acyclic conjugated polyene has zero ring bonds, so
    /// `seed_types` must never fire -- confirms MANCUDE detection doesn't key off
    /// conjugation/alternation alone, only ring membership + the full type pattern.
    #[test]
    fn context_acyclic_polyene_is_not_misdetected_as_mancude() {
        let mut b = MoleculeBuilder::new();
        // hexatriene: C1=C2-C3=C4-C5=C6, all non-aromatic, no ring closure.
        let atoms: Vec<_> = (0..6).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        for i in 0..5 {
            let order = if i % 2 == 0 {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
            b.add_bond(atoms[i], atoms[i + 1], order).unwrap();
        }
        let mol = b.build();
        let ctx = MancudeContext::compute(&mol);
        for i in 0..mol.atom_count() {
            let idx = AtomIdx(i as u32);
            assert_eq!(ctx.fractional_atomic_number(idx), None, "atom {i}");
            assert_eq!(ctx.component_id(idx), None, "atom {i}");
        }
    }
}
