// ---------------------------------------------------------------------------
// Electron donor type (ported from RDKit's `ElectronDonorType`)
// ---------------------------------------------------------------------------

/// Per-atom pi-electron donor classification, computed once per molecule
/// (not per candidate ring). Direct port of RDKit's `ElectronDonorType`.
///
/// Crate-internal: not part of the public API. The only supported entry
/// points are [`assign_aromaticity_rdkit_parity_experimental`] and
/// [`apply_aromaticity_rdkit_parity_experimental`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ElectronDonorType {
    /// No electrons to spare, but an empty p-orbital (e.g. tropylium-type carbocation).
    Vacant,
    /// Exactly 1 electron (a normal sp2 atom with one endocyclic pi bond).
    OneElectron,
    /// Exactly 2 electrons (a lone pair, unconditionally).
    TwoElectron,
    /// Either 1 or 2, ambiguous until a specific candidate ring/subset is evaluated
    /// (RDKit tries every value in this range when checking Hückel's rule).
    ///
    /// Kept for shape-fidelity with RDKit's own `ElectronDonorType` enum
    /// (this port's `get_atom_electron_donor_type` doesn't currently
    /// construct this variant on any input the calibration battery or the
    /// 5,000-molecule corpus exercises -- was previously masked by this
    /// enum being `pub`, which suppresses rustc's dead-code analysis for
    /// externally-constructible items). Not a behavior change to fix here.
    #[allow(dead_code)]
    OneOrTwo,
    /// Dummy-atom wildcard (1 or 2, but at most one such atom per evaluated ring).
    Any,
    /// Not eligible to donate at all (disqualifies any ring it's part of).
    None,
}

/// RDKit's main-group "number of outer-shell (valence) electrons" per
/// element, used by `count_atom_pi_electrons` exactly as
/// `PeriodicTable::getNouterElecs` is used in the source. Small, stable
/// chemistry fact table — not exposed from `chematic-core` since this is the
/// only consumer.
fn outer_shell_electrons(atomic_number: u8) -> Option<u8> {
    match atomic_number {
        1 => Some(1),  // H
        5 => Some(3),  // B
        6 => Some(4),  // C
        7 => Some(5),  // N
        8 => Some(6),  // O
        9 => Some(7),  // F
        14 => Some(4), // Si
        15 => Some(5), // P
        16 => Some(6), // S
        17 => Some(7), // Cl
        33 => Some(5), // As
        34 => Some(6), // Se
        35 => Some(7), // Br
        52 => Some(6), // Te
        53 => Some(7), // I
        _ => None,
    }
}

fn default_valence(atomic_number: u8) -> Option<u8> {
    chematic_core::Element::from_atomic_number(atomic_number)
        .and_then(|e| e.normal_valences().first().copied())
}

fn bond_order_contrib(order: BondOrder) -> f32 {
    match order {
        BondOrder::Single | BondOrder::Up | BondOrder::Down => 1.0,
        BondOrder::Double => 2.0,
        BondOrder::Triple => 3.0,
        BondOrder::Quadruple => 4.0,
        // None of these should occur on pre-kekulized organic input (this
        // module's precondition) -- fall back to a single-bond-equivalent
        // rather than panicking.
        BondOrder::Aromatic
        | BondOrder::Zero
        | BondOrder::Dative
        | BondOrder::QueryAny
        | BondOrder::QuerySingleOrDouble
        | BondOrder::QuerySingleOrAromatic
        | BondOrder::QueryDoubleOrAromatic => 1.0,
    }
}

/// Port of `countAtomElec`: pi electrons available for donation into an
/// aromatic system, from generic valence-shell arithmetic — NOT
/// element-specific branching (RDKit's model is deliberately generic here).
/// Returns `None` for atoms that can never be aromatic (univalent elements,
/// degree > 3, multiple unsaturations already ruled out upstream).
/// The bond orders and aromatic flag the typing rules read for one atom:
/// the molecule's own ([`AsIs`]), or a hypothetical Kekulé form of it
/// ([`KekuleAtom`]). Only the typed atom's own flag and incident bonds are
/// ever queried.
pub(crate) trait OrderView {
    fn order(&self, mol: &Molecule, bond: BondIdx) -> BondOrder;
    fn flagged(&self, mol: &Molecule, atom: AtomIdx) -> bool;
}

/// The molecule as it is.
pub(crate) struct AsIs;

impl OrderView for AsIs {
    #[inline]
    fn order(&self, mol: &Molecule, bond: BondIdx) -> BondOrder {
        mol.bond(bond).order
    }
    #[inline]
    fn flagged(&self, mol: &Molecule, atom: AtomIdx) -> bool {
        mol.atom(atom).aromatic
    }
}

/// One aromatic atom as it is in a Kekulé form of the molecule: flag
/// cleared, every aromatic bond single except `double` (the bond that
/// received the atom's double bond, if it is matched).
pub(crate) struct KekuleAtom {
    double: Option<BondIdx>,
}

impl OrderView for KekuleAtom {
    #[inline]
    fn order(&self, mol: &Molecule, bond: BondIdx) -> BondOrder {
        match mol.bond(bond).order {
            BondOrder::Aromatic if Some(bond) == self.double => BondOrder::Double,
            BondOrder::Aromatic => BondOrder::Single,
            order => order,
        }
    }
    #[inline]
    fn flagged(&self, _mol: &Molecule, _atom: AtomIdx) -> bool {
        false
    }
}

fn count_atom_pi_electrons_in(
    mol: &Molecule,
    atom_idx: AtomIdx,
    view: &impl OrderView,
) -> Option<i32> {
    let atom = mol.atom(atom_idx);
    let an = atom.element.atomic_number();
    let dv = default_valence(an)?;
    if dv <= 1 {
        return None; // univalent elements can't be aromatic or conjugated
    }

    let implicit_h =
        chematic_core::implicit_hcount_with(mol, atom_idx, view.flagged(mol, atom_idx), |b| {
            view.order(mol, b)
        });
    let degree = mol.degree(atom_idx) + implicit_h as usize;
    // A neutral group-16 atom can donate a lone pair only while it remains
    // divalent.  Counting hydrogens matters here: thiophene-1-ol has two
    // ring bonds, one O substituent and one implicit H.  Without this guard,
    // both that atom and four-coordinate S(O)(O) are incorrectly classified
    // as aromatic two-electron donors (#767).
    if atom.charge == 0 && matches!(an, 8 | 16 | 34 | 52) && degree > 2 {
        return None;
    }
    if degree > 3 {
        return None;
    }

    let nlp_raw = outer_shell_electrons(an)? as i32 - dv as i32;
    let nlp = (nlp_raw - atom.charge as i32).max(0);
    let n_radicals = 0i32; // radicals aren't modeled in chematic-core's Atom

    let mut res = (dv as i32 - degree as i32) + nlp - n_radicals;

    if res > 1 {
        let explicit_valence: f32 = mol
            .neighbors(atom_idx)
            .map(|(_, bidx)| bond_order_contrib(view.order(mol, bidx)))
            .sum();
        let n_unsaturations = explicit_valence - mol.degree(atom_idx) as f32;
        if n_unsaturations > 1.0 {
            res = 1;
        }
    }

    Some(res)
}

/// Cyclic-bond membership, either as a bond set or as index-aligned flags
/// (the two agree: an SSSR's bonds are exactly the cyclic bonds).
pub(crate) trait RingBondLookup {
    fn contains(&self, bond: &BondIdx) -> bool;
}

impl RingBondLookup for FxHashSet<BondIdx> {
    fn contains(&self, bond: &BondIdx) -> bool {
        FxHashSet::contains(self, bond)
    }
}

impl RingBondLookup for [bool] {
    fn contains(&self, bond: &BondIdx) -> bool {
        self.get(bond.0 as usize).copied().unwrap_or(false)
    }
}

fn incident_non_cyclic_multiple_bond(
    mol: &Molecule,
    atom_idx: AtomIdx,
    ring_bonds: &(impl RingBondLookup + ?Sized),
    view: &impl OrderView,
) -> Option<AtomIdx> {
    mol.neighbors(atom_idx)
        .find(|&(_, bidx)| {
            !ring_bonds.contains(&bidx) && bond_order_contrib(view.order(mol, bidx)) >= 2.0
        })
        .map(|(nb, _)| nb)
}

fn incident_cyclic_multiple_bond(
    mol: &Molecule,
    atom_idx: AtomIdx,
    ring_bonds: &(impl RingBondLookup + ?Sized),
    view: &impl OrderView,
) -> bool {
    mol.neighbors(atom_idx).any(|(_, bidx)| {
        ring_bonds.contains(&bidx) && bond_order_contrib(view.order(mol, bidx)) >= 2.0
    })
}

fn incident_multiple_bond(mol: &Molecule, atom_idx: AtomIdx, view: &impl OrderView) -> bool {
    let explicit_valence: f32 = mol
        .neighbors(atom_idx)
        .map(|(_, bidx)| bond_order_contrib(view.order(mol, bidx)))
        .sum();
    (explicit_valence - mol.degree(atom_idx) as f32).abs() > 1e-6
}

fn more_electronegative(a: u8, b: u8) -> bool {
    // RDKit's PeriodicTable::moreElectroNegative is Pauling-scale; restricted
    // here to the elements this model's callers actually compare against
    // (the exocyclic-multiple-bond partner check), which are always N/O/S
    // relative to C -- matches every case in `isAtomCandForArom`'s callers.
    fn electronegativity(an: u8) -> f32 {
        match an {
            1 => 2.20,
            5 => 2.04,
            6 => 2.55,
            7 => 3.04,
            8 => 3.44,
            9 => 3.98,
            14 => 1.90,
            15 => 2.19,
            16 => 2.58,
            17 => 3.16,
            34 => 2.55,
            35 => 2.96,
            52 => 2.10,
            53 => 2.66,
            _ => 2.20,
        }
    }
    electronegativity(a) > electronegativity(b)
}

/// Port of `getAtomDonorTypeArom` (default params: `exocyclicBondsStealElectrons = true`).
/// `ring_bonds` = the set of bond indices that are part of *any* SSSR ring in
/// the whole molecule (global, not scoped to one candidate ring/subset).
pub(crate) fn get_atom_electron_donor_type(
    mol: &Molecule,
    atom_idx: AtomIdx,
    ring_bonds: &(impl RingBondLookup + ?Sized),
) -> ElectronDonorType {
    get_atom_electron_donor_type_in(mol, atom_idx, ring_bonds, &AsIs)
}

/// [`get_atom_electron_donor_type`] of `atom_idx` under `view`.
pub(crate) fn get_atom_electron_donor_type_in(
    mol: &Molecule,
    atom_idx: AtomIdx,
    ring_bonds: &(impl RingBondLookup + ?Sized),
    view: &impl OrderView,
) -> ElectronDonorType {
    let atom = mol.atom(atom_idx);
    let an = atom.element.atomic_number();

    let Some(nelec) = count_atom_pi_electrons_in(mol, atom_idx, view) else {
        return ElectronDonorType::None;
    };

    if nelec < 0 {
        ElectronDonorType::None
    } else if nelec == 0 {
        if let Some(_who) = incident_non_cyclic_multiple_bond(mol, atom_idx, ring_bonds, view) {
            ElectronDonorType::Vacant
        } else if incident_cyclic_multiple_bond(mol, atom_idx, ring_bonds, view) {
            ElectronDonorType::OneElectron
        } else {
            ElectronDonorType::None
        }
    } else if nelec == 1 {
        if let Some(who) = incident_non_cyclic_multiple_bond(mol, atom_idx, ring_bonds, view) {
            let other_an = mol.atom(who).element.atomic_number();
            if more_electronegative(other_an, an) {
                ElectronDonorType::Vacant
            } else {
                ElectronDonorType::OneElectron
            }
        } else if incident_multiple_bond(mol, atom_idx, view) {
            ElectronDonorType::OneElectron
        } else if atom.charge == 1 {
            // tropylium / cyclopropenyl cation
            ElectronDonorType::Vacant
        } else {
            ElectronDonorType::None
        }
    } else {
        let mut nelec = nelec;
        if let Some(who) = incident_non_cyclic_multiple_bond(mol, atom_idx, ring_bonds, view) {
            let other_an = mol.atom(who).element.atomic_number();
            if more_electronegative(other_an, an) {
                nelec -= 1;
            }
        }
        if nelec % 2 == 1 {
            ElectronDonorType::OneElectron
        } else {
            ElectronDonorType::TwoElectron
        }
    }
}

/// Port of `isAtomCandForArom` with the DEFAULT model's parameters
/// (`allowThirdRow=true, allowTripleBonds=true, allowHigherExceptions=true,
/// onlyCorN=false, allowExocyclicMultipleBonds=true`).
pub(crate) fn is_atom_candidate_for_aromaticity(
    mol: &Molecule,
    atom_idx: AtomIdx,
    donor_type: ElectronDonorType,
) -> bool {
    is_atom_candidate_for_aromaticity_in(mol, atom_idx, donor_type, &AsIs)
}

/// [`is_atom_candidate_for_aromaticity`] of `atom_idx` under `view`.
pub(crate) fn is_atom_candidate_for_aromaticity_in(
    mol: &Molecule,
    atom_idx: AtomIdx,
    donor_type: ElectronDonorType,
    view: &impl OrderView,
) -> bool {
    let atom = mol.atom(atom_idx);
    let an = atom.element.atomic_number();

    // First two rows, plus Se/Te (allowHigherExceptions).
    if an > 18 && an != 34 && an != 52 {
        return false;
    }

    if matches!(donor_type, ElectronDonorType::None) {
        return false;
    }

    // Atoms not in their default valence state are shut out.
    if let Some(dv) = default_valence(an) {
        let total_valence: f32 = mol
            .neighbors(atom_idx)
            .map(|(_, bidx)| bond_order_contrib(view.order(mol, bidx)))
            .sum::<f32>()
            + chematic_core::implicit_hcount_with(mol, atom_idx, view.flagged(mol, atom_idx), |b| {
                view.order(mol, b)
            }) as f32;
        let an_neutral = (an as i32 - atom.charge as i32).max(0) as u8;
        if let Some(dv_neutral) = default_valence(an_neutral)
            && total_valence.round() as i32 > dv_neutral as i32
        {
            return false;
        }
        let _ = dv;
    }

    // No more than one double/triple bond (rules out cumulated dienes like C=C=N).
    let explicit_valence: f32 = mol
        .neighbors(atom_idx)
        .map(|(_, bidx)| bond_order_contrib(view.order(mol, bidx)))
        .sum();
    let n_unsaturations = explicit_valence - mol.degree(atom_idx) as f32;
    if n_unsaturations > 1.0 {
        let n_mult = mol
            .neighbors(atom_idx)
            .filter(|(_, bidx)| {
                matches!(
                    view.order(mol, *bidx),
                    BondOrder::Double | BondOrder::Triple
                )
            })
            .count();
        if n_mult > 1 {
            return false;
        }
    }

    true
}

// ---------------------------------------------------------------------------
// Hückel evaluation (ported from `applyHuckel` / `applyHuckelToFused`)
// ---------------------------------------------------------------------------

fn min_max_atom_electrons(dtype: ElectronDonorType) -> (i32, i32) {
    match dtype {
        ElectronDonorType::Any | ElectronDonorType::OneOrTwo => (1, 2),
        ElectronDonorType::OneElectron => (1, 1),
        ElectronDonorType::TwoElectron => (2, 2),
        ElectronDonorType::None | ElectronDonorType::Vacant => (0, 0),
    }
}

/// Port of `applyHuckel`: given a candidate atom union, checks whether ANY
/// electron count in `[sum_of_lower_bounds, sum_of_upper_bounds]` satisfies
/// 4n+2 -- or the `rup == 2` special case for tiny rings (e.g. cyclopropenyl
/// cation).
pub(crate) fn apply_huckel(mol: &Molecule, atoms: &[AtomIdx], donor: &[ElectronDonorType]) -> bool {
    let _ = mol;
    let mut rlw = 0i32;
    let mut rup = 0i32;
    let mut n_any = 0u32;
    for &a in atoms {
        let dtype = donor[a.0 as usize];
        if dtype == ElectronDonorType::Any {
            n_any += 1;
            if n_any > 1 {
                return false;
            }
        }
        let (lo, hi) = min_max_atom_electrons(dtype);
        rlw += lo;
        rup += hi;
    }

    if rup >= 6 {
        (rlw..=rup).any(|rie| (rie - 2).rem_euclid(4) == 0)
    } else {
        rup == 2
    }
}

/// RDKit's default cap for rings considered in a fused-aromatic subsystem.
/// `makeRingNeighborMap()` keeps larger rings as isolated candidates.
const MAX_FUSED_AROMATIC_RING_SIZE: usize = 24;

/// One connected group of candidate rings, adjacent only when they share
/// exactly one bond and both rings fit RDKit's fused-aromatic size bound.
/// This is the `maxSize=24, maxOverlapSize=1` call to RDKit's
/// `makeRingNeighborMap`, not a generic "shares any bond" relation: rings
/// with a two-bond overlap must remain independent candidates.
fn fused_ring_groups(ring_bond_ids: &[Vec<BondIdx>]) -> Vec<Vec<usize>> {
    let n = ring_bond_ids.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut [usize], x: usize) -> usize {
        if parent[x] != x {
            parent[x] = find(parent, parent[x]);
        }
        parent[x]
    }
    for i in 0..n {
        for j in (i + 1)..n {
            let overlap = ring_bond_ids[i]
                .iter()
                .filter(|b| ring_bond_ids[j].contains(b))
                .count();
            if ring_bond_ids[i].len() <= MAX_FUSED_AROMATIC_RING_SIZE
                && ring_bond_ids[j].len() <= MAX_FUSED_AROMATIC_RING_SIZE
                && overlap == 1
            {
                let (pi, pj) = (find(&mut parent, i), find(&mut parent, j));
                if pi != pj {
                    parent[pi] = pj;
                }
            }
        }
    }
    let mut groups: FxHashMap<usize, Vec<usize>> = FxHashMap::default();
    for i in 0..n {
        groups.entry(find(&mut parent, i)).or_default().push(i);
    }
    let mut out: Vec<Vec<usize>> = groups.into_values().collect();
    out.sort_by_key(|g| g[0]);
    out
}

/// All `k`-combinations of `0..n`, in RDKit's `nextCombination` order
/// (ascending indices, lexicographic).
fn combinations(n: usize, k: usize) -> Vec<Vec<usize>> {
    if k == 0 || k > n {
        return vec![];
    }
    let mut result = Vec::new();
    let mut combo: Vec<usize> = (0..k).collect();
    loop {
        result.push(combo.clone());
        let mut i = k;
        loop {
            if i == 0 {
                return result;
            }
            i -= 1;
            if combo[i] != i + n - k {
                break;
            }
        }
        combo[i] += 1;
        for j in (i + 1)..k {
            combo[j] = combo[j - 1] + 1;
        }
    }
}

/// Port of `applyHuckelToFused`: within one fused ring group, tries every
/// connected subset of rings (size 1, then 2, ... up to `max_num_fused_rings`),
/// unions each subset's atoms (RDKit's #2895 rule: an atom counts only if it
/// appears in exactly 1 or 2 of the subset's rings), and marks the subset's
/// *outer perimeter* bonds/atoms aromatic if `apply_huckel` accepts. Stops
/// once every bond in the fused group has been assigned a verdict.
/// Candidate rings, bundled so `apply_huckel_to_fused` stays under clippy's
/// too-many-arguments limit -- `atoms[i]`/`bonds[i]` describe the same ring.
struct CandidateRings<'a> {
    atoms: &'a [&'a [AtomIdx]],
    bonds: &'a [Vec<BondIdx>],
}

fn apply_huckel_to_fused(
    mol: &Molecule,
    rings: &CandidateRings<'_>,
    group: &[usize],
    donor: &[ElectronDonorType],
    max_num_fused_rings: usize,
    aromatic_atoms: &mut FxHashSet<AtomIdx>,
    aromatic_bonds: &mut FxHashSet<BondIdx>,
) {
    let mut state = FusedHuckelState::new(mol, rings, group);
    for size in 1..=group.len().min(max_num_fused_rings) {
        if state.complete() {
            break;
        }
        for combo in combinations(group.len(), size) {
            let subset: Vec<usize> = combo.iter().map(|&index| group[index]).collect();
            state.apply_subset(mol, rings, &subset, donor, aromatic_atoms, aromatic_bonds);
        }
    }
}

struct FusedHuckelState {
    atom_count: Vec<u32>,
    bond_count: Vec<u32>,
    done: Vec<bool>,
    done_count: usize,
    target_count: usize,
    union: Vec<AtomIdx>,
}

impl FusedHuckelState {
    fn new(mol: &Molecule, rings: &CandidateRings<'_>, group: &[usize]) -> Self {
        let mut seen = vec![false; mol.bond_count()];
        let mut target_count = 0;
        for &ring in group {
            for &bond in &rings.bonds[ring] {
                if !seen[bond.0 as usize] {
                    seen[bond.0 as usize] = true;
                    target_count += 1;
                }
            }
        }
        Self {
            atom_count: vec![0; mol.atom_count()],
            bond_count: vec![0; mol.bond_count()],
            done: vec![false; mol.bond_count()],
            done_count: 0,
            target_count,
            union: Vec::new(),
        }
    }

    fn complete(&self) -> bool {
        self.done_count >= self.target_count
    }

    fn apply_subset(
        &mut self,
        mol: &Molecule,
        rings: &CandidateRings<'_>,
        subset: &[usize],
        donor: &[ElectronDonorType],
        aromatic_atoms: &mut FxHashSet<AtomIdx>,
        aromatic_bonds: &mut FxHashSet<BondIdx>,
    ) {
        if !connected_ring_subset(rings, subset) {
            return;
        }
        self.collect_union_atoms(rings, subset);
        if !apply_huckel(mol, &self.union, donor) {
            return;
        }
        self.mark_outer_perimeter(mol, rings, subset, aromatic_atoms, aromatic_bonds);
    }

    fn collect_union_atoms(&mut self, rings: &CandidateRings<'_>, subset: &[usize]) {
        for &ring in subset {
            for &atom in rings.atoms[ring] {
                self.atom_count[atom.0 as usize] += 1;
            }
        }
        self.union.clear();
        for &ring in subset {
            for &atom in rings.atoms[ring] {
                let count = &mut self.atom_count[atom.0 as usize];
                if matches!(*count, 1 | 2) {
                    self.union.push(atom);
                }
                *count = u32::MAX;
            }
        }
        for &ring in subset {
            for &atom in rings.atoms[ring] {
                self.atom_count[atom.0 as usize] = 0;
            }
        }
    }

    fn mark_outer_perimeter(
        &mut self,
        mol: &Molecule,
        rings: &CandidateRings<'_>,
        subset: &[usize],
        aromatic_atoms: &mut FxHashSet<AtomIdx>,
        aromatic_bonds: &mut FxHashSet<BondIdx>,
    ) {
        for &ring in subset {
            for &bond in &rings.bonds[ring] {
                self.bond_count[bond.0 as usize] += 1;
            }
        }
        for &ring in subset {
            for &bond in &rings.bonds[ring] {
                if self.bond_count[bond.0 as usize] != 1 {
                    continue;
                }
                aromatic_bonds.insert(bond);
                let endpoints = mol.bond(bond);
                aromatic_atoms.insert(endpoints.atom1);
                aromatic_atoms.insert(endpoints.atom2);
                if !self.done[bond.0 as usize] {
                    self.done[bond.0 as usize] = true;
                    self.done_count += 1;
                }
            }
        }
        for &ring in subset {
            for &bond in &rings.bonds[ring] {
                self.bond_count[bond.0 as usize] = 0;
            }
        }
    }
}

fn connected_ring_subset(rings: &CandidateRings<'_>, subset: &[usize]) -> bool {
    if subset.len() <= 1 {
        return true;
    }
    let bonds: Vec<Vec<BondIdx>> = subset
        .iter()
        .map(|&ring| rings.bonds[ring].clone())
        .collect();
    fused_ring_groups(&bonds).len() == 1
}

/// Top-level driver, matching `aromaticityHelper(mol, srings, 0, 0,
/// includeFused=true)` -- the exact function `AROMATICITY_RDKIT`/
/// `AROMATICITY_DEFAULT` call. `maxNumFusedRings` is RDKit's own hardcoded
/// default (`6`), left as a parameter for the calibration battery.
///
/// Requires pre-kekulized `mol` (see module doc comment).
pub fn rdkit_parity_aromaticity(mol: &Molecule) -> (FxHashSet<AtomIdx>, FxHashSet<BondIdx>) {
    rdkit_parity_aromaticity_ex(mol, 6)
}

pub(crate) fn rdkit_parity_aromaticity_ex(
    mol: &Molecule,
    max_num_fused_rings: usize,
) -> (FxHashSet<AtomIdx>, FxHashSet<BondIdx>) {
    verdict_with_typing(
        mol,
        max_num_fused_rings,
        |idx, cyclic| {
            let d = get_atom_electron_donor_type(mol, idx, cyclic);
            (d, is_atom_candidate_for_aromaticity(mol, idx, d))
        },
        true,
    )
    .expect("ranked ring selection always completes")
}

/// Body of [`rdkit_parity_aromaticity_ex`] with each ring atom's donor type
/// and candidacy supplied by `typing`, so the verdict of a hypothetical
/// Kekulé form can be computed on the explicit molecule (same graph). With
/// `allow_ranks == false` the SSSR ring set must be decidable from the graph
/// alone (`None` otherwise), since canonical ranks read flags and orders.
fn verdict_with_typing(
    mol: &Molecule,
    max_num_fused_rings: usize,
    typing: impl Fn(AtomIdx, &[bool]) -> (ElectronDonorType, bool),
    allow_ranks: bool,
) -> Option<(FxHashSet<AtomIdx>, FxHashSet<BondIdx>)> {
    let typing = ParityTyping::new(mol, typing);
    if !typing.keep.iter().any(|&keep| keep) {
        return Some((FxHashSet::default(), FxHashSet::default()));
    }
    let rings = parity_ring_set(mol, &typing.keep, allow_ranks)?;
    let ring_atoms: Vec<&[AtomIdx]> = rings
        .iter()
        .filter(|ring| ring.iter().all(|atom| typing.candidate[atom.0 as usize]))
        .map(|ring| ring.as_slice())
        .collect();
    let ring_bond_ids = candidate_ring_bonds(mol, &ring_atoms);
    let mut aromatic_atoms: FxHashSet<AtomIdx> = FxHashSet::default();
    let mut aromatic_bonds: FxHashSet<BondIdx> = FxHashSet::default();
    let rings = CandidateRings {
        atoms: &ring_atoms,
        bonds: &ring_bond_ids,
    };

    for group in fused_ring_groups(&ring_bond_ids) {
        apply_huckel_to_fused(
            mol,
            &rings,
            &group,
            &typing.donor,
            max_num_fused_rings,
            &mut aromatic_atoms,
            &mut aromatic_bonds,
        );
    }

    Some((aromatic_atoms, aromatic_bonds))
}

struct ParityTyping {
    donor: Vec<ElectronDonorType>,
    candidate: Vec<bool>,
    keep: Vec<bool>,
}

impl ParityTyping {
    fn new(
        mol: &Molecule,
        typing: impl Fn(AtomIdx, &[bool]) -> (ElectronDonorType, bool),
    ) -> Self {
        let atom_count = mol.atom_count();
        let cyclic = crate::sssr::ring_bond_flags_shared(mol);
        let mut ring_atom = vec![false; atom_count];
        let mut components = UnionFind::new(atom_count);
        for (index, bond) in mol.bonds() {
            if cyclic[index.0 as usize] {
                ring_atom[bond.atom1.0 as usize] = true;
                ring_atom[bond.atom2.0 as usize] = true;
                components.union(bond.atom1.0, bond.atom2.0);
            }
        }
        let mut donor = vec![ElectronDonorType::None; atom_count];
        let mut candidate = vec![false; atom_count];
        for atom in 0..atom_count {
            if ring_atom[atom] {
                (donor[atom], candidate[atom]) = typing(AtomIdx(atom as u32), cyclic.as_slice());
            }
        }
        let keep = candidate_component_mask(mol, cyclic.as_slice(), &ring_atom, &candidate, components);
        Self {
            donor,
            candidate,
            keep,
        }
    }
}

fn candidate_component_mask(
    mol: &Molecule,
    cyclic: &[bool],
    ring_atom: &[bool],
    candidate: &[bool],
    mut components: UnionFind,
) -> Vec<bool> {
    let mut candidate_components = UnionFind::new(mol.atom_count());
    let mut keep_root = vec![false; mol.atom_count()];
    for (index, bond) in mol.bonds() {
        let (left, right) = (bond.atom1.0 as usize, bond.atom2.0 as usize);
        if cyclic[index.0 as usize]
            && candidate[left]
            && candidate[right]
            && !candidate_components.union(left as u32, right as u32)
        {
            keep_root[components.find(left as u32) as usize] = true;
        }
    }
    (0..mol.atom_count())
        .map(|atom| {
            ring_atom[atom] && keep_root[components.find(atom as u32) as usize]
        })
        .collect()
}

fn parity_ring_set(
    mol: &Molecule,
    keep: &[bool],
    allow_ranks: bool,
) -> Option<Vec<Vec<AtomIdx>>> {
    let sssr = if allow_ranks {
        crate::sssr::find_sssr_ring_set_in_components(mol, keep)
    } else {
        crate::sssr::find_sssr_ring_set_in_components_unranked(mol, keep)?
    };
    let mut membership = vec![0u8; mol.atom_count()];
    for ring in &sssr {
        for atom in ring {
            membership[atom.0 as usize] = membership[atom.0 as usize].saturating_add(1);
        }
    }
    if membership.iter().any(|&count| count >= 3)
        && let Some(rings) = crate::rdkit_sssr_ring_order(mol)
    {
        return Some(
            rings
                .into_iter()
                .filter(|ring| keep[ring[0].0 as usize])
                .collect(),
        );
    }
    Some(sssr)
}

fn candidate_ring_bonds(mol: &Molecule, rings: &[&[AtomIdx]]) -> Vec<Vec<BondIdx>> {
    rings
        .iter()
        .map(|ring| {
            (0..ring.len())
                .filter_map(|index| {
                    mol.bond_between(ring[index], ring[(index + 1) % ring.len()])
                        .map(|(bond, _)| bond)
                })
                .collect()
        })
        .collect()
}
