//! Aromaticity-A1-1b-0: a faithful, independent reproduction of RDKit's
//! *default* aromaticity model (`AROMATICITY_RDKIT`/`AROMATICITY_DEFAULT`),
//! ported directly from RDKit's own source
//! (`Code/GraphMol/Aromaticity.cpp`, functions `getAtomDonorTypeArom`,
//! `countAtomElec`, `isAtomCandForArom`, `applyHuckel`, `applyHuckelToFused`,
//! `aromaticityHelper`'s `includeFused` branch — the exact path
//! `setAromaticity(mol, AROMATICITY_RDKIT, ...)` calls).
//!
//! `AromaticityAlgorithm::RdkitLike` now routes through this engine when the
//! input can be kekulized; the historical per-ring implementation remains as
//! an infallible fallback for inputs that cannot satisfy the parity engine's
//! precondition. `Huckel` remains the default and is unchanged. This module
//! also backs a separate, explicitly opt-in, fallible production API:
//! [`assign_aromaticity_rdkit_parity_experimental`] and
//! [`apply_aromaticity_rdkit_parity_experimental`], re-exported from the
//! crate root. Every other item in this module (the low-level donor-type/
//! Hückel machinery) is crate-private; the only way to reach it from
//! outside the crate is through those two functions, or, for diagnostics,
//! `diagnostics::rdkit_parity_aromaticity` behind the `diagnostics` feature.
//! See `docs/rfcs/aromaticity_a1_rfc.md`'s "A1-1b-0"/"A1-1b-1" sections for the
//! full design writeup, the calibration battery, and the corpus gate.
//!
//! Ported from RDKit commit `e89c9f656a694fab4105139844cba88d2e013354`, an
//! ancestor of release tag `Release_2026_03_4` (which resolves to
//! `8afba32ec539dcb2369bc84549d802aca3f7eb39`, independently verified via
//! the GitHub tags API during Morgan M4-A0). `Code/GraphMol/Aromaticity.cpp`
//! is byte-identical between the two commits (independently diffed during
//! M4-A0's provenance audit — the 130 commits between them never touch this
//! file), so the 5 functions cited above are unaffected either way. See
//! `THIRD_PARTY_NOTICES.md` at the repo root for the required BSD 3-Clause
//! attribution and license text.
//!
//! Unlike this crate's own `ring_pi_electrons`/`evaluate_atom_pi_contribution`
//! (which evaluate an atom's contribution *per candidate ring/component*),
//! RDKit computes each atom's [`ElectronDonorType`] **once, globally, per
//! molecule** — whether a multiple bond "counts" for aromaticity purposes
//! depends on whether that bond is part of *any* SSSR ring in the whole
//! molecule (`RingInfo::numBondRings(bond) > 0`), not on whether it's inside
//! the *specific* candidate ring currently being evaluated. This is the
//! precise, source-verified point where this crate's own `ring_pi_electrons`
//! diverges from RDKit for the SMARTS-A0/PR #86 false-positive family: its
//! `CarbonExocyclicHeteroatomDouble` rule checks "is the double-bond partner
//! outside *this ring's* atom set" where RDKit checks "is this bond outside
//! *every* ring in the molecule" — an exocyclic-to-the-candidate-ring double
//! bond whose partner is itself a *different* ring's atom (e.g. this crate's
//! reproducer's atom 8, `C=N` where the N is in a second fused ring) still
//! counts as a normal one-electron donor under RDKit's rule, not a
//! zero-electron "spent on the exocyclic bond" donor.
//!
//! Requires pre-kekulized input (no `BondOrder::Aromatic`), matching RDKit's
//! own pipeline (`Kekulize` always runs before `setAromaticity`).

use rustc_hash::{FxHashMap, FxHashSet};

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule};

use crate::aromaticity::AromaticityModel;
use crate::sssr::find_sssr;

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
fn count_atom_pi_electrons(mol: &Molecule, atom_idx: AtomIdx) -> Option<i32> {
    let atom = mol.atom(atom_idx);
    let an = atom.element.atomic_number();
    let dv = default_valence(an)?;
    if dv <= 1 {
        return None; // univalent elements can't be aromatic or conjugated
    }

    let implicit_h = chematic_core::implicit_hcount(mol, atom_idx);
    let degree = mol.degree(atom_idx) + implicit_h as usize;
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
            .map(|(_, bidx)| bond_order_contrib(mol.bond(bidx).order))
            .sum();
        let n_unsaturations = explicit_valence - mol.degree(atom_idx) as f32;
        if n_unsaturations > 1.0 {
            res = 1;
        }
    }

    Some(res)
}

fn incident_non_cyclic_multiple_bond(
    mol: &Molecule,
    atom_idx: AtomIdx,
    ring_bonds: &FxHashSet<BondIdx>,
) -> Option<AtomIdx> {
    mol.neighbors(atom_idx)
        .find(|&(_, bidx)| {
            !ring_bonds.contains(&bidx) && bond_order_contrib(mol.bond(bidx).order) >= 2.0
        })
        .map(|(nb, _)| nb)
}

fn incident_cyclic_multiple_bond(
    mol: &Molecule,
    atom_idx: AtomIdx,
    ring_bonds: &FxHashSet<BondIdx>,
) -> bool {
    mol.neighbors(atom_idx).any(|(_, bidx)| {
        ring_bonds.contains(&bidx) && bond_order_contrib(mol.bond(bidx).order) >= 2.0
    })
}

fn incident_multiple_bond(mol: &Molecule, atom_idx: AtomIdx) -> bool {
    let explicit_valence: f32 = mol
        .neighbors(atom_idx)
        .map(|(_, bidx)| bond_order_contrib(mol.bond(bidx).order))
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
    ring_bonds: &FxHashSet<BondIdx>,
) -> ElectronDonorType {
    let atom = mol.atom(atom_idx);
    let an = atom.element.atomic_number();

    let Some(nelec) = count_atom_pi_electrons(mol, atom_idx) else {
        return ElectronDonorType::None;
    };

    if nelec < 0 {
        ElectronDonorType::None
    } else if nelec == 0 {
        if let Some(_who) = incident_non_cyclic_multiple_bond(mol, atom_idx, ring_bonds) {
            ElectronDonorType::Vacant
        } else if incident_cyclic_multiple_bond(mol, atom_idx, ring_bonds) {
            ElectronDonorType::OneElectron
        } else {
            ElectronDonorType::None
        }
    } else if nelec == 1 {
        if let Some(who) = incident_non_cyclic_multiple_bond(mol, atom_idx, ring_bonds) {
            let other_an = mol.atom(who).element.atomic_number();
            if more_electronegative(other_an, an) {
                ElectronDonorType::Vacant
            } else {
                ElectronDonorType::OneElectron
            }
        } else if incident_multiple_bond(mol, atom_idx) {
            ElectronDonorType::OneElectron
        } else if atom.charge == 1 {
            // tropylium / cyclopropenyl cation
            ElectronDonorType::Vacant
        } else {
            ElectronDonorType::None
        }
    } else {
        let mut nelec = nelec;
        if let Some(who) = incident_non_cyclic_multiple_bond(mol, atom_idx, ring_bonds) {
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
            .map(|(_, bidx)| bond_order_contrib(mol.bond(bidx).order))
            .sum::<f32>()
            + chematic_core::implicit_hcount(mol, atom_idx) as f32;
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
        .map(|(_, bidx)| bond_order_contrib(mol.bond(bidx).order))
        .sum();
    let n_unsaturations = explicit_valence - mol.degree(atom_idx) as f32;
    if n_unsaturations > 1.0 {
        let n_mult = mol
            .neighbors(atom_idx)
            .filter(|(_, bidx)| {
                matches!(mol.bond(*bidx).order, BondOrder::Double | BondOrder::Triple)
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
pub(crate) fn apply_huckel(
    mol: &Molecule,
    atoms: &[AtomIdx],
    donor: &FxHashMap<AtomIdx, ElectronDonorType>,
) -> bool {
    let _ = mol;
    let mut rlw = 0i32;
    let mut rup = 0i32;
    let mut n_any = 0u32;
    for &a in atoms {
        let dtype = donor[&a];
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

/// One connected group of candidate rings, adjacent when they share ≥1 bond
/// (RDKit's `makeRingNeighborMap`).
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
            if ring_bond_ids[i]
                .iter()
                .any(|b| ring_bond_ids[j].contains(b))
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
    atoms: &'a [Vec<AtomIdx>],
    bonds: &'a [Vec<BondIdx>],
}

fn apply_huckel_to_fused(
    mol: &Molecule,
    rings: &CandidateRings<'_>,
    group: &[usize],
    donor: &FxHashMap<AtomIdx, ElectronDonorType>,
    max_num_fused_rings: usize,
    aromatic_atoms: &mut FxHashSet<AtomIdx>,
    aromatic_bonds: &mut FxHashSet<BondIdx>,
) {
    let ring_atoms = rings.atoms;
    let ring_bond_ids = rings.bonds;
    let n_ring_bonds: usize = {
        let mut all: FxHashSet<BondIdx> = FxHashSet::default();
        for &ri in group {
            all.extend(ring_bond_ids[ri].iter().copied());
        }
        all.len()
    };
    let mut done_bonds: FxHashSet<BondIdx> = FxHashSet::default();

    for size in 1..=group.len().min(max_num_fused_rings) {
        if done_bonds.len() >= n_ring_bonds {
            break;
        }
        for combo in combinations(group.len(), size) {
            let cur_rings: Vec<usize> = combo.iter().map(|&i| group[i]).collect();

            // Subset must itself be connected (share bonds pairwise-reachable).
            if size > 1 {
                let sub_bond_ids: Vec<Vec<BondIdx>> = cur_rings
                    .iter()
                    .map(|&ri| ring_bond_ids[ri].clone())
                    .collect();
                if fused_ring_groups(&sub_bond_ids).len() != 1 {
                    continue;
                }
            }

            let mut membership_count: FxHashMap<AtomIdx, u32> = FxHashMap::default();
            for &ri in &cur_rings {
                for &a in &ring_atoms[ri] {
                    *membership_count.entry(a).or_insert(0) += 1;
                }
            }
            let union: Vec<AtomIdx> = membership_count
                .iter()
                .filter(|&(_, &c)| c == 1 || c == 2)
                .map(|(&a, _)| a)
                .collect();

            if apply_huckel(mol, &union, donor) {
                // Mark only the outer-perimeter bonds (appear in exactly one
                // of this subset's rings), matching `markAtomsBondsArom`.
                let mut bond_count: FxHashMap<BondIdx, u32> = FxHashMap::default();
                for &ri in &cur_rings {
                    for &b in &ring_bond_ids[ri] {
                        *bond_count.entry(b).or_insert(0) += 1;
                    }
                }
                for (&b, &c) in &bond_count {
                    if c == 1 {
                        aromatic_bonds.insert(b);
                        let bond = mol.bond(b);
                        aromatic_atoms.insert(bond.atom1);
                        aromatic_atoms.insert(bond.atom2);
                        done_bonds.insert(b);
                    }
                }
            }
        }
    }
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
    let sssr = find_sssr(mol);
    let srings = sssr.rings();

    let all_ring_bonds: FxHashSet<BondIdx> = srings
        .iter()
        .flat_map(|ring| {
            (0..ring.len()).filter_map(move |i| {
                mol.bond_between(ring[i], ring[(i + 1) % ring.len()])
                    .map(|(bidx, _)| bidx)
            })
        })
        .collect();

    let mut donor: FxHashMap<AtomIdx, ElectronDonorType> = FxHashMap::default();
    let mut candidate: FxHashMap<AtomIdx, bool> = FxHashMap::default();
    for ring in srings {
        for &a in ring {
            donor
                .entry(a)
                .or_insert_with(|| get_atom_electron_donor_type(mol, a, &all_ring_bonds));
            let d = donor[&a];
            candidate
                .entry(a)
                .or_insert_with(|| is_atom_candidate_for_aromaticity(mol, a, d));
        }
    }

    let candidate_rings: Vec<&Vec<AtomIdx>> = srings
        .iter()
        .filter(|ring| {
            ring.iter()
                .all(|a| candidate.get(a).copied().unwrap_or(false))
        })
        .collect();

    let ring_atoms: Vec<Vec<AtomIdx>> = candidate_rings.iter().map(|r| (*r).clone()).collect();
    let ring_bond_ids: Vec<Vec<BondIdx>> = ring_atoms
        .iter()
        .map(|ring| {
            (0..ring.len())
                .filter_map(|i| {
                    mol.bond_between(ring[i], ring[(i + 1) % ring.len()])
                        .map(|(bidx, _)| bidx)
                })
                .collect()
        })
        .collect();

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
            &donor,
            max_num_fused_rings,
            &mut aromatic_atoms,
            &mut aromatic_bonds,
        );
    }

    (aromatic_atoms, aromatic_bonds)
}

// ---------------------------------------------------------------------------
// Production entry points (A1-1b-1): fallible opt-in API
// ---------------------------------------------------------------------------

/// Error from the RDKit-parity experimental aromaticity API.
///
/// Unlike [`assign_aromaticity_ex`](crate::assign_aromaticity_ex)/
/// [`apply_aromaticity_ex`](crate::apply_aromaticity_ex) (infallible, and
/// unchanged by this addition), this engine requires an explicit
/// kekulization step it does not control the success of, so its entry
/// points return `Result` rather than silently falling back to another
/// algorithm or panicking.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AromaticityError {
    /// The input could not be reduced to a Kekulé form (no `BondOrder::Aromatic`
    /// bonds), which this engine requires as a precondition -- mirrors RDKit's
    /// own pipeline, where `Kekulize` always runs before `setAromaticity`.
    KekulizationFailed {
        /// Human-readable detail from the underlying `chematic_core::KekuleError`.
        reason: String,
    },
    /// A post-computation sanity check failed (e.g. an aromatic bond with a
    /// non-aromatic endpoint atom) -- should never happen for chemically
    /// valid input; surfaced as an error rather than a panic or a silently
    /// wrong result.
    InternalInvariantViolation {
        /// Human-readable detail of which invariant failed.
        reason: String,
    },
}

impl std::fmt::Display for AromaticityError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AromaticityError::KekulizationFailed { reason } => {
                write!(f, "rdkit-parity aromaticity: kekulization failed: {reason}")
            }
            AromaticityError::InternalInvariantViolation { reason } => {
                write!(
                    f,
                    "rdkit-parity aromaticity: internal invariant violation: {reason}"
                )
            }
        }
    }
}

impl std::error::Error for AromaticityError {}

/// Clone `mol` with every atom's `aromatic` flag reset to `false`.
///
/// Bond orders (including any `BondOrder::Aromatic`) are copied unchanged --
/// only the atom-level annotation is cleared. This engine derives
/// aromaticity purely from element/charge/bond-order structure (never reads
/// `atom.aromatic`), so clearing stale flags here has no effect on the
/// computation itself; it only ensures the *output* molecule's flags come
/// entirely from this engine's own verdict, never from whatever annotation
/// the caller's input happened to carry in.
fn clear_aromatic_flags(mol: &Molecule) -> Molecule {
    use chematic_core::MoleculeBuilder;
    let mut builder = MoleculeBuilder::new();
    for (_, atom) in mol.atoms() {
        let mut a = atom.clone();
        a.aromatic = false;
        builder.add_atom(a);
    }
    for (_, bond) in mol.bonds() {
        let _ = builder.add_bond(bond.atom1, bond.atom2, bond.order);
    }
    builder.copy_stereo_groups_from(mol);
    builder.copy_stereo_from(mol);
    builder.copy_bond_directions_from(mol);
    builder.build()
}

/// Normalize `mol` to pure Kekulé form (this engine's precondition): clear
/// stale aromatic flags, then kekulize. Returns an explicit error instead of
/// falling back to another algorithm or leaving a partially-rewritten
/// molecule behind -- `mol` itself is never mutated, only read.
fn kekulize_for_rdkit_parity(mol: &Molecule) -> Result<Molecule, AromaticityError> {
    let cleared = clear_aromatic_flags(mol);
    match chematic_core::kekulize(&cleared) {
        Ok(k) => Ok(chematic_core::apply_kekule(&cleared, &k)),
        Err(e) => Err(AromaticityError::KekulizationFailed { reason: e.detail }),
    }
}

/// Every aromatic bond's two endpoint atoms must themselves be in the
/// aromatic atom set -- a basic well-formedness property of any Hückel
/// verdict. Cheap to check, catches a class of bug that would otherwise
/// surface downstream as a confusing SMILES/valence inconsistency instead
/// of a clear error at the source.
fn validate_aromaticity_invariants(
    mol: &Molecule,
    atoms: &FxHashSet<AtomIdx>,
    bonds: &FxHashSet<BondIdx>,
) -> Result<(), AromaticityError> {
    for &bidx in bonds {
        let bond = mol.bond(bidx);
        if !atoms.contains(&bond.atom1) || !atoms.contains(&bond.atom2) {
            return Err(AromaticityError::InternalInvariantViolation {
                reason: format!(
                    "aromatic bond {bidx:?} ({:?}-{:?}) has a non-aromatic endpoint atom",
                    bond.atom1, bond.atom2
                ),
            });
        }
    }
    Ok(())
}

fn assign_from_kekulized(kekulized: &Molecule) -> Result<AromaticityModel, AromaticityError> {
    let (atoms, bonds) = rdkit_parity_aromaticity(kekulized);
    validate_aromaticity_invariants(kekulized, &atoms, &bonds)?;
    Ok(AromaticityModel::from_atom_bond_sets(atoms, bonds))
}

/// Assign aromaticity using the RDKit-parity reference engine
/// (`rdkit_parity_aromaticity`, see the module doc comment).
///
/// Explicitly opt-in and separate from [`assign_aromaticity_ex`]/
/// [`AromaticityAlgorithm`] -- those remain infallible and unchanged. This
/// function is fallible because it performs its own kekulization
/// internally (this engine requires pre-kekulized input); on failure,
/// `mol` is never touched and no partial result is produced.
///
/// The returned model's atom/bond indices correspond 1:1 with `mol`'s own
/// indices (kekulization here is index-preserving: it only clears stale
/// aromatic flags and normalizes bond orders, never adds/removes/reorders
/// atoms or bonds).
///
/// [`ring_classifications`](AromaticityModel::ring_classifications) and
/// [`antiaromatic_rings`](AromaticityModel::antiaromatic_rings) are always
/// empty on the returned model -- this engine (like RDKit's own) determines
/// only the aromatic atom/bond sets, not a per-ring classification or
/// antiaromaticity verdict.
///
/// [`assign_aromaticity_ex`]: crate::assign_aromaticity_ex
/// [`AromaticityAlgorithm`]: crate::AromaticityAlgorithm
pub fn assign_aromaticity_rdkit_parity_experimental(
    mol: &Molecule,
) -> Result<AromaticityModel, AromaticityError> {
    let kekulized = kekulize_for_rdkit_parity(mol)?;
    assign_from_kekulized(&kekulized)
}

/// Apply aromaticity using the RDKit-parity reference engine, returning a
/// new [`Molecule`] with atom/bond flags set according to the computed
/// model.
///
/// See [`assign_aromaticity_rdkit_parity_experimental`] for the fallibility
/// contract (kekulization failure is reported, never silently substituted
/// or partially applied) and the index-correspondence guarantee.
pub fn apply_aromaticity_rdkit_parity_experimental(
    mol: &Molecule,
) -> Result<Molecule, AromaticityError> {
    let kekulized = kekulize_for_rdkit_parity(mol)?;
    let model = assign_from_kekulized(&kekulized)?;
    Ok(crate::aromaticity::build_molecule_from_model(
        &kekulized, &model,
    ))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn mol_kekulized(smiles: &str) -> Molecule {
        let mol = chematic_smiles::parse(smiles).expect("valid SMILES");
        let k = chematic_core::kekulize(&mol).expect("kekulizable");
        chematic_core::apply_kekule(&mol, &k)
    }

    // Calibration battery, RDKit-atom-index-verified (not guessed): every
    // entry here was checked against a live `rdkit.Chem.MolFromSmiles(...)`
    // atom-aromaticity dump before being pinned. Covers the exact cases that
    // motivated this module: simple monocyclics (benzene/pyrrole/furan/
    // thiophene), the exocyclic-carbonyl-in-ring rule (tropone/2-pyridone/
    // 4-pyranone), a genuine bridgehead spanning two valid rings
    // (indolizine), a non-alternant fused bicyclic needing the whole-perimeter
    // candidate (azulene), plain fused benzenoids (naphthalene/anthracene),
    // fused heteroaromatics (indole/quinoline/purine), and both open findings
    // from Aromaticity-A1-1a (the false-positive reproducer, purine).
    #[test]
    fn calibration_battery_matches_rdkit() {
        let cases: &[(&str, &str, &[u32])] = &[
            ("benzene", "c1ccccc1", &[0, 1, 2, 3, 4, 5]),
            ("pyrrole", "c1cc[nH]c1", &[0, 1, 2, 3, 4]),
            ("furan", "c1ccoc1", &[0, 1, 2, 3, 4]),
            ("thiophene", "c1ccsc1", &[0, 1, 2, 3, 4]),
            (
                "selenophene (Se analog control, pre-Kekulized input)",
                "C1=C[Se]C=C1",
                &[0, 1, 2, 3, 4],
            ),
            (
                "tellurophene (pre-Kekulized input; regression test for the Te \
                 normal_valences() gap fixed in chematic-core's element.rs — Te previously \
                 had no valence-table entry, so default_valence/count_atom_pi_electrons/\
                 get_atom_electron_donor_type all returned None and \
                 is_atom_candidate_for_aromaticity rejected it outright)",
                "C1=C[Te]C=C1",
                &[0, 1, 2, 3, 4],
            ),
            ("tropone", "O=c1cccccc1", &[1, 2, 3, 4, 5, 6, 7]),
            ("2-pyridone", "O=c1cccc[nH]1", &[1, 2, 3, 4, 5, 6]),
            ("4-pyranone", "O=c1ccocc1", &[1, 2, 3, 4, 5, 6]),
            (
                "indolizine (true bridgehead, both rings valid)",
                "c1ccn2ccccc12",
                &[0, 1, 2, 3, 4, 5, 6, 7, 8],
            ),
            (
                "azulene (non-alternant, needs whole-perimeter candidate)",
                "C1=CC2=CC=CC=CC2=C1",
                &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            ),
            (
                "naphthalene",
                "c1ccc2ccccc2c1",
                &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            ),
            (
                "anthracene",
                "c1ccc2cc3ccccc3cc2c1",
                &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
            ),
            ("indole", "c1ccc2[nH]ccc2c1", &[0, 1, 2, 3, 4, 5, 6, 7, 8]),
            (
                "quinoline",
                "c1ccc2ncccc2c1",
                &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            ),
            (
                "purine (Aromaticity-A1-1a open finding, fixed here)",
                "c1cnc2[nH]cnc2n1",
                &[0, 1, 2, 3, 4, 5, 6, 7, 8],
            ),
            (
                "PR #86 false-positive reproducer (Aromaticity-A1-1a open finding, fixed here)",
                "C1=Cc2ccccc2C2=NCCCN12",
                &[2, 3, 4, 5, 6, 7],
            ),
        ];

        for (name, smi, expected) in cases {
            let mol = mol_kekulized(smi);
            let (atoms, _bonds) = rdkit_parity_aromaticity(&mol);
            let mut got: Vec<u32> = atoms.iter().map(|a| a.0).collect();
            got.sort();
            assert_eq!(&got, expected, "{name} ({smi}): should match RDKit exactly");
        }
    }

    #[test]
    fn te_default_valence_and_aromatic_bond_flags() {
        // Direct check of the dependent-path fix: default_valence(52) must now
        // resolve to Some(2) (was None before the chematic-core element.rs
        // Te normal_valences() fix).
        assert_eq!(default_valence(52), Some(2));

        // The calibration battery above only checks the aromatic ATOM set for
        // tellurophene. RDKit also reports tellurophene's ring BONDS as
        // aromatic (bond order 1.5 on every ring bond, oracle-verified via
        // rdkit.Chem.MolFromSmiles("C1=C[Te]C=C1")); confirm chematic's
        // aromatic bond set matches too, not just the atom set.
        let mol = mol_kekulized("C1=C[Te]C=C1");
        let (atoms, bonds) = rdkit_parity_aromaticity(&mol);
        assert_eq!(
            atoms.len(),
            5,
            "all 5 tellurophene ring atoms must be aromatic"
        );
        assert_eq!(
            bonds.len(),
            5,
            "all 5 tellurophene ring bonds must be aromatic"
        );
    }

    // Purine's Aromaticity-A1-0 finding was that production's answer depends
    // on whether the input was Kekulized before `apply_aromaticity` ran.
    // rdkit_parity_aromaticity must NOT reintroduce that: both a raw
    // aromatic-lowercase parse (kekulized here identically to every other
    // corpus entry, so this mostly re-confirms `mol_kekulized`'s own
    // determinism) and chematic's own `kekulize()` choice must agree with
    // each other and with RDKit.
    #[test]
    fn purine_representation_stable() {
        let smi = "c1cnc2[nH]cnc2n1";
        let raw = chematic_smiles::parse(smi).expect("valid SMILES");
        let k = chematic_core::kekulize(&raw).expect("purine should kekulize");
        let via_own_kekulize = chematic_core::apply_kekule(&raw, &k);

        let (atoms_a, _) = rdkit_parity_aromaticity(&mol_kekulized(smi));
        let (atoms_b, _) = rdkit_parity_aromaticity(&via_own_kekulize);

        let mut a: Vec<u32> = atoms_a.iter().map(|x| x.0).collect();
        let mut b: Vec<u32> = atoms_b.iter().map(|x| x.0).collect();
        a.sort();
        b.sort();
        assert_eq!(a, b, "purine: two Kekulization paths disagree");
        assert_eq!(
            a,
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8],
            "purine: should match RDKit (all 9 atoms aromatic)"
        );
    }

    #[test]
    fn production_api_assign_matches_engine_on_benzene() {
        let mol = chematic_smiles::parse("c1ccccc1").expect("valid SMILES");
        let model = assign_aromaticity_rdkit_parity_experimental(&mol).expect("benzene kekulizes");
        assert_eq!(model.aromatic_atom_count(), 6);
        for (idx, _) in mol.atoms() {
            assert!(
                model.is_atom_aromatic(idx),
                "atom {idx:?} should be aromatic"
            );
        }
        // This engine only determines atom/bond sets, not per-ring
        // classification or antiaromaticity -- both must be empty.
        assert!(model.ring_classifications().is_empty());
        assert!(model.antiaromatic_rings().is_empty());
    }

    #[test]
    fn production_api_apply_sets_aromatic_flags_and_bond_orders() {
        let mol = chematic_smiles::parse("C1=CC=CC=C1").expect("valid SMILES"); // Kekule benzene
        let applied = apply_aromaticity_rdkit_parity_experimental(&mol).expect("benzene kekulizes");
        assert_eq!(applied.atom_count(), mol.atom_count());
        for (_, atom) in applied.atoms() {
            assert!(atom.aromatic, "every benzene atom should end up aromatic");
        }
        for (_, bond) in applied.bonds() {
            assert_eq!(bond.order, BondOrder::Aromatic);
        }
    }

    #[test]
    fn production_api_reports_kekulize_failure_not_panic() {
        // The one known-gap molecule from the full-corpus gate: RDKit itself
        // parses this fine, but chematic's own `kekulize()` rejects a
        // bridgehead N in this fused purine-like system. Must surface as
        // `AromaticityError::KekulizationFailed`, not a panic and not a
        // silent fallback to another algorithm.
        let smi = "Cc1cn2c(=O)c3ncn(COCCO)c3nc2n1C";
        let mol = chematic_smiles::parse(smi).expect("valid SMILES");

        match assign_aromaticity_rdkit_parity_experimental(&mol) {
            Err(AromaticityError::KekulizationFailed { .. }) => {}
            other => panic!("expected KekulizationFailed, got {other:?}"),
        }
        // `Molecule` has no `Debug` impl, so match on the error shape only
        // (discarding the `Ok(Molecule)` payload) rather than formatting
        // the whole `Result` on failure.
        match apply_aromaticity_rdkit_parity_experimental(&mol).map(|_| ()) {
            Err(AromaticityError::KekulizationFailed { .. }) => {}
            other => panic!("expected KekulizationFailed, got {other:?}"),
        }
    }

    #[test]
    fn production_api_does_not_mutate_input_on_failure() {
        // "元の分子を途中まで書き換えてから失敗する経路は作らないでください" --
        // `mol` is only ever taken by `&Molecule` throughout the fallible
        // path (`clear_aromatic_flags`/`kekulize`/`apply_kekule` all build
        // new molecules rather than mutating in place), so this is enforced
        // by the type system. Pin it as an explicit regression: the input's
        // own atom/bond flags and counts are unchanged after a failed call.
        let smi = "Cc1cn2c(=O)c3ncn(COCCO)c3nc2n1C";
        let mol = chematic_smiles::parse(smi).expect("valid SMILES");
        let atom_count_before = mol.atom_count();
        let bond_count_before = mol.bond_count();
        let aromatic_before: Vec<bool> = mol.atoms().map(|(_, a)| a.aromatic).collect();

        let result = assign_aromaticity_rdkit_parity_experimental(&mol);
        assert!(
            result.is_err(),
            "this molecule is a known kekulize-gap case"
        );

        assert_eq!(mol.atom_count(), atom_count_before);
        assert_eq!(mol.bond_count(), bond_count_before);
        let aromatic_after: Vec<bool> = mol.atoms().map(|(_, a)| a.aromatic).collect();
        assert_eq!(aromatic_before, aromatic_after);
    }

    /// Was a known kekulize-gap case (found during Morgan M4-A0, `chematic-fp`'s
    /// `rdkit_morgan_hash.rs`, full-corpus validation): pyridinium's protonated
    /// `[nH+]` used to make `chematic_core::kekulize()` hard-fail because
    /// `atom_must_be_matched`'s N-with-H lone-pair-donor rule was charge-blind
    /// (docs/rfcs/aromaticity_rdkit_parity_rfc.md §1, root cause A). Fixed by
    /// `fix/kekulize-charge-aware-k1` (see
    /// docs/rfcs/kekulize_charge_aware_rdkit_parity.md): the rule now requires
    /// `atom.charge <= 0`, so a protonated ring N routes back to "must be
    /// matched" -- same as neutral pyridine's bare N -- instead of being
    /// wrongly treated like neutral pyrrole's `[nH]`. Kept as a regression
    /// test (not deleted) so a future re-introduction of the charge-blind rule
    /// is caught here, not just in the 40-fixture diagnosis corpus.
    ///
    /// NOTE for whoever picks up the companion fp-side fix: `chematic-fp`'s
    /// `rdkit_morgan_ecfp4.rs` has two tests
    /// (`kekule_pyridinium_reports_kekulization_failed_not_a_fallback_result`,
    /// `hueckel_fallback_would_be_detectable_if_silently_reintroduced`) that
    /// also used this exact SMILES as their "kekulize fails" positive control
    /// for a *different* invariant (the fallible ECFP4 path must not silently
    /// fall back to Hueckel) -- those are out of scope for K1 (chematic-fp is
    /// off limits for this fix) and still fail as of this commit. They need
    /// the same swap this test got: same invariant, a still-failing molecule
    /// as the example (`Cc1cn2c(=O)c3ncn(COCCO)c3nc2n1C`, see
    /// `production_api_reports_kekulize_failure_not_panic` above), not a
    /// deletion. `validation/README.md`'s Morgan M4-A0 section and several
    /// `validation/*.json` artifacts also cite a now-stale "2 of 5,048"
    /// preprocessing-failure count (only the purine molecule remains).
    #[test]
    fn kekulize_charge_aware_k1_fixes_protonated_pyridinium() {
        let smi = "c1cc[nH+]cc1";
        let mol = chematic_smiles::parse(smi).expect("valid SMILES");
        let result = apply_aromaticity_rdkit_parity_experimental(&mol);
        match result {
            Ok(applied) => {
                for (_, atom) in applied.atoms() {
                    assert!(
                        atom.aromatic,
                        "pyridinium's ring must end up fully aromatic"
                    );
                }
            }
            Err(other) => {
                panic!("expected pyridinium to kekulize successfully post-K1, got {other:?}")
            }
        }
    }

    #[test]
    fn production_api_stale_aromatic_flag_is_overridden_not_leaked() {
        // A non-aromatic atom that happens to carry a stale `aromatic=true`
        // flag on input must not leak that flag into the output -- the
        // engine's own verdict is the sole source of truth for the result.
        use chematic_core::MoleculeBuilder;

        let mut base = chematic_smiles::parse("CC").expect("valid SMILES"); // ethane, acyclic
        let mut builder = MoleculeBuilder::new();
        for (_, atom) in base.atoms() {
            let mut a = atom.clone();
            a.aromatic = true; // stale/bogus annotation
            builder.add_atom(a);
        }
        for (_, bond) in base.bonds() {
            let _ = builder.add_bond(bond.atom1, bond.atom2, bond.order);
        }
        base = builder.build();
        assert!(base.atoms().all(|(_, a)| a.aromatic), "test setup sanity");

        let applied = apply_aromaticity_rdkit_parity_experimental(&base)
            .expect("acyclic molecule kekulizes trivially (no-op)");
        assert!(
            applied.atoms().all(|(_, a)| !a.aromatic),
            "stale aromatic=true on an acyclic atom must not survive"
        );
    }
}
