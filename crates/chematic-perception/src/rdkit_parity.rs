//! Aromaticity-A1-1b-0: a faithful, independent reproduction of RDKit's
//! *default* aromaticity model (`AROMATICITY_RDKIT`/`AROMATICITY_DEFAULT`),
//! ported directly from RDKit's own source
//! (`Code/GraphMol/Aromaticity.cpp`, functions `getAtomDonorTypeArom`,
//! `countAtomElec`, `isAtomCandForArom`, `applyHuckel`, `applyHuckelToFused`,
//! `aromaticityHelper`'s `includeFused` branch — the exact path
//! `setAromaticity(mol, AROMATICITY_RDKIT, ...)` calls).
//!
//! **Test/diagnostic-only. Not wired into `assign_aromaticity_ex`,
//! `apply_aromaticity_ex`, `ring_pi_electrons`, or any other production
//! decision path.** See `docs/aromaticity_a1_rfc.md`'s "A1-1b-0" section for
//! the full design writeup, the calibration battery, and the corpus gate.
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

use crate::sssr::find_sssr;

// ---------------------------------------------------------------------------
// Electron donor type (ported from RDKit's `ElectronDonorType`)
// ---------------------------------------------------------------------------

/// Per-atom pi-electron donor classification, computed once per molecule
/// (not per candidate ring). Direct port of RDKit's `ElectronDonorType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ElectronDonorType {
    /// No electrons to spare, but an empty p-orbital (e.g. tropylium-type carbocation).
    Vacant,
    /// Exactly 1 electron (a normal sp2 atom with one endocyclic pi bond).
    OneElectron,
    /// Exactly 2 electrons (a lone pair, unconditionally).
    TwoElectron,
    /// Either 1 or 2, ambiguous until a specific candidate ring/subset is evaluated
    /// (RDKit tries every value in this range when checking Hückel's rule).
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
pub fn get_atom_electron_donor_type(
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
pub fn is_atom_candidate_for_aromaticity(
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
pub fn apply_huckel(
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

pub fn rdkit_parity_aromaticity_ex(
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
}
