//! Hückel aromaticity perception with antiaromaticity detection.
//!
//! Works on kekulized molecules (no `Aromatic` bond orders) **or** on molecules
//! that retain `Aromatic` bond orders from the SMILES parser (pre-kekulization).
//! Call `kekulize` + `apply_kekule` from `chematic-core` before calling
//! `assign_aromaticity` if you need the explicit double-bond form.
//!
//! Algorithm:
//! 1. Find all SSSR rings via `find_sssr`.
//! 2. **Pass 1**: evaluate each ring independently using Hückel electron counting.
//!    Aromatic (`BondOrder::Aromatic`) bonds are treated equivalently to double bonds
//!    so that pre-kekulization input is handled correctly.
//!    A special "bridgehead N" rule covers fused-ring N atoms whose entire valence
//!    is satisfied by single σ-bonds (like indolizine's junction nitrogen).
//! 3. **Pass 2**: iterative propagation. Rings that were `NonAromatic` or
//!    indeterminate in Pass 1 are re-evaluated using the already-aromatic atom set
//!    as context: confirmed-aromatic atoms contribute 1π unconditionally, allowing
//!    fused rings to be recognised bottom-up (e.g. the 6-ring of indolizine).
//! 4. Classify rings by electron count:
//!    - 4n+2 electrons (n >= 0): aromatic (favorable)
//!    - 4n electrons (n > 0): antiaromatic (unfavorable, strongly disfavored)
//!    - Other: non-aromatic
//! 5. Record all aromatic atoms, bonds, and antiaromatic rings in an `AromaticityModel`.

use std::collections::HashSet;

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule, implicit_hcount};

use crate::sssr::find_sssr;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Ring aromaticity classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RingAromaticity {
    /// 4n+2 electrons: aromatic (favorable)
    Aromatic,
    /// 4n electrons (n > 0): antiaromatic (unfavorable)
    Antiaromatic,
    /// Any other electron count: non-aromatic
    NonAromatic,
}

/// Aromaticity assignment for a molecule.
///
/// Records which atoms and bonds belong to aromatic rings according to
/// the Hückel 4n+2 rule applied to SSSR rings (with fused-ring propagation).
/// Also tracks antiaromatic rings (4n electrons) for chemical accuracy.
#[derive(Debug, Clone)]
pub struct AromaticityModel {
    aromatic_atoms: HashSet<AtomIdx>,
    aromatic_bonds: HashSet<BondIdx>,
    antiaromatic_rings: Vec<Vec<AtomIdx>>,
    ring_classifications: Vec<(Vec<AtomIdx>, RingAromaticity, u32)>,
}

impl AromaticityModel {
    /// Whether atom `idx` is part of an aromatic ring.
    pub fn is_atom_aromatic(&self, idx: AtomIdx) -> bool {
        self.aromatic_atoms.contains(&idx)
    }

    /// Whether bond `idx` is part of an aromatic ring.
    pub fn is_bond_aromatic(&self, idx: BondIdx) -> bool {
        self.aromatic_bonds.contains(&idx)
    }

    /// Total number of atoms flagged as aromatic.
    pub fn aromatic_atom_count(&self) -> usize {
        self.aromatic_atoms.len()
    }

    /// Get all rings and their classification with electron counts.
    ///
    /// Each entry is `(ring_atoms, classification, π_electron_count)`.
    /// Rings that could not be evaluated (sp3 atoms, unsupported elements) are omitted.
    pub fn ring_classifications(&self) -> &[(Vec<AtomIdx>, RingAromaticity, u32)] {
        &self.ring_classifications
    }

    /// Get all antiaromatic rings (4n electrons, n > 0).
    pub fn antiaromatic_rings(&self) -> &[Vec<AtomIdx>] {
        &self.antiaromatic_rings
    }

    /// Check if any atom belongs to an antiaromatic ring.
    pub fn has_antiaromaticity(&self) -> bool {
        !self.antiaromatic_rings.is_empty()
    }
}

// ---------------------------------------------------------------------------
// Main entry points
// ---------------------------------------------------------------------------

/// Classify a ring by its pi electron count using Hückel and antiaromaticity rules.
#[allow(clippy::manual_is_multiple_of)]
fn classify_ring_aromaticity(pi_electrons: u32) -> (RingAromaticity, u32) {
    if pi_electrons >= 2 && (pi_electrons - 2) % 4 == 0 {
        (RingAromaticity::Aromatic, pi_electrons)
    } else if pi_electrons > 0 && pi_electrons % 4 == 0 {
        (RingAromaticity::Antiaromatic, pi_electrons)
    } else {
        (RingAromaticity::NonAromatic, pi_electrons)
    }
}

/// Mark all atoms and bonds in `ring` as aromatic in the provided sets.
fn mark_ring_aromatic(
    mol: &Molecule,
    ring: &[AtomIdx],
    aromatic_atoms: &mut HashSet<AtomIdx>,
    aromatic_bonds: &mut HashSet<BondIdx>,
) {
    for &atom in ring {
        aromatic_atoms.insert(atom);
    }
    for i in 0..ring.len() {
        let a = ring[i];
        let b = ring[(i + 1) % ring.len()];
        if let Some((bidx, _)) = mol.bond_between(a, b) {
            aromatic_bonds.insert(bidx);
        }
    }
}

/// Assign aromaticity to a molecule using the Hückel 4n+2 rule with fused-ring
/// propagation (Pass 2) and antiaromaticity detection (4n electrons).
///
/// The molecule may be kekulized (`Single`/`Double` bonds) **or** may retain
/// `BondOrder::Aromatic` bonds from the SMILES parser.  In the latter case,
/// aromatic bonds are treated as equivalent to double bonds for electron
/// counting, allowing correct detection without an explicit kekulization step.
///
/// For kekulized input from aromatic SMILES, call `chematic_core::kekulize`
/// then `chematic_core::apply_kekule` first.
pub fn assign_aromaticity(mol: &Molecule) -> AromaticityModel {
    let ring_set = find_sssr(mol);
    let sssr_rings = ring_set.rings();

    // Augment SSSR rings with smaller XOR sub-rings (GF(2) differences between pairs).
    // This corrects the case where the SSSR algorithm stores a large fundamental cycle
    // instead of its smaller GF(2)-reduced equivalent (e.g. the 5-ring of indolizine).
    let rings: Vec<Vec<AtomIdx>> = augmented_ring_set(mol, sssr_rings);

    let mut aromatic_atoms: HashSet<AtomIdx> = HashSet::new();
    let mut aromatic_bonds: HashSet<BondIdx> = HashSet::new();
    let mut antiaromatic_rings: Vec<Vec<AtomIdx>> = Vec::new();

    // Per-ring classification: None means "not yet evaluated / indeterminate".
    let mut classifications: Vec<Option<(RingAromaticity, u32)>> = vec![None; rings.len()];

    // Indices of rings that are candidates for Pass 2 re-evaluation
    // (returned None or NonAromatic in Pass 1).
    let mut pass2_candidates: Vec<usize> = Vec::new();

    // ----- Pass 1: independent Hückel per ring -----
    let empty_context = HashSet::new();
    for (ring_idx, ring) in rings.iter().enumerate() {
        match ring_pi_electrons(mol, ring, &empty_context) {
            Some(pi) => {
                let (cls, count) = classify_ring_aromaticity(pi);
                classifications[ring_idx] = Some((cls, count));
                match cls {
                    RingAromaticity::Aromatic => {
                        mark_ring_aromatic(mol, ring, &mut aromatic_atoms, &mut aromatic_bonds);
                    }
                    RingAromaticity::Antiaromatic => {
                        antiaromatic_rings.push(ring.to_vec());
                        // Antiaromatic is definitive — do not retry in Pass 2.
                    }
                    RingAromaticity::NonAromatic => {
                        pass2_candidates.push(ring_idx);
                    }
                }
            }
            None => {
                // Indeterminate (sp3 atoms, unsupported elements, etc.).
                pass2_candidates.push(ring_idx);
            }
        }
    }

    // ----- Pass 2: propagate through fused ring systems -----
    // Re-evaluate rings adjacent to already-aromatic rings.  Repeat until
    // convergence (no newly aromatic ring found in the last iteration).
    loop {
        let mut any_new = false;
        let mut still_pending: Vec<usize> = Vec::new();

        for ring_idx in pass2_candidates {
            let ring = &rings[ring_idx];
            // Only rings that share an atom with an already-aromatic ring qualify.
            if !ring.iter().any(|a| aromatic_atoms.contains(a)) {
                still_pending.push(ring_idx);
                continue;
            }
            match ring_pi_electrons(mol, ring, &aromatic_atoms) {
                Some(pi) => {
                    let (cls, count) = classify_ring_aromaticity(pi);
                    classifications[ring_idx] = Some((cls, count));
                    if matches!(cls, RingAromaticity::Aromatic) {
                        mark_ring_aromatic(mol, ring, &mut aromatic_atoms, &mut aromatic_bonds);
                        any_new = true;
                    }
                    // NonAromatic even in Pass 2 context: do not retry further.
                }
                None => {
                    still_pending.push(ring_idx);
                }
            }
        }

        pass2_candidates = still_pending;
        if !any_new {
            break;
        }
    }

    // Build the public ring_classifications list (SSSR rings only, omitting augmented/indeterminate).
    let ring_classifications: Vec<(Vec<AtomIdx>, RingAromaticity, u32)> = rings
        .iter()
        .take(sssr_rings.len()) // only expose SSSR rings in the public API
        .enumerate()
        .filter_map(|(i, ring)| {
            classifications[i].map(|(cls, count)| (ring.to_vec(), cls, count))
        })
        .collect();

    AromaticityModel {
        aromatic_atoms,
        aromatic_bonds,
        antiaromatic_rings,
        ring_classifications,
    }
}

/// Apply aromaticity perception to a molecule.
///
/// Returns a new [`Molecule`] where atoms in Hückel-aromatic rings have
/// `atom.aromatic = true` and their bonds carry [`BondOrder::Aromatic`].
/// Non-aromatic atoms and bonds are unchanged.
///
/// The input may be kekulized (no `Aromatic` bond orders) or may retain
/// aromatic bond orders from the SMILES parser.
pub fn apply_aromaticity(mol: &Molecule) -> Molecule {
    use chematic_core::{BondOrder, MoleculeBuilder};

    let model = assign_aromaticity(mol);
    let mut builder = MoleculeBuilder::new();

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
        let _ = builder.add_bond(bond.atom1, bond.atom2, order);
    }
    builder.build()
}

// ---------------------------------------------------------------------------
// Ring augmentation (XOR sub-rings)
// ---------------------------------------------------------------------------

/// Return the sorted set of bond indices that form `ring`.
fn ring_bond_set(mol: &Molecule, ring: &[AtomIdx]) -> Vec<BondIdx> {
    let n = ring.len();
    let mut bonds: Vec<BondIdx> = (0..n)
        .filter_map(|i| {
            let a = ring[i];
            let b = ring[(i + 1) % n];
            mol.bond_between(a, b).map(|(bidx, _)| bidx)
        })
        .collect();
    bonds.sort();
    bonds
}

/// Sorted symmetric difference of two sorted slices.
fn bond_sym_diff(a: &[BondIdx], b: &[BondIdx]) -> Vec<BondIdx> {
    let mut result: Vec<BondIdx> = Vec::new();
    let mut i = 0;
    let mut j = 0;
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => { result.push(a[i]); i += 1; }
            std::cmp::Ordering::Greater => { result.push(b[j]); j += 1; }
            std::cmp::Ordering::Equal => { i += 1; j += 1; }
        }
    }
    result.extend_from_slice(&a[i..]);
    result.extend_from_slice(&b[j..]);
    result
}

/// Reconstruct an ordered atom sequence from a set of bond indices forming a simple cycle.
/// Returns `None` if the bonds do not form a valid simple cycle.
fn ring_atoms_from_bond_set(mol: &Molecule, bonds: &[BondIdx]) -> Option<Vec<AtomIdx>> {
    if bonds.is_empty() {
        return None;
    }
    let mut adj: std::collections::HashMap<AtomIdx, [Option<AtomIdx>; 2]> =
        std::collections::HashMap::new();
    for &bidx in bonds {
        let bond = mol.bond(bidx);
        for (a, b) in [(bond.atom1, bond.atom2), (bond.atom2, bond.atom1)] {
            let e = adj.entry(a).or_insert([None; 2]);
            if e[0].is_none() {
                e[0] = Some(b);
            } else if e[1].is_none() {
                e[1] = Some(b);
            } else {
                return None; // degree > 2 — not a simple ring
            }
        }
    }
    // All atoms must have exactly 2 neighbours.
    if adj.values().any(|e| e[1].is_none()) {
        return None;
    }
    let start = *adj.keys().next()?;
    let mut path = vec![start];
    let mut prev = start;
    let mut current = adj[&start][0]?;
    while current != start {
        path.push(current);
        let [n0, n1] = adj[&current];
        let next = if n0 == Some(prev) { n1? } else { n0? };
        prev = current;
        current = next;
    }
    if path.len() != bonds.len() {
        return None;
    }
    Some(path)
}

/// Augment the SSSR ring list with smaller XOR sub-rings found by pairwise GF(2)
/// differences between SSSR rings that share atoms.
///
/// The standard SSSR algorithm sometimes stores a large fundamental cycle rather
/// than its smaller GF(2)-reduced equivalent (e.g. the 5-ring of indolizine is
/// the XOR of the 6-ring and the 9-ring the algorithm reports).
/// This augmentation adds such missing smaller rings so that aromaticity
/// perception works on the correct smallest rings without modifying the SSSR.
///
/// The returned `Vec` starts with all SSSR rings in their original order; any
/// additional sub-rings derived by GF(2) pairwise XOR follow.  The function
/// only adds a ring if it is strictly smaller than *both* parents, ensuring
/// that envelope rings (e.g. the 10-membered perimeter of naphthalene) are
/// never introduced.
pub fn augmented_ring_set(mol: &Molecule, sssr_rings: &[Vec<AtomIdx>]) -> Vec<Vec<AtomIdx>> {
    let mut rings: Vec<Vec<AtomIdx>> = sssr_rings.to_vec();

    // Build bond sets for all SSSR rings.
    let bond_sets: Vec<Vec<BondIdx>> = sssr_rings
        .iter()
        .map(|r| ring_bond_set(mol, r))
        .collect();

    let n = sssr_rings.len();
    // Track which atom-sets we already have (as sorted atom lists).
    let mut known: std::collections::HashSet<Vec<AtomIdx>> = sssr_rings
        .iter()
        .map(|r| { let mut s = r.clone(); s.sort(); s })
        .collect();

    for i in 0..n {
        for j in (i + 1)..n {
            // Only consider pairs that share atoms (fused rings).
            let shares_atom = sssr_rings[i].iter().any(|a| sssr_rings[j].contains(a));
            if !shares_atom {
                continue;
            }
            let xor_bonds = bond_sym_diff(&bond_sets[i], &bond_sets[j]);
            // Only interesting if the XOR ring is smaller than both parents.
            if xor_bonds.len() >= sssr_rings[i].len().min(sssr_rings[j].len()) {
                continue;
            }
            if let Some(new_ring) = ring_atoms_from_bond_set(mol, &xor_bonds) {
                let mut key = new_ring.clone();
                key.sort();
                if known.insert(key) {
                    rings.push(new_ring);
                }
            }
        }
    }

    rings
}

/// Count aromatic rings in `mol`.
///
/// Uses the augmented ring set so that fused systems where the SSSR stores a
/// large fundamental cycle (e.g. a 9-ring for indolizine) still report the
/// correct small-ring count.  After building the augmented set, any "envelope"
/// ring — one that equals the bond-symmetric-difference of two smaller aromatic
/// rings — is excluded, preventing double-counting.
pub fn count_aromatic_rings(mol: &Molecule) -> usize {
    let sssr = crate::sssr::find_sssr(mol);
    let aug = augmented_ring_set(mol, sssr.rings());

    // Keep only rings where every atom carries the aromatic flag.
    let aromatic: Vec<Vec<AtomIdx>> = aug
        .into_iter()
        .filter(|ring| ring.iter().all(|&idx| mol.atom(idx).aromatic))
        .collect();

    if aromatic.len() <= 1 {
        return aromatic.len();
    }

    // Build sorted bond-index sets for each aromatic ring.
    let bond_sets: Vec<Vec<BondIdx>> = aromatic
        .iter()
        .map(|r| ring_bond_set(mol, r))
        .collect();

    // Mark rings that are the XOR of two strictly smaller aromatic rings.
    // Such rings are "envelope" cycles introduced when the SSSR chose a large
    // fundamental cycle instead of its two smaller GF(2) components.
    let n = aromatic.len();
    let mut is_envelope = vec![false; n];
    for i in 0..n {
        let si = aromatic[i].len();
        'jk: for j in 0..n {
            if j == i || aromatic[j].len() >= si {
                continue;
            }
            for k in (j + 1)..n {
                if k == i || aromatic[k].len() >= si {
                    continue;
                }
                let xor = bond_sym_diff(&bond_sets[j], &bond_sets[k]);
                if xor == bond_sets[i] {
                    is_envelope[i] = true;
                    break 'jk;
                }
            }
        }
    }

    is_envelope.iter().filter(|&&e| !e).count()
}

// ---------------------------------------------------------------------------
// Per-ring pi electron count
// ---------------------------------------------------------------------------

/// Count pi electrons for a ring atom, returning `None` if the atom is
/// incompatible with aromaticity (e.g. sp3 carbon).
///
/// `aromatic_context`: atoms already confirmed aromatic (from Pass 1 or a
/// previous Pass 2 iteration).  Such atoms contribute 1π unconditionally,
/// without requiring an explicit double bond.
///
/// Rules:
/// - **C**: `has_double_any` (Double or Aromatic bond anywhere) → 1π, else None.
///   If already in `aromatic_context` → 1π (confirmed sp2).
/// - **N**:
///   1. Has H → 2π (pyrrole-type lone pair).
///   2. Has an explicit `Double` bond → 1π (pyridine-type).
///   3. Bridgehead N: total_degree == 3 AND ring_degree < total_degree AND no
///      explicit double bond → 2π (lone pair in p orbital, like indolizine N).
///   4. Has in-ring `Aromatic` bond → 1π (pyridine-like aromatic N).
///   5. Already in `aromatic_context` → 1π.
///   6. Otherwise → None.
/// - **O/S**: ring_degree must be 2; contributes 2π (lone pair).
/// - **Other elements**: None (unsupported).
fn ring_pi_electrons(
    mol: &Molecule,
    ring: &[AtomIdx],
    aromatic_context: &HashSet<AtomIdx>,
) -> Option<u32> {
    let ring_atom_set: HashSet<AtomIdx> = ring.iter().copied().collect();
    let mut total_pi: u32 = 0;

    for &atom_idx in ring {
        // Atoms already confirmed aromatic in an adjacent ring contribute 1π.
        if aromatic_context.contains(&atom_idx) {
            total_pi += 1;
            continue;
        }

        let atom = mol.atom(atom_idx);
        let an = atom.element.atomic_number();

        let ring_degree = mol
            .neighbors(atom_idx)
            .filter(|(nb, _)| ring_atom_set.contains(nb))
            .count();

        let total_degree = mol.degree(atom_idx);

        // Explicit Double bond anywhere (not counting Aromatic).
        let has_explicit_double = mol
            .neighbors(atom_idx)
            .any(|(_, bidx)| mol.bond(bidx).order == BondOrder::Double);

        // Double OR Aromatic bond anywhere (for C sp2 check).
        let has_double_any = has_explicit_double
            || mol
                .neighbors(atom_idx)
                .any(|(_, bidx)| mol.bond(bidx).order == BondOrder::Aromatic);

        // Aromatic bond within the ring (for pyridine-like N in aromatic SMILES).
        let has_aromatic_in_ring = mol
            .neighbors(atom_idx)
            .filter(|(nb, _)| ring_atom_set.contains(nb))
            .any(|(_, bidx)| mol.bond(bidx).order == BondOrder::Aromatic);

        let pi = match an {
            // Carbon: must be sp2 (has a double or aromatic bond somewhere).
            6 => {
                if !has_double_any {
                    return None; // sp3 carbon — ring cannot be aromatic
                }
                1
            }

            // Nitrogen
            7 => {
                if implicit_hcount(mol, atom_idx) > 0 {
                    // Pyrrole-type N with H: lone pair → 2π.
                    2
                } else if has_explicit_double {
                    // Pyridine-type N with explicit double bond → 1π.
                    1
                } else if total_degree == 3 && ring_degree < total_degree {
                    // Bridgehead N (e.g. indolizine): no H, no explicit double bond,
                    // and all three σ-bonds exactly fill the N valence (3).
                    // The lone pair occupies the p orbital → 2π (pyrrole-analogue).
                    2
                } else if has_aromatic_in_ring {
                    // N in an aromatic ring (pre-kekulization input) without an
                    // explicit double bond and not a bridgehead → pyridine-like → 1π.
                    1
                } else {
                    // Cannot determine pi contribution.
                    return None;
                }
            }

            // Oxygen / sulfur: lone-pair donor, must be 2-connected in the ring.
            8 | 16 => {
                if ring_degree != 2 {
                    return None;
                }
                2
            }

            // Unsupported element.
            _ => return None,
        };

        total_pi += pi;
    }

    Some(total_pi)
}


// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};

    // =========================================================================
    // Molecule builder helpers (kekulized, manually constructed)
    // =========================================================================

    fn benzene_kekule() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..6).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        for i in 0..6 {
            let order = if i % 2 == 0 {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
            b.add_bond(atoms[i], atoms[(i + 1) % 6], order).unwrap();
        }
        b.build()
    }

    fn cyclohexane() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..6).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        for i in 0..6 {
            b.add_bond(atoms[i], atoms[(i + 1) % 6], BondOrder::Single)
                .unwrap();
        }
        b.build()
    }

    fn pyridine_kekule() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let n = b.add_atom(Atom::new(Element::N));
        let atoms_c: Vec<_> = (0..5).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        let ring = [
            n, atoms_c[0], atoms_c[1], atoms_c[2], atoms_c[3], atoms_c[4],
        ];
        for i in 0..6 {
            let order = if i % 2 == 0 {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
            b.add_bond(ring[i], ring[(i + 1) % 6], order).unwrap();
        }
        b.build()
    }

    fn furan_kekule() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let o = b.add_atom(Atom::new(Element::O));
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let c3 = b.add_atom(Atom::new(Element::C));
        let c4 = b.add_atom(Atom::new(Element::C));
        let ring = [o, c1, c2, c3, c4];
        b.add_bond(ring[0], ring[1], BondOrder::Single).unwrap();
        b.add_bond(ring[1], ring[2], BondOrder::Double).unwrap();
        b.add_bond(ring[2], ring[3], BondOrder::Single).unwrap();
        b.add_bond(ring[3], ring[4], BondOrder::Double).unwrap();
        b.add_bond(ring[4], ring[0], BondOrder::Single).unwrap();
        b.build()
    }

    fn pyrrole_kekule() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let mut n_atom = Atom::new(Element::N);
        n_atom.hydrogen_count = Some(1);
        let n = b.add_atom(n_atom);
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let c3 = b.add_atom(Atom::new(Element::C));
        let c4 = b.add_atom(Atom::new(Element::C));
        let ring = [n, c1, c2, c3, c4];
        b.add_bond(ring[0], ring[1], BondOrder::Single).unwrap();
        b.add_bond(ring[1], ring[2], BondOrder::Double).unwrap();
        b.add_bond(ring[2], ring[3], BondOrder::Single).unwrap();
        b.add_bond(ring[3], ring[4], BondOrder::Double).unwrap();
        b.add_bond(ring[4], ring[0], BondOrder::Single).unwrap();
        b.build()
    }

    fn naphthalene_kekule() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..10).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        let ring1 = [0usize, 1, 2, 3, 4, 9];
        let orders1 = [
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
        ];
        for i in 0..6 {
            b.add_bond(atoms[ring1[i]], atoms[ring1[(i + 1) % 6]], orders1[i])
                .unwrap();
        }
        let ring2_extra = [(4, 5), (5, 6), (6, 7), (7, 8), (8, 9)];
        let orders2 = [
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
        ];
        for (i, &(a, bb)) in ring2_extra.iter().enumerate() {
            b.add_bond(atoms[a], atoms[bb], orders2[i]).unwrap();
        }
        b.build()
    }

    fn cyclobutadiene_kekule() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..4).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        for i in 0..4 {
            let order = if i % 2 == 0 {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
            b.add_bond(atoms[i], atoms[(i + 1) % 4], order).unwrap();
        }
        b.build()
    }

    fn cyclooctatetraene_kekule() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..8).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        for i in 0..8 {
            let order = if i % 2 == 0 {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
            b.add_bond(atoms[i], atoms[(i + 1) % 8], order).unwrap();
        }
        b.build()
    }

    /// Helper: parse an aromatic SMILES and return the molecule with aromatic bonds
    /// (no kekulization).  Use for compounds where kekulization is unsupported.
    #[cfg(test)]
    fn mol_aromatic(smiles: &str) -> chematic_core::Molecule {
        chematic_smiles::parse(smiles).expect("valid SMILES")
    }

    /// Helper: parse SMILES and kekulize.  Panics if kekulization fails.
    #[cfg(test)]
    fn mol_kekulized(smiles: &str) -> chematic_core::Molecule {
        let mol = chematic_smiles::parse(smiles).expect("valid SMILES");
        let k = chematic_core::kekulize(&mol).expect("kekulizable");
        chematic_core::apply_kekule(&mol, &k)
    }

    // =========================================================================
    // Regression: kekulized single-ring aromatics (Pass 1 only, no context)
    // =========================================================================

    #[test]
    fn test_benzene_is_aromatic() {
        let mol = benzene_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 6, "all 6 benzene atoms aromatic");
        for i in 0..6u32 {
            assert!(model.is_atom_aromatic(AtomIdx(i)));
        }
    }

    #[test]
    fn test_cyclohexane_not_aromatic() {
        let mol = cyclohexane();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 0, "cyclohexane not aromatic");
    }

    #[test]
    fn test_pyridine_is_aromatic() {
        let mol = pyridine_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 6);
    }

    #[test]
    fn test_furan_is_aromatic() {
        let mol = furan_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 5);
    }

    #[test]
    fn test_pyrrole_is_aromatic() {
        let mol = pyrrole_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 5);
    }

    #[test]
    fn test_naphthalene_both_rings_aromatic() {
        let mol = naphthalene_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 10, "all 10 naphthalene atoms aromatic");
    }

    #[test]
    fn test_bond_aromaticity_benzene() {
        let mol = benzene_kekule();
        let model = assign_aromaticity(&mol);
        let count = mol.bonds().filter(|(b, _)| model.is_bond_aromatic(*b)).count();
        assert_eq!(count, 6);
    }

    #[test]
    fn test_apply_aromaticity_benzene() {
        let mol = benzene_kekule();
        let aromatic = apply_aromaticity(&mol);
        for (_, atom) in aromatic.atoms() {
            assert!(atom.aromatic, "every benzene carbon should be aromatic");
        }
        let aromatic_bond_count = aromatic
            .bonds()
            .filter(|(_, b)| b.order == BondOrder::Aromatic)
            .count();
        assert_eq!(aromatic_bond_count, 6);
    }

    #[test]
    fn test_apply_aromaticity_cyclohexane_unchanged() {
        let mol = cyclohexane();
        let result = apply_aromaticity(&mol);
        for (_, atom) in result.atoms() {
            assert!(!atom.aromatic);
        }
        for (_, bond) in result.bonds() {
            assert_ne!(bond.order, BondOrder::Aromatic);
        }
    }

    // =========================================================================
    // Antiaromaticity
    // =========================================================================

    #[test]
    fn test_cyclobutadiene_antiaromatic() {
        let mol = cyclobutadiene_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 0, "cyclobutadiene not aromatic");
        assert!(model.has_antiaromaticity(), "cyclobutadiene antiaromatic");
        assert_eq!(model.antiaromatic_rings().len(), 1);
        let classifications = model.ring_classifications();
        assert_eq!(classifications.len(), 1);
        assert_eq!(classifications[0].1, RingAromaticity::Antiaromatic);
        assert_eq!(classifications[0].2, 4);
    }

    #[test]
    fn test_cyclooctatetraene_antiaromatic() {
        let mol = cyclooctatetraene_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 0, "COT not aromatic");
        assert!(model.has_antiaromaticity(), "COT antiaromatic");
        assert_eq!(model.antiaromatic_rings().len(), 1);
        let cls = &model.ring_classifications()[0];
        assert_eq!(cls.1, RingAromaticity::Antiaromatic);
        assert_eq!(cls.2, 8);
    }

    // =========================================================================
    // Ring classifications
    // =========================================================================

    #[test]
    fn test_ring_classifications_benzene() {
        let mol = benzene_kekule();
        let model = assign_aromaticity(&mol);
        let classifications = model.ring_classifications();
        assert_eq!(classifications.len(), 1);
        assert_eq!(classifications[0].1, RingAromaticity::Aromatic);
        assert_eq!(classifications[0].2, 6);
    }

    #[test]
    fn test_ring_classifications_naphthalene() {
        let mol = naphthalene_kekule();
        let model = assign_aromaticity(&mol);
        let classifications = model.ring_classifications();
        assert_eq!(classifications.len(), 2, "naphthalene has two rings");
        for (_, classification, count) in classifications {
            assert_eq!(*classification, RingAromaticity::Aromatic);
            assert_eq!(*count, 6);
        }
    }

    #[test]
    fn test_non_aromatic_cyclohexane() {
        let mol = cyclohexane();
        let model = assign_aromaticity(&mol);
        for (_, classification, _) in model.ring_classifications() {
            assert_ne!(*classification, RingAromaticity::Aromatic);
            assert_ne!(*classification, RingAromaticity::Antiaromatic);
        }
    }

    // =========================================================================
    // Electron distribution
    // =========================================================================

    #[test]
    fn test_thiophene_aromatic() {
        let mut b = MoleculeBuilder::new();
        let s = b.add_atom(Atom::new(Element::S));
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let c3 = b.add_atom(Atom::new(Element::C));
        let c4 = b.add_atom(Atom::new(Element::C));
        let ring = [s, c1, c2, c3, c4];
        b.add_bond(ring[0], ring[1], BondOrder::Single).unwrap();
        b.add_bond(ring[1], ring[2], BondOrder::Double).unwrap();
        b.add_bond(ring[2], ring[3], BondOrder::Single).unwrap();
        b.add_bond(ring[3], ring[4], BondOrder::Double).unwrap();
        b.add_bond(ring[4], ring[0], BondOrder::Single).unwrap();
        let mol = b.build();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 5);
        assert_eq!(model.ring_classifications()[0].2, 6);
    }

    #[test]
    fn test_electron_distribution_tracking() {
        let mol = benzene_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.ring_classifications()[0].2, 6, "benzene: 6 × 1π = 6");

        let mol = pyrrole_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.ring_classifications()[0].2,
            6,
            "pyrrole: N(2π) + 4C(1π) = 6"
        );

        let mol = furan_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.ring_classifications()[0].2,
            6,
            "furan: O(2π) + 4C(1π) = 6"
        );
    }

    // =========================================================================
    // Aromatic-SMILES input (BondOrder::Aromatic, no kekulization)
    // Verifies that assign_aromaticity works on pre-kekulization molecules.
    // =========================================================================

    #[test]
    fn test_benzene_aromatic_smiles() {
        // c1ccccc1 — parsed with BondOrder::Aromatic bonds
        let mol = mol_aromatic("c1ccccc1");
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 6, "benzene from aromatic SMILES");
    }

    #[test]
    fn test_naphthalene_aromatic_smiles() {
        let mol = mol_aromatic("c1ccc2ccccc2c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            10,
            "naphthalene from aromatic SMILES"
        );
    }

    #[test]
    fn test_pyridine_aromatic_smiles() {
        let mol = mol_aromatic("c1ccncc1");
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 6, "pyridine from aromatic SMILES");
    }

    #[test]
    fn test_furan_aromatic_smiles() {
        let mol = mol_aromatic("c1ccoc1");
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 5, "furan from aromatic SMILES");
    }

    #[test]
    fn test_pyrrole_aromatic_smiles() {
        // [nH] bracket atom: hydrogen_count = Some(1)
        let mol = mol_aromatic("c1cc[nH]c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 5, "pyrrole from aromatic SMILES");
    }

    #[test]
    fn test_thiophene_aromatic_smiles() {
        let mol = mol_aromatic("c1ccsc1");
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 5, "thiophene from aromatic SMILES");
    }

    // =========================================================================
    // Fused-ring kekulized systems (Pass 2 propagation)
    // =========================================================================

    #[test]
    fn test_indole_aromatic() {
        // c1ccc2[nH]ccc2c1 — indole (9 atoms, 5-ring + 6-ring fused)
        let mol = mol_kekulized("c1ccc2[nH]ccc2c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            9,
            "all 9 indole atoms aromatic"
        );
    }

    #[test]
    fn test_benzimidazole_aromatic() {
        // Two N atoms in fused 5+6 ring system
        let mol = mol_kekulized("c1ccc2[nH]cnc2c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 9, "all 9 benzimidazole atoms");
    }

    #[test]
    fn test_quinoline_aromatic() {
        let mol = mol_kekulized("c1ccc2ncccc2c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 10, "all 10 quinoline atoms");
    }

    #[test]
    fn test_acridine_aromatic() {
        // 3 fused 6-membered rings, central N: 13 atoms
        let mol = mol_kekulized("c1ccc2nc3ccccc3cc2c1");
        let model = assign_aromaticity(&mol);
        // acridine is C13H9N → 14 heavy atoms (13 C + 1 N), all aromatic
        assert_eq!(model.aromatic_atom_count(), 14, "all 14 acridine atoms");
    }

    // =========================================================================
    // Fused-ring aromatic-SMILES input (BondOrder::Aromatic, kekulize fails)
    // =========================================================================

    #[test]
    fn test_indolizine_aromatic() {
        // c1ccn2cccc2c1 — indolizine: bridgehead N, kekulization unsupported.
        // The SSSR finds a 6-ring and a 9-ring; the 5-ring is recovered via
        // augmentation (XOR of 6- and 9-ring).
        // Pass 1: 5-ring (augmented) detected via bridgehead-N rule → 6π.
        // Pass 2: 6-ring detected using N already aromatic from 5-ring → 6π.
        // The 9-ring (SSSR artifact) is NonAromatic (9π ≠ 4n+2), but all
        // 9 atoms are correctly flagged aromatic via the 5- and 6-ring.
        let mol = mol_aromatic("c1ccn2cccc2c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            9,
            "all 9 indolizine atoms aromatic"
        );
        // At least the 6-ring should be classified as Aromatic in the SSSR set.
        let has_aromatic_ring = model
            .ring_classifications()
            .iter()
            .any(|(_, cls, _)| *cls == RingAromaticity::Aromatic);
        assert!(has_aromatic_ring, "at least one SSSR ring aromatic");
    }

    #[test]
    fn test_purine_aromatic() {
        // c1cnc2[nH]cnc2n1 — purine: 9 atoms, kekulizable
        let mol = mol_kekulized("c1cnc2[nH]cnc2n1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            9,
            "all 9 purine atoms aromatic"
        );
    }

    #[test]
    fn test_purine_aromatic_from_aromatic_smiles() {
        let mol = mol_aromatic("c1cnc2[nH]cnc2n1");
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 9, "purine from aromatic SMILES");
    }

    #[test]
    fn test_2_pyridinone_aromatic() {
        // O=c1ccncc1 — 2-pyridinone (aromatic SMILES, N without H, exo C=O).
        // Kekulization fails; tested on the aromatic-bond form directly.
        // The exo C=O gives the C atom has_double_any=true → 1π.
        // N has Aromatic bonds in ring → 1π (pyridine-like).
        // Total: 6 × 1π = 6π → aromatic.
        let mol = mol_aromatic("O=c1ccncc1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            6,
            "all 6 ring atoms of 2-pyridinone aromatic"
        );
    }

    #[test]
    fn test_quinolone_aromatic() {
        // O=c1ccc2ncccc2c1 — quinolone: fused 6+6 with exo C=O, kekulize fails
        let mol = mol_aromatic("O=c1ccc2ncccc2c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            10,
            "all 10 quinolone ring atoms aromatic"
        );
        assert_eq!(
            model.ring_classifications().len(),
            2,
            "two rings classified"
        );
    }

    #[test]
    fn test_indole_aromatic_smiles() {
        let mol = mol_aromatic("c1ccc2[nH]ccc2c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 9, "indole from aromatic SMILES");
    }

    // =========================================================================
    // Bridgehead N rule: specifically test that the rule fires correctly
    // =========================================================================

    #[test]
    fn test_bridgehead_n_contributes_lone_pair() {
        // Indolizine: the bridgehead N (degree 3, no H, no explicit double bond)
        // must be detected as a 2π contributor for the 5-membered ring.
        // We verify by checking the 5-ring classification (if accessible).
        let mol = mol_aromatic("c1ccn2cccc2c1");
        let model = assign_aromaticity(&mol);
        // All 9 atoms aromatic: both rings must be aromatic.
        assert_eq!(model.aromatic_atom_count(), 9);
        // The bridgehead N itself must be in the aromatic set.
        // In the SMILES c1ccn2cccc2c1, n is atom index 3.
        assert!(model.is_atom_aromatic(AtomIdx(3)), "bridgehead N must be aromatic");
    }

    #[test]
    fn test_non_bridgehead_n_no_false_positive() {
        // Pyrimidine: two N atoms in a 6-membered ring, no bridgehead.
        // Both N have ring_degree == total_degree == 2.
        // Should be detected as aromatic via has_aromatic_in_ring (Aromatic bonds).
        let mol = mol_aromatic("c1ccncn1");
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 6, "pyrimidine is aromatic");
    }

    #[test]
    fn test_imidazole_aromatic() {
        // c1cn[nH]c1 / c1c[nH]cn1 — imidazole: one pyridine-type N, one pyrrole-type N
        let mol = mol_aromatic("c1cn[nH]c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 5, "imidazole is aromatic");
    }

    // =========================================================================
    // Pass 2 specifically: rings that need fused-ring context
    // =========================================================================

    #[test]
    fn test_pass2_needed_for_indolizine_6ring() {
        // The augmented 5-ring (XOR of SSSR 6-ring and 9-ring) is detected aromatic in Pass 1.
        // The SSSR 6-ring is then detected aromatic in Pass 2 (N already aromatic → 1π).
        // The SSSR 9-ring (9π) remains NonAromatic per Hückel.
        // Key assertion: all 9 atoms are aromatic (correct overall perception).
        let mol = mol_aromatic("c1ccn2cccc2c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 9, "all 9 indolizine atoms aromatic");
        // The bridgehead N must be aromatic.
        assert!(model.is_atom_aromatic(AtomIdx(3)), "bridgehead N is aromatic");
        // The 6-ring (SSSR ring, improved by Pass 2) should be classified Aromatic.
        let aromatic_count = model
            .ring_classifications()
            .iter()
            .filter(|(_, cls, _)| *cls == RingAromaticity::Aromatic)
            .count();
        assert!(aromatic_count >= 1, "at least one SSSR ring is aromatic");
    }

    #[test]
    fn test_no_pass2_needed_for_naphthalene() {
        // Naphthalene: both rings pass independently in Pass 1.
        // Verifies Pass 2 doesn't break things that already work.
        let mol = naphthalene_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 10);
        let classes = model.ring_classifications();
        assert_eq!(classes.len(), 2);
        for (_, cls, _) in classes {
            assert_eq!(*cls, RingAromaticity::Aromatic);
        }
    }

    #[test]
    fn test_anthracene_aromatic() {
        // c1ccc2cc3ccccc3cc2c1 — anthracene: 3 linearly fused 6-rings, 14 atoms
        let mol = mol_kekulized("c1ccc2cc3ccccc3cc2c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 14, "all 14 anthracene atoms");
    }

    // =========================================================================
    // Regression: aromatic-bond path must not perturb kekulized correctness
    // =========================================================================

    #[test]
    fn test_kekulized_path_unaffected_by_aromatic_bond_changes() {
        // Kekulized benzene: bonds are Double/Single, not Aromatic.
        // The new Aromatic-bond branches must stay dormant.
        let mol = benzene_kekule();
        // Verify no aromatic bonds in input.
        for (_, bond) in mol.bonds() {
            assert_ne!(bond.order, BondOrder::Aromatic, "input must be kekulized");
        }
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 6);
        // All 6 bonds in benzene ring should be aromatic.
        let aromatic_bonds = mol
            .bonds()
            .filter(|(b, _)| model.is_bond_aromatic(*b))
            .count();
        assert_eq!(aromatic_bonds, 6);
    }

    #[test]
    fn test_keto_pyridinone_not_huckel_aromatic() {
        // O=C1NC=CC=C1 — 2-pyridinone keto form with N-H.
        // π count: C(=O)(1π) + N-H(2π) + 3×C(1π each) + C6(1π) = 7π → not 4n+2.
        // This is a known scope boundary: keto pyridinone has partial aromatic
        // character by resonance but is NOT Hückel 4n+2 aromatic.
        // RDKit classifies it aromatic using an extended model; our strict Hückel
        // implementation correctly returns 0 aromatic atoms.
        let mol = mol_kekulized("O=C1NC=CC=C1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            0,
            "keto pyridinone is not Hückel aromatic (7π ≠ 4n+2)"
        );
    }

    #[test]
    fn test_cyclopentadienyl_not_aromatic_kekulized() {
        // C1=CC=CC1 — cyclopentadiene (4 C with doubles + 1 sp3 CH2): not aromatic.
        let mut b = MoleculeBuilder::new();
        let c0 = b.add_atom(Atom::new(Element::C)); // sp3
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let c3 = b.add_atom(Atom::new(Element::C));
        let c4 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c0, c1, BondOrder::Single).unwrap();
        b.add_bond(c1, c2, BondOrder::Double).unwrap();
        b.add_bond(c2, c3, BondOrder::Single).unwrap();
        b.add_bond(c3, c4, BondOrder::Double).unwrap();
        b.add_bond(c4, c0, BondOrder::Single).unwrap();
        let mol = b.build();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 0, "cyclopentadiene not aromatic");
    }
}

