//! Hückel aromaticity perception with antiaromaticity detection.
//!
//! Works on a **kekulized** molecule (no `Aromatic` bond orders).
//! Call `kekulize` + `apply_kekule` from `chematic-core` before calling
//! `assign_aromaticity` if the input contains aromatic bonds.
//!
//! Algorithm:
//! 1. Find all SSSR rings via `find_sssr`.
//! 2. For each ring, check whether it could participate in a planar pi system.
//! 3. Count pi electrons contributed by each ring atom with explicit electron distribution tracking.
//! 4. Classify rings by electron count:
//!    - 4n+2 electrons (n >= 0): aromatic (favorable)
//!    - 4n electrons (n > 0): antiaromatic (unfavorable, strongly disfavored)
//!    - Other: non-aromatic
//! 5. Record all aromatic atoms, bonds, and antiaromatic rings in an `AromaticityModel`.

use std::collections::HashSet;

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule};

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
/// the Hückel 4n+2 rule applied to SSSR rings.
/// Also tracks antiaromatic rings (4n electrons) for chemical accuracy.
#[derive(Debug, Clone)]
pub struct AromaticityModel {
    aromatic_atoms: HashSet<AtomIdx>,
    aromatic_bonds: HashSet<BondIdx>,
    antiaromatic_rings: Vec<Vec<AtomIdx>>, // Rings with 4n electrons
    ring_classifications: Vec<(Vec<AtomIdx>, RingAromaticity, u32)>, // (ring, classification, electron_count)
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
// Main entry point
// ---------------------------------------------------------------------------

/// Classify a ring by its pi electron count using Hückel and antiaromaticity rules.
///
/// Returns (classification, electron_count) tuple.
fn classify_ring_aromaticity(pi_electrons: u32) -> (RingAromaticity, u32) {
    // Hückel rule: 4n+2 electrons → aromatic
    if pi_electrons >= 2 && (pi_electrons - 2) % 4 == 0 {
        (RingAromaticity::Aromatic, pi_electrons)
    }
    // Antiaromaticity: 4n electrons (n > 0) → antiaromatic
    else if pi_electrons > 0 && pi_electrons % 4 == 0 {
        (RingAromaticity::Antiaromatic, pi_electrons)
    }
    // Everything else: non-aromatic
    else {
        (RingAromaticity::NonAromatic, pi_electrons)
    }
}

/// Assign aromaticity to a kekulized molecule using the Hückel 4n+2 rule
/// and antiaromaticity detection (4n electrons).
///
/// The molecule must use `Single` / `Double` bond orders (no `Aromatic` bonds).
/// For input parsed from aromatic SMILES call `chematic_core::kekulize` then
/// `chematic_core::apply_kekule` before passing the molecule here.
pub fn assign_aromaticity(mol: &Molecule) -> AromaticityModel {
    let ring_set = find_sssr(mol);

    let mut aromatic_atoms: HashSet<AtomIdx> = HashSet::new();
    let mut aromatic_bonds: HashSet<BondIdx> = HashSet::new();
    let mut antiaromatic_rings: Vec<Vec<AtomIdx>> = Vec::new();
    let mut ring_classifications: Vec<(Vec<AtomIdx>, RingAromaticity, u32)> = Vec::new();

    for ring in ring_set.rings() {
        if let Some(pi_electrons) = ring_pi_electrons(mol, ring) {
            let (classification, count) = classify_ring_aromaticity(pi_electrons);

            ring_classifications.push((ring.to_vec(), classification, count));

            match classification {
                RingAromaticity::Aromatic => {
                    // Mark all atoms in this ring as aromatic.
                    for &atom in ring {
                        aromatic_atoms.insert(atom);
                    }
                    // Mark all bonds in this ring as aromatic.
                    for i in 0..ring.len() {
                        let a = ring[i];
                        let b = ring[(i + 1) % ring.len()];
                        if let Some((bidx, _)) = mol.bond_between(a, b) {
                            aromatic_bonds.insert(bidx);
                        }
                    }
                }
                RingAromaticity::Antiaromatic => {
                    // Record antiaromatic ring (do NOT mark atoms as aromatic).
                    antiaromatic_rings.push(ring.to_vec());
                }
                RingAromaticity::NonAromatic => {
                    // Neither aromatic nor antiaromatic; do nothing.
                }
            }
        }
    }

    AromaticityModel {
        aromatic_atoms,
        aromatic_bonds,
        antiaromatic_rings,
        ring_classifications,
    }
}

/// Apply aromaticity perception to a kekulized molecule.
///
/// Returns a new [`Molecule`] where atoms in Hückel-aromatic rings have
/// `atom.aromatic = true` and their bonds carry [`BondOrder::Aromatic`].
/// Non-aromatic atoms and bonds are unchanged.
///
/// The input must be kekulized (no `Aromatic` bond orders).  See
/// [`chematic_core::kekulize`] / [`chematic_core::apply_kekule`] if you
/// need to de-aromatize an aromatic-SMILES input first.
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
// Per-ring pi electron count
// ---------------------------------------------------------------------------

/// Detailed electron contribution from a ring atom.
///
/// Tracks how many pi electrons are contributed and why.
/// This struct is kept for documentation purposes of the algorithm.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
struct AtomElectronContribution {
    pi_electrons: u32,
    is_lone_pair: bool,      // true if from lone pair (N-H, O, S)
    is_double_bond: bool,    // true if from double bond
    is_sp2_compatible: bool, // true if atom can participate in pi system
}

/// Try to count pi electrons for a ring with detailed distribution analysis.
///
/// Returns `None` if any atom in the ring is not sp2-compatible (i.e. the ring
/// cannot be aromatic regardless of electron count).
/// Returns `Some(count)` otherwise.
///
/// Algorithm for electron distribution:
/// 1. **Carbon**: Contributes 1 π-electron from sp2 hybridization.
///    - Must have a double bond somewhere (ring or exocyclic).
///    - Without a double bond, it's sp3 → ring cannot be aromatic.
/// 2. **Nitrogen (N)**:
///    - With explicit H (pyrrole-type): Lone pair → 2 π-electrons.
///    - Without H (pyridine-type): 1 π-electron from sp2 (accepts electrons).
///    - Ambiguous cases (no H, no double) → skip ring.
/// 3. **Oxygen/Sulfur (O, S)**: Lone pair donor → 2 π-electrons.
///    - Must be 2-connected in the ring (not linked to exocyclic groups).
/// 4. **Other elements**: Unsupported → ring cannot be aromatic.
fn ring_pi_electrons(mol: &Molecule, ring: &[AtomIdx]) -> Option<u32> {
    let ring_atom_set: HashSet<AtomIdx> = ring.iter().copied().collect();
    let mut total_pi: u32 = 0;

    for &atom_idx in ring {
        let atom = mol.atom(atom_idx);
        let an = atom.element.atomic_number();

        // Count heavy-atom bonds to ring neighbors (ring degree).
        let ring_degree = mol
            .neighbors(atom_idx)
            .filter(|(nb, _)| ring_atom_set.contains(nb))
            .count();

        // Determine whether the atom has a double bond to a ring neighbor.
        let has_double_in_ring = mol
            .neighbors(atom_idx)
            .filter(|(nb, _)| ring_atom_set.contains(nb))
            .any(|(_, bidx)| mol.bond(bidx).order == BondOrder::Double);

        // Determine whether the atom has any double bond (inside or outside the ring).
        let has_double_any = mol
            .neighbors(atom_idx)
            .any(|(_, bidx)| mol.bond(bidx).order == BondOrder::Double);

        let pi = match an {
            // Carbon: must have a double bond somewhere to be sp2.
            6 => {
                if !has_double_any {
                    // sp3 carbon — ring cannot be aromatic.
                    return None;
                }
                1 // One pi electron from the C=C double bond (shared between C atoms).
            }
            // Nitrogen (atomic number 7)
            7 => {
                // Determine H count: use explicit hydrogen_count if set (bracket atom),
                // otherwise derive from valence.
                if hydrogen_count(mol, atom_idx) > 0 {
                    // Pyrrole-type N with H: contributes a lone pair → 2 pi electrons.
                    2
                } else if has_double_in_ring || has_double_any {
                    // Pyridine-type N (or N with exocyclic C=N): contributes 1.
                    1
                } else {
                    // N in ring with no double bond and no H — cannot determine pi system.
                    return None;
                }
            }
            // Oxygen / sulfur: lone-pair donor, contributes 2 (must be 2-connected in ring).
            8 | 16 => {
                if ring_degree != 2 {
                    return None;
                }
                2
            }
            // Unsupported element for Hückel aromaticity.
            _ => return None,
        };

        total_pi += pi;
    }

    Some(total_pi)
}

// ---------------------------------------------------------------------------
// Hydrogen count helper
// ---------------------------------------------------------------------------

/// Return the total hydrogen count for `atom_idx`.
///
/// - Bracket atoms: use the stored `hydrogen_count`.
/// - Organic-subset atoms: derive from the standard valence table.
fn hydrogen_count(mol: &Molecule, atom_idx: AtomIdx) -> u8 {
    let atom = mol.atom(atom_idx);

    // Bracket atoms store the count directly.
    if let Some(h) = atom.hydrogen_count {
        return h;
    }

    // Organic-subset atoms: compute from valence.
    let normal_valences = atom.element.normal_valences();
    if normal_valences.is_empty() {
        return 0;
    }

    let bond_sum: i32 = mol
        .neighbors(atom_idx)
        .map(|(_, bidx)| mol.bond(bidx).order.order_int() as i32)
        .sum();

    let charge = atom.charge as i32;

    for &v in normal_valences {
        let target = v as i32 + charge;
        if target >= bond_sum {
            return (target - bond_sum) as u8;
        }
    }

    0
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};

    // Build a kekulized benzene ring (alternating single/double bonds).
    fn benzene_kekule() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..6).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        for i in 0..6 {
            let order = if i % 2 == 0 { BondOrder::Double } else { BondOrder::Single };
            b.add_bond(atoms[i], atoms[(i + 1) % 6], order).unwrap();
        }
        b.build()
    }

    // Build cyclohexane (6 single bonds — no pi system).
    fn cyclohexane() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..6).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        for i in 0..6 {
            b.add_bond(atoms[i], atoms[(i + 1) % 6], BondOrder::Single).unwrap();
        }
        b.build()
    }

    // Build kekulized pyridine: 5 C + 1 N (pyridine-type, no H).
    // Atoms: N(0)-C(1)=C(2)-C(3)=C(4)-C(5)=N(0) ring
    fn pyridine_kekule() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let n = b.add_atom(Atom::new(Element::N));
        let atoms_c: Vec<_> = (0..5).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        let ring = [n, atoms_c[0], atoms_c[1], atoms_c[2], atoms_c[3], atoms_c[4]];
        // N=C-C=C-C=N  alternating, starting with double
        for i in 0..6 {
            let order = if i % 2 == 0 { BondOrder::Double } else { BondOrder::Single };
            b.add_bond(ring[i], ring[(i + 1) % 6], order).unwrap();
        }
        b.build()
    }

    // Build kekulized furan: O + 4 C, 5-membered ring.
    // O-C=C-C=C-O  (O contributes lone pair, 2 double bonds in ring)
    fn furan_kekule() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let o = b.add_atom(Atom::new(Element::O));
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let c3 = b.add_atom(Atom::new(Element::C));
        let c4 = b.add_atom(Atom::new(Element::C));
        let ring = [o, c1, c2, c3, c4];
        // O-C=C-C=C ring; O-C single, C=C double, C-C single, C=C double, C-O single (back)
        b.add_bond(ring[0], ring[1], BondOrder::Single).unwrap();
        b.add_bond(ring[1], ring[2], BondOrder::Double).unwrap();
        b.add_bond(ring[2], ring[3], BondOrder::Single).unwrap();
        b.add_bond(ring[3], ring[4], BondOrder::Double).unwrap();
        b.add_bond(ring[4], ring[0], BondOrder::Single).unwrap();
        b.build()
    }

    // Build kekulized pyrrole: [NH] + 4 C, 5-membered ring.
    // N has an explicit H (hydrogen_count = Some(1)).
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

    // Build kekulized naphthalene: 10 C, 11 bonds, 2 fused 6-membered rings.
    fn naphthalene_kekule() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..10).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        // Ring 1: 0-1-2-3-4-9
        // Ring 2: 4-5-6-7-8-9 sharing bond 4-9
        // Kekulé pattern: 0=1-2=3-4=9-0, then 4-5=6-7=8-9 (4-9 already single)
        // Bond orders for a valid kekulé:
        //   0-1: double, 1-2: single, 2-3: double, 3-4: single, 4-9: double, 9-0: single
        //   4-5: single, 5-6: double, 6-7: single, 7-8: double, 8-9: single
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
            b.add_bond(atoms[ring1[i]], atoms[ring1[(i + 1) % 6]], orders1[i]).unwrap();
        }
        // Ring 2 extra bonds (4-9 already added):
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

    #[test]
    fn test_benzene_is_aromatic() {
        let mol = benzene_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 6, "all 6 benzene atoms are aromatic");
        for i in 0..6u32 {
            assert!(model.is_atom_aromatic(AtomIdx(i)), "atom {} should be aromatic", i);
        }
    }

    #[test]
    fn test_cyclohexane_not_aromatic() {
        let mol = cyclohexane();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 0, "cyclohexane has no aromatic atoms");
    }

    #[test]
    fn test_pyridine_is_aromatic() {
        let mol = pyridine_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 6, "all 6 pyridine atoms are aromatic");
    }

    #[test]
    fn test_furan_is_aromatic() {
        let mol = furan_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 5, "all 5 furan atoms are aromatic");
    }

    #[test]
    fn test_pyrrole_is_aromatic() {
        let mol = pyrrole_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 5, "all 5 pyrrole atoms are aromatic");
    }

    #[test]
    fn test_naphthalene_both_rings_aromatic() {
        let mol = naphthalene_kekule();
        let model = assign_aromaticity(&mol);
        // Both 6-membered rings are aromatic; all 10 atoms should be flagged.
        assert_eq!(
            model.aromatic_atom_count(),
            10,
            "all 10 naphthalene atoms are aromatic"
        );
    }

    #[test]
    fn test_bond_aromaticity_benzene() {
        let mol = benzene_kekule();
        let model = assign_aromaticity(&mol);
        // All 6 ring bonds should be aromatic.
        let mut count = 0;
        for (bidx, _) in mol.bonds() {
            if model.is_bond_aromatic(bidx) {
                count += 1;
            }
        }
        assert_eq!(count, 6, "benzene has 6 aromatic bonds");
    }

    #[test]
    fn test_apply_aromaticity_benzene() {
        use chematic_core::BondOrder;

        let mol = benzene_kekule();
        let aromatic = apply_aromaticity(&mol);

        // All 6 atoms must have aromatic=true
        for (_, atom) in aromatic.atoms() {
            assert!(atom.aromatic, "every benzene carbon should be aromatic");
        }

        // All 6 bonds must have BondOrder::Aromatic
        let mut aromatic_bond_count = 0;
        for (_, bond) in aromatic.bonds() {
            if bond.order == BondOrder::Aromatic {
                aromatic_bond_count += 1;
            }
        }
        assert_eq!(aromatic_bond_count, 6);
    }

    #[test]
    fn test_apply_aromaticity_cyclohexane_unchanged() {
        use chematic_core::BondOrder;

        let mol = cyclohexane();
        let result = apply_aromaticity(&mol);

        // No atom should gain aromatic flag
        for (_, atom) in result.atoms() {
            assert!(!atom.aromatic);
        }
        // No bond should become Aromatic
        for (_, bond) in result.bonds() {
            assert_ne!(bond.order, BondOrder::Aromatic);
        }
    }

    // =========================================================================
    // Antiaromaticity Tests (4n electrons, n > 0)
    // =========================================================================

    /// Build kekulized cyclobutadiene: 4-membered ring with 4 electrons.
    /// C=C-C=C pattern (4 π electrons → 4n antiaromatic)
    fn cyclobutadiene_kekule() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..4).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        for i in 0..4 {
            let order = if i % 2 == 0 { BondOrder::Double } else { BondOrder::Single };
            b.add_bond(atoms[i], atoms[(i + 1) % 4], order).unwrap();
        }
        b.build()
    }

    #[test]
    fn test_cyclobutadiene_antiaromatic() {
        let mol = cyclobutadiene_kekule();
        let model = assign_aromaticity(&mol);

        // Cyclobutadiene should NOT be aromatic (4 electrons = 4n).
        assert_eq!(model.aromatic_atom_count(), 0, "cyclobutadiene should not be aromatic");

        // But it SHOULD be detected as antiaromatic.
        assert!(model.has_antiaromaticity(), "cyclobutadiene should be antiaromatic");
        assert_eq!(model.antiaromatic_rings().len(), 1, "one antiaromatic ring");

        // Verify ring classification
        let classifications = model.ring_classifications();
        assert_eq!(classifications.len(), 1);
        assert_eq!(classifications[0].1, RingAromaticity::Antiaromatic);
        assert_eq!(classifications[0].2, 4, "cyclobutadiene has 4 π electrons");
    }

    /// Build cyclooctatetraene (COT): 8-membered ring with 8 electrons.
    /// C=C-C=C-C=C-C=C pattern (8 π electrons → 4n antiaromatic)
    fn cyclooctatetraene_kekule() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..8).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        for i in 0..8 {
            let order = if i % 2 == 0 { BondOrder::Double } else { BondOrder::Single };
            b.add_bond(atoms[i], atoms[(i + 1) % 8], order).unwrap();
        }
        b.build()
    }

    #[test]
    fn test_cyclooctatetraene_antiaromatic() {
        let mol = cyclooctatetraene_kekule();
        let model = assign_aromaticity(&mol);

        // COT should NOT be aromatic (8 electrons = 4n).
        assert_eq!(model.aromatic_atom_count(), 0, "cyclooctatetraene should not be aromatic");

        // But it SHOULD be detected as antiaromatic.
        assert!(model.has_antiaromaticity(), "cyclooctatetraene should be antiaromatic");
        assert_eq!(model.antiaromatic_rings().len(), 1, "one antiaromatic ring");

        let classifications = model.ring_classifications();
        assert_eq!(classifications[0].1, RingAromaticity::Antiaromatic);
        assert_eq!(classifications[0].2, 8, "cyclooctatetraene has 8 π electrons");
    }

    // =========================================================================
    // Ring Classification Tests
    // =========================================================================

    #[test]
    fn test_ring_classifications_benzene() {
        let mol = benzene_kekule();
        let model = assign_aromaticity(&mol);
        let classifications = model.ring_classifications();

        assert_eq!(classifications.len(), 1, "benzene has one ring");
        assert_eq!(classifications[0].1, RingAromaticity::Aromatic);
        assert_eq!(classifications[0].2, 6, "benzene has 6 π electrons");
    }

    #[test]
    fn test_ring_classifications_mixed() {
        // Naphthalene: both rings should be aromatic (6 π electrons each).
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
        let classifications = model.ring_classifications();

        // Cyclohexane may or may not be detected by SSSR depending on ring perception.
        // What matters is that it should never be classified as aromatic or antiaromatic.
        for (_, classification, _count) in classifications {
            assert_ne!(*classification, RingAromaticity::Aromatic, "cyclohexane should not be aromatic");
            assert_ne!(*classification, RingAromaticity::Antiaromatic, "cyclohexane should not be antiaromatic");
        }
    }

    // =========================================================================
    // Electron Distribution Edge Cases
    // =========================================================================

    #[test]
    fn test_thiophene_aromatic() {
        // Thiophene: 5-membered ring with S contributing a lone pair.
        // S-C=C-C=C (2 + 1 + 1 + 1 + 1 = 6 π electrons → aromatic)
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
        assert_eq!(model.aromatic_atom_count(), 5, "thiophene is aromatic");

        let classifications = model.ring_classifications();
        assert_eq!(classifications[0].1, RingAromaticity::Aromatic);
        assert_eq!(classifications[0].2, 6);
    }

    #[test]
    fn test_electron_distribution_tracking() {
        // Benzene: all carbons contribute 1 π electron each = 6 total.
        let mol = benzene_kekule();
        let model = assign_aromaticity(&mol);
        let classifications = model.ring_classifications();

        assert_eq!(
            classifications[0].2, 6,
            "benzene electron distribution: 6 × 1π (from C) = 6 total"
        );

        // Pyrrole: N contributes 2 (lone pair), 4 C contribute 1 each = 6 total.
        let mol = pyrrole_kekule();
        let model = assign_aromaticity(&mol);
        let classifications = model.ring_classifications();

        assert_eq!(
            classifications[0].2, 6,
            "pyrrole electron distribution: 2π (N lone pair) + 4 × 1π (C) = 6 total"
        );

        // Furan: O contributes 2 (lone pair), 4 C contribute 1 each = 6 total.
        let mol = furan_kekule();
        let model = assign_aromaticity(&mol);
        let classifications = model.ring_classifications();

        assert_eq!(
            classifications[0].2, 6,
            "furan electron distribution: 2π (O lone pair) + 4 × 1π (C) = 6 total"
        );
    }
}
