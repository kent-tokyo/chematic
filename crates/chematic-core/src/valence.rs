//! Valence model: implicit hydrogen count for organic-subset atoms.
//!
//! Reference: OpenSMILES specification, section 3.4 (Implicit hydrogen)
//! <http://opensmiles.org/opensmiles-spec.html>

use crate::bond::BondOrder;
use crate::molecule::{AtomIdx, Molecule};

/// Compute the implicit hydrogen count for atom `idx`.
///
/// - Bracket atoms (`hydrogen_count.is_some()`): return the stored value directly.
/// - Wildcard atoms: return 0.
/// - Organic-subset atoms: derive from the normal-valence table.
/// - All other atoms: return 0 (no implicit H rule defined).
///
/// # Algorithm
/// 1. Sum the integer bond orders of all bonds on the atom.
/// 2. Find the smallest normal valence >= bond_sum (adjusted for formal charge).
/// 3. implicit_H = adjusted_valence - bond_sum.
///
/// Charge adjustment:
/// - Positive charge: increases target valence (e.g. [NH4]+ has 4 bonds → valence 4).
/// - Negative charge: decreases target valence.
pub fn implicit_hcount(mol: &Molecule, idx: AtomIdx) -> u8 {
    let atom = mol.atom(idx);

    // Wildcards have no defined implicit H.
    if atom.wildcard {
        return 0;
    }

    // Bracket atoms store the explicit H count.
    if let Some(h) = atom.hydrogen_count {
        return h;
    }

    // Only the organic subset gets implicit H.
    if !atom.element.is_organic_subset() {
        return 0;
    }

    let normal_valences = atom.element.normal_valences();
    if normal_valences.is_empty() {
        return 0;
    }

    let charge = atom.charge as i32;

    // Separate aromatic bonds from non-aromatic bonds.
    let mut aromatic_count: usize = 0;
    let mut non_aromatic_sum: i32 = 0;
    for (_, bidx) in mol.neighbors(idx) {
        let order = mol.bond(bidx).order;
        if order == BondOrder::Aromatic {
            aromatic_count += 1;
        } else {
            non_aromatic_sum += order.order_int() as i32;
        }
    }

    if aromatic_count > 0 {
        // Aromatic molecule (pre-Kekulization): each aromatic bond contributes 1.5
        // to the effective bond order (OpenSMILES convention).
        //
        // floor(1.5 × n) gives the contribution from n aromatic bonds:
        //   n=2 → 3  benzene CH:   4−3=1H ✓   pyridine N: 3−3=0H ✓
        //   n=3 → 4  junction C:   4−4=0H ✓
        //
        // Combined with non-aromatic substituents (e.g. N−CH₃) this correctly yields
        // 0 H for all substituted aromatic atoms without needing Kekulization.
        // Always use the lowest normal valence; aromatic atoms cannot be hypervalent.
        let effective_sum =
            (aromatic_count as f64 * 1.5).floor() as i32 + non_aromatic_sum;
        let v = normal_valences[0] as i32 + charge;
        if v <= 0 || effective_sum >= v {
            return 0;
        }
        return (v - effective_sum) as u8;
    }

    // Non-aromatic path (or post-Kekulization molecule where all bonds are explicit).
    let bond_sum = non_aromatic_sum;

    // For atoms that carry the aromatic flag but reside in a kekulized molecule
    // (bonds are Single/Double, not Aromatic), use only the lowest normal valence.
    // Rationale: after Kekulization, a substituted aromatic N (e.g. N−CH₃ in caffeine
    // with one ring double bond) has bond_sum=4, which would select valence 5 and
    // give 1 implicit H.  Capping at the primary valence (3) returns 0 H instead.
    let valences_to_check: &[u8] = if atom.aromatic {
        &normal_valences[..1]
    } else {
        normal_valences
    };

    // Iterate through valences (ascending) and pick the smallest ≥ bond_sum.
    for &v in valences_to_check {
        let target = v as i32 + charge;
        if target < 0 {
            continue;
        }
        if target >= bond_sum {
            return (target - bond_sum) as u8;
        }
    }

    // bond_sum exceeds all consulted valences → 0 implicit H.
    0
}

/// Compute the *total* hydrogen count (explicit implicit H + explicit bracket H).
/// For bracket atoms the explicit H is already returned by `implicit_hcount`.
pub fn total_hcount(mol: &Molecule, idx: AtomIdx) -> u8 {
    implicit_hcount(mol, idx)
}

/// Sum of integer bond orders for heavy-atom bonds on `idx`.
/// Aromatic bonds count as 1 (pre-Kekulization representation).
pub fn bond_order_sum(mol: &Molecule, idx: AtomIdx) -> u8 {
    mol.neighbors(idx)
        .map(|(_, bidx)| mol.bond(bidx).order.order_int())
        .fold(0u8, |acc, x| acc.saturating_add(x))
}

/// Returns true if the bond is counted as a "double bond equivalent" in valence sums.
pub fn is_pi_bond(order: BondOrder) -> bool {
    matches!(order, BondOrder::Double | BondOrder::Triple | BondOrder::Quadruple)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::bond::BondOrder;
    use crate::element::Element;
    use crate::molecule::MoleculeBuilder;

    fn single_atom(elem: Element) -> Molecule {
        let mut b = MoleculeBuilder::new();
        b.add_atom(Atom::organic(elem));
        b.build()
    }

    fn two_atoms(e1: Element, e2: Element, order: BondOrder) -> Molecule {
        let mut b = MoleculeBuilder::new();
        let a = b.add_atom(Atom::organic(e1));
        let c = b.add_atom(Atom::organic(e2));
        b.add_bond(a, c, order).unwrap();
        b.build()
    }

    #[test]
    fn test_methane() {
        // C alone: 0 bonds, valence 4 → 4 implicit H
        let mol = single_atom(Element::C);
        assert_eq!(implicit_hcount(&mol, AtomIdx(0)), 4);
    }

    #[test]
    fn test_ethane_c() {
        // CC: each C has 1 single bond → valence 4 → 3 implicit H
        let mol = two_atoms(Element::C, Element::C, BondOrder::Single);
        assert_eq!(implicit_hcount(&mol, AtomIdx(0)), 3);
        assert_eq!(implicit_hcount(&mol, AtomIdx(1)), 3);
    }

    #[test]
    fn test_ethylene_c() {
        // C=C: double bond → bond_sum=2 → 4-2=2 implicit H
        let mol = two_atoms(Element::C, Element::C, BondOrder::Double);
        assert_eq!(implicit_hcount(&mol, AtomIdx(0)), 2);
    }

    #[test]
    fn test_acetylene_c() {
        // C#C: triple bond → bond_sum=3 → 4-3=1 implicit H
        let mol = two_atoms(Element::C, Element::C, BondOrder::Triple);
        assert_eq!(implicit_hcount(&mol, AtomIdx(0)), 1);
    }

    #[test]
    fn test_nitrogen_amine() {
        // N alone: 0 bonds, first normal valence=3 → 3 implicit H (NH3)
        let mol = single_atom(Element::N);
        assert_eq!(implicit_hcount(&mol, AtomIdx(0)), 3);
    }

    #[test]
    fn test_nitrogen_triple() {
        // N#C: N has triple bond → bond_sum=3 → 3-3=0 (nitrile N)
        let mol = two_atoms(Element::N, Element::C, BondOrder::Triple);
        assert_eq!(implicit_hcount(&mol, AtomIdx(0)), 0);
    }

    #[test]
    fn test_oxygen_ether() {
        // O alone: 0 bonds, valence 2 → 2 implicit H (water)
        let mol = single_atom(Element::O);
        assert_eq!(implicit_hcount(&mol, AtomIdx(0)), 2);
    }

    #[test]
    fn test_fluorine() {
        // F alone: valence 1 → 1 implicit H (HF)
        let mol = single_atom(Element::F);
        assert_eq!(implicit_hcount(&mol, AtomIdx(0)), 1);
    }

    #[test]
    fn test_bracket_atom_explicit_h() {
        // [NH4+] — bracket atom: explicit H=4 returned directly
        let mut b = MoleculeBuilder::new();
        let atom = Atom::bracket(Element::N, None, Default::default(), 4, 1, None);
        b.add_atom(atom);
        let mol = b.build();
        assert_eq!(implicit_hcount(&mol, AtomIdx(0)), 4);
    }

    #[test]
    fn test_hypervalent_sulfur() {
        // S with four single bonds: bond_sum=4, S valences=[2,4,6] → target=4 → 0 H
        let mut b = MoleculeBuilder::new();
        let s = b.add_atom(Atom::organic(Element::S));
        for _ in 0..4 {
            let c = b.add_atom(Atom::organic(Element::C));
            b.add_bond(s, c, BondOrder::Single).unwrap();
        }
        let mol = b.build();
        assert_eq!(implicit_hcount(&mol, AtomIdx(0)), 0);
    }
}
