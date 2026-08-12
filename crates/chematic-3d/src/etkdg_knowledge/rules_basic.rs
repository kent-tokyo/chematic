//! Basic chemical knowledge --
//! [`TorsionKnowledgeSource::BasicChemicalKnowledge`].
//!
//! Structural rules not tied to a specific fitted SMARTS data table --
//! either textbook VSEPR facts (sp centers are linear) or a direct
//! translation of RDKit's own `useBasicKnowledge` block
//! (`TorsionPreferences.cpp` lines ~271-362, fetched and hashed in the
//! sources manifest): flat, all-sp2, 4-6-membered rings get a 2-fold
//! planarity-enforcing term.
//!
//! RDKit's own comment on the flat-ring force constant is internally
//! inconsistent (`// 7.0 is MMFF force constants for aromatic rings` next to
//! `fconsts[1] = 100.0;` -- the comment cites 7.0, the code uses 100.0).
//! This translation uses the value the code actually executes (100.0), not
//! the stale comment, and flags the discrepancy here rather than silently
//! picking one without saying so.

use chematic_core::{AtomIdx, BondOrder, Molecule};

use super::types::FourierTorsionTerm;

/// RDKit's `useBasicKnowledge` flat-ring rule applies to rings of exactly
/// these sizes (`rSize < 4 || rSize > 6` is skipped in the original --
/// i.e. 4, 5, or 6 qualify). Translated verbatim as a range.
pub const FLAT_RING_SIZES: std::ops::RangeInclusive<usize> = 4..=6;

/// `true` when `atom` is sp2 by this crate's existing convention (aromatic,
/// or has at least one incident double bond) -- matching how
/// `classify_atom_type` (legacy) and `classify_bond` (v2) both already
/// approximate hybridization from bond order rather than a dedicated
/// hybridization field (`chematic_core::Atom` has no explicit hybridization
/// slot).
pub fn is_sp2(mol: &Molecule, atom: AtomIdx) -> bool {
    mol.atom(atom).aromatic
        || mol.neighbors(atom).any(|(n, _)| {
            mol.bond_between(atom, n)
                .map(|(_, b)| b.order == BondOrder::Double)
                .unwrap_or(false)
        })
}

/// The 2-fold planarity term RDKit's `useBasicKnowledge` applies to every
/// bond of a 4-6-membered ring whose surrounding 4-atom tetrad is entirely
/// sp2 (translated from `TorsionPreferences.cpp`'s
/// `signs[1] = -1; fconsts[1] = 100.0;`, i.e. `n=2, s=-1, V=100.0`).
pub fn flat_ring_term() -> FourierTorsionTerm {
    FourierTorsionTerm::from_rdkit(2, -1, 100.0)
}

/// `true` when a ring of size `ring_size` containing the tetrad
/// `(a, b, c, d)` (with `b`-`c` the central bond) qualifies for
/// [`flat_ring_term`]: ring size in [`FLAT_RING_SIZES`] and all four atoms
/// sp2 (per [`is_sp2`]).
pub fn flat_ring_applies(mol: &Molecule, ring_size: usize, atoms: [AtomIdx; 4]) -> bool {
    FLAT_RING_SIZES.contains(&ring_size) && atoms.iter().all(|&a| is_sp2(mol, a))
}

/// `true` when `atom` sits at the center of a linear (sp-hybridized)
/// system: a triple bond, or a cumulated double-bond system (e.g. the
/// central carbon of an isocyanate `N=C=O` or a carbodiimide `N=C=N`).
/// This crate's own structural rule (not translated from RDKit -- see
/// module docs), reflecting the same VSEPR fact the legacy heuristic's
/// L306/L719/L745/L751 branches encode (see
/// `docs/rfcs/3d_torsion_knowledge_audit.md`), reimplemented independently here
/// so the v2 architecture does not depend on the legacy cascade at all.
pub fn is_linear_sp_center(mol: &Molecule, atom: AtomIdx) -> bool {
    let has_triple = mol.neighbors(atom).any(|(n, _)| {
        mol.bond_between(atom, n)
            .map(|(_, b)| b.order == BondOrder::Triple)
            .unwrap_or(false)
    });
    if has_triple {
        return true;
    }
    let double_bond_count = mol
        .neighbors(atom)
        .filter(|&(n, _)| {
            mol.bond_between(atom, n)
                .map(|(_, b)| b.order == BondOrder::Double)
                .unwrap_or(false)
        })
        .count();
    double_bond_count >= 2
}

/// This crate's own strong-linearity term for a bond adjacent to a linear
/// sp center: a stiff single-well 180 degree preference. Not RDKit-sourced
/// (see module docs) -- classified `BasicChemicalKnowledge`, never
/// `StandardExperimental`.
pub fn linear_sp_term() -> FourierTorsionTerm {
    FourierTorsionTerm::from_rdkit(1, -1, 50.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn benzene_atoms_are_all_sp2() {
        let mol = parse("c1ccccc1").unwrap();
        for i in 0..mol.atom_count() {
            assert!(is_sp2(&mol, AtomIdx(i as u32)));
        }
    }

    #[test]
    fn butane_carbons_are_not_sp2() {
        let mol = parse("CCCC").unwrap();
        for i in 0..mol.atom_count() {
            assert!(!is_sp2(&mol, AtomIdx(i as u32)));
        }
    }

    #[test]
    fn flat_ring_applies_for_benzene_tetrad() {
        let mol = parse("c1ccccc1").unwrap();
        let atoms = [AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3)];
        assert!(flat_ring_applies(&mol, 6, atoms));
    }

    #[test]
    fn flat_ring_does_not_apply_for_cyclohexane_tetrad() {
        let mol = parse("C1CCCCC1").unwrap();
        let atoms = [AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3)];
        assert!(!flat_ring_applies(&mol, 6, atoms));
    }

    #[test]
    fn flat_ring_does_not_apply_outside_4_6_range() {
        // A hypothetical all-sp2 7-membered ring tetrad should not qualify
        // (RDKit's own `rSize < 4 || rSize > 6` cutoff, translated verbatim).
        let mol = parse("c1ccccc1").unwrap();
        let atoms = [AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3)];
        assert!(!flat_ring_applies(&mol, 7, atoms));
        assert!(!flat_ring_applies(&mol, 3, atoms));
    }

    #[test]
    fn nitrile_carbon_is_linear_sp_center() {
        let mol = parse("CCC#N").unwrap();
        assert!(is_linear_sp_center(&mol, AtomIdx(2))); // triple-bonded C
    }

    #[test]
    fn alkane_carbon_is_not_linear_sp_center() {
        let mol = parse("CCCC").unwrap();
        assert!(!is_linear_sp_center(&mol, AtomIdx(1)));
    }
}
