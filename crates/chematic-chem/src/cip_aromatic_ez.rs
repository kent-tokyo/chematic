//! Use hierarchical CIP priority on marked E/Z bonds beside aromatic ligands.

use chematic_core::{AtomIdx, BondIdx, CipCode, Molecule};

/// The legacy sphere walk counts explicit double-bond duplicates, but not
/// MANCUDE duplicates of aromatic bonds. At an aromatic imine junction that
/// can reverse the two substituent priorities. Preserve the fast legacy path
/// for other bonds and fall back to it if the accurate ranker cannot resolve
/// both ends.
pub(super) fn try_assign(mol: &Molecule, bond_idx: BondIdx) -> Option<(AtomIdx, CipCode)> {
    let bond = mol.bond(bond_idx);
    let aromatic_substituent = [bond.atom1, bond.atom2].into_iter().any(|end| {
        mol.neighbors(end)
            .any(|(neighbor, _)| mol.atom(neighbor).aromatic)
    });
    if !aromatic_substituent {
        return None;
    }
    let ranker = super::LazyEzRanker::new(mol, chematic_cip::CipBudget::default_budget());
    super::assign_ez_accurate(mol, bond_idx, &ranker)
}
