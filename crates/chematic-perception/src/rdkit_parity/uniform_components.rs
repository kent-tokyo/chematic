//! Exact complete-preperceived case without temporary component facts.
use super::{AtomIdx, BondOrder, Molecule};

/// The caller has rejected RDKit-model disagreement. Every cyclic component
/// here is either fully explicit aromatic or saturated with no incident
/// multiple bond. Thus the original complete-preperceived test returns
/// Identity; mixed components, cyclic multiples and exocyclic donor cases
/// keep that original test. Ring cache preparation remains unchanged.
pub(super) fn preperceived(mol: &Molecule) -> bool {
    let cyclic = crate::sssr::ring_bond_flags_shared(mol);
    let non_aromatic_ring_atom = |atom: AtomIdx| {
        !mol.atom(atom).aromatic
            && mol.neighbor_slice(atom).iter().any(|&(_, b)| cyclic[b.0 as usize])
    };
    for (idx, bond) in mol.bonds() {
        let (a, b) = (bond.atom1, bond.atom2);
        let aromatic_a = mol.atom(a).aromatic;
        let aromatic_b = mol.atom(b).aromatic;
        match bond.order {
            BondOrder::Aromatic => {
                if !cyclic[idx.0 as usize] || !aromatic_a || !aromatic_b {
                    return false;
                }
            }
            BondOrder::Single | BondOrder::Up | BondOrder::Down => {
                if cyclic[idx.0 as usize] && (aromatic_a || aromatic_b) {
                    return false;
                }
            }
            BondOrder::Double | BondOrder::Triple => {
                if cyclic[idx.0 as usize]
                    || non_aromatic_ring_atom(a)
                    || non_aromatic_ring_atom(b)
                {
                    return false;
                }
            }
            _ => return false,
        }
    }
    // Preserve ComponentInfo's explicit-representation consistency check.
    !mol.atoms().any(|(idx, atom)| {
        atom.aromatic
            && atom.charge == 0
            && matches!(atom.element.atomic_number(), 8 | 16 | 34 | 52)
            && mol.degree(idx) + chematic_core::implicit_hcount(mol, idx) as usize > 2
    })
}
