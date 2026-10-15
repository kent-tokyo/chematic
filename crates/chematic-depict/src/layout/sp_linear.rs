use chematic_core::{AtomIdx, BondOrder, Molecule};

/// Continue through a degree-2 sp centre without the ordinary zigzag turn.
/// A fork's compensating `back` turn must also be skipped at that centre.
pub(super) fn continuation_direction(
    mol: &Molecule,
    center: AtomIdx,
    incoming: f64,
    back: f64,
    sign: f64,
    deflection: f64,
) -> f64 {
    if mol.degree(center) == 2 && !mol.atom(center).aromatic {
        let bonds = mol.neighbor_slice(center);
        let first = mol.bond(bonds[0].1).order;
        let second = mol.bond(bonds[1].1).order;
        if matches!(
            (first, second),
            (BondOrder::Triple, _)
                | (_, BondOrder::Triple)
                | (BondOrder::Double, BondOrder::Double)
        ) {
            return incoming;
        }
    }
    incoming + back + sign * deflection
}
