//! Exact unspecified E/Z eligibility with bounded substituent storage.

use super::{AtomIdx, BondIdx, BondOrder, HashMap, Molecule, in_ring_smaller_than};

/// Store only the two substituents that a qualifying end can have, while
/// counting all neighbors exactly as the original unbounded Vec did.
fn end_subs(
    mol: &Molecule,
    end: AtomIdx,
    other: AtomIdx,
    db: BondIdx,
) -> Option<([AtomIdx; 2], usize)> {
    let mut subs = [AtomIdx(0); 2];
    let mut len = 0;
    for (nb, bi) in mol.neighbors(end) {
        if nb != other && mol.bond(bi).order != BondOrder::Dative {
            if len < subs.len() {
                subs[len] = nb;
            }
            len += 1;
        }
    }
    let cumulated = mol
        .neighbors(end)
        .any(|(_, bi)| bi != db && mol.bond(bi).order == BondOrder::Double);
    let h = chematic_core::implicit_hcount(mol, end) as usize;
    (!cumulated && len != 0 && len + h <= 2).then_some((subs, len))
}

/// Double bonds outside small rings whose ends each have one or two
/// distinct substituents (an implicit H or a lone pair making up the
/// rest), with no E/Z in `declared`.
pub(super) fn unspecified_double_bonds(
    mol: &Molecule,
    declared: &HashMap<BondIdx, chematic_core::atom::CipCode>,
) -> Vec<BondIdx> {
    let classes = std::cell::OnceCell::new();
    let distinct = |subs: &[AtomIdx]| {
        subs.len() < 2 || {
            let classes =
                classes.get_or_init(|| chematic_smiles::topological_equivalence_classes(mol));
            classes[subs[0].0 as usize] != classes[subs[1].0 as usize]
        }
    };
    let mut out: Vec<BondIdx> = mol
        .bonds()
        .filter(|&(bi, b)| {
            if b.order != BondOrder::Double
                || declared.contains_key(&bi)
                || (mol.atom(b.atom1).aromatic && mol.atom(b.atom2).aromatic)
            {
                return false;
            }
            let (Some((s1, n1)), Some((s2, n2))) = (
                end_subs(mol, b.atom1, b.atom2, bi),
                end_subs(mol, b.atom2, b.atom1, bi),
            ) else {
                return false;
            };
            !in_ring_smaller_than(mol, bi, 8) && distinct(&s1[..n1]) && distinct(&s2[..n2])
        })
        .map(|(bi, _)| bi)
        .collect();
    out.sort_by_key(|b| b.0);
    out
}

#[cfg(test)]
#[path = "unspecified_reference_tests.rs"]
mod reference_tests;
