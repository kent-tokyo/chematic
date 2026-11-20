//! Conservative remaining-atom bounds under the selected atom comparator.
use chematic_core::BondOrder;
use rustc_hash::FxHashMap;

use super::{
    AtomCompare, AtomIdx, BondCompare, McsConfig, McsObjective, Molecule, normalize_bond,
    rdkit_heavy_atom_class,
};
use crate::match_vf2::effective_atomic_number;

const CLASS_COUNT: usize = u8::MAX as usize + 1;

pub(super) fn atom_class(
    atom: &chematic_core::Atom,
    compare: &AtomCompare,
    rdkit_heavy: bool,
) -> Option<u8> {
    if rdkit_heavy {
        return Some(rdkit_heavy_atom_class(atom));
    }
    match compare {
        AtomCompare::Elements => Some(effective_atomic_number(atom)),
        AtomCompare::Any => Some(0),
        AtomCompare::AnyHeavyAtom => (effective_atomic_number(atom) > 1).then_some(0),
    }
}

pub(super) fn objective_limits(
    mols: &[&Molecule],
    config: &McsConfig,
    rdkit_heavy: bool,
    objective: McsObjective,
    atoms: &AtomUpperBound,
) -> (usize, usize) {
    match objective {
        McsObjective::AtomCount => (
            mols.iter().map(|mol| mol.atom_count()).min().unwrap_or(0),
            mols.iter().map(|mol| mol.bond_count()).min().unwrap_or(0),
        ),
        McsObjective::BondCount => {
            let mut common = bond_classes(mols[0], config, rdkit_heavy);
            for mol in &mols[1..] {
                let classes = bond_classes(mol, config, rdkit_heavy);
                common.retain(|class, count| {
                    *count = (*count).min(classes.get(class).copied().unwrap_or(0));
                    *count != 0
                });
            }
            let bonds: usize = common.values().sum();
            // A connected graph with at least one vertex has at least V - 1 edges.
            (atoms.total().min(bonds + 1), bonds)
        }
    }
}

/// Every shared bond uses the same unordered pair of atom classes and,
/// when required, the same normalized bond order in each input.
fn bond_classes(
    mol: &Molecule,
    config: &McsConfig,
    rdkit_heavy: bool,
) -> FxHashMap<(u8, u8, BondOrder), usize> {
    let mut counts = FxHashMap::default();
    let compare_order = config.match_bonds && config.bond_compare == BondCompare::OrderOrAromatic;
    for (_, bond) in mol.bonds() {
        let class = |idx| atom_class(mol.atom(idx), &config.atom_compare, rdkit_heavy);
        let (Some(a), Some(b)) = (class(bond.atom1), class(bond.atom2)) else {
            continue;
        };
        let order = if compare_order {
            normalize_bond(bond.order)
        } else {
            BondOrder::QueryAny
        };
        *counts.entry((a.min(b), a.max(b), order)).or_insert(0) += 1;
    }
    counts
}

/// Upper atom bound from per-class minima across all inputs. All atoms in a
/// mapped row have the same class, so subtracting mapped atoms would reduce
/// the bound by exactly the size of the mapping. Only excluding a source atom
/// changes the bound on the final total.
pub(super) struct AtomUpperBound {
    source_classes: Vec<Option<u8>>,
    excluded: Vec<bool>,
    source_counts: [usize; CLASS_COUNT],
    target_minima: [usize; CLASS_COUNT],
    total: usize,
}

impl AtomUpperBound {
    pub(super) fn new(mols: &[&Molecule], compare: &AtomCompare, rdkit_heavy: bool) -> Self {
        let source = mols[0];
        let source_classes: Vec<_> = source
            .atoms()
            .map(|(_, atom)| atom_class(atom, compare, rdkit_heavy))
            .collect();
        let mut source_counts = [0; CLASS_COUNT];
        for &class in source_classes.iter().flatten() {
            source_counts[class as usize] += 1;
        }
        let mut target_minima = [usize::MAX; CLASS_COUNT];
        for mol in &mols[1..] {
            let mut counts = [0; CLASS_COUNT];
            for (_, atom) in mol.atoms() {
                if let Some(class) = atom_class(atom, compare, rdkit_heavy) {
                    counts[class as usize] += 1;
                }
            }
            for (minimum, count) in target_minima.iter_mut().zip(counts) {
                *minimum = (*minimum).min(count);
            }
        }
        let total = source_counts
            .iter()
            .zip(target_minima.iter())
            .map(|(source, target)| (*source).min(*target))
            .sum();
        Self {
            source_classes,
            excluded: vec![false; source.atom_count()],
            source_counts,
            target_minima,
            total,
        }
    }

    pub(super) fn total(&self) -> usize {
        self.total
    }

    pub(super) fn is_excluded(&self, atom: AtomIdx) -> bool {
        self.excluded[atom.0 as usize]
    }

    pub(super) fn set_excluded(&mut self, atom: AtomIdx, excluded: bool) {
        let idx = atom.0 as usize;
        debug_assert_ne!(self.excluded[idx], excluded);
        self.excluded[idx] = excluded;
        let Some(class) = self.source_classes[idx] else {
            return;
        };
        let count = &mut self.source_counts[class as usize];
        let minimum = self.target_minima[class as usize];
        if excluded {
            self.total -= usize::from(*count <= minimum);
            *count -= 1;
        } else {
            self.total += usize::from(*count < minimum);
            *count += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exclusions_change_only_the_limiting_atom_class() {
        let source = chematic_smiles::parse("CC[*]O").unwrap();
        let left = chematic_smiles::parse("CCO").unwrap();
        let right = chematic_smiles::parse("[*]OC").unwrap();
        let mols = [&source, &left, &right];
        for (compare, rdkit_heavy, expected) in [
            (AtomCompare::Elements, false, [2, 2, 1, 1]),
            (AtomCompare::Any, false, [3, 3, 2, 1]),
            (AtomCompare::AnyHeavyAtom, false, [2, 2, 1, 1]),
            (AtomCompare::Any, true, [2, 2, 1, 1]),
        ] {
            let mut bound = AtomUpperBound::new(&mols, &compare, rdkit_heavy);
            assert_eq!(bound.total(), expected[0]);
            for (i, &size) in expected[1..].iter().enumerate() {
                let atom = AtomIdx(i as u32);
                bound.set_excluded(atom, true);
                assert!(bound.is_excluded(atom));
                assert_eq!(bound.total(), size);
            }
            for i in (0..3).rev() {
                let atom = AtomIdx(i as u32);
                bound.set_excluded(atom, false);
                assert!(!bound.is_excluded(atom));
                assert_eq!(bound.total(), expected[i]);
            }
        }
    }
}
