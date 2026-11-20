//! Conservative remaining-atom bounds under the selected atom comparator.
use chematic_core::{BondEntry, BondIdx, BondOrder};
use rustc_hash::FxHashMap;

use super::{
    AtomCompare, AtomIdx, BondCompare, McsConfig, Molecule, normalize_bond, rdkit_heavy_atom_class,
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

type BondClass = (u8, u8, BondOrder);

/// Upper bond bound from compatible classes across inputs. The include/exclude
/// search removes source bonds, so only that input's class counts need updating.
pub(super) struct BondUpperBound {
    source_classes: Vec<Option<usize>>,
    excluded: Vec<bool>,
    source_counts: Vec<usize>,
    target_minima: Vec<usize>,
    total: usize,
}

/// Every shared bond uses the same unordered pair of atom classes and,
/// when required, the same normalized bond order in each input.
fn bond_class(
    mol: &Molecule,
    bond: &BondEntry,
    config: &McsConfig,
    rdkit_heavy: bool,
) -> Option<BondClass> {
    let class = |idx| atom_class(mol.atom(idx), &config.atom_compare, rdkit_heavy);
    let (a, b) = (class(bond.atom1)?, class(bond.atom2)?);
    let compare_order = config.match_bonds && config.bond_compare == BondCompare::OrderOrAromatic;
    let order = if compare_order {
        normalize_bond(bond.order)
    } else {
        BondOrder::QueryAny
    };
    Some((a.min(b), a.max(b), order))
}

impl BondUpperBound {
    pub(super) fn new(mols: &[&Molecule], config: &McsConfig, rdkit_heavy: bool) -> Self {
        let source = mols[0];
        let mut classes = FxHashMap::default();
        let mut source_counts = Vec::new();
        let source_classes: Vec<_> = source
            .bonds()
            .map(|(_, bond)| {
                let class = bond_class(source, bond, config, rdkit_heavy)?;
                let index = *classes.entry(class).or_insert_with(|| {
                    source_counts.push(0);
                    source_counts.len() - 1
                });
                source_counts[index] += 1;
                Some(index)
            })
            .collect();
        let mut target_minima = vec![usize::MAX; source_counts.len()];
        for &mol in &mols[1..] {
            let mut counts = vec![0; source_counts.len()];
            for (_, bond) in mol.bonds() {
                if let Some(class) = bond_class(mol, bond, config, rdkit_heavy)
                    && let Some(&index) = classes.get(&class)
                {
                    counts[index] += 1;
                }
            }
            for (minimum, count) in target_minima.iter_mut().zip(counts) {
                *minimum = (*minimum).min(count);
            }
        }
        let total = source_counts
            .iter()
            .zip(&target_minima)
            .map(|(&source, &target)| source.min(target))
            .sum();
        let excluded = source_classes
            .iter()
            .map(|class| class.is_none_or(|index| target_minima[index] == 0))
            .collect();
        Self {
            source_classes,
            excluded,
            source_counts,
            target_minima,
            total,
        }
    }

    pub(super) fn total(&self) -> usize {
        self.total
    }

    pub(super) fn is_allowed(&self, bond: BondIdx) -> bool {
        !self.excluded[bond.0 as usize]
    }

    pub(super) fn set_excluded(&mut self, bond: BondIdx, excluded: bool) {
        let idx = bond.0 as usize;
        debug_assert_ne!(self.excluded[idx], excluded);
        self.excluded[idx] = excluded;
        let index = self.source_classes[idx].expect("only compatible bonds can be toggled");
        update_bound(
            &mut self.total,
            &mut self.source_counts[index],
            self.target_minima[index],
            excluded,
        );
    }
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
        update_bound(
            &mut self.total,
            &mut self.source_counts[class as usize],
            self.target_minima[class as usize],
            excluded,
        );
    }
}

fn update_bound(total: &mut usize, count: &mut usize, minimum: usize, excluded: bool) {
    if excluded {
        *total -= usize::from(*count <= minimum);
        *count -= 1;
    } else {
        *total += usize::from(*count < minimum);
        *count += 1;
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
