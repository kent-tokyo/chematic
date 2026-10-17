//! Conservative remaining-atom bounds under the selected atom comparator.
use std::collections::HashMap;

use super::{AtomCompare, AtomIdx, Molecule, PartialMapping};

fn atom_class(atom: &chematic_core::Atom, compare: &AtomCompare) -> Option<u8> {
    match compare {
        AtomCompare::Elements => Some(atom.element.atomic_number()),
        AtomCompare::Any => Some(0),
        AtomCompare::AnyHeavyAtom => (atom.element.atomic_number() > 1).then_some(0),
    }
}

/// Ignore extra charge/isotope/ring constraints: this must overestimate, never
/// underestimate. Element buckets are valid only for the Elements comparator.
pub(super) fn upper_bound_for_seeds(
    mols: &[&Molecule],
    a0: AtomIdx,
    seed: &[AtomIdx],
    compare: &AtomCompare,
) -> usize {
    let mut counts: HashMap<u8, Vec<usize>> = HashMap::new();
    for (mi, mol) in mols.iter().enumerate() {
        let exclude = if mi == 0 { a0 } else { seed[mi - 1] };
        for (ai, atom) in mol.atoms() {
            if ai == exclude {
                continue;
            }
            if let Some(class) = atom_class(atom, compare) {
                counts.entry(class).or_insert_with(|| vec![0; mols.len()])[mi] += 1;
            }
        }
    }
    1 + counts
        .values()
        .map(|v| *v.iter().min().unwrap_or(&0))
        .sum::<usize>()
}

/// Minimum compatible-class counts across unmapped, non-excluded atoms.
pub(super) fn upper_bound_additional(
    mols: &[&Molecule],
    mapping: &PartialMapping,
    excluded: &[bool],
    compare: &AtomCompare,
) -> usize {
    let mut counts: HashMap<u8, Vec<usize>> = HashMap::new();
    for (mi, mol) in mols.iter().enumerate() {
        for (ai, atom) in mol.atoms() {
            if mapping.is_mapped(mi, ai) || (mi == 0 && excluded[ai.0 as usize]) {
                continue;
            }
            if let Some(class) = atom_class(atom, compare) {
                counts.entry(class).or_insert_with(|| vec![0; mols.len()])[mi] += 1;
            }
        }
    }
    counts.values().map(|v| *v.iter().min().unwrap_or(&0)).sum()
}
