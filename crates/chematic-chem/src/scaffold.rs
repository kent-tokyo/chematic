//! Murcko scaffold decomposition.
//!
//! Provides functions to extract the Murcko scaffold from a molecule:
//! - `murcko_scaffold`: ring atoms plus atoms on paths connecting ring systems.
//! - `generic_murcko_scaffold`: scaffold with all atoms replaced by C and all bonds by Single.

#![forbid(unsafe_code)]

use std::collections::{HashMap, HashSet};

use chematic_core::{Atom, AtomIdx, BondIdx, BondOrder, Element, Molecule, MoleculeBuilder};
use chematic_perception::find_sssr;

/// Extract the Murcko scaffold from `mol`.
///
/// The scaffold consists of:
/// - All ring atoms (atoms participating in at least one ring).
/// - Linker atoms: a non-ring atom is included iff it has >=2 heavy-atom
///   neighbors that are already in the scaffold. The expansion repeats until
///   no new linkers are added.
///
/// Returns an empty `Molecule` if `mol` contains no rings.
pub fn murcko_scaffold(mol: &Molecule) -> Molecule {
    let rings = find_sssr(mol);
    if rings.ring_count() == 0 {
        return MoleculeBuilder::new().build();
    }

    let mut scaffold_atoms: HashSet<AtomIdx> = rings
        .rings()
        .iter()
        .flat_map(|r| r.iter().copied())
        .collect();

    // Iteratively pull in linker atoms until stable.
    loop {
        let mut changed = false;
        for i in 0..mol.atom_count() {
            let idx = AtomIdx(i as u32);
            if scaffold_atoms.contains(&idx) {
                continue;
            }
            let scaffold_neighbors = mol
                .neighbors(idx)
                .filter(|(nb, _)| scaffold_atoms.contains(nb))
                .count();
            if scaffold_neighbors >= 2 {
                scaffold_atoms.insert(idx);
                changed = true;
            }
        }
        if !changed { break; }
    }

    build_subgraph(mol, &scaffold_atoms)
}

/// Generic Murcko scaffold: every atom becomes C and every bond becomes Single.
///
/// Returns an empty `Molecule` if `mol` has no rings.
pub fn generic_murcko_scaffold(mol: &Molecule) -> Molecule {
    let scaffold = murcko_scaffold(mol);
    if scaffold.atom_count() == 0 {
        return scaffold;
    }

    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();
    for i in 0..scaffold.atom_count() {
        let new_idx = builder.add_atom(Atom::organic(Element::C));
        remap.insert(AtomIdx(i as u32), new_idx);
    }
    for i in 0..scaffold.bond_count() {
        let bond = scaffold.bond(BondIdx(i as u32));
        if let (Some(&new_a), Some(&new_b)) = (remap.get(&bond.atom1), remap.get(&bond.atom2)) {
            let _ = builder.add_bond(new_a, new_b, BondOrder::Single);
        }
    }
    builder.build()
}

/// Build a subgraph of `mol` containing only the atoms in `atom_set`,
/// preserving bonds whose both endpoints are in the set.
fn build_subgraph(mol: &Molecule, atom_set: &HashSet<AtomIdx>) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();
    for i in 0..mol.atom_count() {
        let old_idx = AtomIdx(i as u32);
        if atom_set.contains(&old_idx) {
            let new_idx = builder.add_atom(mol.atom(old_idx).clone());
            remap.insert(old_idx, new_idx);
        }
    }
    for i in 0..mol.bond_count() {
        let bond = mol.bond(BondIdx(i as u32));
        if let (Some(&new_a), Some(&new_b)) = (remap.get(&bond.atom1), remap.get(&bond.atom2)) {
            let _ = builder.add_bond(new_a, new_b, bond.order);
        }
    }
    builder.build()
}

// ---------------------------------------------------------------------------
// Schuffenhauer scaffold network
// ---------------------------------------------------------------------------

/// Return the full scaffold network of `mol` as a `Vec<Molecule>`.
///
/// Starting from the Murcko scaffold, rings are iteratively removed one at a
/// time following Schuffenhauer's priority rules (Schuffenhauer et al. 2007)
/// until a single-ring scaffold remains.  Each intermediate is included.
///
/// The returned vector starts with the Murcko scaffold and ends with the
/// smallest core ring.  Returns an empty `Vec` if the molecule has no rings.
pub fn scaffold_network(mol: &Molecule) -> Vec<Molecule> {
    let start = murcko_scaffold(mol);
    if start.atom_count() == 0 {
        return Vec::new();
    }

    let mut network: Vec<Molecule> = Vec::new();
    let mut current = start;

    loop {
        // Store a copy of current before mutating.
        let snapshot: HashSet<AtomIdx> = (0..current.atom_count())
            .map(|i| AtomIdx(i as u32))
            .collect();
        network.push(build_subgraph(&current, &snapshot));

        let rings = find_sssr(&current);
        if rings.ring_count() <= 1 {
            break; // single ring or no ring — stop
        }

        match schuffenhauer_remove_ring(&current, &rings) {
            Some(next) => {
                current = next;
            }
            None => break,
        }
    }

    network
}

/// Return the direct Schuffenhauer parent scaffolds of `mol` (one ring removed).
///
/// In most cases this returns a single scaffold.  Returns an empty `Vec` if
/// the Murcko scaffold has 0 or 1 rings.
pub fn schuffenhauer_parents(mol: &Molecule) -> Vec<Molecule> {
    let start = murcko_scaffold(mol);
    let rings = find_sssr(&start);
    if rings.ring_count() <= 1 {
        return Vec::new();
    }
    schuffenhauer_remove_ring(&start, &rings)
        .into_iter()
        .collect()
}

/// Remove one ring from `mol` according to Schuffenhauer priority rules.
///
/// Rules (applied in order until a winner is found):
/// 1. Remove an outermost ring (one that, when removed, reduces ring count by 1).
/// 2. Among candidates, prefer all-carbon rings over heteroaromatic rings.
/// 3. Among candidates with same heteroatom content, prefer the smallest ring.
/// 4. Break remaining ties by preferring the ring with the lowest-priority
///    attachment point (smallest atom index in the ring).
fn schuffenhauer_remove_ring(mol: &Molecule, rings: &chematic_perception::RingSet) -> Option<Molecule> {
    let n_rings = rings.ring_count();
    if n_rings == 0 {
        return None;
    }

    // Find "removable" rings: removing all atoms in the ring (that don't belong
    // to other rings) still leaves a connected scaffold or is the last ring.
    let all_rings: Vec<Vec<AtomIdx>> = rings.rings().to_vec();

    // For each ring, count how many of its atoms appear in exactly 1 ring.
    // Those atoms can be safely deleted when we remove the ring.
    let mut candidates: Vec<usize> = (0..n_rings)
        .filter(|&ri| {
            // A ring is outermost if it shares at most one ring member with
            // any other ring (i.e. it has atoms exclusive to it).
            let exclusive: usize = all_rings[ri].iter()
                .filter(|&&atom| {
                    all_rings.iter().enumerate()
                        .filter(|(j, _)| *j != ri)
                        .all(|(_, other)| !other.contains(&atom))
                })
                .count();
            exclusive > 0
        })
        .collect();

    if candidates.is_empty() {
        // Fallback: all rings are fused. Pick the one with the most shared atoms.
        candidates = (0..n_rings).collect();
    }

    // Rule 2: prefer all-carbon rings (no heteroatoms).
    let carbon_only: Vec<usize> = candidates.iter().copied()
        .filter(|&ri| all_rings[ri].iter().all(|&a| mol.atom(a).element.atomic_number() == 6))
        .collect();
    if !carbon_only.is_empty() {
        candidates = carbon_only;
    }

    // Rule 3: prefer smallest ring.
    let min_size = candidates.iter().map(|&ri| all_rings[ri].len()).min().unwrap();
    candidates.retain(|&ri| all_rings[ri].len() == min_size);

    // Rule 4: tie-break by smallest atom index in the ring.
    candidates.sort_by_key(|&ri| all_rings[ri].iter().map(|a| a.0).min().unwrap_or(0));
    let chosen_ring = candidates[0];

    // Build the set of atoms to DELETE: atoms exclusive to the chosen ring.
    let to_delete: HashSet<AtomIdx> = all_rings[chosen_ring].iter().copied()
        .filter(|&atom| {
            all_rings.iter().enumerate()
                .filter(|(j, _)| *j != chosen_ring)
                .all(|(_, other)| !other.contains(&atom))
        })
        .collect();

    if to_delete.is_empty() {
        return None; // fused ring with no exclusive atoms — can't remove
    }

    // Build new scaffold: keep all atoms NOT in to_delete.
    let keep: HashSet<AtomIdx> = (0..mol.atom_count())
        .map(|i| AtomIdx(i as u32))
        .filter(|a| !to_delete.contains(a))
        .collect();

    if keep.is_empty() {
        return None;
    }

    Some(build_subgraph(mol, &keep))
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn murcko_benzene_preserves_all_atoms() {
        // Benzene: all 6 atoms are ring atoms; scaffold == entire molecule.
        let mol = parse("c1ccccc1").unwrap();
        let scaffold = murcko_scaffold(&mol);
        assert_eq!(scaffold.atom_count(), 6);
    }

    #[test]
    fn murcko_toluene_removes_methyl() {
        // Toluene: benzene ring (6 atoms) + one methyl group (1 C); methyl is a side chain.
        let mol = parse("Cc1ccccc1").unwrap();
        let scaffold = murcko_scaffold(&mol);
        assert_eq!(scaffold.atom_count(), 6, "methyl group should be removed");
    }

    #[test]
    fn murcko_ethylbenzene_removes_chain() {
        // Ethylbenzene: benzene ring + two-atom chain; chain is a side chain.
        let mol = parse("CCc1ccccc1").unwrap();
        let scaffold = murcko_scaffold(&mol);
        assert_eq!(scaffold.atom_count(), 6, "ethyl chain should be removed");
    }

    #[test]
    fn murcko_acyclic_returns_empty() {
        // Ethane has no rings; scaffold should be empty.
        let mol = parse("CC").unwrap();
        let scaffold = murcko_scaffold(&mol);
        assert_eq!(scaffold.atom_count(), 0);
    }

    #[test]
    fn generic_murcko_benzene_all_carbon_single() {
        // Generic scaffold of benzene: 6 C atoms, all bonds Single.
        let mol = parse("c1ccccc1").unwrap();
        let generic = generic_murcko_scaffold(&mol);
        assert_eq!(generic.atom_count(), 6);
        for i in 0..generic.atom_count() {
            let atom = generic.atom(AtomIdx(i as u32));
            assert_eq!(atom.element, Element::C, "all atoms should be carbon");
        }
        for i in 0..generic.bond_count() {
            let bond = generic.bond(chematic_core::BondIdx(i as u32));
            assert_eq!(
                bond.order,
                BondOrder::Single,
                "all bonds should be Single"
            );
        }
    }

    #[test]
    fn murcko_biphenyl_keeps_all_ring_atoms() {
        // Biphenyl: two fused/connected phenyl rings, all atoms are ring atoms.
        let mol = parse("c1ccccc1c1ccccc1").unwrap();
        let scaffold = murcko_scaffold(&mol);
        assert!(
            scaffold.atom_count() >= 12,
            "biphenyl scaffold should have at least 12 atoms, got {}",
            scaffold.atom_count()
        );
    }
}
