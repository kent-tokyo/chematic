//! Murcko scaffold decomposition.
//!
//! Provides functions to extract the Murcko scaffold from a molecule:
//! - `murcko_scaffold`: ring atoms plus atoms on paths connecting ring systems.
//! - `generic_murcko_scaffold`: scaffold with all atoms replaced by C and all bonds by Single.
//! - `scaffold_network_with_counts`: compute scaffold network with molecule frequencies across a library.

#![forbid(unsafe_code)]

use std::collections::{HashMap, HashSet};

use chematic_core::{Atom, AtomIdx, BondIdx, BondOrder, Element, Molecule, MoleculeBuilder};
use chematic_perception::find_sssr;
use chematic_smiles::canonical_smiles;

/// Scaffold network of a molecule library: hierarchical scaffolds with occurrence counts.
///
/// Represents a directed acyclic graph (DAG) where each node is a scaffold and edges
/// represent parent-child relationships (parent = one ring removed). Each scaffold
/// tracks how many molecules in the input library had that scaffold in their hierarchy.
pub struct ScaffoldNetwork {
    /// Unique scaffolds in canonical order (by SMILES).
    pub scaffolds: Vec<Molecule>,
    /// Occurrence count for each scaffold across the input library.
    pub counts: Vec<usize>,
    /// Parent scaffold index for each scaffold (if any). `None` indicates root (single ring).
    pub parents: Vec<Option<usize>>,
}

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
        if !changed {
            break;
        }
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
/// 4. Prefer rings with fewer nitrogen atoms (nitrogen-containing rings deprioritized).
/// 5. Prefer rings with fewer heavy-atom substituents (attachments to other rings).
/// 6. Prefer 5-membered rings over 6-membered and larger (5-ring priority).
/// 7. Prefer rings with fewer linker bonds (fewest connections to other rings).
/// 8. Break remaining ties by smallest atom index in the ring.
fn schuffenhauer_remove_ring(
    mol: &Molecule,
    rings: &chematic_perception::RingSet,
) -> Option<Molecule> {
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
            let exclusive: usize = all_rings[ri]
                .iter()
                .filter(|&&atom| {
                    all_rings
                        .iter()
                        .enumerate()
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
    let carbon_only: Vec<usize> = candidates
        .iter()
        .copied()
        .filter(|&ri| {
            all_rings[ri]
                .iter()
                .all(|&a| mol.atom(a).element.atomic_number() == 6)
        })
        .collect();
    if !carbon_only.is_empty() {
        candidates = carbon_only;
    }

    // Rule 3: prefer smallest ring.
    let min_size = candidates
        .iter()
        .map(|&ri| all_rings[ri].len())
        .min()
        .unwrap();
    candidates.retain(|&ri| all_rings[ri].len() == min_size);

    // Rule 4: prefer rings with fewer nitrogen atoms.
    let count_heteroatoms = |ri: usize| {
        all_rings[ri]
            .iter()
            .filter(|&&a| mol.atom(a).element.atomic_number() == 7)
            .count()
    };
    let min_nitrogen = candidates
        .iter()
        .map(|&ri| count_heteroatoms(ri))
        .min()
        .unwrap();
    candidates.retain(|&ri| count_heteroatoms(ri) == min_nitrogen);

    // Rule 5: prefer rings with fewer heavy-atom substituents (linker attachments).
    let count_substituents = |ri: usize| {
        all_rings[ri]
            .iter()
            .filter(|&&atom| {
                mol.neighbors(atom).any(|(nb, _)| {
                    // A substituent is a neighbor not in this ring and not hydrogen
                    !all_rings[ri].contains(&nb) && mol.atom(nb).element.atomic_number() != 1
                })
            })
            .count()
    };
    let min_subs = candidates
        .iter()
        .map(|&ri| count_substituents(ri))
        .min()
        .unwrap();
    candidates.retain(|&ri| count_substituents(ri) == min_subs);

    // Rule 6: prefer 5-membered rings over larger rings (5-ring priority).
    let prefer_5ring: Vec<usize> = candidates
        .iter()
        .copied()
        .filter(|&ri| all_rings[ri].len() == 5)
        .collect();
    if !prefer_5ring.is_empty() {
        candidates = prefer_5ring;
    }

    // Rule 7: prefer rings with fewer linker bonds (fewest inter-ring connections).
    let count_linker_bonds = |ri: usize| {
        all_rings[ri]
            .iter()
            .filter(|&&atom| {
                all_rings
                    .iter()
                    .enumerate()
                    .any(|(j, other)| j != ri && other.contains(&atom))
            })
            .count()
    };
    let min_linker = candidates
        .iter()
        .map(|&ri| count_linker_bonds(ri))
        .min()
        .unwrap();
    candidates.retain(|&ri| count_linker_bonds(ri) == min_linker);

    // Rule 8: tie-break by smallest atom index in the ring.
    //
    // This sorts by raw AtomIdx, which depends on SMILES parse order, not a
    // canonical rank -- same non-canonical-tie-break shape as the
    // is_fused_ring bug fixed in chematic-perception's ring_family.rs.
    // Flagged, not fixed: a direct discriminator (comparing scaffold()
    // atom counts across disagreeing worst-of-10 traversals, full 5000-mol
    // corpus) found 0/40 topology-selection mismatches -- every measured
    // scaffold instability was a canonical-SMILES writer issue, not this
    // tie-break choosing a different ring. Worth auditing on its own before
    // changing it; don't assume this is broken just because it's
    // non-canonical.
    candidates.sort_by_key(|&ri| all_rings[ri].iter().map(|a| a.0).min().unwrap_or(0));
    let chosen_ring = candidates[0];

    // Build the set of atoms to DELETE: atoms exclusive to the chosen ring.
    let to_delete: HashSet<AtomIdx> = all_rings[chosen_ring]
        .iter()
        .copied()
        .filter(|&atom| {
            all_rings
                .iter()
                .enumerate()
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

/// Compute the scaffold network from a molecule library with occurrence counts.
///
/// For each molecule in `mols`, extracts its scaffold hierarchy via `scaffold_network()`,
/// then aggregates all unique scaffolds by canonical SMILES, counting occurrences.
///
/// Returns a `ScaffoldNetwork` with:
/// - `scaffolds`: unique scaffolds (deduplicated by canonical SMILES)
/// - `counts`: occurrence frequency of each scaffold across the library
/// - `parents`: parent scaffold index (one ring removed in the Schuffenhauer hierarchy)
///
/// Returns an empty network if no scaffolds are found.
pub fn scaffold_network_with_counts(mols: &[Molecule]) -> ScaffoldNetwork {
    // Collect scaffolds and metadata: SMILES → (count, parent_smiles_opt)
    let mut smi_counts: HashMap<String, usize> = HashMap::new();
    let mut smi_parents: HashMap<String, Option<String>> = HashMap::new();
    let mut scaffolds_list: Vec<(String, Molecule)> = Vec::new(); // Ordered by first encounter
    let mut seen_smiles: HashSet<String> = HashSet::new();

    for mol in mols {
        let network = scaffold_network(mol);
        for (i, scaffold) in network.iter().enumerate() {
            let smiles = canonical_smiles(scaffold);

            // Track count
            smi_counts
                .entry(smiles.clone())
                .and_modify(|c| *c += 1)
                .or_insert(1);

            // Track parent (by computing SMILES of network[i-1])
            if i > 0 {
                let parent_smiles = canonical_smiles(&network[i - 1]);
                smi_parents
                    .entry(smiles.clone())
                    .or_insert(Some(parent_smiles));
            } else {
                smi_parents.entry(smiles.clone()).or_insert(None);
            }

            // Store first occurrence by reconstructing the scaffold
            if !seen_smiles.contains(&smiles) {
                seen_smiles.insert(smiles.clone());
                // Reconstruct the scaffold since we're borrowing from network
                let atom_set: HashSet<AtomIdx> = (0..scaffold.atom_count())
                    .map(|i| AtomIdx(i as u32))
                    .collect();
                let rebuilt = build_subgraph(scaffold, &atom_set);
                scaffolds_list.push((smiles, rebuilt));
            }
        }
    }

    // Sort by SMILES for stable output
    scaffolds_list.sort_by(|a, b| a.0.cmp(&b.0));

    // Build SMILES → index mapping before consuming scaffolds_list
    let smiles_vec: Vec<String> = scaffolds_list.iter().map(|(smi, _)| smi.clone()).collect();
    let smi_to_idx: HashMap<String, usize> = smiles_vec
        .iter()
        .enumerate()
        .map(|(idx, smi)| (smi.clone(), idx))
        .collect();

    // Extract molecules and build count/parent arrays
    let (_, scaffolds_vec): (Vec<_>, Vec<_>) = scaffolds_list.into_iter().unzip();
    let counts: Vec<usize> = smiles_vec
        .iter()
        .map(|smi| smi_counts.get(smi).copied().unwrap_or(0))
        .collect();
    let parents: Vec<Option<usize>> = smiles_vec
        .iter()
        .map(|smi| {
            smi_parents
                .get(smi)
                .cloned()
                .flatten()
                .and_then(|parent_smi| smi_to_idx.get(&parent_smi).copied())
        })
        .collect();

    ScaffoldNetwork {
        scaffolds: scaffolds_vec,
        counts,
        parents,
    }
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
            assert_eq!(bond.order, BondOrder::Single, "all bonds should be Single");
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

    // scaffold_network_with_counts tests

    #[test]
    fn scaffold_network_with_counts_empty_input() {
        let network = scaffold_network_with_counts(&[]);
        assert!(
            network.scaffolds.is_empty(),
            "empty input should yield empty network"
        );
    }

    #[test]
    fn scaffold_network_with_counts_single_molecule() {
        // Single molecule with multiple scaffold layers (e.g. naphthalene)
        let mol = parse("c1ccc2ccccc2c1").unwrap();
        let network = scaffold_network_with_counts(&[mol]);
        assert!(
            !network.scaffolds.is_empty(),
            "naphthalene should yield at least one scaffold"
        );
        assert_eq!(
            network.scaffolds.len(),
            network.counts.len(),
            "scaffolds and counts must have same length"
        );
        assert_eq!(
            network.scaffolds.len(),
            network.parents.len(),
            "scaffolds and parents must have same length"
        );
    }

    #[test]
    fn scaffold_network_with_counts_duplicate_scaffolds() {
        // Two molecules sharing a common scaffold should increment its count
        let mol1 = parse("c1ccc2c(c1)ccc1ccccc12").unwrap(); // Anthracene
        let mol2 = parse("c1ccc2c(c1)ccc1ccccc12").unwrap(); // Same anthracene
        let network = scaffold_network_with_counts(&[mol1, mol2]);

        // Both molecules have the same scaffold hierarchy,
        // so each scaffold should be counted twice
        for count in &network.counts {
            assert!(
                *count >= 1,
                "each scaffold should appear at least once (but may appear multiple times)"
            );
        }

        // At least the root scaffold should appear twice
        let has_count_2_or_more = network.counts.iter().any(|&c| c >= 2);
        assert!(
            has_count_2_or_more,
            "at least one scaffold should have count >= 2 from duplicate molecules"
        );
    }

    #[test]
    fn scaffold_network_with_counts_parent_relationships() {
        // Scaffold network should have parent relationships
        let mol = parse("c1ccc2c(c1)ccc1ccccc12").unwrap(); // Anthracene (3 fused rings)
        let network = scaffold_network_with_counts(&[mol]);

        // Should have multiple scaffolds (one for each layer)
        assert!(
            network.scaffolds.len() > 1,
            "Anthracene should yield multiple scaffold layers"
        );

        // Parent relationships should form a chain (root has None, intermediate have Some)
        let root_count = network.parents.iter().filter(|p| p.is_none()).count();
        assert_eq!(root_count, 1, "should have exactly one root scaffold");
    }

    #[test]
    fn scaffold_network_with_counts_multiple_molecules() {
        // Multiple molecules with different scaffolds
        let benzene = parse("c1ccccc1").unwrap();
        let naphthalene = parse("c1ccc2ccccc2c1").unwrap();
        let network = scaffold_network_with_counts(&[benzene, naphthalene]);

        // Should have scaffolds from both molecules
        assert!(
            network.scaffolds.len() >= 2,
            "network should include scaffolds from both benzene and naphthalene"
        );

        // Benzene appears once, naphthalene scaffolds appear once each
        assert_eq!(
            network.scaffolds.len(),
            network.counts.len(),
            "lengths must match"
        );
    }
}
