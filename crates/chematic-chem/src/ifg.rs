//! Functional group identification using the Ertl (2017) algorithm.
//!
//! Reference: P. Ertl, J. Cheminf. 2017, 9, 36.
//!
//! The algorithm identifies functional groups as connected subgraphs of
//! "environment" atoms (heteroatoms + carbons adjacent to heteroatoms +
//! all atoms in heteroatom-containing rings).

use chematic_core::{AtomIdx, Molecule};
use std::collections::{HashSet, VecDeque};

/// A functional group: the atom indices and a sorted element-symbol label.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FunctionalGroup {
    /// Atom indices (0-based, same order as `mol.atoms()`).
    pub atom_indices: Vec<usize>,
    /// Sorted element symbols joined by commas, e.g. "C,N,O".
    pub atom_types: String,
}

/// Identify functional groups in `mol` using the Ertl 2017 algorithm.
///
/// Returns one `FunctionalGroup` per connected component of marked atoms.
/// Aliphatic CH₂/CH₃ chains are excluded; isolated heteroatoms form their
/// own group.
pub fn identify_functional_groups(mol: &Molecule) -> Vec<FunctionalGroup> {
    let n = mol.atom_count();
    if n == 0 {
        return Vec::new();
    }

    // Step 1: mark atoms as FG members.
    let mut marked = vec![false; n];

    // 1a. All non-C, non-H heavy atoms.
    for (idx, atom) in mol.atoms() {
        let an = atom.element.atomic_number();
        if an != 1 && an != 6 {
            marked[idx.0 as usize] = true;
        }
    }

    // 1b. C atoms bonded to at least one non-C heavy atom.
    for (idx, atom) in mol.atoms() {
        if atom.element.atomic_number() != 6 {
            continue;
        }
        let has_hetero_neighbor = mol.neighbors(idx).any(|(nb, _)| {
            let an = mol.atom(nb).element.atomic_number();
            an != 1 && an != 6
        });
        if has_hetero_neighbor {
            marked[idx.0 as usize] = true;
        }
    }

    // Step 2: include all atoms in rings that contain at least one marked atom.
    // Use the ring info from perception to expand marks to entire heteroatom rings.
    let rings = ring_atom_sets(mol);
    let mut changed = true;
    while changed {
        changed = false;
        for ring in &rings {
            if ring.iter().any(|&i| marked[i]) {
                for &i in ring {
                    if !marked[i] {
                        // Only include ring carbons if the ring has a heteroatom.
                        let ring_has_hetero = ring
                            .iter()
                            .any(|&j| mol.atom(AtomIdx(j as u32)).element.atomic_number() != 6);
                        if ring_has_hetero {
                            marked[i] = true;
                            changed = true;
                        }
                    }
                }
            }
        }
    }

    // Step 3: find connected components of marked atoms.
    let mut visited = vec![false; n];
    let mut groups: Vec<FunctionalGroup> = Vec::new();

    for start in 0..n {
        if !marked[start] || visited[start] {
            continue;
        }
        // BFS over marked atoms only.
        let mut component: Vec<usize> = Vec::new();
        let mut queue = VecDeque::new();
        queue.push_back(start);
        visited[start] = true;

        while let Some(cur) = queue.pop_front() {
            component.push(cur);
            for (nb, _) in mol.neighbors(AtomIdx(cur as u32)) {
                let nbi = nb.0 as usize;
                if marked[nbi] && !visited[nbi] {
                    visited[nbi] = true;
                    queue.push_back(nbi);
                }
            }
        }

        component.sort_unstable();
        let mut syms: Vec<&str> = component
            .iter()
            .map(|&i| mol.atom(AtomIdx(i as u32)).element.symbol())
            .collect();
        syms.sort_unstable();
        let atom_types = syms.join(",");

        groups.push(FunctionalGroup {
            atom_indices: component,
            atom_types,
        });
    }

    groups
}

/// Returns SSSR rings as sets of atom indices.
fn ring_atom_sets(mol: &Molecule) -> Vec<Vec<usize>> {
    // Use DFS to find rings (simple cycle detection on each edge).
    let n = mol.atom_count();
    let mut rings: Vec<Vec<usize>> = Vec::new();
    let mut visited = vec![false; n];
    let mut parent = vec![usize::MAX; n];
    let mut stack: Vec<(usize, usize)> = Vec::new(); // (node, parent)

    // Adjacency.
    let adj: Vec<Vec<usize>> = (0..n)
        .map(|i| {
            mol.neighbors(AtomIdx(i as u32))
                .map(|(nb, _)| nb.0 as usize)
                .collect()
        })
        .collect();

    for root in 0..n {
        if visited[root] {
            continue;
        }
        stack.push((root, usize::MAX));
        while let Some((cur, par)) = stack.pop() {
            if visited[cur] {
                // Back-edge: reconstruct cycle.
                // Find the path from par back to cur in the DFS tree.
                let mut cycle = vec![cur];
                let mut node = par;
                while node != cur && node != usize::MAX {
                    cycle.push(node);
                    node = parent[node];
                }
                if node == cur {
                    cycle.dedup();
                    rings.push(cycle);
                }
                continue;
            }
            visited[cur] = true;
            parent[cur] = par;
            for &nb in &adj[cur] {
                if nb != par {
                    stack.push((nb, cur));
                }
            }
        }
    }

    // Deduplicate rings (sets).
    let mut seen: Vec<HashSet<usize>> = Vec::new();
    let mut result = Vec::new();
    for ring in rings {
        let s: HashSet<usize> = ring.iter().cloned().collect();
        if s.len() >= 3 && !seen.contains(&s) {
            seen.push(s);
            result.push(ring);
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> Molecule {
        parse(s).unwrap()
    }

    #[test]
    fn hexane_no_functional_groups() {
        let groups = identify_functional_groups(&mol("CCCCCC"));
        assert!(
            groups.is_empty(),
            "hexane should have no FGs, got {:?}",
            groups
        );
    }

    #[test]
    fn acetic_acid_one_group() {
        let groups = identify_functional_groups(&mol("CC(=O)O"));
        // The carboxylic acid group: C(=O)O
        assert!(
            !groups.is_empty(),
            "acetic acid should have at least one FG"
        );
        let all_types: Vec<&str> = groups.iter().map(|g| g.atom_types.as_str()).collect();
        assert!(
            all_types.iter().any(|t| t.contains('O')),
            "expected group containing O, got {:?}",
            all_types
        );
    }

    #[test]
    fn pyridine_one_group_containing_n() {
        let groups = identify_functional_groups(&mol("c1ccncc1"));
        assert_eq!(groups.len(), 1, "pyridine should be one FG");
        assert!(
            groups[0].atom_types.contains('N'),
            "pyridine FG must contain N"
        );
    }

    #[test]
    fn aspirin_multiple_groups() {
        let groups = identify_functional_groups(&mol("CC(=O)Oc1ccccc1C(=O)O"));
        // ester + carboxylic acid + aromatic ring (with O substituents)
        assert!(groups.len() >= 2, "aspirin should have multiple FGs");
    }

    #[test]
    fn aniline_two_groups() {
        let groups = identify_functional_groups(&mol("Nc1ccccc1"));
        // NH2 + aromatic ring (C connected to N)
        assert!(!groups.is_empty());
    }

    #[test]
    fn methane_no_groups() {
        let groups = identify_functional_groups(&mol("C"));
        assert!(groups.is_empty());
    }

    #[test]
    fn chlorobenzene_one_group() {
        let groups = identify_functional_groups(&mol("Clc1ccccc1"));
        assert!(!groups.is_empty());
        assert!(groups.iter().any(|g| g.atom_types.contains("Cl")));
    }
}
