use chematic_core::{Molecule, AtomIdx};
use chematic_smiles::canonical::canonical_atom_order;
use std::collections::{HashMap, HashSet};

/// Generate connectivity layer (/c) for InChI.
/// Returns None if molecule has no non-H heavy atoms.
pub fn connectivity_layer(mol: &Molecule) -> Option<String> {
    // Get canonical atom ordering
    let canonical_order = canonical_atom_order(mol);

    // Filter out hydrogen atoms and create mapping to InChI indices (1-indexed)
    let mut inchi_index: HashMap<AtomIdx, usize> = HashMap::new();
    let mut inchi_num = 0;
    for &canon_idx in &canonical_order {
        let atom_idx = AtomIdx(canon_idx as u32);
        let atom = mol.atom(atom_idx);
        if atom.element.atomic_number() != 1 {
            inchi_num += 1;
            inchi_index.insert(atom_idx, inchi_num);
        }
    }

    if inchi_num == 0 {
        return None;
    }

    // If only one heavy atom, return just "1"
    if inchi_num == 1 {
        return Some("1".to_string());
    }

    // Build connection string via DFS from first atom
    let first_atom = canonical_order
        .iter()
        .find_map(|&idx| {
            let atom_idx = AtomIdx(idx as u32);
            if inchi_index.contains_key(&atom_idx) {
                Some(atom_idx)
            } else {
                None
            }
        })
        .expect("at least one heavy atom");

    let mut visited = HashSet::new();
    let mut result = String::new();
    dfs_connection(&first_atom, None, mol, &inchi_index, &mut visited, &mut result);

    Some(result)
}

fn dfs_connection(
    atom: &AtomIdx,
    parent: Option<AtomIdx>,
    mol: &Molecule,
    inchi_index: &HashMap<AtomIdx, usize>,
    visited: &mut HashSet<AtomIdx>,
    result: &mut String,
) {
    if visited.contains(atom) {
        return;
    }
    visited.insert(*atom);

    let my_index = inchi_index[atom];

    // Add current atom number (only if not already added by parent)
    if result.is_empty() {
        result.push_str(&my_index.to_string());
    } else if let Some(_parent_idx) = parent {
        result.push('-');
        result.push_str(&my_index.to_string());
    }

    // Get neighbors (heavy atoms only)
    let mut neighbors: Vec<AtomIdx> = mol
        .neighbors(*atom)
        .filter_map(|(n_idx, _bond_idx)| {
            let n_atom = mol.atom(n_idx);
            if n_atom.element.atomic_number() != 1 && inchi_index.contains_key(&n_idx) {
                Some(n_idx)
            } else {
                None
            }
        })
        .collect();

    // Sort neighbors to ensure deterministic order (by InChI index)
    neighbors.sort_by_key(|n| inchi_index[n]);

    let mut first = true;
    for &neighbor in &neighbors {
        if !visited.contains(&neighbor) && parent != Some(neighbor) {
            if first {
                first = false;
                dfs_connection(&neighbor, Some(*atom), mol, inchi_index, visited, result);
            } else {
                // Branch: wrap in parentheses
                let mut branch = String::new();
                dfs_connection(&neighbor, Some(*atom), mol, inchi_index, visited, &mut branch);
                if !branch.is_empty() {
                    result.push('(');
                    result.push_str(&branch);
                    result.push(')');
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_connectivity_methane() {
        let mol = parse("C").expect("methane");
        let c_layer = connectivity_layer(&mol);
        assert_eq!(c_layer, Some("1".to_string()));
    }

    #[test]
    fn test_connectivity_ethane() {
        let mol = parse("CC").expect("ethane");
        let c_layer = connectivity_layer(&mol);
        assert_eq!(c_layer, Some("1-2".to_string()));
    }

    #[test]
    fn test_connectivity_benzene() {
        let mol = parse("c1ccccc1").expect("benzene");
        let c_layer = connectivity_layer(&mol);
        assert!(c_layer.is_some());
        let c_str = c_layer.unwrap();
        // Just check that connectivity is generated for benzene
        assert!(!c_str.is_empty());
        assert!(c_str.contains("1"));
    }
}
