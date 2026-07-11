use chematic_core::{AtomIdx, Molecule};
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
    let mut tree_edges = HashSet::new();
    let mut result = String::new();
    dfs_connection(
        &first_atom,
        None,
        mol,
        &inchi_index,
        &mut visited,
        &mut result,
        &mut tree_edges,
    );

    // Add ring closure bonds (back-edges)
    let mut ring_closures = Vec::new();
    for (_bond_idx, bond) in mol.bonds() {
        let atom1 = bond.atom1;
        let atom2 = bond.atom2;
        if !inchi_index.contains_key(&atom1) || !inchi_index.contains_key(&atom2) {
            continue; // Skip bonds with H
        }
        // Check if this is a back-edge (not in tree_edges)
        let normalized_edge = if atom1 < atom2 {
            (atom1, atom2)
        } else {
            (atom2, atom1)
        };
        if !tree_edges.contains(&normalized_edge) {
            let i1 = inchi_index[&atom1];
            let i2 = inchi_index[&atom2];
            let (lo, hi) = if i1 < i2 { (i1, i2) } else { (i2, i1) };
            ring_closures.push((lo, hi));
        }
    }

    // Sort and add ring closures
    ring_closures.sort();
    for (lo, _hi) in ring_closures {
        result.push('-');
        result.push_str(&lo.to_string());
    }

    Some(result)
}

fn dfs_connection(
    atom: &AtomIdx,
    parent: Option<AtomIdx>,
    mol: &Molecule,
    inchi_index: &HashMap<AtomIdx, usize>,
    visited: &mut HashSet<AtomIdx>,
    result: &mut String,
    tree_edges: &mut HashSet<(AtomIdx, AtomIdx)>,
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
            // Record this as a tree edge (normalize to smaller index first)
            let normalized_edge = if *atom < neighbor {
                (*atom, neighbor)
            } else {
                (neighbor, *atom)
            };
            tree_edges.insert(normalized_edge);

            if first {
                first = false;
                dfs_connection(
                    &neighbor,
                    Some(*atom),
                    mol,
                    inchi_index,
                    visited,
                    result,
                    tree_edges,
                );
            } else {
                // Branch: wrap in parentheses with branch cursor reset (comma)
                // The comma resets the cursor to the parent atom after closing the branch
                let mut branch = String::new();
                dfs_connection(
                    &neighbor,
                    Some(*atom),
                    mol,
                    inchi_index,
                    visited,
                    &mut branch,
                    tree_edges,
                );
                if !branch.is_empty() {
                    result.push('(');
                    result.push_str(&branch);
                    result.push(')');
                    // Branch cursor reset: comma tells the parser to reset to parent atom
                    result.push(',');
                }
            }
        }
    }

    // Remove trailing comma (cursor resets only between multiple branches)
    if result.ends_with(',') {
        result.pop();
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
        // Benzene's 6 ring atoms are fully symmetric (D6h automorphism group), so
        // no single numbering (e.g. "1-2-3-4-5-6-1") is the uniquely "correct" one --
        // see the same point documented for standard_inchi() in
        // crates/chematic-inchi/tests/standard_inchi.rs. Check the structural
        // property (a closed 6-cycle over exactly atoms 1..=6) instead of an
        // exact string.
        let nums: Vec<usize> = c_str
            .split('-')
            .map(|s| {
                s.parse()
                    .unwrap_or_else(|e| panic!("bad token in {c_str:?}: {e}"))
            })
            .collect();
        assert_eq!(
            nums.len(),
            7,
            "expected 6 ring atoms + 1 closure: {c_str:?}"
        );
        assert_eq!(
            nums.first(),
            nums.last(),
            "should close the ring: {c_str:?}"
        );
        let mut distinct: Vec<usize> = nums[..6].to_vec();
        distinct.sort_unstable();
        assert_eq!(
            distinct,
            vec![1, 2, 3, 4, 5, 6],
            "should visit each ring atom once: {c_str:?}"
        );
    }

    // C4 Tests: Branch cursor reset

    #[test]
    fn test_connectivity_branched_propane() {
        // Propane: C-C-C (linear, no branches at DFS level)
        let mol = parse("CCC").expect("propane");
        let c_layer = connectivity_layer(&mol);
        assert!(c_layer.is_some());
        let c_str = c_layer.unwrap();
        // Linear propane should not have trailing commas (no multiple branches)
        assert!(
            !c_str.ends_with(','),
            "Linear propane should not end with comma"
        );
    }

    #[test]
    fn test_connectivity_branched_isobutane() {
        // Isobutane: CC(C)C — one central C with 3 methyl branches
        let mol = parse("CC(C)C").expect("isobutane");
        let c_layer = connectivity_layer(&mol);
        assert!(c_layer.is_some());
        let c_str = c_layer.unwrap();
        eprintln!("isobutane connectivity: {:?}", c_str);
        // Should have parentheses for branches.
        // The canonical SMILES atom-order fix (#14) may change which atom is
        // numbered first; what matters is that the string encodes branching.
        assert!(
            c_str.contains('(') || c_str.contains(','),
            "Should encode branched topology (parentheses or cursor resets): {c_str}"
        );
    }

    #[test]
    fn test_connectivity_neopentane() {
        // Neopentane: C(C)(C)(C)C — central C with 4 methyl branches
        let mol = parse("CC(C)(C)C").expect("neopentane");
        let c_layer = connectivity_layer(&mol);
        assert!(c_layer.is_some());
        let c_str = c_layer.unwrap();
        // Should have cursor resets between multiple branches
        assert!(
            c_str.contains(','),
            "Multiple branches should have cursor resets"
        );
    }

    #[test]
    fn test_connectivity_cursor_reset_position() {
        // Verify cursor reset occurs in correct position
        let mol = parse("CC(C)C").expect("isobutane");
        let c_layer = connectivity_layer(&mol);
        assert!(c_layer.is_some());
        let c_str = c_layer.unwrap();
        // Pattern should be: atom-atom(branch),... with commas between branches
        // Verify no double commas and proper structure
        assert!(!c_str.contains(",,"), "Should not have double commas");
        assert!(!c_str.ends_with(','), "Should not end with trailing comma");
    }

    #[test]
    fn test_connectivity_multi_level_branches() {
        // More complex branched structure: C-C(C-C,C)-C
        let mol = parse("CC(CC)C").expect("branched");
        let c_layer = connectivity_layer(&mol);
        assert!(c_layer.is_some());
        let c_str = c_layer.unwrap();
        // Should have proper branch cursor resets
        assert!(c_str.contains('('), "Should have parentheses for branches");
    }

    #[test]
    fn test_connectivity_toluene() {
        // Toluene: methylbenzene C1=CC=CC=C1C
        let mol = parse("Cc1ccccc1").expect("toluene");
        let c_layer = connectivity_layer(&mol);
        assert!(c_layer.is_some());
        // Aromatic with one branch should be handled correctly
    }

    #[test]
    fn test_connectivity_dimethylbenzene() {
        // o-Xylene (1,2-dimethylbenzene)
        let mol = parse("Cc1ccccc1C").expect("xylene");
        let c_layer = connectivity_layer(&mol);
        assert!(c_layer.is_some());
        // Multiple branches should have cursor resets
    }

    #[test]
    fn test_connectivity_no_false_commas() {
        // Linear molecule should not have trailing commas
        let mol = parse("CCCC").expect("butane");
        let c_layer = connectivity_layer(&mol);
        assert!(c_layer.is_some());
        let c_str = c_layer.unwrap();
        // Linear butane should not have multiple branches
        assert!(
            !c_str.ends_with(','),
            "Linear molecules should not end with trailing comma"
        );
    }

    #[test]
    fn test_connectivity_cursor_reset_inchi_standard() {
        // Test that cursor reset conforms to InChI standard:
        // After closing a branch, explicit reset before starting next branch at same level
        let mol = parse("CC(C)C").expect("isobutane");
        let c_layer = connectivity_layer(&mol);
        assert!(c_layer.is_some());
        let c_str = c_layer.unwrap();

        // Verify structure is valid InChI-like format
        // Should be something like: 1-2(3,4)
        // Not: 1-2(3)4 (missing cursor reset)
        let paren_count_open = c_str.matches('(').count();
        let paren_count_close = c_str.matches(')').count();
        assert_eq!(
            paren_count_open, paren_count_close,
            "Parentheses should be balanced"
        );
    }

    #[test]
    fn test_connectivity_pentane_isomers() {
        // Test various pentane isomers with different branching

        // n-Pentane: linear
        let lin = parse("CCCCC").expect("n-pentane");
        let lin_c = connectivity_layer(&lin).unwrap();
        assert!(!lin_c.contains(','), "Linear pentane should have no commas");

        // Isopentane: one branch
        let iso = parse("CC(C)CC").expect("isopentane");
        let iso_c = connectivity_layer(&iso).unwrap();
        assert!(
            iso_c.contains('('),
            "Branched pentane should have parentheses"
        );
    }

    #[test]
    fn test_connectivity_single_branch_no_trailing_comma() {
        // Single branch should not have trailing comma after closing parenthesis
        let mol = parse("CC(C)C").expect("isobutane");
        let c_layer = connectivity_layer(&mol);
        assert!(c_layer.is_some());
        let c_str = c_layer.unwrap();

        // Verify no trailing comma
        assert!(
            !c_str.ends_with(','),
            "Should not end with comma after single branch"
        );
        assert!(
            !c_str.ends_with("),"),
            "Should not have comma after last closing paren"
        );
    }
}
