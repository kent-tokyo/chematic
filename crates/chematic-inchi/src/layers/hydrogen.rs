use chematic_core::{AtomIdx, Molecule};
use chematic_smiles::canonical::canonical_atom_order;
use std::collections::HashMap;

/// Generate hydrogen layer (/h) for InChI.
/// Returns None if no non-H atoms or if all atoms have no implicit hydrogens and no explicit H.
pub fn hydrogen_layer(mol: &Molecule) -> Option<String> {
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

    // Collect H counts for each atom in canonical order
    let h_counts: Vec<(usize, u8)> = canonical_order
        .iter()
        .filter_map(|&canon_idx| {
            let atom_idx = AtomIdx(canon_idx as u32);
            if let Some(&inchi_idx) = inchi_index.get(&atom_idx) {
                let h_count = mol.implicit_hydrogen_count(atom_idx);
                Some((inchi_idx, h_count))
            } else {
                None
            }
        })
        .collect();

    // Check if all atoms have 0 hydrogens
    if h_counts.iter().all(|(_, h)| *h == 0) {
        return None;
    }

    // Compress consecutive atoms with same H count into range notation
    let mut result = String::new();
    let mut i = 0;
    while i < h_counts.len() {
        let (start_idx, h_count) = h_counts[i];
        let mut end_idx = start_idx;

        // Find consecutive atoms with same H count
        while i + 1 < h_counts.len()
            && h_counts[i + 1].0 == end_idx + 1
            && h_counts[i + 1].1 == h_count
        {
            i += 1;
            end_idx = h_counts[i].0;
        }

        // Append to result
        if !result.is_empty() {
            result.push(',');
        }

        if start_idx == end_idx {
            // Single atom
            if h_count == 1 {
                result.push_str(&format!("{}H", start_idx));
            } else {
                result.push_str(&format!("{}H{}", start_idx, h_count));
            }
        } else {
            // Range of atoms
            if h_count == 1 {
                result.push_str(&format!("{}-{}H", start_idx, end_idx));
            } else {
                result.push_str(&format!("{}-{}H{}", start_idx, end_idx, h_count));
            }
        }

        i += 1;
    }

    if result.is_empty() {
        None
    } else {
        Some(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_hydrogen_methane() {
        let mol = parse("C").expect("methane");
        let h_layer = hydrogen_layer(&mol);
        assert_eq!(h_layer, Some("1H4".to_string()));
    }

    #[test]
    fn test_hydrogen_ethane() {
        let mol = parse("CC").expect("ethane");
        let h_layer = hydrogen_layer(&mol);
        assert_eq!(h_layer, Some("1-2H3".to_string()));
    }

    #[test]
    fn test_hydrogen_benzene() {
        let mol = parse("c1ccccc1").expect("benzene");
        let h_layer = hydrogen_layer(&mol);
        assert_eq!(h_layer, Some("1-6H".to_string()));
    }

    #[test]
    fn test_hydrogen_ethanol() {
        let mol = parse("CCO").expect("ethanol");
        let h_layer = hydrogen_layer(&mol);
        assert!(h_layer.is_some());
        let h_str = h_layer.unwrap();
        assert!(h_str.contains('H'));
    }
}
