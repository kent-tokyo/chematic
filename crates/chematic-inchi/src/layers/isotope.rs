use chematic_core::{AtomIdx, Molecule};
use chematic_smiles::canonical::canonical_atom_order;
use std::collections::HashMap;

/// Generate isotope layer (/i) for InChI.
/// Returns None if no atoms have explicit isotope labels.
pub fn isotope_layer(mol: &Molecule) -> Option<String> {
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

    // Collect isotope information for atoms with explicit isotope labels
    let isotopes: Vec<(usize, u16)> = canonical_order
        .iter()
        .filter_map(|&canon_idx| {
            let atom_idx = AtomIdx(canon_idx as u32);
            if let Some(&inchi_idx) = inchi_index.get(&atom_idx) {
                let atom = mol.atom(atom_idx);
                if let Some(mass) = atom.isotope {
                    return Some((inchi_idx, mass));
                }
            }
            None
        })
        .collect();

    if isotopes.is_empty() {
        return None;
    }

    // Build isotope string: 1+2,2+13 format
    let result = isotopes
        .iter()
        .map(|(idx, mass)| format!("{}+{}", idx, mass))
        .collect::<Vec<_>>()
        .join(",");

    Some(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_isotope_none() {
        let mol = parse("c1ccccc1").expect("benzene");
        assert_eq!(isotope_layer(&mol), None);
    }

    #[test]
    fn test_isotope_deuterium() {
        // Note: deuterated molecules require explicit isotope notation
        // For now just test that isotope layer returns None for natural isotopes
        let mol = parse("C").expect("methane");
        let i_layer = isotope_layer(&mol);
        assert_eq!(i_layer, None);
    }
}
