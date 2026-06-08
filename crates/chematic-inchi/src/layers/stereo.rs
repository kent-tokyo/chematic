use chematic_chem::assign_cip;
use chematic_core::{CipCode, Molecule, AtomIdx, BondOrder};
use std::collections::HashMap;

/// Generate /t layer (tetrahedral stereo) from R/S assignment.
///
/// Returns `atom_index-parity` pairs, e.g., "1-,2+"
pub fn tetrahedral_stereo_layer(
    mol: &Molecule,
    inchi_index: &HashMap<AtomIdx, usize>,
) -> Option<String> {
    let cip = assign_cip(mol);
    let mut t_list: Vec<(usize, char)> = vec![];

    for (atom_idx, cip_code) in &cip.assignments {
        let sign = match cip_code {
            CipCode::R => '+',
            CipCode::S => '-',
            _ => continue, // Skip E/Z codes
        };
        if let Some(&inchi_num) = inchi_index.get(atom_idx) {
            t_list.push((inchi_num, sign));
        }
    }

    if t_list.is_empty() {
        return None;
    }

    t_list.sort_by_key(|(n, _)| *n);
    let parts: Vec<String> = t_list
        .iter()
        .map(|(n, s)| format!("{}{}", n, s))
        .collect();
    Some(parts.join(","))
}

/// Generate /b layer (E/Z double bond stereo) from CIP assignment.
///
/// Returns `atom1-atom2+/-` pairs, e.g., "2-3+,5-6-"
pub fn ez_stereo_layer(
    mol: &Molecule,
    inchi_index: &HashMap<AtomIdx, usize>,
) -> Option<String> {
    let cip = assign_cip(mol);
    let mut b_list: Vec<String> = vec![];

    for (atom_idx, cip_code) in &cip.assignments {
        let sign = match cip_code {
            CipCode::Z => '+',
            CipCode::E => '-',
            _ => continue, // Skip R/S codes
        };

        // Find the double bond where atom_idx is one endpoint
        for (_bond_idx, bond) in mol.bonds() {
            if (bond.atom1 == *atom_idx || bond.atom2 == *atom_idx)
                && bond.order == BondOrder::Double
            {
                if let (Some(&i1), Some(&i2)) = (
                    inchi_index.get(&bond.atom1),
                    inchi_index.get(&bond.atom2),
                ) {
                    let (lo, hi) = if i1 < i2 { (i1, i2) } else { (i2, i1) };
                    b_list.push(format!("{}-{}{}", lo, hi, sign));
                    break;
                }
            }
        }
    }

    if b_list.is_empty() {
        return None;
    }

    b_list.sort();
    b_list.dedup();
    Some(b_list.join(","))
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;
    use chematic_smiles::canonical::canonical_atom_order;

    fn build_inchi_index(mol: &Molecule) -> HashMap<AtomIdx, usize> {
        let canonical_order = canonical_atom_order(mol);
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
        inchi_index
    }

    #[test]
    fn test_tetrahedral_l_alanine() {
        // L-alanine is S configuration
        let mol = parse("N[C@@H](C)C(=O)O").expect("L-alanine");
        let idx = build_inchi_index(&mol);
        let t = tetrahedral_stereo_layer(&mol, &idx);
        assert!(t.is_some());
        let t_str = t.unwrap();
        assert!(t_str.contains('-'), "L-alanine (S) should have minus parity");
    }

    #[test]
    fn test_tetrahedral_d_alanine() {
        // D-alanine is R configuration
        let mol = parse("N[C@H](C)C(=O)O").expect("D-alanine");
        let idx = build_inchi_index(&mol);
        let t = tetrahedral_stereo_layer(&mol, &idx);
        assert!(t.is_some());
        let t_str = t.unwrap();
        assert!(t_str.contains('+'), "D-alanine (R) should have plus parity");
    }

    #[test]
    fn test_tetrahedral_none_for_achiral() {
        let mol = parse("CC").expect("ethane");
        let idx = build_inchi_index(&mol);
        assert_eq!(
            tetrahedral_stereo_layer(&mol, &idx),
            None,
            "Ethane has no stereo centers"
        );
    }

    #[test]
    fn test_ez_none_for_no_double_bonds() {
        let mol = parse("CC").expect("ethane");
        let idx = build_inchi_index(&mol);
        assert_eq!(
            ez_stereo_layer(&mol, &idx),
            None,
            "Ethane has no double bonds"
        );
    }
}
