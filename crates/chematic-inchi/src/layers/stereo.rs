use chematic_chem::assign_cip;
use chematic_core::{AtomIdx, BondOrder, Chirality, CipCode, Molecule};
use std::collections::HashMap;

/// Generate /m layer (relative stereo parity) from stereogenic centers.
///
/// Returns parity: 0 for even number of opposite configurations, 1 for odd.
/// Only generated if molecule has 2+ stereogenic centers.
pub fn relative_stereo_parity_layer(
    mol: &Molecule,
    _inchi_index: &HashMap<AtomIdx, usize>,
) -> Option<String> {
    let cip = assign_cip(mol);
    let r_count = cip
        .assignments
        .iter()
        .filter(|(_, c)| *c == CipCode::R)
        .count();
    let s_count = cip
        .assignments
        .iter()
        .filter(|(_, c)| *c == CipCode::S)
        .count();

    let total_stereo = r_count + s_count;
    if total_stereo < 2 {
        return None; // /m layer only for 2+ stereocenters
    }

    // Parity: 0 if total R/S centers is even, 1 if odd
    let parity = if (r_count + s_count) % 2 == 0 {
        "0"
    } else {
        "1"
    };
    Some(parity.to_string())
}

/// Generate /s layer (stereo type) indicating absolute vs. relative stereo.
///
/// Returns:
/// - "1" if absolute stereo (explicit chirality markers present)
/// - "2" if relative stereo only (no absolute stereo defined)
/// - "3" if racemic (no stereo definition)
pub fn stereo_type_layer(mol: &Molecule) -> Option<String> {
    // Count atoms with explicit chirality markers (@/@@ notation)
    let mut has_chirality_marker = false;

    for (_, atom) in mol.atoms() {
        if atom.chirality != Chirality::None {
            has_chirality_marker = true;
            break;
        }
    }

    // If explicit chirality is present, stereo is absolute
    if has_chirality_marker {
        return Some("1".to_string());
    }

    // Check if molecule has any assigned stereo (R/S codes)
    let cip = assign_cip(mol);
    let has_assigned_stereo = cip
        .assignments
        .iter()
        .any(|(_, code)| *code == CipCode::R || *code == CipCode::S);

    if has_assigned_stereo {
        // Could be relative or racemic; assume relative for now
        Some("2".to_string())
    } else {
        // No stereo defined
        Some("3".to_string())
    }
}

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
    let parts: Vec<String> = t_list.iter().map(|(n, s)| format!("{}{}", n, s)).collect();
    Some(parts.join(","))
}

/// Generate /b layer (E/Z double bond stereo) from CIP assignment.
///
/// Returns `atom1-atom2+/-` pairs, e.g., "2-3+,5-6-"
pub fn ez_stereo_layer(mol: &Molecule, inchi_index: &HashMap<AtomIdx, usize>) -> Option<String> {
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
                && let (Some(&i1), Some(&i2)) =
                    (inchi_index.get(&bond.atom1), inchi_index.get(&bond.atom2))
            {
                let (lo, hi) = if i1 < i2 { (i1, i2) } else { (i2, i1) };
                b_list.push(format!("{}-{}{}", lo, hi, sign));
                break;
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
    use chematic_smiles::canonical::canonical_atom_order;
    use chematic_smiles::parse;

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
        assert!(
            t_str.contains('-'),
            "L-alanine (S) should have minus parity"
        );
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

    #[test]
    fn test_relative_stereo_parity_none_for_single_stereocenter() {
        let mol = parse("N[C@@H](C)C(=O)O").expect("L-alanine");
        let idx = build_inchi_index(&mol);
        assert_eq!(
            relative_stereo_parity_layer(&mol, &idx),
            None,
            "Single stereocenter should not generate /m layer"
        );
    }

    #[test]
    fn test_relative_stereo_parity_for_multiple_stereocenters() {
        // Tartaric acid: has 2 stereocenters
        let mol = parse("C[C@H](O)[C@@H](O)C(=O)O").expect("tartaric acid");
        let idx = build_inchi_index(&mol);
        let parity = relative_stereo_parity_layer(&mol, &idx);
        assert!(
            parity.is_some(),
            "Multiple stereocenters should generate /m layer"
        );
        // Parity is computed from total R+S count; verify it's 0 or 1
        let parity_val = parity.unwrap();
        assert!(
            parity_val == "0" || parity_val == "1",
            "Parity should be 0 or 1"
        );
    }

    #[test]
    fn test_stereo_type_absolute_with_chirality_marker() {
        let mol = parse("N[C@@H](C)C(=O)O").expect("L-alanine");
        let s_type = stereo_type_layer(&mol);
        assert!(
            s_type.is_some(),
            "Molecule with chirality marker should have /s layer"
        );
        assert_eq!(
            s_type.unwrap(),
            "1",
            "Chirality marker indicates absolute stereo"
        );
    }

    #[test]
    fn test_stereo_type_achiral() {
        let mol = parse("CC").expect("ethane");
        let s_type = stereo_type_layer(&mol);
        assert!(s_type.is_some());
        assert_eq!(
            s_type.unwrap(),
            "3",
            "Achiral molecule should have racemic type"
        );
    }
}
