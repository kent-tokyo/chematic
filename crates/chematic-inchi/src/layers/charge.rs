use chematic_core::{AtomIdx, Molecule};

/// Generate charge layer (/q) for InChI.
/// Returns None if net formal charge is zero.
pub fn charge_layer(mol: &Molecule) -> Option<String> {
    let mut total_charge: i32 = 0;
    for i in 0..mol.atom_count() {
        let atom_idx = AtomIdx(i as u32);
        let atom = mol.atom(atom_idx);
        total_charge += atom.charge as i32;
    }

    if total_charge == 0 {
        None
    } else if total_charge > 0 {
        Some(format!("+{}", total_charge))
    } else {
        Some(format!("{}", total_charge))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_charge_neutral() {
        let mol = parse("c1ccccc1").expect("benzene");
        assert_eq!(charge_layer(&mol), None);
    }

    #[test]
    fn test_charge_positive() {
        let mol = parse("[NH4+]").expect("ammonium");
        let q_layer = charge_layer(&mol);
        assert_eq!(q_layer, Some("+1".to_string()));
    }

    #[test]
    fn test_charge_negative() {
        let mol = parse("[O-]").expect("hydroxide");
        let q_layer = charge_layer(&mol);
        assert_eq!(q_layer, Some("-1".to_string()));
    }
}
