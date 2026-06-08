use chematic_core::Molecule;

/// Generate the formula prefix for InChI.
/// Uses Hill convention: C first, H second, then alphabetical.
pub fn formula_layer(mol: &Molecule) -> String {
    mol.total_formula()
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_formula_methane() {
        let mol = parse("C").expect("methane");
        assert_eq!(formula_layer(&mol), "CH4");
    }

    #[test]
    fn test_formula_ethane() {
        let mol = parse("CC").expect("ethane");
        assert_eq!(formula_layer(&mol), "C2H6");
    }

    #[test]
    fn test_formula_benzene() {
        let mol = parse("c1ccccc1").expect("benzene");
        assert_eq!(formula_layer(&mol), "C6H6");
    }

    #[test]
    fn test_formula_ethanol() {
        let mol = parse("CCO").expect("ethanol");
        assert_eq!(formula_layer(&mol), "C2H6O");
    }
}
