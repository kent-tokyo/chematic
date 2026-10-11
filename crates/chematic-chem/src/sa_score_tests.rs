#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> Molecule {
        parse(s).unwrap()
    }

    #[test]
    fn sa_score_in_range_1_to_10() {
        for s in &["C", "CC(=O)O", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"] {
            let score = sa_score(&mol(s));
            assert!(
                (1.0..=10.0).contains(&score),
                "SA score for {s} = {score:.2} is out of [1, 10]"
            );
        }
    }

    #[test]
    fn simple_molecule_lower_than_complex() {
        let simple = sa_score(&mol("CC"));
        let complex = sa_score(&mol("C1CC2(CC1)CCCC2")); // spiro bicyclic
        assert!(
            simple <= complex,
            "simple {simple:.2} should be ≤ complex {complex:.2}"
        );
    }

    #[test]
    fn methane_is_easy() {
        let score = sa_score(&mol("C"));
        assert!(
            score < 6.0,
            "methane SA score should be < 6, got {score:.2}"
        );
    }

    #[test]
    fn taxol_harder_than_aspirin() {
        let aspirin = sa_score(&mol("CC(=O)Oc1ccccc1C(=O)O"));
        let taxol = sa_score(&mol(
            "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(=C2OC1=O)C)(C)C)OC(=O)C5=CC=CC=C5)OC(=O)C4)O)(C)C",
        ));
        assert!(
            aspirin < taxol,
            "aspirin ({aspirin:.2}) should be easier than taxol ({taxol:.2})"
        );
    }
}
