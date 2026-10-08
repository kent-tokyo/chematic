//! #769: a neutral oxygen or fluorine with more bonds than its element
//! permits is not a molecule. The parser used to accept `O=O1C=CC=C1`, and
//! aromaticity perception then aromatized the ring through the invalid
//! oxygen. It is now a typed `InvalidValence` error; charged and hypervalent
//! controls keep parsing.

use chematic_smiles::{SmilesError, parse};

#[test]
fn over_valent_neutral_oxygen_and_fluorine_are_rejected() {
    for (smiles, element, atom, valence, max) in [
        ("O=O1C=CC=C1", "O", 1, 4, 2),
        ("CO(C)C", "O", 1, 3, 2),
        ("C=O=C", "O", 1, 4, 2),
        ("[O](C)(C)C", "O", 0, 3, 2),
        ("[OH3]", "O", 0, 3, 2),
        ("[H]O([H])[H]", "O", 1, 3, 2),
        ("c1cco(C)c1", "O", 3, 3, 2),
        ("CC(C)=O(C)", "O", 3, 3, 2),
        ("F(C)C", "F", 0, 2, 1),
        // The one such row in the NCI first-5k corpus; RDKit rejects it too.
        ("CCO1=O=C1C2=CC=C(O2)C(C)(C)C", "O", 2, 4, 2),
    ] {
        assert_eq!(
            parse(smiles).err(),
            Some(SmilesError::InvalidValence {
                element,
                atom,
                valence,
                max
            }),
            "{smiles}"
        );
    }
    let message = parse("O=O1C=CC=C1").err().unwrap().to_string();
    assert!(message.contains("greater than permitted"), "{message}");
}

#[test]
fn charged_hypervalent_and_dative_controls_still_parse() {
    for smiles in [
        "C[O+](C)C",
        "c1cc[o+]cc1",
        "O=O",
        "c1ccoc1",
        "O",
        "[OH2]",
        "[O]",
        "[O-]C",
        "[2H]O[2H]",
        "CC(=O)[O-].[Na+]",
        "CS(=O)(=O)C",
        "CS(=O)(=O)(=O)C",
        "CN(=O)=O",
        "c1ccn(=O)cc1",
        "FC(F)(F)F",
        "[F-]",
        "ClF",
        "FCl(F)F",
        "[Fe]<-O(C)C",
        "C[O]~[Fe]",
    ] {
        assert!(parse(smiles).is_ok(), "{smiles}");
    }
}
