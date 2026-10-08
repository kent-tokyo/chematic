//! A lowercase letter in a bracket SMARTS atom starts a two-letter symbol
//! only for the aromatic `se`, `as`, `te` and `si`, as in RDKit. Otherwise
//! it is the one-letter aromatic atom and the next letter is the next
//! primitive: `[ca]` is aromatic carbon AND aromatic, not calcium; `[cr6]`
//! is aromatic carbon in a six-membered ring, not chromium; `[na]` is not
//! sodium.
//! Found by the property-based SMIRKS fuzz lane (#754), which wrote
//! `[ca:3]` and `[na:1]` templates that RDKit matched and chematic did not.
//! Expected counts are RDKit 2026.03.6's `GetSubstructMatches`.

use chematic_smarts::{find_matches, parse_smarts};

fn count(smiles: &str, smarts: &str) -> usize {
    let mol = chematic_smiles::parse(smiles).unwrap();
    let query = parse_smarts(smarts).unwrap_or_else(|e| panic!("{smarts}: {e:?}"));
    find_matches(&query, &mol).len()
}

#[test]
fn lowercase_letters_after_an_aromatic_atom_are_primitives() {
    for (smiles, smarts, expected) in [
        ("Clc1ccccn1", "[ca]", 5),
        ("Clc1ccccn1", "[na]", 1),
        ("Clc1ccccn1", "[cr6]", 5),
        ("Clc1ccccn1", "[nr6]", 1),
        ("Clc1ccccn1", "[cH1a]", 4),
        ("Clc1ccccn1", "[cD2]", 4),
        ("Clc1ccccn1", "[Cl]", 1),
        ("c1ccsc1", "[sa]", 1),
        ("c1ccsc1", "[ca]", 4),
        ("c1cc[se]c1", "[se]", 1),
        ("c1cc[se]c1", "[cr6]", 0),
        ("C[Se]C", "[se]", 0),
        ("[Na+].[Cl-]", "[Na]", 1),
    ] {
        assert_eq!(count(smiles, smarts), expected, "{smarts} on {smiles}");
    }
}

#[test]
fn a_lowercase_pair_that_is_no_aromatic_symbol_is_not_an_element() {
    // RDKit rejects `[cl]`; it is not chlorine.
    assert!(parse_smarts("[cl]").is_err());
}
