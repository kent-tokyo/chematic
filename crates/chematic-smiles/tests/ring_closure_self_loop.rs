//! A ring closure whose opening and closing atom are the same (`C11`) would
//! create a self-loop bond. OpenSMILES forbids it and RDKit 2026.09 rejects
//! every input below with "duplicated ring closure N bonds atom A to itself";
//! the parser must fail closed instead of building a looped graph.

use chematic_smiles::{SmilesError, parse};

#[test]
fn self_loop_ring_closures_are_rejected() {
    for (smiles, ring_num, atom) in [
        ("C11", 1, 0),
        ("n11n11", 1, 0),
        ("CC11C", 1, 1),
        ("C%11%11", 11, 0),
        ("C%(300)%(300)", 300, 0),
        ("[CH4]11", 1, 0),
        ("C=1=1", 1, 0),
        ("C11CC1", 1, 0),
        ("c1ccccc1C22", 2, 6),
    ] {
        match parse(smiles) {
            Err(SmilesError::InvalidBracketAtom { detail, .. }) => assert_eq!(
                detail,
                format!("ring closure {ring_num} bonds atom {atom} to itself"),
                "{smiles}"
            ),
            Err(other) => panic!("{smiles}: expected self-loop rejection, got {other:?}"),
            Ok(_) => panic!("{smiles}: self-loop ring closure was accepted"),
        }
    }
}

#[test]
fn self_loop_error_message_names_ring_and_atom() {
    let Err(err) = parse("CC11C") else {
        panic!("CC11C must be rejected");
    };
    assert_eq!(
        err.to_string(),
        "invalid bracket atom at position 4: ring closure 1 bonds atom 1 to itself"
    );
}

#[test]
fn reopening_a_ring_digit_on_the_closing_atom_is_still_unmatched() {
    // `C1CC11`: the third atom closes ring 1 and then reopens it, which is
    // an unclosed ring (RDKit: "unclosed ring"), not a self-loop.
    assert!(matches!(
        parse("C1CC11"),
        Err(SmilesError::UnmatchedRingClosure { ring_num: 1, .. })
    ));
}

#[test]
fn legitimate_ring_closures_still_parse() {
    for smiles in ["C1CC1", "C1CCCCC1", "c1ccccc1", "C12CC1C2", "C1CC1C1CC1"] {
        let mol = parse(smiles).unwrap_or_else(|e| panic!("{smiles}: {e}"));
        for (_, bond) in mol.bonds() {
            assert_ne!(bond.atom1, bond.atom2, "{smiles}: self-loop bond");
        }
    }
}
