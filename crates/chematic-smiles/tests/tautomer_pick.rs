use chematic_smiles::{RdkitTautomerPickError as Error, parse, rdkit_pick_tautomer};

#[test]
fn score_precedes_lexical_order_without_enumeration() {
    let enol = parse("C=C(C)O").unwrap();
    let keto = parse("CC(C)=O").unwrap();
    let pick = rdkit_pick_tautomer([&enol, &keto], 2).unwrap();
    assert_eq!(
        (pick.index, pick.score, pick.smiles.as_str()),
        (1, 5, "CC(C)=O")
    );
    let singleton = rdkit_pick_tautomer([&enol], 1).unwrap();
    assert_eq!(singleton.smiles, "C=C(C)O");
}

#[test]
fn lexical_ties_and_exact_ties() {
    let a = parse("CCC").unwrap();
    let b = parse("CC(C)O").unwrap();
    for inputs in [[&a, &b], [&b, &a]] {
        let pick = rdkit_pick_tautomer(inputs, 2).unwrap();
        assert_eq!(pick.score, 2);
        assert_eq!(pick.smiles, "CC(C)O");
    }
    assert_eq!(rdkit_pick_tautomer([&a, &a], 2).unwrap().index, 0);
}

#[test]
fn empty_inputs_and_invalid_limits() {
    assert_eq!(rdkit_pick_tautomer([], 1), Err(Error::Empty));
    assert_eq!(rdkit_pick_tautomer([], 0), Err(Error::InvalidLimit));
    assert_eq!(rdkit_pick_tautomer([], 4097), Err(Error::InvalidLimit));
    // An empty molecule is a legitimate candidate, distinct from no candidates.
    let empty = chematic_core::MoleculeBuilder::new().build();
    assert_eq!(rdkit_pick_tautomer([&empty], 1).unwrap().smiles, "");
}

#[test]
fn overflow_consumes_one_extra_and_never_returns_a_partial_winner() {
    let mol = parse("CCO").unwrap();
    let consumed = std::cell::Cell::new(0);
    let repeated = std::iter::repeat(&mol).inspect(|_| consumed.set(consumed.get() + 1));
    assert_eq!(
        rdkit_pick_tautomer(repeated, 3),
        Err(Error::TooManyCandidates { limit: 3 })
    );
    assert_eq!(consumed.get(), 4);
}

#[test]
fn invalid_candidate_is_indexed_even_when_it_cannot_win() {
    let valid = parse("c1ccccc1").unwrap();
    let invalid = parse("[CH5]").unwrap();
    assert!(matches!(
        rdkit_pick_tautomer([&valid, &invalid], 2),
        Err(Error::Candidate { index: 1, .. })
    ));
    let large = parse(&"C".repeat(1025)).unwrap();
    assert_eq!(
        rdkit_pick_tautomer([&large], 1),
        Err(Error::CandidateTooLarge { index: 0 })
    );
}
