//! Exercise the bounded local namer over the checked-in special-chemistry set.

use chematic_iupac::{IupacError, name};
use chematic_smiles::parse;

const STRESS: &str = include_str!("../../../validation/cosmolkit_comparison/corpora/stress-v1.smi");

#[test]
fn local_namer_classifies_the_special_chemistry_corpus_without_panicking() {
    let mut parsed = 0usize;
    let mut named = 0usize;
    for line in STRESS.lines() {
        let smiles = line.split_whitespace().next().unwrap_or("");
        let Ok(mol) = parse(smiles) else { continue };
        parsed += 1;
        match name(&mol) {
            Ok(value) => {
                named += 1;
                assert!(!value.is_empty());
            }
            Err(IupacError::Empty | IupacError::NotSupported) => {}
        }
    }
    assert!(parsed > 3_500);
    assert!(named > 100);
}

#[test]
fn supported_family_matrix_has_stable_nonempty_names() {
    let families = [
        "C",
        "CC",
        "CCC",
        "CCCC",
        "C=C",
        "CC=C",
        "C#C",
        "CC#C",
        "CO",
        "CCO",
        "CC(C)O",
        "C=O",
        "CC=O",
        "CC(=O)C",
        "CC(=O)O",
        "CC(=O)OC",
        "CC(=O)N",
        "CN",
        "CCN",
        "C#N",
        "CC#N",
        "CS",
        "CCS",
        "CSC",
        "CF",
        "CCCl",
        "CCBr",
        "CCI",
        "C1CCCCC1",
        "CC1CCCCC1",
        "C1CCC2CCCCC2C1",
        "C1CCC2(CC1)CCCC2",
        "c1ccccc1",
        "n1ccccc1",
        "o1cccc1",
        "s1cccc1",
        "[nH]1cccc1",
        "c1ncc[nH]1",
        "n1ccnc(n1)",
        "C1COCC1",
        "C1CNCC1",
        "O1CCNCC1",
    ];
    let mut named = 0usize;
    for smiles in families {
        let mol = parse(smiles).unwrap_or_else(|error| panic!("{smiles}: {error}"));
        if let Ok(value) = name(&mol) {
            named += 1;
            assert!(!value.is_empty(), "{smiles}");
        }
    }
    assert!(named >= 30, "only {named} curated families were named");
}
