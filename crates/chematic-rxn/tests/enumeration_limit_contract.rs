use chematic_rxn::{LibraryConfig, LibraryError, enumerate_library};
use chematic_smiles::{canonical_smiles, parse};

#[test]
fn enumeration_obeys_product_caps_for_multi_product_reactions() {
    let template = "[C:1][C:2]>>[C:1].[C:2]";
    let ethane = parse("CC").unwrap();
    let config = |max_size| LibraryConfig {
        skip_failures: true,
        max_size,
    };
    let products = enumerate_library(template, vec![vec![ethane.clone()]], &config(None)).unwrap();
    assert_eq!(products.len(), 4);
    assert!(products.iter().all(|mol| canonical_smiles(mol) == "C"));
    assert!(matches!(
        enumerate_library(template, vec![vec![ethane]], &config(Some(2))),
        Err(LibraryError::EnumerationLimitExceeded)
    ));
}

#[test]
fn zero_capacity_rejects_a_nonempty_fragment_set() {
    let config = LibraryConfig {
        skip_failures: true,
        max_size: Some(0),
    };
    let result = enumerate_library("[N:1]>>[N:1]", vec![vec![parse("C").unwrap()]], &config);
    assert!(matches!(
        result,
        Err(LibraryError::EnumerationTooLarge(1, 0))
    ));
}
