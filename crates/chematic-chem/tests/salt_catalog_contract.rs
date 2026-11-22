use chematic_chem::standardize::{SaltCatalog, prefer_organic, remove_salts};
use chematic_core::MoleculeBuilder;
use chematic_smiles::parse;

#[test]
fn custom_salt_pattern_is_applied_after_a_malformed_entry() {
    let xenon = parse("[Xe]").unwrap();
    let mut catalog = SaltCatalog::new();
    assert!(!catalog.is_salt(&xenon));
    catalog.add("malformed", "[");
    assert!(!catalog.is_salt(&xenon));
    catalog.add("xenon", "[#54]");
    assert!(catalog.is_salt(&xenon));
    assert!(!catalog.is_salt(&parse("[He]").unwrap()));
}

#[test]
fn fragment_selection_accepts_empty_molecules() {
    let empty = MoleculeBuilder::new().build();
    assert_eq!(remove_salts(&empty).atom_count(), 0);
    assert_eq!(prefer_organic(&empty).atom_count(), 0);
}
