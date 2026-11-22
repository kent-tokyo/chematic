use chematic_core::MoleculeBuilder;
use chematic_iupac::{IupacError, name};

#[test]
fn naming_errors_distinguish_empty_and_out_of_scope_structures() {
    let empty = name(&MoleculeBuilder::new().build()).unwrap_err();
    assert_eq!(empty, IupacError::Empty);
    assert_eq!(empty.to_string(), "empty molecule");
    let mixture = chematic_smiles::parse("C.O").unwrap();
    let unsupported = name(&mixture).unwrap_err();
    assert_eq!(unsupported, IupacError::NotSupported);
    assert!(unsupported.to_string().contains("not supported"));
}
