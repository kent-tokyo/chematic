use super::*;

#[test]
fn graph_verification_rejects_element_and_hydrogen_loss() {
    let source = chematic_smiles::parse("C").unwrap();
    assert!(verify_graph(&source, &source).is_ok());
    let different_element = chematic_smiles::parse("N").unwrap();
    assert!(
        verify_graph(&source, &different_element)
            .unwrap_err()
            .to_string()
            .contains("indexed graph")
    );
    let mut different_hydrogens = source.clone();
    different_hydrogens.set_hydrogen_count(chematic_core::AtomIdx(0), Some(0));
    assert!(
        verify_graph(&source, &different_hydrogens)
            .unwrap_err()
            .to_string()
            .contains("hydrogen count")
    );
}

#[test]
fn unsupported_cml_bond_order_cannot_enter_projection_profile() {
    let cml = r#"<molecule><atomArray><atom id="a1" elementType="C"/><atom id="a2" elementType="C"/></atomArray><bondArray><bond atomRefs2="a1 a2" order="1.5"/></bondArray></molecule>"#;
    assert!(
        profile_name(cml)
            .unwrap_err()
            .to_string()
            .contains("outside single/double/triple")
    );
}
