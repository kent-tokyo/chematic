use crate::stereo_geometry::StereoConfiguration;
use crate::*;

#[test]
fn stereo_errors_identify_invalid_ligands_without_changing_the_configuration() {
    let duplicate =
        StereoConfiguration::new(StereoGeometry::Tetrahedral, [7, 11, 7, 23]).unwrap_err();
    assert!(duplicate.to_string().contains("duplicate ligand id 7"));
    let original = StereoConfiguration::new(StereoGeometry::Tetrahedral, [7, 11, 19, 23]).unwrap();
    let missing = original
        .renumber(|id| (id != 19).then_some(id))
        .unwrap_err();
    assert!(missing.to_string().contains("ligand id 19"));
    assert_eq!(original.renumber(Some).unwrap(), original);
    let mismatched = remap_tetrahedral_parity([7, 11, 19, 23], [7, 11, 19, 29]).unwrap_err();
    assert!(mismatched.to_string().contains("same ligand ids"));
}

#[test]
fn valence_diagnostics_preserve_atom_index_observed_valence_and_allowed_values() {
    let mut builder = MoleculeBuilder::new();
    builder.add_atom(Atom::new(Element::HE));
    let mut carbon = Atom::new(Element::C);
    carbon.hydrogen_count = Some(5);
    builder.add_atom(carbon);
    let errors = validate_valence(&builder.build());
    assert_eq!(errors.len(), 1);
    assert_eq!(errors[0].atom, AtomIdx(1));
    assert_eq!(errors[0].actual, 5);
    assert!(errors[0].to_string().contains("atom 1 has valence 5"));
    assert!(errors[0].to_string().contains("[4]"));
}

struct RequiresAtoms;
impl MoleculeExtension for RequiresAtoms {
    fn id(&self) -> &str {
        "test.requires-atoms.v1"
    }
    fn version(&self) -> u32 {
        1
    }
    fn run(&self, mol: &Molecule) -> Result<ExtensionValue, ExtensionError> {
        if mol.atom_count() == 0 {
            Err(ExtensionError::InvalidInput(
                "requires at least one atom".into(),
            ))
        } else {
            Ok(ExtensionValue::Scalar(mol.atom_count() as f64))
        }
    }
}

#[test]
fn registry_errors_keep_extension_identity_and_original_failure_context() {
    let mol = MoleculeBuilder::new().build();
    let mut registry = ExtensionRegistry::new();
    let missing = registry.run("missing.extension.v2", &mol).unwrap_err();
    assert!(
        missing
            .to_string()
            .contains("unknown extension id: missing.extension.v2")
    );
    registry.register(RequiresAtoms).unwrap();
    let duplicate = registry.register(RequiresAtoms).unwrap_err();
    assert!(
        duplicate
            .to_string()
            .contains("duplicate extension id: test.requires-atoms.v1")
    );
    let failed = registry.run("test.requires-atoms.v1", &mol).unwrap_err();
    assert!(matches!(failed, ExtensionError::InvalidInput(_)));
    assert!(
        failed
            .to_string()
            .contains("invalid extension input: requires at least one atom")
    );
    assert_eq!(registry.extensions().count(), 1);
}
