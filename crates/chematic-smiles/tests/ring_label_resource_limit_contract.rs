use chematic_smiles::{SmilesError, parse};

#[test]
fn too_many_open_ring_labels_are_rejected_without_constructing_a_molecule() {
    let mut smiles = String::from("C");
    for label in 100..=100_100 {
        smiles.push_str(&format!("%({label})"));
    }
    let error = parse(&smiles).err().expect("open ring label limit");
    assert!(matches!(
        error,
        SmilesError::ResourceLimit {
            resource: "open ring labels",
            actual: 100_001,
            limit: 100_000
        }
    ));
}

#[test]
fn oversized_ring_label_is_a_typed_resource_error() {
    let error = parse("C%(4294967296)").err().expect("ring label overflow");
    assert!(matches!(
        error,
        SmilesError::ResourceLimit {
            resource: "ring label value",
            ..
        }
    ));
}
