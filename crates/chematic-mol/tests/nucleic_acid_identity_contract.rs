use chematic_mol::{
    NucleicAcidDocument, NucleicAcidError, NucleicAcidLimits, validate_nucleic_acid_json,
};
use serde_json::{Value, json};

fn dna() -> Value {
    json!({
        "schema":"chematic.nucleic-acid.v1", "atom_ids":["a"],
        "strands":[{"id":"dna", "kind":"dna", "residues":[{"id":"r1", "base":"A", "sugar":"deoxyribose", "atom_refs":["a"]}]}],
        "linkages":[]
    })
}

#[test]
fn dna_rejects_incompatible_base_and_modification_without_rewriting() {
    for (base, modification) in [("U", None), ("other", None), ("A", Some("5mC"))] {
        let mut document = dna();
        document["strands"][0]["residues"][0]["base"] = json!(base);
        if let Some(modification) = modification {
            document["strands"][0]["residues"][0]["modification"] = json!(modification);
        }
        let report =
            validate_nucleic_acid_json(&document.to_string(), &NucleicAcidLimits::default());
        assert_eq!(report["ok"], false, "{report}");
        assert_eq!(
            report["error"]["code"], "invalid_residue_identity",
            "{report}"
        );
    }
}

#[test]
fn document_object_rejects_wrong_field_types_as_invalid_json() {
    let mut document = dna();
    document["strands"][0]["residues"] = json!("not an array");
    assert!(matches!(
        NucleicAcidDocument::from_json(&document),
        Err(NucleicAcidError::InvalidJson { .. })
    ));
}

#[test]
fn pseudouridine_requires_rna_and_a_compatible_base() {
    let mut document = dna();
    document["strands"][0]["kind"] = json!("rna");
    document["strands"][0]["residues"][0]["sugar"] = json!("ribose");
    document["strands"][0]["residues"][0]["modification"] = json!("pseudouridine");
    for base in ["U", "other"] {
        document["strands"][0]["residues"][0]["base"] = json!(base);
        assert!(NucleicAcidDocument::from_json(&document).is_ok(), "{base}");
    }
    document["strands"][0]["residues"][0]["base"] = json!("A");
    assert!(matches!(
        NucleicAcidDocument::from_json(&document),
        Err(NucleicAcidError::InvalidResidueIdentity { .. })
    ));
    document["strands"][0]["kind"] = json!("dna");
    document["strands"][0]["residues"][0]["sugar"] = json!("deoxyribose");
    document["strands"][0]["residues"][0]["base"] = json!("other");
    assert!(matches!(
        NucleicAcidDocument::from_json(&document),
        Err(NucleicAcidError::InvalidResidueIdentity { .. })
    ));
}
