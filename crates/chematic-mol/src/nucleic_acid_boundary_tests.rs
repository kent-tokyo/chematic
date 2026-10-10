use crate::nucleic_acid::*;
use serde_json::{Value, json};

fn document() -> Value {
    let fixture: Value = serde_json::from_str(include_str!(
        "../../../validation/nucleic_acid_document_contract.json"
    ))
    .unwrap();
    fixture["valid_cases"][0]["document"].clone()
}

#[test]
fn invalid_nucleic_acid_topology_reports_the_source_path() {
    let original = document();
    for (path, bad, diagnostic) in [
        ("/schema", json!("bad"), "schema"),
        ("/atom_ids/1", json!("d1-5"), "duplicate"),
        ("/atom_ids/0", json!(""), "id"),
        ("/strands", json!([]), "strand"),
        ("/strands/0/id", json!(""), "id"),
        ("/strands/0/residues", json!([]), "residue"),
        ("/strands/0/residues/1/id", json!("d1"), "duplicate"),
        ("/strands/0/residues/0/atom_refs", json!([]), "atom_refs"),
        ("/strands/0/residues/0/atom_refs/0", json!("gone"), "gone"),
        (
            "/strands/0/residues/0/atom_refs/1",
            json!("d1-5"),
            "ambiguous",
        ),
        ("/strands/0/residues/0/sugar", json!("ribose"), "sugar"),
        ("/strands/0/residues/0/base", json!("U"), "base"),
        ("/strands/0/residues/0/base", json!("other"), "modification"),
        (
            "/strands/0/residues/0/modification",
            json!("5mC"),
            "modification",
        ),
        (
            "/strands/0/residues/0/modification",
            json!("unknown"),
            "unknown",
        ),
        ("/linkages/0/id", json!("dna-1"), "duplicate"),
        ("/linkages/0/kind", json!("other"), "topology"),
        ("/linkages/0/three_prime_residue_id", json!("gone"), "gone"),
        ("/linkages/0/five_prime_residue_id", json!("gone"), "gone"),
        ("/linkages/0/three_prime_atom_ref", json!("d2-5"), "atom"),
        ("/linkages/0/five_prime_atom_ref", json!("d1-5"), "atom"),
        ("/linkages", json!([]), "linkages"),
        (
            "/strands/0/residues/0/five_prime_linkage",
            json!("dl1"),
            "terminus",
        ),
        (
            "/strands/0/residues/1/three_prime_linkage",
            json!("dl1"),
            "terminus",
        ),
        (
            "/strands/0/residues/0/three_prime_linkage",
            Value::Null,
            "linkage",
        ),
        (
            "/strands/0/residues/1/five_prime_linkage",
            json!("gone"),
            "linkage",
        ),
    ] {
        let mut bad_document = original.clone();
        if let Some(slot) = bad_document.pointer_mut(path) {
            *slot = bad;
        } else {
            let (parent, key) = path.rsplit_once('/').unwrap();
            bad_document
                .pointer_mut(parent)
                .unwrap()
                .as_object_mut()
                .unwrap()
                .insert(key.into(), bad);
        }
        let error = NucleicAcidDocument::from_json(&bad_document).unwrap_err();
        assert!(error.to_string().contains(diagnostic), "{path}: {error}");
        assert_eq!(document(), original);
    }
}

#[test]
fn each_nucleic_acid_resource_limit_is_enforced() {
    let original = NucleicAcidDocument::from_json(&document()).unwrap();
    let defaults = NucleicAcidLimits::default();
    for limits in [
        NucleicAcidLimits {
            max_strands: 0,
            ..defaults.clone()
        },
        NucleicAcidLimits {
            max_residues: 1,
            ..defaults.clone()
        },
        NucleicAcidLimits {
            max_atom_ids: 3,
            ..defaults.clone()
        },
        NucleicAcidLimits {
            max_atom_refs_per_residue: 1,
            ..defaults.clone()
        },
        NucleicAcidLimits {
            max_linkages: 0,
            ..defaults.clone()
        },
        NucleicAcidLimits {
            max_annotations: 0,
            ..defaults.clone()
        },
        NucleicAcidLimits {
            max_annotation_bytes: 0,
            ..defaults.clone()
        },
        NucleicAcidLimits {
            max_id_bytes: 1,
            ..defaults.clone()
        },
    ] {
        let error = original.validate_with_limits(&limits).unwrap_err();
        assert!(matches!(
            error,
            NucleicAcidError::ResourceLimit { .. } | NucleicAcidError::InvalidId { .. }
        ));
        assert!(!error.to_string().is_empty());
    }
}
