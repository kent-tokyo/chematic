use chematic_mol::{NucleicAcidLimits, validate_nucleic_acid_json};
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
