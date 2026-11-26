use chematic_mol::SemanticModel;
use serde_json::{Value, json};

fn model() -> Value {
    json!({
        "schema": "chematic.semantic.v1",
        "atom_ids": ["a1", "a2"], "bond_ids": [], "r_groups": [], "polymer_units": [],
        "s_groups": [{
            "id": "sg1", "kind": "SRU", "member_atoms": ["a1", "a2"],
            "linkage": {"id": "link1", "left_atom": "a1", "right_atom": "a2", "bond_order": 1}
        }]
    })
}

#[test]
fn s_group_linkage_requires_a_byte_sized_integer_order() {
    let mut document = model();
    for value in [json!(-1), json!(256), json!(1.5), json!("1"), Value::Null] {
        document["s_groups"][0]["linkage"]["bond_order"] = value;
        let error = SemanticModel::from_json(&document).unwrap_err();
        assert!(
            error.to_string().contains("bond_order must be a u8"),
            "{error}"
        );
    }
}
