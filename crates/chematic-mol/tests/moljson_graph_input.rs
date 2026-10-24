use chematic_core::{AtomIdx, BondOrder, implicit_hcount};
use chematic_mol::{parse_moljson, write_moljson};
use serde_json::Value;

fn packet() -> Value {
    serde_json::from_str(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../validation/moljson_graph_input/cases.json"
    )))
    .unwrap()
}

#[test]
fn supported_graph_fields_survive_three_read_write_cycles() {
    for case in packet()["valid"].as_array().unwrap() {
        let expected = &case["expected"];
        let mut text = case["input"].as_str().unwrap().to_owned();
        for _ in 0..3 {
            let mol = parse_moljson(&text).unwrap();
            assert_eq!(
                mol.atom_count(),
                expected["atoms"].as_array().unwrap().len()
            );
            assert_eq!(
                mol.bond_count(),
                expected["bonds"].as_array().unwrap().len()
            );
            for (i, e) in expected["atoms"].as_array().unwrap().iter().enumerate() {
                let idx = AtomIdx(i as u32);
                let a = mol.atom(idx);
                assert_eq!(a.element.symbol(), e["element"].as_str().unwrap());
                assert_eq!(i64::from(a.charge), e["charge"].as_i64().unwrap());
                assert_eq!(a.isotope.map(u64::from), e["isotope"].as_u64());
                assert_eq!(
                    u64::from(implicit_hcount(&mol, idx)),
                    e["hydrogens"].as_u64().unwrap()
                );
                assert_eq!(a.aromatic, e["aromatic"].as_bool().unwrap());
            }
            for ((_, b), e) in mol.bonds().zip(expected["bonds"].as_array().unwrap()) {
                let order = match e["order"].as_f64().unwrap() {
                    1.0 => BondOrder::Single,
                    1.5 => BondOrder::Aromatic,
                    2.0 => BondOrder::Double,
                    3.0 => BondOrder::Triple,
                    4.0 => BondOrder::Quadruple,
                    _ => panic!("fixture order"),
                };
                assert_eq!(b.order, order, "{}", case["id"]);
                assert_eq!(u64::from(b.atom1.0), e["source"].as_u64().unwrap());
                assert_eq!(u64::from(b.atom2.0), e["target"].as_u64().unwrap());
            }
            text = write_moljson(&mol);
        }
    }
}

#[test]
fn malformed_graph_fields_and_ambiguous_references_are_refused() {
    for case in packet()["invalid"].as_array().unwrap() {
        assert!(
            parse_moljson(case["input"].as_str().unwrap()).is_err(),
            "{}",
            case["id"]
        );
    }
}
