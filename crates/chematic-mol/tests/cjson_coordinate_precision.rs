//! Finite f64 coordinates must survive numeric JSON output and re-reading.
use chematic_core::{Atom, Element, MoleculeBuilder};
use chematic_mol::{parse_cjson, parse_cjson_record, write_cjson};

#[test]
fn finite_coordinate_roundtrips_preserve_all_bits() {
    let mut builder = MoleculeBuilder::new();
    builder.add_atom(Atom::new(Element::C));
    let mol = builder.build();
    let mut bits = 0x243f_6a88_85a3_08d3_u64;
    for _ in 0..1024 {
        // Fixed finite samples across signs, exponents and mantissas. The
        // f64 source, not another JSON reader, defines the expected value.
        bits ^= bits << 13;
        bits ^= bits >> 7;
        bits ^= bits << 17;
        for value in [
            f64::from_bits(bits),
            f64::from_bits((bits & 0x800f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000),
        ] {
            if !value.is_finite() {
                continue;
            }
            let json = write_cjson(&mol, &[(value, -0.0, -value)]);
            let (_, coords) = parse_cjson(&json).unwrap();
            assert_eq!(coords[0].0.to_bits(), value.to_bits(), "x={value}");
            assert_eq!(coords[0].1.to_bits(), (-0.0_f64).to_bits());
            assert_eq!(coords[0].2.to_bits(), (-value).to_bits(), "z={}", -value);
            let record = parse_cjson_record(&json).unwrap();
            let strict = record.coordinates()[0];
            assert_eq!(strict.0.to_bits(), value.to_bits());
            assert_eq!(strict.1.to_bits(), (-0.0_f64).to_bits());
            assert_eq!(strict.2.to_bits(), (-value).to_bits());
        }
    }
}

#[test]
fn integral_bond_orders_are_json_integers_and_aromatic_orders_remain_fractional() {
    for (smiles, expected) in [("CC", 1), ("C=C", 2), ("C#C", 3)] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let text = write_cjson(&mol, &[]);
        let json: serde_json::Value = serde_json::from_str(&text).unwrap();
        assert_eq!(json["bonds"]["order"][0].as_u64(), Some(expected));
        let (after, _) = parse_cjson(&text).unwrap();
        assert_eq!(
            mol.bonds().next().unwrap().1.order,
            after.bonds().next().unwrap().1.order
        );
    }
    let aromatic = chematic_smiles::parse("c1ccncc1").unwrap();
    let json: serde_json::Value = serde_json::from_str(&write_cjson(&aromatic, &[])).unwrap();
    assert!(
        json["bonds"]["order"]
            .as_array()
            .unwrap()
            .iter()
            .all(|n| n.as_f64() == Some(1.5) && n.as_u64().is_none())
    );
}
