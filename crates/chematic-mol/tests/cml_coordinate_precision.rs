//! Original two-atom serialization cases; MIT OR Apache-2.0.
//! Oracle: supplied binary64 values, compared bitwise after both CML readers.
//! Representation limits below do not assert physically meaningful geometries.
use chematic_mol::{parse_cml, parse_cml_strict, write_cml};

const SOURCE: &str = "<molecule>\n<atom id=\"a\" elementType=\"C\"/>\n<atom id=\"b\" elementType=\"O\"/>\n<bond atomRefs2=\"a b\" order=\"1\"/>\n</molecule>";

fn assert_roundtrip(values: &[f64]) {
    let (mol, _) = parse_cml_strict(SOURCE).unwrap();
    for &x in values {
        let coords = [(x, -x), (1.523456789012345, -0.334567890123456)];
        let text = write_cml(&mol, Some(&coords));
        for parser in [parse_cml, parse_cml_strict] {
            let (after, got) = parser(&text).unwrap();
            assert_eq!((after.atom_count(), after.bond_count()), (2, 1));
            for (a, b) in coords.iter().zip(&got) {
                assert_eq!(a.0.to_bits(), b.0.to_bits(), "x changed: {x}");
                assert_eq!(a.1.to_bits(), b.1.to_bits(), "y changed: {x}");
            }
            assert_eq!(write_cml(&after, Some(&got)), text);
        }
    }
}

#[test]
fn supplied_layout_coordinates_survive_without_decimal_quantization() {
    assert_roundtrip(&[
        0.0,
        -0.0,
        0.1,
        -0.2,
        0.123456789012345,
        -1.23456789012345,
        1.23456789012345e-8,
        1000000.123456789,
        0.000049999999999,
        -0.000050000000001,
        0.3333333333333333,
    ]);
}

#[test]
fn binary64_boundaries_and_neighboring_values_remain_distinguishable() {
    assert_roundtrip(&[
        f64::from_bits(1),
        f64::from_bits(2),
        f64::MIN_POSITIVE,
        f64::from_bits(f64::MIN_POSITIVE.to_bits() - 1),
        f64::MAX,
        f64::MIN,
        1e-100,
        1e100,
        f64::from_bits(1.0f64.to_bits() - 1),
        1.0,
        f64::from_bits(1.0f64.to_bits() + 1),
    ]);
}
