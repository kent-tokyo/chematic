use super::*;

#[test]
fn radical_and_unknown_charge_codes_decode_as_neutral() {
    assert_eq!(decode_charge(4), 0, "doublet radical");
    assert_eq!(decode_charge(0), 0);
    assert_eq!(decode_charge(9), 0);
    assert_eq!(decode_charge(1), 3);
    assert_eq!(decode_charge(7), -3);
}

#[test]
fn zero_mass_difference_in_any_spelling_is_no_isotope() {
    assert_eq!(decode_mass_difference(Element::C, " 0"), None);
    assert_eq!(decode_mass_difference(Element::C, "00"), None);
    assert_eq!(decode_mass_difference(Element::C, "x"), None);
    assert_eq!(decode_mass_difference(Element::C, " 1"), Some(13));
    assert_eq!(decode_mass_difference(Element::H, "-3"), None);
}
