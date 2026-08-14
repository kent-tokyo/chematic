//! Serde round-trip tests. Only compiled when the `serde` feature is on
//! (see `Cargo.toml`'s `[[test]] required-features`).

use chematic_core::Element;
use chematic_crystal::{
    CartesianCoord, FractionalCoord, Lattice, Occupancy, PeriodicSite, PeriodicStructure,
    SiteSpecies,
};

fn nacl_structure() -> PeriodicStructure {
    let lattice = Lattice::from_parameters(5.64, 5.64, 5.64, 90.0, 90.0, 90.0).unwrap();
    let sites = vec![
        PeriodicSite::new(
            vec![SiteSpecies::full(Element::NA)],
            FractionalCoord::new([0.0, 0.0, 0.0]),
            Some("Na1".to_string()),
        )
        .unwrap(),
        PeriodicSite::new(
            vec![SiteSpecies::full(Element::CL)],
            FractionalCoord::new([0.5, 0.5, 0.5]),
            Some("Cl1".to_string()),
        )
        .unwrap(),
    ];
    PeriodicStructure::new(lattice, sites).unwrap()
}

#[test]
fn lattice_json_round_trip() {
    let l = Lattice::from_parameters(5.0, 6.0, 7.0, 80.0, 95.0, 110.0).unwrap();
    let json = serde_json::to_string(&l).unwrap();
    let back: Lattice = serde_json::from_str(&json).unwrap();
    assert_eq!(l, back);
}

#[test]
fn periodic_structure_json_round_trip() {
    let s = nacl_structure();
    let json = serde_json::to_string(&s).unwrap();
    let back: PeriodicStructure = serde_json::from_str(&json).unwrap();
    assert_eq!(s, back);
}

#[test]
fn fractional_and_cartesian_coord_round_trip() {
    let f = FractionalCoord::new([0.1, 0.2, 0.3]);
    let json = serde_json::to_string(&f).unwrap();
    let back: FractionalCoord = serde_json::from_str(&json).unwrap();
    assert_eq!(f, back);

    let c = CartesianCoord::new([1.5, -2.5, 3.5]);
    let json = serde_json::to_string(&c).unwrap();
    let back: CartesianCoord = serde_json::from_str(&json).unwrap();
    assert_eq!(c, back);
}

#[test]
fn occupancy_round_trip() {
    let o = Occupancy::new(0.75).unwrap();
    let json = serde_json::to_string(&o).unwrap();
    assert_eq!(json, "0.75");
    let back: Occupancy = serde_json::from_str(&json).unwrap();
    assert_eq!(o.value(), back.value());
}

#[test]
fn site_species_round_trips_through_element_symbol() {
    let s = SiteSpecies::full(Element::PT);
    let json = serde_json::to_string(&s).unwrap();
    assert!(
        json.contains("\"Pt\""),
        "expected element symbol in JSON, got {json}"
    );
    let back: SiteSpecies = serde_json::from_str(&json).unwrap();
    assert_eq!(s, back);
}

// -- invalid numeric rejection ------------------------------------------

#[test]
fn lattice_deserialize_rejects_nan() {
    let json = r#"{"matrix":[[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,"NaN"]]}"#;
    // serde_json has no NaN literal token in valid JSON at all, so a raw
    // NaN can't even parse as a number here -- the more realistic "invalid
    // numeric value" case is a matrix that parses fine but fails Lattice's
    // own singular/near-singular checks, exercised below. This test
    // confirms serde_json itself rejects the malformed literal.
    assert!(serde_json::from_str::<Lattice>(json).is_err());
}

#[test]
fn lattice_deserialize_rejects_singular_matrix() {
    let json = r#"{"matrix":[[1.0,0.0,0.0],[0.0,1.0,0.0],[1.0,1.0,0.0]]}"#;
    let err = serde_json::from_str::<Lattice>(json).unwrap_err();
    assert!(err.to_string().contains("singular"));
}

#[test]
fn occupancy_deserialize_rejects_negative() {
    let err = serde_json::from_str::<Occupancy>("-0.5").unwrap_err();
    assert!(err.to_string().contains("occupancy"));
}

#[test]
fn fractional_coord_deserialize_rejects_infinity_via_serde_json_null() {
    // serde_json serializes non-finite f64 as `null`; deserializing that
    // back must fail (not silently produce Infinity/NaN), whether it fails
    // at the array-of-f64 stage or a hypothetical finite check.
    let json = "[1.0, null, 2.0]";
    assert!(serde_json::from_str::<FractionalCoord>(json).is_err());
}

#[test]
fn site_species_deserialize_rejects_unknown_element_symbol() {
    let json = r#"{"element":"Zz","occupancy":1.0}"#;
    let err = serde_json::from_str::<SiteSpecies>(json).unwrap_err();
    assert!(err.to_string().contains("Zz"));
}

#[test]
fn periodic_site_deserialize_rejects_occupancy_sum_over_tolerance() {
    let json = r#"{"species":[{"element":"Fe","occupancy":0.7},{"element":"Ni","occupancy":0.5}],"fractional":[0.0,0.0,0.0],"label":null}"#;
    let err = serde_json::from_str::<PeriodicSite>(json).unwrap_err();
    assert!(err.to_string().contains("occupanc"));
}

// -- field-name stability -------------------------------------------------

#[test]
fn periodic_structure_field_names_are_lattice_and_sites() {
    let s = nacl_structure();
    let json = serde_json::to_value(&s).unwrap();
    let obj = json.as_object().unwrap();
    assert!(obj.contains_key("lattice"));
    assert!(obj.contains_key("sites"));
    assert_eq!(obj.len(), 2);
}

#[test]
fn periodic_site_field_names_are_species_fractional_label() {
    let site = PeriodicSite::new(
        vec![SiteSpecies::full(Element::NA)],
        FractionalCoord::new([0.0, 0.0, 0.0]),
        Some("Na1".to_string()),
    )
    .unwrap();
    let json = serde_json::to_value(&site).unwrap();
    let obj = json.as_object().unwrap();
    assert!(obj.contains_key("species"));
    assert!(obj.contains_key("fractional"));
    assert!(obj.contains_key("label"));
    assert_eq!(obj.len(), 3);
}

#[test]
fn lattice_field_name_is_matrix_only_no_inverse() {
    let l = Lattice::cubic(4.0).unwrap();
    let json = serde_json::to_value(&l).unwrap();
    let obj = json.as_object().unwrap();
    assert_eq!(obj.len(), 1);
    assert!(obj.contains_key("matrix"));
    assert!(!obj.contains_key("inverse"));
}
