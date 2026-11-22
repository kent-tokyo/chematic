use chematic_chem::cns_mpo_from_parts;
use chematic_smiles::parse;

#[test]
fn cns_mpo_polar_surface_area_uses_both_shoulders_and_plateau() {
    let methane = parse("C").unwrap();
    // Neutral methane's LogD equals the supplied LogP; all other terms are optimal.
    for (tpsa, expected) in [
        (-1.0, 5.0),
        (0.0, 5.0),
        (20.0, 5.5),
        (40.0, 6.0),
        (90.0, 6.0),
        (105.0, 5.5),
        (120.0, 5.0),
        (121.0, 5.0),
    ] {
        let score = cns_mpo_from_parts(&methane, 0.0, tpsa, 100.0, 0, 0.0);
        assert!((score - expected).abs() < 1e-10, "TPSA {tpsa}: {score}");
    }
}

#[test]
fn cns_mpo_applies_linear_cutoffs_to_supplied_properties() {
    let methane = parse("C").unwrap();
    let optimal = cns_mpo_from_parts(&methane, 0.0, 60.0, 360.0, 0, 8.0);
    let boundary = cns_mpo_from_parts(&methane, 5.0, 60.0, 500.0, 2, 10.0);
    assert!((optimal - 6.0).abs() < 1e-10);
    assert!((boundary - 1.0).abs() < 1e-10);
}
