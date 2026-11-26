use chematic_mol::bgf::parse_bgf;

#[test]
fn malformed_metadata_and_fixed_width_connectivity_are_rejected() {
    let valid = include_str!("../../../validation/bgf/carbonyl.bgf");
    for (input, expected) in [
        (
            valid.replace("FORCEFIELD CUSTOM", "FORCEF INVALID"),
            "invalid metadata record",
        ),
        (
            valid.replace("CONECT     1     2", "CONECT     1 2"),
            "connection/order fields",
        ),
    ] {
        let error = parse_bgf(&input).err().expect("invalid BGF must fail");
        assert!(error.to_string().contains(expected), "{error}");
    }
}
