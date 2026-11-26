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

#[test]
fn malformed_atom_keyword_and_excess_connectivity_are_rejected() {
    use chematic_mol::bgf::{BgfLimits, parse_bgf_with_limits};
    let valid = include_str!("../../../validation/bgf/carbonyl.bgf");
    let bad_atom = valid.replacen("\nATOM  ", "\n  ATOM", 1);
    let error = parse_bgf(&bad_atom).err().expect("misaligned atom keyword");
    assert!(
        error.to_string().contains("invalid atom keyword alignment"),
        "{error}"
    );

    let excess = valid.replace(
        "ORDER      2     2\nEND",
        "ORDER      2     2\nCONECT     1     2\nEND",
    );
    let limits = BgfLimits {
        max_atoms: 2,
        ..Default::default()
    };
    let error = parse_bgf_with_limits(&excess, &limits)
        .err()
        .expect("too many connectivity records");
    assert!(
        error.to_string().contains("connectivity count limit"),
        "{error}"
    );
}
