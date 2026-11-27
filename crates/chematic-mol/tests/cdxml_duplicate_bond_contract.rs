use chematic_mol::parse_cdxml;

#[test]
fn duplicate_bonds_are_rejected_even_when_endpoints_are_reversed() {
    for second in ["<b B=\"1\" E=\"2\"/>", "<b B=\"2\" E=\"1\"/>"] {
        let xml = format!(
            "<CDXML><page><fragment><n id=\"1\" Element=\"6\"/><n id=\"2\" Element=\"6\"/><b B=\"1\" E=\"2\"/>{second}</fragment></page></CDXML>"
        );
        assert!(parse_cdxml(&xml).is_err());
    }
}
