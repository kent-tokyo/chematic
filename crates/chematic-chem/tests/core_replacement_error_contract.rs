use chematic_chem::core_replacement::replace_core_checked;
use chematic_core::AtomIdx;
use chematic_smiles::parse;

#[test]
fn core_replacement_explains_invalid_source_mappings() {
    let source = parse("CC").unwrap();
    let core = parse("CC").unwrap();
    for (mapping, expected) in [
        (
            vec![AtomIdx(0), AtomIdx(0)],
            "matched source atoms must be distinct",
        ),
        (
            vec![AtomIdx(0), AtomIdx(2)],
            "matched source atom is out of range",
        ),
    ] {
        let error = replace_core_checked(&source, &core, &mapping)
            .err()
            .expect("invalid mapping");
        assert_eq!(error.to_string(), expected);
    }
}
