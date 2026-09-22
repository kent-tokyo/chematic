use chematic_core::{Atom, AtomIdx, Element, Molecule, MoleculeBuilder, RGroupLabel};
use chematic_mol::{
    MolMetadata, parse_cdxml, parse_mol, parse_mol_v3000, write_cdxml, write_mol, write_mol_v3000,
    write_sdf,
};

fn pseudoatom_fixture() -> Molecule {
    let mut builder = MoleculeBuilder::new();
    builder.add_atom(Atom::new(Element::C));
    builder.add_atom(Atom::wildcard());
    let r = builder.add_atom(Atom::wildcard());
    builder.set_r_group(r, RGroupLabel::unnumbered());
    let r1 = builder.add_atom(Atom::wildcard());
    builder.set_r_group(r1, RGroupLabel::numbered(1).unwrap());
    let r2 = builder.add_atom(Atom::wildcard());
    builder.set_r_group(r2, RGroupLabel::numbered(2).unwrap());
    builder.build()
}

fn assert_fixture(mol: &Molecule) {
    assert_eq!(mol.atom_count(), 5);
    let carbon = mol.atom(AtomIdx(0));
    assert!(!carbon.wildcard);
    assert_eq!(carbon.element, Element::C);

    let wildcard = mol.atom(AtomIdx(1));
    assert!(wildcard.wildcard);
    assert_eq!(mol.r_group_label(AtomIdx(1)), None);

    let r = mol.atom(AtomIdx(2));
    assert!(r.wildcard);
    assert_eq!(mol.r_group_label(AtomIdx(2)).unwrap().number(), None);

    let r1 = mol.atom(AtomIdx(3));
    assert!(r1.wildcard);
    assert_eq!(mol.r_group_label(AtomIdx(3)).unwrap().number(), Some(1));

    let r2 = mol.atom(AtomIdx(4));
    assert!(r2.wildcard);
    assert_eq!(mol.r_group_label(AtomIdx(4)).unwrap().number(), Some(2));
}

fn assert_real_r_elements(mol: &Molecule) {
    let expected = ["Ru", "Rh", "Re", "Rn"];
    assert_eq!(mol.atom_count(), expected.len());
    for (index, symbol) in expected.into_iter().enumerate() {
        let idx = AtomIdx(index as u32);
        let atom = mol.atom(idx);
        assert!(!atom.wildcard, "{symbol} was misclassified as a wildcard");
        assert_eq!(atom.element.symbol(), symbol);
        assert_eq!(mol.r_group_label(idx), None);
    }
}

fn real_r_element_fixture() -> Molecule {
    let mut builder = MoleculeBuilder::new();
    for symbol in ["Ru", "Rh", "Re", "Rn"] {
        builder.add_atom(Atom::new(Element::from_symbol(symbol).unwrap()));
    }
    builder.build()
}

#[test]
fn v2000_preserves_carbon_wildcard_and_r_groups() {
    let text = write_mol(&pseudoatom_fixture(), &MolMetadata::default());
    assert!(text.contains(" * "));
    assert!(text.contains(" R "));
    assert!(text.contains(" R#"));
    assert!(text.contains("M  RGP  2"));
    let (parsed, _) = parse_mol(&text).unwrap();
    assert_fixture(&parsed);
}

#[test]
fn v3000_preserves_carbon_wildcard_and_r_groups() {
    let text = write_mol_v3000(&pseudoatom_fixture(), &MolMetadata::default(), &[]);
    assert!(text.contains(" 2 * "));
    assert!(text.contains(" 3 R "));
    assert!(text.contains(" 4 R# "));
    assert!(text.contains("RGROUPS=(1 1)"));
    assert!(text.contains("RGROUPS=(1 2)"));
    let (parsed, _) = parse_mol_v3000(&text).unwrap();
    assert_fixture(&parsed);
}

#[test]
fn sdf_preserves_carbon_wildcard_and_r_groups() {
    let fixture = pseudoatom_fixture();
    let metadata = MolMetadata::default();
    let text = write_sdf(&[(&fixture, &metadata, &[])]);
    let parsed = chematic_mol::sdf::parse_sdf(&text).unwrap();
    assert_eq!(parsed.len(), 1);
    assert_fixture(&parsed[0].0);
}

#[test]
fn cdxml_preserves_carbon_wildcard_and_r_groups() {
    let text = write_cdxml(&pseudoatom_fixture(), &[]);
    assert!(text.contains("Element=\"6\""));
    assert!(text.contains("GenericNickname=\"*\""));
    assert!(text.contains("GenericNickname=\"R\""));
    assert!(text.contains("GenericNickname=\"R1\""));
    assert!(text.contains("GenericNickname=\"R2\""));
    let (parsed, _) = parse_cdxml(&text).unwrap();
    assert_fixture(&parsed);
}

#[test]
fn unsupported_cdxml_pseudoatom_fails_closed() {
    let input = r#"<CDXML>
<page>
<fragment>
<n id="1" p="0 0" NodeType="GenericNickname" GenericNickname="Foo"/>
</fragment>
</page>
</CDXML>"#;
    let err = match parse_cdxml(input) {
        Ok(_) => panic!("unsupported pseudoatom must fail closed"),
        Err(err) => err,
    };
    assert!(
        err.to_string()
            .contains("unsupported CDXML pseudoatom label")
    );
}

#[test]
fn r_prefixed_elements_remain_real_elements_in_mol_formats() {
    let fixture = real_r_element_fixture();

    let v2000 = write_mol(&fixture, &MolMetadata::default());
    let (parsed_v2000, _) = parse_mol(&v2000).unwrap();
    assert_real_r_elements(&parsed_v2000);

    let v3000 = write_mol_v3000(&fixture, &MolMetadata::default(), &[]);
    let (parsed_v3000, _) = parse_mol_v3000(&v3000).unwrap();
    assert_real_r_elements(&parsed_v3000);
}
