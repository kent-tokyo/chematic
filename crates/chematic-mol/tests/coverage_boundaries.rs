use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};
use chematic_mol::{
    Mol2Error,
    mol2_source::{parse_mol2_source, write_mol2_source},
    parse_cjson_record, parse_cml_strict, write_cml,
};

const MOL2: &str = "@<TRIPOS>MOLECULE\nfixture\n2 1\nSMALL\nNO_CHARGES\n\
@<TRIPOS>ATOM\n1 C1 0 0 0 C.3 1 LIG 0\n2 O1 1 0 0 O.2 1 LIG 0\n\
@<TRIPOS>BOND\n1 1 2 1\n";

#[test]
fn malformed_mol2_fields_are_rejected_at_their_boundary() {
    let cases = [
        MOL2.replace("1 C1 0 0 0 C.3 1 LIG 0", "1 C1 0 0 0"),
        MOL2.replace("1 C1", "atom C1"),
        MOL2.replace("1 C1 0 0 0", "1 C1 bad 0 0"),
        MOL2.replace("1 C1 0 0 0", "1 C1 0 bad 0"),
        MOL2.replace("1 C1 0 0 0", "1 C1 0 0 bad"),
        MOL2.replace("C.3 1 LIG", "C.3 residue LIG"),
        MOL2.replace("C.3 1 LIG", ".3 1 LIG"),
        MOL2.replace("C.3 1 LIG", "Xx.3 1 LIG"),
        format!("{MOL2}@<TRIPOS>ATOM\n3 C3 2 0 0 C.3 1 LIG 0\n"),
        MOL2.replace("1 1 2 1", "1 1 2"),
        MOL2.replace("1 1 2 1", "bond 1 2 1"),
        MOL2.replace("1 1 2 1", "1 from 2 1"),
        MOL2.replace("1 1 2 1", "1 1 to 1"),
        MOL2.replace("1 1 2 1", "1 9 2 1"),
    ];
    for input in cases {
        assert!(parse_mol2_source(&input).is_err(), "accepted:\n{input}");
    }
}

#[test]
fn mol2_error_details_remain_typed_and_actionable() {
    let errors = [
        Mol2Error::MissingSection("ATOM".into()),
        Mol2Error::MultipleMolecules { count: 2 },
        Mol2Error::InvalidAtomLine {
            line: 7,
            detail: "bad atom".into(),
        },
        Mol2Error::InvalidBondLine {
            line: 9,
            detail: "bad bond".into(),
        },
        Mol2Error::UnknownElement {
            symbol: "Xx".into(),
            line: 7,
        },
        Mol2Error::ResourceLimit {
            resource: "atoms",
            actual: 2,
            limit: 1,
        },
    ];
    for error in errors {
        assert!(error.to_string().starts_with("MOL2:"));
    }
}

#[test]
fn graph_writer_covers_lossy_bond_order_spellings_explicitly() {
    for (order, token) in [
        (BondOrder::Zero, "nc"),
        (BondOrder::Triple, "3"),
        (BondOrder::Aromatic, "ar"),
        (BondOrder::Quadruple, "4"),
        (BondOrder::QueryAny, "un"),
        (BondOrder::QuerySingleOrDouble, "un"),
        (BondOrder::QuerySingleOrAromatic, "un"),
        (BondOrder::QueryDoubleOrAromatic, "un"),
    ] {
        let mut builder = MoleculeBuilder::new();
        let a = builder.add_atom(Atom::new(Element::C));
        let b = builder.add_atom(Atom::new(Element::C));
        builder.add_bond(a, b, order).unwrap();
        let output = write_mol2_source(&builder.build(), &[]);
        let bond_line = output.lines().last().unwrap();
        assert_eq!(bond_line.split_whitespace().last(), Some(token));
    }
}

#[test]
fn cjson_projection_rejects_container_type_confusion() {
    for source in [
        "[]",
        r#"{"atoms":[],"vendor":[true,"text",null]}"#,
        r#"{"atoms":{"elements":[]}}"#,
        r#"{"atoms":{"elements":{"number":[6]}},"bonds":{"connections":[]}}"#,
        r#"{"atoms":{"elements":{"number":[6]}},"bonds":{"connections":{"index":{}}}}"#,
    ] {
        assert!(parse_cjson_record(source).is_err(), "accepted: {source}");
    }
}

#[test]
fn cjson_unique_value_visitor_accepts_opaque_json_scalar_kinds() {
    let source = r#"{
        "atoms":{"elements":{"number":[6]}},
        "vendor":{"bool":true,"string":"plain","escaped":"\u0061","null":null,"array":[false,1]}
    }"#;
    let record = parse_cjson_record(source).unwrap();
    assert_eq!(record.opaque_fields(), ["/vendor"]);
}

#[test]
fn strict_cml_rejects_empty_and_unbalanced_documents() {
    assert!(parse_cml_strict("<molecule></molecule>").is_err());
    assert!(parse_cml_strict("</molecule>").is_err());
    assert!(parse_cml_strict("<molecule><></molecule>").is_err());
}

#[test]
fn strict_cml_keeps_special_xml_sections_outside_the_molecule_graph() {
    let input = concat!(
        "text<!-- comment -->text<![CDATA[content]]>text<?processing instruction?>",
        "<molecule><atomArray><atom id=\"a1\" elementType=\"C\"/></atomArray></molecule>"
    );
    let (molecule, _) = parse_cml_strict(input).unwrap();
    assert_eq!(molecule.atom_count(), 1);
}

#[test]
fn cml_writer_spells_quadruple_bonds_without_downgrading_them() {
    let mut builder = MoleculeBuilder::new();
    let a = builder.add_atom(Atom::new(Element::C));
    let b = builder.add_atom(Atom::new(Element::C));
    builder.add_bond(a, b, BondOrder::Quadruple).unwrap();
    assert!(write_cml(&builder.build(), None).contains("order=\"4\""));
}
