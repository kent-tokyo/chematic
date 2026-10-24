use chematic_core::{BondOrder, Molecule, implicit_hcount};
use chematic_mol::{
    CdxmlError, CdxmlParseOptions, parse_cdxml, parse_cdxml_all, parse_cdxml_with_options,
    write_cdxml,
};
use serde_json::Value;

fn packet() -> Value {
    serde_json::from_str(include_str!(
        "../../../validation/cdxml_bond_meaning/cases.json"
    ))
    .unwrap()
}

fn check(mol: &Molecule, expected: &Value) {
    let atoms = expected["atoms"].as_array().unwrap();
    assert_eq!(mol.atom_count(), atoms.len());
    let position = |index: chematic_core::AtomIdx| {
        atoms
            .iter()
            .position(|a| {
                a["number"].as_u64().unwrap() == u64::from(mol.atom(index).element.atomic_number())
            })
            .unwrap()
    };
    for (index, atom) in mol.atoms() {
        let e = &atoms[position(index)];
        assert_eq!(i64::from(atom.charge), e["charge"].as_i64().unwrap());
        assert_eq!(
            u64::from(atom.isotope.unwrap_or(0)),
            e["isotope"].as_u64().unwrap()
        );
        assert_eq!(
            u64::from(implicit_hcount(mol, index)),
            e["hydrogens"].as_u64().unwrap()
        );
    }
    assert_eq!(mol.bond_count(), 1);
    let (_, bond) = mol.bonds().next().unwrap();
    let kind = expected["bond"]["kind"].as_str().unwrap();
    let order = match kind {
        "single" => BondOrder::Single,
        "double" => BondOrder::Double,
        "triple" => BondOrder::Triple,
        "quadruple" => BondOrder::Quadruple,
        "dative" => BondOrder::Dative,
        _ => unreachable!(),
    };
    assert_eq!(bond.order, order);
    let actual = [position(bond.atom1), position(bond.atom2)];
    let expected = [
        expected["bond"]["begin"].as_u64().unwrap() as usize,
        expected["bond"]["end"].as_u64().unwrap() as usize,
    ];
    if kind == "dative" {
        assert_eq!(actual, expected);
    } else {
        let mut a = actual;
        let mut e = expected;
        a.sort();
        e.sort();
        assert_eq!(a, e);
    }
}

#[test]
fn supported_classes_direction_and_native_writer_preserve_authored_graphs() {
    for c in packet()["cases"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|c| c["group"] != "refusal")
    {
        let input = c["input"].as_str().unwrap();
        for infer_nondirectional_stereo in [false, true] {
            let (mol, coords) = parse_cdxml_with_options(
                input,
                &CdxmlParseOptions {
                    infer_nondirectional_stereo,
                },
            )
            .unwrap();
            check(&mol, &c["expected"]);
            for (i, (index, _)) in mol.atoms().enumerate() {
                let slot = c["expected"]["atoms"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .position(|a| {
                        a["number"].as_u64().unwrap()
                            == u64::from(mol.atom(index).element.atomic_number())
                    })
                    .unwrap();
                assert_eq!(coords[i], (slot as f64 * 1.5, slot as f64 * 0.25));
            }
            check(
                &parse_cdxml(&write_cdxml(&mol, &coords)).unwrap().0,
                &c["expected"],
            );
        }
        check(&parse_cdxml_all(input).unwrap()[0].0, &c["expected"]);
        let oracle = chematic_smiles::parse(c["expected"]["smiles"].as_str().unwrap()).unwrap();
        check(
            &parse_cdxml(&write_cdxml(&oracle, &[])).unwrap().0,
            &c["expected"],
        );
    }
}

#[test]
fn declared_unsupported_orders_refuse_instead_of_becoming_single() {
    for c in packet()["cases"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|c| c["group"] == "refusal")
    {
        let input = c["input"].as_str().unwrap();
        for infer_nondirectional_stereo in [false, true] {
            assert!(
                matches!(parse_cdxml_with_options(input,&CdxmlParseOptions {infer_nondirectional_stereo}),
                Err(CdxmlError::UnsupportedBondOrder(value)) if value==c["declared_order"].as_str().unwrap())
            );
        }
        assert!(matches!(
            parse_cdxml_all(input),
            Err(CdxmlError::UnsupportedBondOrder(_))
        ));
    }
}
