//! Open Babel #3000: HIN `a` is aromatic, never MDL query type 5.
use chematic_core::{AtomIdx, BondOrder};
use chematic_mol::{
    HinError, HinParseLimits, parse_hin_record, parse_hin_record_with_limits, parse_mol,
    write_hin_record, write_mol,
};

const SOURCE: &str = include_str!("../../../validation/openbabel_issues/hin/aromatic-ring.hin");
const PAIR: &str = "mol 1 \"pair\"\natom 3 N1 N N3 - -0.75 0 0 0 1 9 s\natom 9 C1 C CT - 0.75 1 2 3 1 3 s\nendmol 1\n";

#[test]
fn issue3000_hin_aromatic_bonds_write_as_mdl_type4() {
    let record = parse_hin_record(SOURCE).unwrap();
    let mol = record.to_molecule();
    assert_eq!((mol.atom_count(), mol.bond_count()), (6, 6));
    assert!(mol.atoms().all(|(_, a)| a.aromatic));
    assert!(mol.bonds().all(|(_, b)| b.order == BondOrder::Aromatic));
    let block = write_mol(&mol, &Default::default());
    for row in block.lines().skip(10).take(6) {
        assert_eq!(
            row[6..9].trim(),
            "4",
            "aromatic bond must not become query type 5"
        );
    }
    let (again, _) = parse_mol(&block).unwrap();
    assert_eq!(again.bond_count(), 6);
    assert!(again.bonds().all(|(_, b)| b.order == BondOrder::Aromatic));
}

#[test]
fn source_identity_charge_and_coordinate_precision_are_retained() {
    let record = parse_hin_record(SOURCE).unwrap();
    assert_eq!(record.molecule_id(), 7);
    assert_eq!(record.title(), "aromatic ring — source fields");
    assert_eq!(
        record
            .atoms()
            .iter()
            .map(|a| a.source_id)
            .collect::<Vec<_>>(),
        [1, 2, 3, 4, 5, 6]
    );
    assert_eq!(
        record.atoms()[0].coordinates[0].to_bits(),
        1.2345678901234567_f64.to_bits()
    );
    assert_eq!(record.atoms()[0].partial_charge, -0.75);
    assert!(record.to_molecule().atoms().all(|(_, a)| a.charge == 0));
    assert_eq!(record.atoms()[0].atom_type, "**");
    assert_eq!(write_hin_record(&record), SOURCE);
    let json = record.to_json();
    assert_eq!(json["schema"], "chematic.hin-record.v1");
    assert_eq!(json["atoms"][0]["name"], "C1");
    assert_eq!(json["bonds"][0]["atom1"], 0);
    assert_eq!(json["bonds"][0]["atom2"], 1);
    assert_eq!(json["bonds"][0]["bond_type"], "a");
    assert_eq!(json["graph_projection"]["stereo"], "not_perceived");
}

#[test]
fn both_line_endings_and_missing_terminal_newline_preserve_exact_bytes() {
    for text in [
        SOURCE.to_owned(),
        SOURCE.replace('\n', "\r\n"),
        SOURCE.trim_end().to_owned(),
    ] {
        assert_eq!(write_hin_record(&parse_hin_record(&text).unwrap()), text);
    }
}

#[test]
fn all_supported_bond_types_and_elements_use_literal_source_order() {
    for (token, order) in [
        ("s", BondOrder::Single),
        ("d", BondOrder::Double),
        ("t", BondOrder::Triple),
        ("a", BondOrder::Aromatic),
    ] {
        let input = PAIR.replace(" s", &format!(" {token}"));
        let record = parse_hin_record(&input).unwrap();
        assert_eq!(
            record
                .atoms()
                .iter()
                .map(|a| a.source_id)
                .collect::<Vec<_>>(),
            [3, 9]
        );
        assert_eq!(record.bonds()[0].order, order);
        assert_eq!(
            record.to_molecule().bond(chematic_core::BondIdx(0)).order,
            order
        );
    }
    let m = parse_hin_record(&PAIR.replace("N1 N N3", "Cl1 Cl CL"))
        .unwrap()
        .to_molecule();
    assert_eq!(m.atom(AtomIdx(0)).element, chematic_core::Element::CL);
}

#[test]
fn disconnected_atoms_and_explicit_hydrogen_nodes_are_not_discarded() {
    let text =
        "mol 1 \"ions\"\natom 10 H1 H ** - 0.1 0 0 0 0\natom 20 Br1 Br ** - -0.1 2 0 0 0\nendmol 1";
    let record = parse_hin_record(text).unwrap();
    assert_eq!(record.to_molecule().atom_count(), 2);
    assert_eq!(record.to_molecule().bond_count(), 0);
}

fn rejects(text: &str) {
    let error = parse_hin_record(text).unwrap_err();
    assert!(!error.to_string().is_empty());
}

#[test]
fn malformed_identity_and_record_boundaries_refuse() {
    for text in [
        "".into(),
        PAIR.replace("endmol 1", ""),
        PAIR.replace("endmol 1", "endmol 2"),
        PAIR.replace("mol 1 \"pair\"", "mol 0 \"pair\""),
        PAIR.replace("\"pair\"", "pair"),
        PAIR.replace("atom 9", "atom 3"),
        PAIR.replace("atom 9", "atom 0"),
        format!("{PAIR}{PAIR}"),
        format!("{PAIR}unknown directive\n"),
        "mol 1 \"empty\"\nendmol 1".into(),
        PAIR.replace("mol 1", "atom 1"),
        PAIR.replace("atom 9", "atom 4294967296"),
    ] {
        rejects(&text);
    }
}

#[test]
fn unknown_or_nonfinite_fields_refuse_with_line_context() {
    for (from, to) in [
        ("N1 N N3", "N1 Xx N3"),
        ("-0.75", "NaN"),
        ("-0.75", "inf"),
        ("1 2 3", "1 -inf 3"),
        ("1 2 3", "1 2 not-a-number"),
        ("1 9 s", "1 9 vendor"),
        ("1 9 s", "2 9 s"),
        ("1 9 s", "degree 9 s"),
        ("1 9 s", "1 9 s ignored"),
        ("1 9 s", "1 9"),
    ] {
        let error = parse_hin_record(&PAIR.replace(from, to)).unwrap_err();
        assert!(matches!(error, HinError::InvalidRecord { line: 2 | 3, .. }));
        assert!(error.to_string().contains("line"));
    }
}

#[test]
fn dangling_self_duplicate_and_nonreciprocal_connectivity_refuse() {
    for (from, to) in [
        ("1 9 s", "1 99 s"),
        ("1 9 s", "1 3 s"),
        ("1 9 s", "2 9 s 9 s"),
        ("1 3 s", "0"),
        ("1 3 s", "1 3 d"),
        ("1 3 s", "1 0 s"),
    ] {
        rejects(&PAIR.replace(from, to));
    }
}

#[test]
fn unsupported_hin_directives_are_not_silently_lost() {
    for directive in ["res 1 ALA 1", "forcefield mm+", "sys 1", "vel 1 0 0 0"] {
        rejects(&PAIR.replace("atom 3", &format!("{directive}\natom 3")));
    }
}

#[test]
fn limits_are_exact_inclusive_boundaries() {
    let limits = HinParseLimits {
        max_bytes: PAIR.len(),
        max_atoms: 2,
        max_bonds: 1,
    };
    assert!(parse_hin_record_with_limits(PAIR, limits).is_ok());
    for limits in [
        HinParseLimits {
            max_bytes: PAIR.len() - 1,
            ..limits
        },
        HinParseLimits {
            max_atoms: 1,
            ..limits
        },
        HinParseLimits {
            max_bonds: 0,
            ..limits
        },
    ] {
        assert_eq!(
            parse_hin_record_with_limits(PAIR, limits).unwrap_err(),
            HinError::ResourceLimit
        );
    }
    let over = PAIR.replace("1 9 s", "18446744073709551615 9 s");
    assert!(parse_hin_record(&over).is_err());
    assert_eq!(
        parse_hin_record(&" ".repeat((1 << 20) + 1)).unwrap_err(),
        HinError::ResourceLimit
    );
}

#[test]
fn numeric_spelling_and_signed_zero_do_not_force_source_reformatting() {
    let text = PAIR
        .replace("-0.75", "-7.5e-1")
        .replace("0 0 0 1 9", "-0.0 1e-300 1.7976931348623157e308 1 9");
    let record = parse_hin_record(&text).unwrap();
    assert_eq!(
        record.atoms()[0].coordinates[0].to_bits(),
        (-0.0_f64).to_bits()
    );
    assert_eq!(record.atoms()[0].coordinates[2], f64::MAX);
    assert_eq!(record.atoms()[0].partial_charge, -0.75);
    assert_eq!(write_hin_record(&record), text);
    rejects(&PAIR.replace("-0.75", "1e309"));
}
