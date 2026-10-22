use chematic_core::Element;
use chematic_mol::{
    PqrAtomRecord, PqrWriteLimits, parse_pqr, write_pqr_checked, write_pqr_checked_with_limits,
};

fn atom() -> PqrAtomRecord {
    PqrAtomRecord {
        group_pdb: "ATOM".into(),
        serial: 17,
        atom_name: "CA".into(),
        res_name: "ALA".into(),
        chain_id: Some("A".into()),
        res_seq: -3,
        icode: Some('B'),
        x: 1.2345678901234567,
        y: -0.0,
        z: 5e-324,
        charge: -0.12345678901234567,
        radius: 1.2345678901234567,
        element: Element::C,
    }
}

fn numeric_bits(a: &PqrAtomRecord) -> [u64; 5] {
    [a.x, a.y, a.z, a.charge, a.radius].map(f64::to_bits)
}

#[test]
fn ordered_fields_survive_checked_rewrite_without_rounding() {
    let mut second = atom();
    second.group_pdb = "HETATM".into();
    second.serial = 29;
    second.atom_name = "ZN".into();
    second.res_name = "ZN".into();
    second.element = Element::ZN;
    second.chain_id = None;
    second.icode = None;
    let input = vec![atom(), second];
    let text = write_pqr_checked(&input).unwrap();
    let actual = parse_pqr(&text).unwrap();
    for (before, after) in input.iter().zip(&actual.atoms) {
        assert_eq!(before, after);
        assert_eq!(numeric_bits(before), numeric_bits(after));
    }
    let legacy = parse_pqr(&chematic_mol::write_pqr(&input)).unwrap();
    assert_ne!(numeric_bits(&input[0]), numeric_bits(&legacy.atoms[0]));
}

#[test]
fn finite_binary64_patterns_round_trip_in_all_five_fields() {
    let mut seed = 0x243f6a8885a308d3u64;
    let mut sample = atom();
    for _ in 0..4096 {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let value = f64::from_bits(seed);
        if !value.is_finite() {
            continue;
        }
        sample.x = value;
        sample.y = -value;
        sample.z = value;
        sample.charge = value;
        sample.radius = value.abs();
        let text = write_pqr_checked(std::slice::from_ref(&sample)).unwrap();
        assert_eq!(
            numeric_bits(&sample),
            numeric_bits(&parse_pqr(&text).unwrap().atoms[0])
        );
    }
    for value in [
        0.,
        -0.,
        f64::from_bits(1),
        -f64::from_bits(1),
        f64::MIN_POSITIVE,
        f64::MAX,
        -f64::MAX,
    ] {
        sample.x = value;
        sample.y = value;
        sample.z = value;
        sample.charge = value;
        sample.radius = value.abs();
        assert_eq!(
            numeric_bits(&sample),
            numeric_bits(
                &parse_pqr(&write_pqr_checked(std::slice::from_ref(&sample)).unwrap())
                    .unwrap()
                    .atoms[0]
            )
        );
    }
    sample.radius = -0.;
    assert_eq!(
        parse_pqr(&write_pqr_checked(&[sample]).unwrap())
            .unwrap()
            .atoms[0]
            .radius
            .to_bits(),
        (-0.0f64).to_bits()
    );
}

#[test]
fn incompatible_metadata_and_nonfinite_fields_return_no_text() {
    let mut bad = Vec::new();
    for value in ["", "C H", "C\nATOM", "é", &"C".repeat(65)] {
        let mut a = atom();
        a.atom_name = value.into();
        bad.push(a);
    }
    for code in ['1', '\n', 'é'] {
        let mut a = atom();
        a.icode = Some(code);
        bad.push(a);
    }
    for value in ["", "A B", "A\rB"] {
        let mut a = atom();
        a.chain_id = Some(value.into());
        bad.push(a);
    }
    let mut a = atom();
    a.element = Element::CA;
    bad.push(a);
    let mut a = atom();
    a.radius = -1.;
    bad.push(a);
    let mut a = atom();
    a.group_pdb = "MODEL".into();
    bad.push(a);
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for index in 0..5 {
            let mut a = atom();
            match index {
                0 => a.x = value,
                1 => a.y = value,
                2 => a.z = value,
                3 => a.charge = value,
                _ => a.radius = value,
            };
            bad.push(a);
        }
    }
    assert!(write_pqr_checked(&[]).is_err());
    for a in bad {
        assert!(write_pqr_checked(&[atom(), a]).is_err());
    }
}

#[test]
fn signed_integer_ids_and_token_boundaries_are_retained() {
    let mut a = atom();
    a.serial = i64::MAX;
    a.res_seq = i64::MIN;
    a.atom_name = format!("C{}", "1".repeat(63));
    a.res_name = "R".repeat(64);
    a.chain_id = Some("A".repeat(64));
    let after = parse_pqr(&write_pqr_checked(std::slice::from_ref(&a)).unwrap()).unwrap();
    assert_eq!(a, after.atoms[0]);
}

#[test]
fn atom_and_output_byte_limits_include_final_newline() {
    let a = atom();
    let text = write_pqr_checked(std::slice::from_ref(&a)).unwrap();
    let exact = PqrWriteLimits {
        max_atoms: 1,
        max_output_bytes: text.len(),
    };
    assert_eq!(
        write_pqr_checked_with_limits(std::slice::from_ref(&a), &exact).unwrap(),
        text
    );
    assert!(
        write_pqr_checked_with_limits(
            std::slice::from_ref(&a),
            &PqrWriteLimits {
                max_output_bytes: text.len() - 1,
                ..exact
            }
        )
        .is_err()
    );
    assert!(write_pqr_checked_with_limits(&[a.clone(), a.clone()], &exact).is_err());
    let atoms = vec![a; 100_000];
    assert!(write_pqr_checked(&atoms).is_ok());
    let mut too_many = atoms;
    too_many.push(atom());
    assert!(write_pqr_checked(&too_many).is_err());
}
