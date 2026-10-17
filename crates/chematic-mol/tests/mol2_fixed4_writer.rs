use chematic_mol::{Mol2Record, parse_mol2_record, write_mol2_record};

fn original_row(record: &Mol2Record) -> String {
    let source = &record.atoms[0];
    let atom_id = source.atom_id;
    let atom_name = &source.atom_name;
    let atom_type = &source.atom_type;
    let (x, y, z) = record.coords[0];
    let subst_id = source.subst_id.unwrap_or(1);
    let subst_name = source.subst_name.as_deref().unwrap_or("LIG");
    let partial_charge = source.partial_charge.unwrap_or(0.0);
    let mut row = format!(
        "{atom_id:>6} {atom_name:<8} {x:>10.4} {y:>10.4} {z:>10.4} {atom_type:<8} {subst_id:>4}  {subst_name:<8} {partial_charge:>10.4}"
    );
    for status in &source.status_bits {
        row.push(' ');
        row.push_str(status);
    }
    row
}

fn check(record: &mut Mol2Record, value: f64) {
    record.coords[0] = (value, -value, value * 0.5);
    record.atoms[0].partial_charge = Some(value);
    let expected = original_row(record);
    let output = write_mol2_record(record);
    let actual = output
        .split("@<TRIPOS>ATOM\n")
        .nth(1)
        .unwrap()
        .lines()
        .next()
        .unwrap();
    assert_eq!(actual, expected, "bits={:016x}", value.to_bits());
}

#[test]
fn atom_rows_match_original_float_formatter_including_ieee_boundaries() {
    let mut record =
        parse_mol2_record(include_str!("../../../benchmarks/fixtures/ethanol.mol2")).unwrap();
    let edges = [
        0.0,
        -0.0,
        f64::MIN_POSITIVE,
        f64::from_bits(1),
        f64::MAX,
        f64::MIN,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
        -0.000049999,
        -0.00005,
        1e11,
        1e11 - 1.0,
        99999.99995,
        1.23445,
    ];
    for value in edges {
        check(&mut record, value);
    }
    for integer in -4096..4096 {
        let tie = (f64::from(integer) + 0.5) / 10000.0;
        for bits in [
            tie.to_bits().wrapping_sub(1),
            tie.to_bits(),
            tie.to_bits().wrapping_add(1),
        ] {
            check(&mut record, f64::from_bits(bits));
        }
    }
    let mut state = 0x71a5_92fe_61ad_be23_u64;
    for _ in 0..8192 {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        check(&mut record, f64::from_bits(state));
        check(&mut record, (state as i64 as f64) / 123_456_789_000.0);
    }
}

#[test]
fn split_fields_keep_unicode_padding_oversize_names_and_optional_metadata() {
    let mut record =
        parse_mol2_record(include_str!("../../../benchmarks/fixtures/ethanol.mol2")).unwrap();
    for name in ["名α", "ordinary", "very-long-atom-name"] {
        record.atoms[0].atom_id = u32::MAX;
        record.atoms[0].atom_name = name.into();
        record.atoms[0].atom_type = "long-atom-type".into();
        record.atoms[0].subst_id = Some(i32::MIN);
        record.atoms[0].subst_name = Some("残基α-long".into());
        record.atoms[0].status_bits = vec!["first".into(), "長い".into()];
        check(&mut record, -12345.67895);
        record.atoms[0].subst_id = None;
        record.atoms[0].subst_name = None;
        record.atoms[0].partial_charge = None;
        let expected = original_row(&record);
        let output = write_mol2_record(&record);
        assert_eq!(
            output
                .split("@<TRIPOS>ATOM\n")
                .nth(1)
                .unwrap()
                .lines()
                .next()
                .unwrap(),
            expected
        );
    }
}
