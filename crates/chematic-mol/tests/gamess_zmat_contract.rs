use chematic_core::Element;
use chematic_mol::{GamessZmatWriteOptions, write_gamess_zmat};

fn fixture() -> (Vec<Element>, Vec<[f64; 3]>, Vec<Vec<usize>>) {
    (
        vec![Element::C, Element::O, Element::H, Element::N],
        vec![[2., -1., 3.], [3.2, -1., 3.], [2., 0., 3.], [2.2, 0.3, 4.]],
        vec![vec![], vec![0], vec![0, 1], vec![0, 1, 2]],
    )
}

#[test]
fn preserves_frame_and_emits_explicit_internal_definitions() {
    let (e, p, r) = fixture();
    let out = write_gamess_zmat(&e, &p, &r, &GamessZmatWriteOptions::new(-1, 2)).unwrap();
    assert!(out.starts_with(" $CONTRL COORD=UNIQUE UNITS=ANGS\n ICHARG=-1 MULT=2 NZVAR=6\n"));
    assert!(
        out.contains(
            "\nC1\nC 6.0 2.0000000000000000E0 -1.0000000000000000E0 3.0000000000000000E0\n"
        )
    );
    assert!(out.ends_with(
        " $ZMAT IZMAT(1)=\n 1,2,1,\n 1,3,1,\n 2,3,1,2,\n 1,4,1,\n 2,4,1,2,\n 3,4,1,2,3\n $END\n"
    ));
    assert!(out.lines().all(|l| l.len() <= 80));
    assert!(!out.contains("$BASIS"));
}

#[test]
fn one_and_two_atom_cases_have_the_correct_nzvar() {
    let o = GamessZmatWriteOptions::new(0, 1);
    let atom = write_gamess_zmat(&[Element::HE], &[[1., 2., 3.]], &[vec![]], &o).unwrap();
    assert!(atom.contains("NZVAR=0"));
    assert!(!atom.contains("$ZMAT"));
    let pair = write_gamess_zmat(
        &[Element::H; 2],
        &[[0.; 3], [0., 0., 0.74]],
        &[vec![], vec![0]],
        &o,
    )
    .unwrap();
    assert!(pair.contains("NZVAR=1"));
    assert!(pair.contains(" $ZMAT IZMAT(1)=\n 1,2,1\n $END\n"));
}

#[test]
fn current_coordinates_and_references_are_validated() {
    let (e, p, r) = fixture();
    let o = GamessZmatWriteOptions::new(0, 1);
    for bad in [
        vec![vec![], vec![0], vec![0, 1], vec![0, 0, 2]],
        vec![vec![], vec![0], vec![0, 1], vec![0, 1, 4]],
    ] {
        assert!(write_gamess_zmat(&e, &p, &bad, &o).is_err());
    }
    assert!(write_gamess_zmat(&e, &p[..3], &r, &o).is_err());
    let mut bad = p.clone();
    bad[3][2] = f64::NAN;
    assert!(write_gamess_zmat(&e, &bad, &r, &o).is_err());
    let mut bad = p.clone();
    bad[2] = [4., -1., 3.];
    assert!(write_gamess_zmat(&e, &bad, &r, &o).is_err());
}

#[test]
fn rejects_invalid_options_and_header_injection() {
    let (e, p, r) = fixture();
    for title in ["", " ", "$END", "hello!comment", "line\nbreak", "café"] {
        let mut o = GamessZmatWriteOptions::new(0, 1);
        o.title = title.into();
        assert!(write_gamess_zmat(&e, &p, &r, &o).is_err());
    }
    for digits in [0, 18, 255] {
        let mut o = GamessZmatWriteOptions::new(0, 1);
        o.significant_digits = digits;
        assert!(write_gamess_zmat(&e, &p, &r, &o).is_err());
    }
    assert!(write_gamess_zmat(&e, &p, &r, &GamessZmatWriteOptions::new(0, 0)).is_err());
    let mut o = GamessZmatWriteOptions::new(0, 1);
    o.title = "x".repeat(81);
    assert!(write_gamess_zmat(&e, &p, &r, &o).is_err());
}

#[test]
fn rounded_cartesian_coordinates_must_keep_reference_frames_defined() {
    let e = [Element::C; 3];
    let p = [[1., 0., 0.], [1.01, 0., 0.], [1., 1., 0.]];
    let r = [vec![], vec![0], vec![0, 1]];
    let mut o = GamessZmatWriteOptions::new(0, 1);
    o.significant_digits = 1;
    assert!(write_gamess_zmat(&e, &p, &r, &o).is_err());
    o.significant_digits = 17;
    assert!(write_gamess_zmat(&e, &p, &r, &o).is_ok());
}

#[test]
fn card_overflow_and_atom_limits_fail_explicitly() {
    let o = GamessZmatWriteOptions::new(0, 1);
    assert!(write_gamess_zmat(&[], &[], &[], &o).is_err());
    assert!(write_gamess_zmat(&vec![Element::C; 1025], &[], &[], &o).is_err());
    let out = write_gamess_zmat(
        &[Element::OG],
        &[[-1.2345678901234567e-300; 3]],
        &[vec![]],
        &o,
    );
    assert!(out.is_err());
    let mut o = o;
    o.significant_digits = 15;
    assert!(
        write_gamess_zmat(
            &[Element::OG],
            &[[-1.2345678901234567e-300; 3]],
            &[vec![]],
            &o
        )
        .is_ok()
    );
}

#[test]
fn exact_cartesian_f64_round_trip_and_signed_zero() {
    let (e, mut p, r) = fixture();
    p[0] = [-0., -1.2345678901234567, 3.123456789012345];
    let out = write_gamess_zmat(&e, &p, &r, &GamessZmatWriteOptions::new(0, 1)).unwrap();
    let data = out.split("C1\n").nth(1).unwrap();
    for (line, expected) in data.lines().take(4).zip(&p) {
        let got: Vec<f64> = line
            .split_whitespace()
            .skip(2)
            .map(|s| s.parse().unwrap())
            .collect();
        assert_eq!(got.as_slice(), expected);
    }
    assert!(!out.contains("-0.0000000000000000E0"));
}

#[test]
fn direct_xyz_conversion_retains_hydrogen_and_rejects_additional_frames() {
    let xyz = "3\nOriginal water\nO 0 0 0\nH .96 0 0\nH -.24 .93 0\n";
    let r = [vec![], vec![0], vec![0, 1]];
    let o = GamessZmatWriteOptions::new(0, 1);
    let out = chematic_mol::xyz_to_gamess_zmat(xyz, &r, &o).unwrap();
    assert_eq!(out.lines().filter(|l| l.starts_with("H 1.0 ")).count(), 2);
    assert!(out.contains("NZVAR=3"));
    for bad in [
        format!("{xyz}{xyz}"),
        format!("{xyz}garbage\n"),
        String::new(),
        " ".repeat((1 << 20) + 1),
    ] {
        assert!(chematic_mol::xyz_to_gamess_zmat(&bad, &r, &o).is_err());
    }
}

#[test]
fn xyz_adapter_refuses_extra_columns_and_key_value_metadata() {
    let o = GamessZmatWriteOptions::new(0, 1);
    for xyz in ["1\nLattice=cell\nHe 0 0 0\n", "1\nplain\nHe 0 0 0 .5\n"] {
        assert!(chematic_mol::xyz_to_gamess_zmat(xyz, &[vec![]], &o).is_err());
    }
}
