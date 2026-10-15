//! Original minimal analytic geometries; MIT OR Apache-2.0.
use chematic_core::Element;
use chematic_mol::{
    GzmatWriteOptions, ZMatrixRow, cartesian_to_zmatrix, write_gzmat, zmatrix_to_cartesian,
};

fn tetra() -> Vec<[f64; 3]> {
    vec![
        [0.0, 0.0, 0.0],
        [1.234567891234567, 0.0, 0.0],
        [1.1, 1.33333333333333, 0.0],
        [0.1, 0.3, 1.444444444444444],
    ]
}
fn refs(n: usize) -> Vec<Vec<usize>> {
    (0..n).map(|i| (0..i).rev().take(3).collect()).collect()
}
fn rows() -> Vec<ZMatrixRow> {
    diagnostic(cartesian_to_zmatrix(
        &[Element::C, Element::N, Element::O, Element::F],
        &tetra(),
        &refs(4),
    ))
    .unwrap()
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-12, "{a} != {b}");
}
fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    (a[0] - b[0]).hypot(a[1] - b[1]).hypot(a[2] - b[2])
}
fn volume(p: &[[f64; 3]]) -> f64 {
    let v: Vec<_> = p[1..]
        .iter()
        .map(|q| [q[0] - p[0][0], q[1] - p[0][1], q[2] - p[0][2]])
        .collect();
    v[0][0] * (v[1][1] * v[2][2] - v[1][2] * v[2][1])
        - v[0][1] * (v[1][0] * v[2][2] - v[1][2] * v[2][0])
        + v[0][2] * (v[1][0] * v[2][1] - v[1][1] * v[2][0])
}
#[test]
fn analytic_lengths_angles_and_signed_torsions() {
    let p = vec![
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 1.0, 1.0],
    ];
    let r = diagnostic(cartesian_to_zmatrix(&[Element::C; 4], &p, &refs(4))).unwrap();
    close(r[1].values[0], 1.0);
    close(r[2].values[0], 1.0);
    close(r[2].values[1], 90.0);
    close(r[3].values[0], 1.0);
    close(r[3].values[1], 90.0);
    close(r[3].values[2], 90.0);
    let mut mirror = p;
    mirror[3][1] = -1.0;
    close(
        diagnostic(cartesian_to_zmatrix(&[Element::C; 4], &mirror, &refs(4))).unwrap()[3].values[2],
        -90.0,
    );
}
#[test]
fn reconstruction_preserves_distances_and_handedness_with_nonchain_refs() {
    let p = tetra();
    for r in [
        refs(4),
        vec![vec![], vec![0], vec![0, 1], vec![0, 2, 1]],
        vec![vec![], vec![0], vec![1, 0], vec![1, 0, 2]],
    ] {
        let internal = diagnostic(cartesian_to_zmatrix(&[Element::C; 4], &p, &r)).unwrap();
        let out = zmatrix_to_cartesian(&internal).unwrap();
        for i in 0..4 {
            for j in 0..i {
                close(distance(p[i], p[j]), distance(out[i], out[j]));
            }
        }
        assert_eq!(
            volume(&p).is_sign_positive(),
            volume(&out).is_sign_positive()
        );
    }
}
#[test]
fn seventeen_digits_round_trip_every_computed_scalar() {
    let r = rows();
    let text = diagnostic(write_gzmat(&r, &GzmatWriteOptions::new(-1, 2))).unwrap();
    let variables: Vec<_> = text
        .lines()
        .filter_map(|line| {
            line.split_once("= ")
                .map(|(_, v)| v.parse::<f64>().unwrap())
        })
        .collect();
    assert_eq!(
        variables,
        r.iter()
            .flat_map(|row| row.values.iter().copied())
            .collect::<Vec<_>>()
    );
    assert!(text.contains("-1 2\nC\nN 1 r2\nO 2 r3 1 a3\nF 3 r4 2 a4 1 d4"));
}
#[test]
fn requested_significant_digits_affect_each_value() {
    let mut o = GzmatWriteOptions::new(0, 1);
    o.significant_digits = 8;
    let text = diagnostic(write_gzmat(&rows(), &o)).unwrap();
    assert!(text.contains("r2= 1.2345679e0"));
    o.significant_digits = 4;
    assert!(
        diagnostic(write_gzmat(&rows(), &o))
            .unwrap()
            .contains("r2= 1.235e0")
    );
}
#[test]
fn monatomic_diatomic_and_tiny_distances() {
    for length in [1.234567891234567, 1e-200, 1e200] {
        let r = diagnostic(cartesian_to_zmatrix(
            &[Element::H; 2],
            &[[0.0; 3], [0.0, 0.0, length]],
            &refs(2),
        ))
        .unwrap();
        assert_eq!(r[1].values[0], length);
        assert!(diagnostic(write_gzmat(&r, &GzmatWriteOptions::new(0, 1))).is_ok());
    }
    let r = diagnostic(cartesian_to_zmatrix(
        &[Element::HE],
        &[[123.0, 4.0, 5.0]],
        &refs(1),
    ))
    .unwrap();
    assert_eq!(zmatrix_to_cartesian(&r).unwrap(), [[0.0; 3]]);
    assert!(
        !diagnostic(write_gzmat(&r, &GzmatWriteOptions::new(0, 1)))
            .unwrap()
            .contains("Variables:")
    );
}
#[test]
fn references_are_checked_before_access() {
    for r in [
        vec![vec![], vec![1]],
        vec![vec![], vec![2]],
        vec![vec![0], vec![0]],
        vec![vec![], vec![]],
    ] {
        assert!(
            diagnostic(cartesian_to_zmatrix(
                &[Element::C; 2],
                &[[0.0; 3], [1.0; 3]],
                &r
            ))
            .is_err()
        );
    }
    let mut r = refs(4);
    r[3] = vec![2, 1, 1];
    assert!(diagnostic(cartesian_to_zmatrix(&[Element::C; 4], &tetra(), &r)).is_err());
}
#[test]
fn missing_nonfinite_and_overflowing_coordinates_fail() {
    assert!(diagnostic(cartesian_to_zmatrix(&[], &[], &[])).is_err());
    assert!(diagnostic(cartesian_to_zmatrix(&[Element::C], &[], &[vec![]])).is_err());
    for v in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(
            diagnostic(cartesian_to_zmatrix(
                &[Element::C],
                &[[v, 0.0, 0.0]],
                &refs(1)
            ))
            .is_err()
        );
    }
    assert!(
        diagnostic(cartesian_to_zmatrix(
            &[Element::C; 2],
            &[[f64::MAX, 0.0, 0.0], [-f64::MAX, 0.0, 0.0]],
            &refs(2)
        ))
        .is_err()
    );
}
#[test]
fn coincident_and_linear_frames_fail_without_invented_torsions() {
    assert!(
        diagnostic(cartesian_to_zmatrix(
            &[Element::C; 2],
            &[[0.0; 3]; 2],
            &refs(2)
        ))
        .is_err()
    );
    for y in [0.0, 1e-12] {
        let p = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, y, 0.0]];
        assert!(diagnostic(cartesian_to_zmatrix(&[Element::C; 3], &p, &refs(3))).is_err());
    }
    let p = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let r = vec![vec![], vec![0], vec![1, 0], vec![2, 1, 0], vec![0, 1, 3]];
    let e = diagnostic(cartesian_to_zmatrix(&[Element::C; 5], &p, &r)).unwrap_err();
    assert_eq!(e.atom, Some(4));
    assert!(e.detail.contains("plane"));
}
#[test]
fn options_validate_precision_headers_and_spin() {
    for digits in [0, 18, 255] {
        let mut o = GzmatWriteOptions::new(0, 1);
        o.significant_digits = digits;
        assert!(diagnostic(write_gzmat(&rows(), &o)).is_err());
    }
    let mut o = GzmatWriteOptions::new(0, 0);
    assert!(diagnostic(write_gzmat(&rows(), &o)).is_err());
    o.multiplicity = 1;
    for title in ["", "x\n0 1", "x\rfoo", "x\0y"] {
        o.title = title.into();
        assert!(diagnostic(write_gzmat(&rows(), &o)).is_err());
    }
    o.title = "valid".into();
    o.route = "HF".into();
    assert!(diagnostic(write_gzmat(&rows(), &o)).is_err());
    o.route = "# HF/STO-3G".into();
    assert!(diagnostic(write_gzmat(&rows(), &o)).is_ok());
}
#[test]
fn mutable_rows_fail_closed() {
    for values in [vec![], vec![0.0], vec![f64::NAN], vec![1.0, 90.0]] {
        let mut r = rows();
        r[1].values = values;
        assert!(diagnostic(write_gzmat(&r, &GzmatWriteOptions::new(0, 1))).is_err());
    }
    for angle in [0.0, 180.0, 1e-12, 181.0] {
        let mut r = rows();
        r[2].values[1] = angle;
        assert!(zmatrix_to_cartesian(&r).is_err());
    }
    let mut r = rows();
    r[3].values[2] = 181.0;
    assert!(diagnostic(write_gzmat(&r, &GzmatWriteOptions::new(0, 1))).is_err());
}
#[test]
fn rounding_must_not_create_a_linear_angle() {
    let mut r = rows();
    r[2].values[1] = 179.99999;
    let mut o = GzmatWriteOptions::new(0, 1);
    o.significant_digits = 4;
    assert!(diagnostic(write_gzmat(&r, &o)).is_err());
    o.significant_digits = 17;
    assert!(diagnostic(write_gzmat(&r, &o)).is_ok());
}
#[test]
fn rounding_must_not_collapse_an_earlier_reference_plane() {
    let p = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [2.0, 1.0, 0.0],
        [2.0, 0.00001, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let r = vec![
        vec![],
        vec![0],
        vec![1, 0],
        vec![2, 1, 0],
        vec![3, 2, 1],
        vec![0, 1, 4],
    ];
    let rows = diagnostic(cartesian_to_zmatrix(&[Element::C; 6], &p, &r)).unwrap();
    let mut o = GzmatWriteOptions::new(0, 1);
    o.significant_digits = 4;
    assert!(diagnostic(write_gzmat(&rows, &o)).is_err());
    o.significant_digits = 17;
    assert!(diagnostic(write_gzmat(&rows, &o)).is_ok());
}
#[test]
fn large_atom_count_is_rejected_before_other_shapes() {
    let e = vec![Element::C; chematic_mol::gzmat::MAX_ZMATRIX_ATOMS + 1];
    assert!(
        diagnostic(cartesian_to_zmatrix(&e, &[], &[]))
            .unwrap_err()
            .detail
            .contains("count")
    );
}

#[path = "support/upstream_diagnostics.rs"]
mod upstream_diagnostics;
use upstream_diagnostics::diagnostic;
