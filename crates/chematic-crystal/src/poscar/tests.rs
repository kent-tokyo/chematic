use super::*;

const NACL_POSCAR: &str = "NaCl\n\
1.0\n\
   5.6400000000000000    0.0000000000000000    0.0000000000000000\n\
   0.0000000000000000    5.6400000000000000    0.0000000000000000\n\
   0.0000000000000000    0.0000000000000000    5.6400000000000000\n\
Na Cl\n\
1 1\n\
Direct\n\
  0.0 0.0 0.0\n\
  0.5 0.5 0.5\n";

#[test]
fn parses_basic_vasp5_direct() {
    let doc = parse_poscar(NACL_POSCAR).unwrap();
    assert_eq!(doc.comment, "NaCl");
    assert_eq!(doc.species_order, vec![Element::NA, Element::CL]);
    assert_eq!(doc.structure.site_count(), 2);
    assert_eq!(doc.structure.sites()[0].species[0].element, Element::NA);
    assert_eq!(doc.structure.sites()[1].species[0].element, Element::CL);
    assert_eq!(doc.structure.sites()[0].fractional.0, [0.0, 0.0, 0.0]);
    assert_eq!(doc.structure.sites()[1].fractional.0, [0.5, 0.5, 0.5]);
    assert!(doc.selective_dynamics.is_none());
    assert!(doc.velocities.is_none());
    assert!(doc.predictor_corrector.is_none());
    // No per-atom labels fabricated (see module docs).
    assert!(doc.structure.sites()[0].label.is_none());
}

#[test]
fn empty_comment_line_is_preserved() {
    let text = NACL_POSCAR.replacen("NaCl\n", "\n", 1);
    let doc = parse_poscar(&text).unwrap();
    assert_eq!(doc.comment, "");
    // The rest of the header must still line up correctly.
    assert_eq!(doc.structure.site_count(), 2);
    assert_eq!(doc.structure.lattice().lengths()[0], 5.64);
}

#[test]
fn parses_cartesian_mode() {
    let text = "cart\n\
1.0\n\
4.0 0.0 0.0\n\
0.0 4.0 0.0\n\
0.0 0.0 4.0\n\
C\n\
1\n\
Cartesian\n\
2.0 2.0 2.0\n";
    let doc = parse_poscar(text).unwrap();
    let f = doc.structure.sites()[0].fractional.0;
    assert!((f[0] - 0.5).abs() < 1e-12);
    assert!((f[1] - 0.5).abs() < 1e-12);
    assert!((f[2] - 0.5).abs() < 1e-12);
}

#[test]
fn cartesian_mode_case_insensitive_first_letter_only() {
    for mode in ["c", "C", "k", "K", "Kartesian", "cartesian"] {
        let text =
            format!("t\n1.0\n4.0 0.0 0.0\n0.0 4.0 0.0\n0.0 0.0 4.0\nC\n1\n{mode}\n2.0 2.0 2.0\n");
        let doc = parse_poscar(&text).unwrap();
        let f = doc.structure.sites()[0].fractional.0;
        assert!(
            (f[0] - 0.5).abs() < 1e-9,
            "mode line '{mode}' not treated as Cartesian"
        );
    }
    for mode in ["d", "D", "Direct", "direct", "Fractional"] {
        let text = format!(
            "t\n1.0\n4.0 0.0 0.0\n0.0 4.0 0.0\n0.0 0.0 4.0\nC\n1\n{mode}\n0.25 0.25 0.25\n"
        );
        let doc = parse_poscar(&text).unwrap();
        let f = doc.structure.sites()[0].fractional.0;
        assert_eq!(
            f,
            [0.25, 0.25, 0.25],
            "mode line '{mode}' not treated as Direct"
        );
    }
}

#[test]
fn negative_scale_factor_is_target_volume() {
    // Raw cubic cell 2x2x2 = volume 8. Request volume 64 -> linear
    // factor 2 -> cell becomes 4x4x4.
    let text = "vol\n\
-64.0\n\
2.0 0.0 0.0\n\
0.0 2.0 0.0\n\
0.0 0.0 2.0\n\
C\n\
1\n\
Direct\n\
0.0 0.0 0.0\n";
    let doc = parse_poscar(text).unwrap();
    assert!((doc.structure.lattice().volume() - 64.0).abs() < 1e-9);
    for len in doc.structure.lattice().lengths() {
        assert!((len - 4.0).abs() < 1e-9);
    }
}

#[test]
fn negative_scale_factor_with_cartesian_coords_scales_consistently() {
    // Same target-volume cell as above, but the ion position is given
    // in Cartesian mode using the *raw* (unscaled) vector's units --
    // must be scaled by the same derived linear factor as the lattice,
    // not left as a raw-magnitude mirror (a naive `scale * coord` with
    // scale still negative would flip the sign here).
    let text = "vol_cart\n\
-64.0\n\
2.0 0.0 0.0\n\
0.0 2.0 0.0\n\
0.0 0.0 2.0\n\
C\n\
1\n\
Cartesian\n\
1.0 1.0 1.0\n";
    let doc = parse_poscar(text).unwrap();
    let f = doc.structure.sites()[0].fractional.0;
    // Raw Cartesian (1,1,1) with raw cell edge 2 -> fractional (0.5,
    // 0.5, 0.5) regardless of the linear factor applied uniformly to
    // both lattice and coordinate -- and must NOT be negative/mirrored.
    for c in f {
        assert!(
            (c - 0.5).abs() < 1e-9,
            "fractional {f:?} should be (0.5,0.5,0.5)"
        );
    }
}

#[test]
fn three_component_scale_factor_scales_columns_not_rows() {
    // Non-diagonal (triclinic-ish) raw matrix with 3 distinct factors:
    // a column-scale and a row-scale convention diverge here, and
    // produce different angles -- pinning down which one this reader
    // implements (column, per the VASP wiki wording).
    let text = "tri3\n\
2.0 3.0 5.0\n\
1.0 0.0 0.0\n\
0.3 1.0 0.0\n\
0.2 0.1 1.0\n\
Fe\n\
1\n\
Direct\n\
0.0 0.0 0.0\n";
    let doc = parse_poscar(text).unwrap();
    let m = doc.structure.lattice().matrix();
    // Column-scale convention: matrix[row][col] = raw[row][col] * scale[col].
    let expected = [
        [1.0 * 2.0, 0.0 * 3.0, 0.0 * 5.0],
        [0.3 * 2.0, 1.0 * 3.0, 0.0 * 5.0],
        [0.2 * 2.0, 0.1 * 3.0, 1.0 * 5.0],
    ];
    for r in 0..3 {
        for c in 0..3 {
            assert!(
                (m[r][c] - expected[r][c]).abs() < 1e-9,
                "m[{r}][{c}] = {}, expected {}",
                m[r][c],
                expected[r][c]
            );
        }
    }
    // Cross-check via lengths and angles too. The entry-by-entry matrix
    // assertion above already fully pins down column-vs-row scaling on
    // its own; lengths alone would not (a row-scale convention with
    // these particular factors happens to produce different lengths
    // too, just not by coincidence the *same* ones), so angles are the
    // more targeted secondary check -- row-scaling would change the
    // triangle each pair of (now differently-scaled) row vectors forms
    // to a different degree than column-scaling does.
    let lengths = doc.structure.lattice().lengths();
    let expected_lengths = [
        (expected[0][0].powi(2) + expected[0][1].powi(2) + expected[0][2].powi(2)).sqrt(),
        (expected[1][0].powi(2) + expected[1][1].powi(2) + expected[1][2].powi(2)).sqrt(),
        (expected[2][0].powi(2) + expected[2][1].powi(2) + expected[2][2].powi(2)).sqrt(),
    ];
    for i in 0..3 {
        assert!((lengths[i] - expected_lengths[i]).abs() < 1e-9);
    }
    let angles = doc.structure.lattice().angles_degrees();
    let dot = |u: [f64; 3], v: [f64; 3]| u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    let norm = |u: [f64; 3]| dot(u, u).sqrt();
    let angle_between = |u: [f64; 3], v: [f64; 3]| {
        (dot(u, v) / (norm(u) * norm(v)))
            .clamp(-1.0, 1.0)
            .acos()
            .to_degrees()
    };
    let expected_angles = [
        angle_between(expected[1], expected[2]),
        angle_between(expected[0], expected[2]),
        angle_between(expected[0], expected[1]),
    ];
    for i in 0..3 {
        assert!(
            (angles[i] - expected_angles[i]).abs() < 1e-6,
            "angle[{i}] = {}, expected {}",
            angles[i],
            expected_angles[i]
        );
    }
}

#[test]
fn three_component_scale_rejects_non_positive() {
    let text = "bad3\n\
1.0 -2.0 3.0\n\
1.0 0.0 0.0\n\
0.0 1.0 0.0\n\
0.0 0.0 1.0\n\
C\n\
1\n\
Direct\n\
0.0 0.0 0.0\n";
    let err = parse_poscar(text).unwrap_err();
    assert!(matches!(err, PoscarError::InvalidLine { line: 2, .. }));
}

#[test]
fn selective_dynamics_parsed() {
    let text = "sd\n\
1.0\n\
4.0 0.0 0.0\n\
0.0 4.0 0.0\n\
0.0 0.0 4.0\n\
C O\n\
1 1\n\
Selective dynamics\n\
Direct\n\
0.0 0.0 0.0 T F T\n\
0.5 0.5 0.5 F F F\n";
    let doc = parse_poscar(text).unwrap();
    let sd = doc.selective_dynamics.unwrap();
    assert_eq!(sd, vec![[true, false, true], [false, false, false]]);
}

#[test]
fn velocities_parsed_and_no_mode_line() {
    let text = "v\n\
1.0\n\
4.0 0.0 0.0\n\
0.0 4.0 0.0\n\
0.0 0.0 4.0\n\
C\n\
2\n\
Direct\n\
0.0 0.0 0.0\n\
0.5 0.5 0.5\n\
\n\
0.1 0.2 0.3\n\
-0.1 0.0 0.05\n";
    let doc = parse_poscar(text).unwrap();
    let v = doc.velocities.unwrap();
    assert_eq!(v, vec![[0.1, 0.2, 0.3], [-0.1, 0.0, 0.05]]);
}

#[test]
fn velocities_bad_line_fails_closed() {
    let text = "v\n\
1.0\n\
4.0 0.0 0.0\n\
0.0 4.0 0.0\n\
0.0 0.0 4.0\n\
C\n\
1\n\
Direct\n\
0.0 0.0 0.0\n\
\n\
not a velocity\n";
    let err = parse_poscar(text).unwrap_err();
    assert!(matches!(err, PoscarError::InvalidLine { .. }));
}

#[test]
fn predictor_corrector_round_trips_verbatim() {
    let text = "pc\n\
1.0\n\
4.0 0.0 0.0\n\
0.0 4.0 0.0\n\
0.0 0.0 4.0\n\
C\n\
1\n\
Direct\n\
0.0 0.0 0.0\n\
\n\
0.1 0.2 0.3\n\
\n\
1\n\
0.5\n\
1.0 2.0\n\
0.01 0.02 0.03\n\
0.04 0.05 0.06\n\
0.07 0.08 0.09\n";
    let doc = parse_poscar(text).unwrap();
    let pc = doc.predictor_corrector.clone().unwrap();
    assert_eq!(
        pc.preamble,
        ["1".to_string(), "0.5".to_string(), "1.0 2.0".to_string()]
    );
    assert_eq!(
        pc.lines,
        vec![
            "0.01 0.02 0.03".to_string(),
            "0.04 0.05 0.06".to_string(),
            "0.07 0.08 0.09".to_string(),
        ]
    );

    let written = write_poscar(&doc).unwrap();
    let reparsed = parse_poscar(&written).unwrap();
    assert_eq!(reparsed.predictor_corrector, doc.predictor_corrector);
    assert_eq!(reparsed.velocities, doc.velocities);
}

#[test]
fn lattice_velocities_block_rejected_as_unsupported() {
    let text = "lv\n\
1.0\n\
4.0 0.0 0.0\n\
0.0 4.0 0.0\n\
0.0 0.0 4.0\n\
C\n\
1\n\
Direct\n\
0.0 0.0 0.0\n\
Lattice velocities and vectors\n\
1\n\
0.0 0.0 0.0\n\
0.0 0.0 0.0\n\
0.0 0.0 0.0\n\
4.0 0.0 0.0\n\
0.0 4.0 0.0\n\
0.0 0.0 4.0\n";
    let err = parse_poscar(text).unwrap_err();
    assert!(matches!(err, PoscarError::UnsupportedSection { .. }));
}

#[test]
fn vasp4_style_rejected_not_silently_misparsed() {
    // No species-name line: line 6 goes straight to atom counts.
    let text = "v4\n\
1.0\n\
4.0 0.0 0.0\n\
0.0 4.0 0.0\n\
0.0 0.0 4.0\n\
1 1\n\
Direct\n\
0.0 0.0 0.0\n\
0.5 0.5 0.5\n";
    let err = parse_poscar(text).unwrap_err();
    assert_eq!(err, PoscarError::Vasp4NotSupported);
}

#[test]
fn unknown_element_rejected() {
    let text = "bad\n1.0\n4.0 0.0 0.0\n0.0 4.0 0.0\n0.0 0.0 4.0\nXx\n1\nDirect\n0.0 0.0 0.0\n";
    let err = parse_poscar(text).unwrap_err();
    assert!(matches!(err, PoscarError::UnknownElement { .. }));
}

#[test]
fn atom_count_mismatch_detected() {
    let text = "short\n1.0\n4.0 0.0 0.0\n0.0 4.0 0.0\n0.0 0.0 4.0\nC\n2\nDirect\n0.0 0.0 0.0\n";
    let err = parse_poscar(text).unwrap_err();
    assert_eq!(
        err,
        PoscarError::AtomCountMismatch {
            declared: 2,
            found: 1
        }
    );
}

#[test]
fn atom_count_mismatch_detected_when_trailing_section_follows() {
    // Coordinate block is short by one, but (unlike the fixture above,
    // which simply ends at EOF) a properly blank-line-separated
    // trailing section follows -- must still report a clear
    // AtomCountMismatch, not consume the trailing section's line as
    // bogus coordinate data or surface a confusing "coordinate line
    // needs 3 fields, got ''" error from the blank separator itself.
    let text = "short-trailer\n1.0\n4.0 0.0 0.0\n0.0 4.0 0.0\n0.0 0.0 4.0\nC\n2\nDirect\n0.0 0.0 0.0\n\n0.1 0.2 0.3\n";
    let err = parse_poscar(text).unwrap_err();
    assert_eq!(
        err,
        PoscarError::AtomCountMismatch {
            declared: 2,
            found: 1
        }
    );
}

#[test]
fn contcar_alias_matches_poscar_parse() {
    assert_eq!(
        parse_contcar(NACL_POSCAR).unwrap(),
        parse_poscar(NACL_POSCAR).unwrap()
    );
}

// -- round-trip -----------------------------------------------------

#[test]
fn roundtrip_triclinic_fixture_text() {
    // A genuine triclinic POSCAR (not built via Lattice::from_parameters
    // -- this is literal file text, per the task's round-trip
    // requirement): a, b, c = 5, 6, 7 with distinct, non-90-degree
    // angles, two species, non-alphabetical group order (O before C)
    // to also exercise "preserve file order, don't alphabetize".
    let text = "triclinic fixture\n\
1.0\n\
5.0000000000000000 0.0000000000000000 0.0000000000000000\n\
1.5529142706151024 5.7996130342078438 0.0000000000000000\n\
-1.4790559332561876 1.4998727130474476 6.5850758980592793\n\
O C\n\
2 1\n\
Direct\n\
0.1 0.2 0.3\n\
0.6 0.1 0.9\n\
0.4 0.4 0.4\n";
    let doc = parse_poscar(text).unwrap();
    assert_eq!(doc.species_order, vec![Element::O, Element::C]);

    let written = write_poscar(&doc).unwrap();
    let reparsed = parse_poscar(&written).unwrap();

    assert_eq!(doc.structure, reparsed.structure);
    assert_eq!(doc.comment, reparsed.comment);
    assert_eq!(doc.selective_dynamics, reparsed.selective_dynamics);
    assert_eq!(doc.velocities, reparsed.velocities);
    assert_eq!(doc.predictor_corrector, reparsed.predictor_corrector);
    // Species header order is re-derived from site order (first
    // appearance), which for this fixture matches the original file's
    // own O-then-C group order -- but note write_poscar does not
    // consult doc.species_order to produce it (see docs).
    assert_eq!(reparsed.species_order, vec![Element::O, Element::C]);

    // Lattice vectors match to float precision (write emits exact
    // Display-formatted f64s, which round-trip losslessly).
    let (m0, m1) = (
        doc.structure.lattice().matrix(),
        reparsed.structure.lattice().matrix(),
    );
    for r in 0..3 {
        for c in 0..3 {
            assert!((m0[r][c] - m1[r][c]).abs() < 1e-12);
        }
    }
}

#[test]
fn roundtrip_triclinic_cartesian_mode_fixture_text() {
    // Same triclinic lattice as `roundtrip_triclinic_fixture_text`, but
    // ion positions given in Cartesian mode -- this is the path
    // write_poscar's "always emit Direct" decision actually transforms
    // (cart_to_frac on read, never converted back), so it needs its own
    // coverage rather than relying on the Direct-mode fixture, and a
    // tolerance comparison rather than `==` since cart_to_frac's matrix
    // solve introduces ordinary floating-point rounding.
    let text = "triclinic cartesian fixture\n\
1.0\n\
5.0000000000000000 0.0000000000000000 0.0000000000000000\n\
1.5529142706151024 5.7996130342078438 0.0000000000000000\n\
-1.4790559332561876 1.4998727130474476 6.5850758980592793\n\
O C\n\
2 1\n\
Cartesian\n\
1.1 0.9 1.3\n\
2.4 3.1 4.2\n\
0.5 0.5 0.5\n";
    let doc = parse_poscar(text).unwrap();
    assert_eq!(doc.structure.site_count(), 3);

    let written = write_poscar(&doc).unwrap();
    let reparsed = parse_poscar(&written).unwrap();

    assert_eq!(doc.structure.site_count(), reparsed.structure.site_count());
    let (f0, f1) = (
        doc.structure.fractional_positions(),
        reparsed.structure.fractional_positions(),
    );
    for (a, b) in f0.iter().zip(f1.iter()) {
        for k in 0..3 {
            assert!(
                (a.0[k] - b.0[k]).abs() < 1e-9,
                "fractional coord mismatch: {a:?} vs {b:?}"
            );
        }
    }
    let (c0, c1) = (
        doc.structure.cartesian_positions(),
        reparsed.structure.cartesian_positions(),
    );
    for (a, b) in c0.iter().zip(c1.iter()) {
        for k in 0..3 {
            assert!(
                (a.0[k] - b.0[k]).abs() < 1e-9,
                "cartesian coord mismatch: {a:?} vs {b:?}"
            );
        }
    }
    let (m0, m1) = (
        doc.structure.lattice().matrix(),
        reparsed.structure.lattice().matrix(),
    );
    for r in 0..3 {
        for c in 0..3 {
            assert!((m0[r][c] - m1[r][c]).abs() < 1e-9);
        }
    }
}

#[test]
fn roundtrip_with_selective_dynamics_and_velocities() {
    let text = "sdv\n\
1.0\n\
4.0 0.0 0.0\n\
0.0 4.0 0.0\n\
0.0 0.0 4.0\n\
Fe Ni\n\
1 1\n\
Selective dynamics\n\
Direct\n\
0.0 0.0 0.0 T T F\n\
0.5 0.5 0.5 F F F\n\
\n\
0.01 0.0 0.0\n\
0.0 0.0 0.0\n";
    let doc = parse_poscar(text).unwrap();
    let written = write_poscar(&doc).unwrap();
    let reparsed = parse_poscar(&written).unwrap();
    assert_eq!(doc, reparsed);
}

#[test]
fn write_species_grouped_by_first_appearance_not_alphabetical() {
    let lattice = Lattice::cubic(4.0).unwrap();
    let sites = vec![
        PeriodicSite::new(
            vec![SiteSpecies::full(Element::ZN)],
            FractionalCoord::new([0.0, 0.0, 0.0]),
            None,
        )
        .unwrap(),
        PeriodicSite::new(
            vec![SiteSpecies::full(Element::AL)],
            FractionalCoord::new([0.25, 0.25, 0.25]),
            None,
        )
        .unwrap(),
        PeriodicSite::new(
            vec![SiteSpecies::full(Element::ZN)],
            FractionalCoord::new([0.5, 0.5, 0.5]),
            None,
        )
        .unwrap(),
    ];
    let structure = PeriodicStructure::new(lattice, sites).unwrap();
    let doc = PoscarDocument {
        structure,
        comment: "zn-al-zn".to_string(),
        species_order: vec![],
        selective_dynamics: None,
        velocities: None,
        predictor_corrector: None,
    };
    let written = write_poscar(&doc).unwrap();
    let lines: Vec<&str> = written.lines().collect();
    // Zn appears before Al in site order (first appearance), and
    // alphabetically Al < Zn -- confirms grouping is by first
    // appearance, not alphabetical.
    assert_eq!(lines[5], "Zn Al");
    assert_eq!(lines[6], "2 1");

    let reparsed = parse_poscar(&written).unwrap();
    assert_eq!(reparsed.structure.site_count(), 3);
    // Original site order was Zn, Al, Zn -- after grouping, Zn's two
    // sites (original fractional (0,0,0) and (0.5,0.5,0.5)) must stay
    // in that relative order, now contiguous, followed by Al.
    let frac: Vec<[f64; 3]> = reparsed
        .structure
        .sites()
        .iter()
        .map(|s| s.fractional.0)
        .collect();
    assert_eq!(frac[0], [0.0, 0.0, 0.0]);
    assert_eq!(frac[1], [0.5, 0.5, 0.5]);
    assert_eq!(frac[2], [0.25, 0.25, 0.25]);
}

#[test]
fn write_rejects_disordered_site() {
    let lattice = Lattice::cubic(4.0).unwrap();
    let sites = vec![
        PeriodicSite::new(
            vec![
                SiteSpecies {
                    element: Element::FE,
                    occupancy: Occupancy::new(0.6).unwrap(),
                },
                SiteSpecies {
                    element: Element::NI,
                    occupancy: Occupancy::new(0.4).unwrap(),
                },
            ],
            FractionalCoord::new([0.0, 0.0, 0.0]),
            None,
        )
        .unwrap(),
    ];
    let structure = PeriodicStructure::new(lattice, sites).unwrap();
    let doc = PoscarDocument {
        structure,
        comment: "disordered".to_string(),
        species_order: vec![],
        selective_dynamics: None,
        velocities: None,
        predictor_corrector: None,
    };
    let err = write_poscar(&doc).unwrap_err();
    assert!(matches!(err, PoscarError::Unwritable { .. }));
}

#[test]
fn write_rejects_partial_occupancy() {
    let lattice = Lattice::cubic(4.0).unwrap();
    let sites = vec![
        PeriodicSite::new(
            vec![SiteSpecies {
                element: Element::FE,
                occupancy: Occupancy::new(0.5).unwrap(),
            }],
            FractionalCoord::new([0.0, 0.0, 0.0]),
            None,
        )
        .unwrap(),
    ];
    let structure = PeriodicStructure::new(lattice, sites).unwrap();
    let doc = PoscarDocument {
        structure,
        comment: "vacancy".to_string(),
        species_order: vec![],
        selective_dynamics: None,
        velocities: None,
        predictor_corrector: None,
    };
    let err = write_poscar(&doc).unwrap_err();
    assert!(matches!(err, PoscarError::Unwritable { .. }));
}

#[test]
fn write_rejects_predictor_corrector_when_reorder_needed() {
    // Site order Al, Zn, Al is NOT contiguous per species -- writing
    // requires reordering, which is unsafe with opaque
    // predictor_corrector data present.
    let lattice = Lattice::cubic(4.0).unwrap();
    let sites = vec![
        PeriodicSite::new(
            vec![SiteSpecies::full(Element::AL)],
            FractionalCoord::new([0.0, 0.0, 0.0]),
            None,
        )
        .unwrap(),
        PeriodicSite::new(
            vec![SiteSpecies::full(Element::ZN)],
            FractionalCoord::new([0.25, 0.25, 0.25]),
            None,
        )
        .unwrap(),
        PeriodicSite::new(
            vec![SiteSpecies::full(Element::AL)],
            FractionalCoord::new([0.5, 0.5, 0.5]),
            None,
        )
        .unwrap(),
    ];
    let structure = PeriodicStructure::new(lattice, sites).unwrap();
    let doc = PoscarDocument {
        structure,
        comment: "needs-reorder".to_string(),
        species_order: vec![],
        selective_dynamics: None,
        velocities: Some(vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
        predictor_corrector: Some(PredictorCorrector {
            preamble: ["1".to_string(), "0.5".to_string(), "1.0 2.0".to_string()],
            lines: vec!["0.0 0.0 0.0".to_string()],
        }),
    };
    let err = write_poscar(&doc).unwrap_err();
    assert!(matches!(err, PoscarError::Unwritable { .. }));
}

#[test]
fn write_rejects_predictor_corrector_without_velocities() {
    let lattice = Lattice::cubic(4.0).unwrap();
    let sites = vec![
        PeriodicSite::new(
            vec![SiteSpecies::full(Element::AL)],
            FractionalCoord::new([0.0, 0.0, 0.0]),
            None,
        )
        .unwrap(),
    ];
    let structure = PeriodicStructure::new(lattice, sites).unwrap();
    let doc = PoscarDocument {
        structure,
        comment: "pc-no-vel".to_string(),
        species_order: vec![],
        selective_dynamics: None,
        velocities: None,
        predictor_corrector: Some(PredictorCorrector {
            preamble: ["1".to_string(), "0.5".to_string(), "1.0 2.0".to_string()],
            lines: vec!["0.0 0.0 0.0".to_string()],
        }),
    };
    let err = write_poscar(&doc).unwrap_err();
    assert!(matches!(err, PoscarError::Unwritable { .. }));
}

#[test]
fn empty_input_rejected() {
    assert_eq!(parse_poscar("").unwrap_err(), PoscarError::Empty);
}
