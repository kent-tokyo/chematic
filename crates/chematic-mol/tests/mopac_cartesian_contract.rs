//! Original synthetic MIT OR Apache-2.0 contracts; see the fixture manifest.
use chematic_mol::{
    MopacParseLimits, parse_mopac_cartesian, parse_mopac_cartesian_with_limits,
    write_mopac_cartesian,
};
const TEXT: &str = "PM7 XYZ CHARGE=1 MINI\noriginal synthetic\nflag and label contract\n14(site alpha) 1.25D-1 -2 -1.25 -1 2.5 0\nCl() -1.5 1 1.25d-1 2 0.5 T\nTv 4 0 1 1 0 0\nTv -1 0 5 0 2 0\nTv 1 0 0 0 -6 0\n\n";

#[test]
fn flags_labels_headers_and_vectors_roundtrip() {
    let g = diagnostic(parse_mopac_cartesian(TEXT)).unwrap();
    assert_eq!(g.source_text, TEXT);
    assert_eq!(g.keywords, "PM7 XYZ CHARGE=1 MINI");
    assert_eq!(g.title, "original synthetic");
    assert_eq!(g.comment, "flag and label contract");
    assert_eq!(g.atoms[0].element.symbol(), "Si");
    assert_eq!(g.atoms[0].label.as_deref(), Some("site alpha"));
    assert_eq!(g.atoms[1].label.as_deref(), Some(""));
    assert_eq!(g.atoms[0].row.position_angstrom, [0.125, -1.25, 2.5]);
    for (atom, codes) in g.atoms.iter().zip([["-2", "-1", "0"], ["1", "2", "T"]]) {
        assert_eq!(
            atom.row.optimization_flags.unwrap().map(|f| f.code()),
            codes
        );
    }
    assert_eq!(g.periodicity(), 3);
    assert_eq!(g.translation_vectors[2].position_angstrom, [1., 0., -6.]);
    let back = diagnostic(parse_mopac_cartesian(
        &diagnostic(write_mopac_cartesian(&g)).unwrap(),
    ))
    .unwrap();
    assert_eq!(g.atoms, back.atoms);
    assert_eq!(g.translation_vectors, back.translation_vectors);
}

#[test]
fn mixed_flags_absence_partial_periodicity_and_no_terminator() {
    for n in 0..=3 {
        let mut text = String::from("\n\n\nC 0 0 0\nH 1 +0 2 +1 3 +2\n");
        for row in ["4 1 0", "-1 5 2", "1 0 -6"].iter().take(n) {
            text.push_str(&format!("Tv {row}\n"));
        }
        let g = diagnostic(parse_mopac_cartesian(&text)).unwrap();
        assert_eq!(g.periodicity(), n);
        assert!(g.atoms[0].row.optimization_flags.is_none());
        assert_eq!(
            g.atoms[1].row.optimization_flags.unwrap().map(|f| f.code()),
            ["0", "1", "2"]
        );
        let back = diagnostic(parse_mopac_cartesian(
            &diagnostic(write_mopac_cartesian(&g)).unwrap(),
        ))
        .unwrap();
        assert_eq!(g.atoms, back.atoms);
        assert_eq!(g.translation_vectors, back.translation_vectors);
    }
}

#[test]
fn precision_and_crlf_are_preserved() {
    let text = "\tPM7 XYZ\r\ntitle \r\ncomment\r\nH -0.0 4.9406564584124654e-324 1.7976931348623157e308\r\n";
    let g = diagnostic(parse_mopac_cartesian(text)).unwrap();
    assert_eq!(g.source_text, text);
    assert_eq!(g.title, "title ");
    let back = diagnostic(parse_mopac_cartesian(
        &diagnostic(write_mopac_cartesian(&g)).unwrap(),
    ))
    .unwrap();
    for (a, b) in g.atoms[0]
        .row
        .position_angstrom
        .iter()
        .zip(back.atoms[0].row.position_angstrom)
    {
        assert_eq!(a.to_bits(), b.to_bits());
    }
}

#[test]
fn unsupported_geometry_cannot_be_silently_reinterpreted() {
    for body in [
        "",
        "C 0 1 0 1 0 1 0 0 0",
        "C 0 1 0 1 0 1 0.5",
        "C 0 3 0 1 0 1",
        "C 0 t 0 1 0 1",
        "C 0.0 0.0",
        "c 0 0 0",
        "C14.0 0 0 0",
        "X 0 0 0",
        "XX 0 0 0",
        "0 0 0 0",
        "99 0 0 0",
        "At 0 0 0",
        "Fr 0 0 0",
        "Es 0 0 0",
        "+ 0 0 0",
        "C(unclosed 0 0 0",
        "C(nested(a)) 0 0 0",
        "C(label)0 0 0",
        "C(λ) 0 0 0",
        "C NaN 0 0",
        "C 1e309 0 0",
        "C 0 0 0\n\nC 1 2 3",
        "Tv 1 0 0\nC 0 0 0",
        "C 0 0 0\nTv(label) 1 0 0",
    ] {
        assert!(
            diagnostic(parse_mopac_cartesian(&format!(
                "PM7 XYZ\noriginal\n\n{body}"
            )))
            .is_err(),
            "{body}"
        );
    }
    for keywords in ["PM7 A0", "a0=1", "++", "PM7 &"] {
        assert_eq!(
            diagnostic(parse_mopac_cartesian(&format!("{keywords}\n\n\nC 0 0 0")))
                .unwrap_err()
                .line,
            1
        );
    }
    assert!(diagnostic(parse_mopac_cartesian("PM7\ntitle")).is_err());
    assert!(
        diagnostic(parse_mopac_cartesian(&format!(
            "\n\n\nC({}) 0 0 0",
            "a".repeat(39)
        )))
        .is_err()
    );
    assert!(
        diagnostic(parse_mopac_cartesian(&format!(
            "\n\n\nC({}) 0 0 0",
            "a".repeat(38)
        )))
        .is_ok()
    );
}

#[test]
fn degenerate_translation_vectors_fail() {
    for rows in [
        "Tv 0 0 0",
        "Tv 1 0 0\nTv 2 0 0",
        "Tv 1 0 0\nTv 1 1e-12 0",
        "Tv 1 0 0\nTv 0 1 0\nTv 1 1 1e-12",
        "Tv 1 0 0\nTv 0 1 0\nTv 0 0 1\nTv 1 1 1",
        "Tv 1.79e308 1.79e308 0",
    ] {
        assert!(
            diagnostic(parse_mopac_cartesian(&format!("\n\n\nC 0 0 0\n{rows}"))).is_err(),
            "{rows}"
        );
    }
}

#[test]
fn limits_and_line_numbers() {
    for limits in [
        MopacParseLimits {
            max_atoms: 1,
            ..Default::default()
        },
        MopacParseLimits {
            max_input_bytes: 5,
            ..Default::default()
        },
        MopacParseLimits {
            max_line_len: 8,
            ..Default::default()
        },
    ] {
        assert!(diagnostic(parse_mopac_cartesian_with_limits(TEXT, &limits)).is_err());
    }
    assert_eq!(
        diagnostic(parse_mopac_cartesian("\n\n\nC 0 0 0\nH bad 0 0"))
            .unwrap_err()
            .line,
        5
    );
}

#[test]
fn writer_validates_current_data_instead_of_replaying_source() {
    let mut g = diagnostic(parse_mopac_cartesian(TEXT)).unwrap();
    g.atoms[0].row.position_angstrom[0] = 0.12345678901234566;
    g.atoms[0].row.optimization_flags = None;
    g.atoms[0].label = None;
    g.keywords = "PM6 XYZ".into();
    let back = diagnostic(parse_mopac_cartesian(
        &diagnostic(write_mopac_cartesian(&g)).unwrap(),
    ))
    .unwrap();
    assert_eq!(back.atoms, g.atoms);
    assert_eq!(back.keywords, g.keywords);
    g.title.push('\n');
    assert!(diagnostic(write_mopac_cartesian(&g)).is_err());
    g.title.clear();
    g.atoms[0].row.position_angstrom[1] = f64::NAN;
    assert!(diagnostic(write_mopac_cartesian(&g)).is_err());
}

#[path = "support/upstream_diagnostics.rs"]
mod upstream_diagnostics;
use upstream_diagnostics::diagnostic;
