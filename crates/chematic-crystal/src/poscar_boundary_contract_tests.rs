use crate::*;

const POSCAR: &str = "contract\n1\n2 0 0\n0 2 0\n0 0 2\nSi O\n1 1\nDirect\n0 0 0\n0.5 0.5 0.5\n";
#[test]
fn poscar_restart_sections_reject_malformed_or_unsupported_data() {
    assert!(matches!(parse_poscar(""), Err(PoscarError::Empty)));
    for suffix in ["\n0 0 0\n", "Lattice velocities and vectors\n"] {
        assert!(matches!(
            parse_poscar(&format!("{POSCAR}{suffix}")),
            Err(PoscarError::UnsupportedSection { .. })
        ));
    }
    for bad in ["1 2", "x 0 0", "NaN 0 0", "0 inf 0"] {
        let error = parse_poscar(&format!("{POSCAR}\n{bad}\n0 0 0\n")).unwrap_err();
        match error {
            PoscarError::InvalidLine { line, .. } | PoscarError::NonFiniteValue { line, .. } => {
                assert_eq!(line, 12)
            }
            other => panic!("unexpected velocity error: {other}"),
        }
    }
    let velocities = format!("{POSCAR}\n1 2 3\n4 5 6\n");
    assert!(matches!(
        parse_poscar(&format!("{velocities}\na\nb\n")),
        Err(PoscarError::InvalidLine { .. })
    ));
    assert!(matches!(
        parse_poscar(&format!("{velocities}\na\nb\nc\n\nextra\n")),
        Err(PoscarError::UnsupportedSection { .. })
    ));
    let mut doc = parse_poscar(&velocities).unwrap();
    doc.velocities = Some(vec![[0.0; 3]]);
    assert!(matches!(
        write_poscar(&doc),
        Err(PoscarError::Unwritable { .. })
    ));
    doc.velocities = None;
    doc.selective_dynamics = Some(vec![[true; 3]]);
    assert!(matches!(
        write_poscar(&doc),
        Err(PoscarError::Unwritable { .. })
    ));
}
#[test]
fn malformed_poscar_fields_report_the_affected_line() {
    for (line, replacement, diagnostic) in [
        (1, "nope", "scale"),
        (1, "1 2", "scale"),
        (1, "NaN", "finite"),
        (1, "1 0 1", "positive"),
        (2, "1 0", "lattice"),
        (2, "1 x 0", "lattice"),
        (2, "NaN 0 0", "finite"),
        (5, "", "species"),
        (5, "1 1", "VASP"),
        (5, "Qq O", "element"),
        (6, "1 nope", "count"),
        (6, "1", "count"),
        (8, "0 0", "coordinate"),
        (8, "0 nope 0", "coordinate"),
        (8, "0 NaN 0", "finite"),
    ] {
        let mut lines = POSCAR.lines().map(str::to_string).collect::<Vec<_>>();
        lines[line] = replacement.into();
        let error = parse_poscar(&(lines.join("\n") + "\n")).err().unwrap();
        assert!(
            error
                .to_string()
                .to_lowercase()
                .contains(&diagnostic.to_lowercase()),
            "line {line}: {error}"
        );
    }
    for n in 0..9 {
        let text = POSCAR.lines().take(n).collect::<Vec<_>>().join("\n");
        assert!(parse_poscar(&text).is_err(), "truncation {n}");
    }
    let selective = POSCAR.replace("Direct", "Selective dynamics\nDirect");
    assert!(
        parse_poscar(&selective)
            .err()
            .unwrap()
            .to_string()
            .contains("6")
    );
}

#[test]
fn poscar_cartesian_and_direct_scaling_describe_the_same_positions() {
    let direct = parse_poscar(POSCAR).unwrap();
    let cartesian = POSCAR
        .replace("Direct", "Cartesian")
        .replace("0.5 0.5 0.5", "1 1 1");
    let cart = parse_poscar(&cartesian).unwrap();
    assert_eq!(direct.structure, cart.structure);
    let volume = POSCAR.replacen("\n1\n", "\n-8\n", 1);
    assert_eq!(direct.structure, parse_poscar(&volume).unwrap().structure);
    let round = write_poscar(&direct).unwrap();
    assert_eq!(parse_contcar(&round).unwrap().structure, direct.structure);
}
