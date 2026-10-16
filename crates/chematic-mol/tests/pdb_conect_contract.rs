use chematic_core::AtomIdx;
use chematic_mol::{
    PdbConectError, PdbConectLimits, PdbConectOptions, parse_pdb_conect,
    parse_pdb_conect_with_options, parse_pdb_conect_with_options_and_limits,
};

const TEXT: &str = include_str!("../../../validation/rdkit_issues/pdb-conect/sulfur.pdb");
fn edges(text: &str) -> Vec<(u32, u32, u8)> {
    parse_pdb_conect(text)
        .unwrap_or_else(|e| panic!("{e}"))
        .molecule
        .bonds()
        .map(|(_, b)| (b.atom1.0, b.atom2.0, b.order.order_int()))
        .collect()
}
fn atoms() -> String {
    TEXT.lines()
        .filter(|s| s.starts_with("HETATM"))
        .map(|s| format!("{s}\n"))
        .collect()
}
fn replace_field(start: usize, end: usize, value: &str) -> String {
    let mut source = TEXT.to_owned();
    source.replace_range(start..end, value);
    source
}
#[test]
fn sulfur_has_order_sum_six_without_proximity_bonds() {
    let f = parse_pdb_conect(TEXT).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(f.serials, [10, 20, 30, 40, 50, 90]);
    assert_eq!(f.molecule.atom_count(), 6);
    assert_eq!(edges(TEXT), [(0, 1, 2), (0, 2, 2), (0, 3, 1), (0, 4, 1)]);
    assert_eq!(
        f.molecule
            .atoms()
            .map(|(_, a)| (a.element.atomic_number(), a.charge))
            .collect::<Vec<_>>(),
        [(16, 0), (8, 0), (8, 0), (8, -1), (8, -1), (1, 0)]
    );
    assert_eq!(f.coords.points.len(), 6);
    assert_eq!(
        (
            f.coords.points[5].x,
            f.coords.points[5].y,
            f.coords.points[5].z
        ),
        (0.0, 0.0, 0.8)
    );
    assert_eq!(f.molecule.neighbors(AtomIdx(5)).count(), 0);
}
#[test]
fn no_conect_means_no_edges_and_coordinates_do_not_create_edges() {
    assert!(edges(&atoms()).is_empty());
    let moved = TEXT
        .lines()
        .map(|line| {
            if line.starts_with("HETATM") {
                format!(
                    "{}{:8.3}{:8.3}{:8.3}{}",
                    &line[..30],
                    9.0,
                    9.0,
                    9.0,
                    &line[54..]
                )
            } else {
                line.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    assert_eq!(edges(&moved), edges(TEXT));
}
#[test]
fn one_way_reciprocal_and_continuation_multiplicities() {
    for order in 1..=3 {
        let targets = "   20".repeat(order);
        let a = format!("{}CONECT   10{targets}\n", atoms());
        let b = format!("{a}CONECT   20{}\n", "   10".repeat(order));
        let c = format!("{}{}", atoms(), "CONECT   10   20\n".repeat(order));
        assert_eq!(edges(&a), [(0, 1, order as u8)]);
        assert_eq!(edges(&b), edges(&a));
        assert_eq!(edges(&c), edges(&a));
    }
}
#[test]
fn extended_targets_require_explicit_opt_in() {
    let s = TEXT.replace(
        "CONECT   10   20   20   30   30\nCONECT   10   40   50",
        "CONECT   10   20   20   30   30   40   50",
    );
    assert!(matches!(
        parse_pdb_conect(&s),
        Err(PdbConectError::UnsupportedRecord { .. })
    ));
    let f = parse_pdb_conect_with_options(
        &s,
        &PdbConectOptions {
            allow_extended_targets: true,
        },
    )
    .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(f.molecule.bond_count(), 4);
    let legacy = format!("{}CONECT   10                       90\n", atoms());
    assert!(parse_pdb_conect(&legacy).is_err());
}
#[test]
fn invalid_references_and_multiplicities_are_errors() {
    for record in [
        "CONECT   10   99",
        "CONECT   99",
        "CONECT    0   20",
        "CONECT   10    0",
        "CONECT   10   10",
        "CONECT   10   20   20   20   20",
        "CONECT   10   20\nCONECT   20   10   10",
        "CONECT   10  20",
        "CONECT   10    λ",
    ] {
        assert!(
            parse_pdb_conect(&format!("{}{record}\n", atoms())).is_err(),
            "{record}"
        );
    }
    assert!(matches!(
        parse_pdb_conect(&replace_field(6, 11, "    0")),
        Err(PdbConectError::InvalidRecord { .. })
    ));
    assert!(matches!(
        parse_pdb_conect(&TEXT.replacen("HETATM   20", "HETATM   10", 1)),
        Err(PdbConectError::DuplicateSerial { .. })
    ));
}
#[test]
fn explicit_elements_charges_and_single_location_profile() {
    for (start, end, value) in [
        (76, 78, "  "),
        (76, 78, " X"),
        (76, 78, " D"),
        (78, 80, "+1"),
        (78, 80, "0+"),
        (16, 17, "A"),
        (30, 38, "     NaN"),
    ] {
        assert!(
            parse_pdb_conect(&replace_field(start, end, value)).is_err(),
            "{value}"
        );
    }
    for prefix in ["MODEL        1\n", "ENDMDL\n", "END   \n"] {
        assert!(parse_pdb_conect(&format!("{prefix}{TEXT}")).is_err());
    }
    assert!(parse_pdb_conect(&format!("{TEXT}CONECT   10   20\n")).is_err());
    for symbol in ["Cl", "CL"] {
        let f = parse_pdb_conect(&replace_field(76, 78, symbol)).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(f.molecule.atom(AtomIdx(0)).element.atomic_number(), 17);
    }
    let f = parse_pdb_conect(&replace_field(78, 80, "2+")).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(f.molecule.atom(AtomIdx(0)).charge, 2);
}
#[test]
fn bounded_input_atoms_lines_bonds_and_targets() {
    let defaults = PdbConectLimits::default();
    let mut cases = Vec::new();
    let mut l = defaults;
    l.atoms.max_input_bytes = TEXT.len() - 1;
    cases.push(l);
    let mut l = defaults;
    l.atoms.max_line_bytes = 79;
    cases.push(l);
    let mut l = defaults;
    l.atoms.max_atoms = 5;
    cases.push(l);
    let mut l = defaults;
    l.max_bonds = 3;
    cases.push(l);
    let mut l = defaults;
    l.max_targets = 11;
    cases.push(l);
    for limits in cases {
        assert!(
            parse_pdb_conect_with_options_and_limits(TEXT, &PdbConectOptions::default(), &limits)
                .is_err()
        );
    }
    let exact = PdbConectLimits {
        max_bonds: 4,
        max_targets: 12,
        ..defaults
    };
    assert!(
        parse_pdb_conect_with_options_and_limits(TEXT, &PdbConectOptions::default(), &exact)
            .is_ok()
    );
}
#[test]
fn empty_and_crlf_are_deterministic() {
    assert!(edges("").is_empty());
    assert_eq!(edges(TEXT), edges(&TEXT.replace('\n', "\r\n")));
}
