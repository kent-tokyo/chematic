use chematic_mol::{
    Mol2KekuleError, kekulize_mol2_record, parse_mol2_source_record, write_mol2_source_record,
};

fn cases() -> Vec<serde_json::Value> {
    serde_json::from_str::<serde_json::Value>(include_str!(
        "../../../validation/openbabel_issues/mol2-kekule/cases.json"
    ))
    .unwrap()["cases"]
        .as_array()
        .unwrap()
        .clone()
}

fn assert_token_only(before: &str, after: &str) -> usize {
    let (a, b) = (
        before.split_inclusive('\n').collect::<Vec<_>>(),
        after.split_inclusive('\n').collect::<Vec<_>>(),
    );
    assert_eq!(a.len(), b.len());
    let mut bonds = false;
    let mut changed = 0;
    for (x, y) in a.into_iter().zip(b) {
        if x.trim().starts_with("@<TRIPOS>") {
            bonds = x.trim() == "@<TRIPOS>BOND";
        }
        if bonds && x.split_whitespace().nth(3) == Some("ar") {
            let old = x.split_whitespace().nth(3).unwrap();
            let new = y.split_whitespace().nth(3).unwrap();
            assert!(matches!(new, "1" | "2"));
            let start = old.as_ptr() as usize - x.as_ptr() as usize;
            assert_eq!(&x[..start], &y[..start]);
            assert_eq!(&x[start + 2..], &y[start + 1..]);
            changed += 1;
        } else {
            assert_eq!(x, y);
        }
    }
    changed
}

#[test]
fn neutral_controls_change_only_ar_bond_tokens_and_keep_source_immutable() {
    for c in cases() {
        let source = c["text"].as_str().unwrap();
        let original = parse_mol2_source_record(source).unwrap();
        let output = kekulize_mol2_record(source).unwrap_or_else(|e| panic!("{}: {e}", c["id"]));
        assert_eq!(
            assert_token_only(source, &output),
            c["aromatic_bonds"].as_u64().unwrap() as usize
        );
        assert_eq!(write_mol2_source_record(&original), source);
        let parsed = parse_mol2_source_record(&output).unwrap();
        assert_eq!(original.atoms(), parsed.atoms());
        assert_eq!(original.coordinates(), parsed.coordinates());
        assert_eq!(original.opaque_sections(), parsed.opaque_sections());
        assert_eq!(kekulize_mol2_record(&output).unwrap(), output);
    }
}

#[test]
fn explicit_assignment_has_expected_double_bonds_for_pyridine() {
    let input = cases().into_iter().find(|c| c["id"] == "pyridine").unwrap();
    let output = kekulize_mol2_record(input["text"].as_str().unwrap()).unwrap();
    let m = parse_mol2_source_record(&output).unwrap();
    assert_eq!(m.bonds().iter().filter(|b| b.bond_type == "2").count(), 3);
    assert!(m.atoms().iter().any(|a| a.atom_type == "N.ar"));
}

#[test]
fn sparse_ids_comments_and_rich_metadata_are_retained_without_retyping() {
    let source = include_str!("../../../validation/openbabel_issues/mol2-kekule/rich-source.mol2");
    let output = kekulize_mol2_record(source).unwrap();
    assert_eq!(assert_token_only(source, &output), 6);
    assert!(output.contains("FLAG_ar"));
    assert!(output.contains("# bond comments stay ar unchanged"));
    assert!(output.contains(" café\r\n"));
}

fn ring(n: usize, element: &str) -> String {
    let mut s = format!("@<TRIPOS>MOLECULE\nring\n{n} {n}\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n");
    for i in 1..=n {
        s.push_str(&format!("{i} {element}{i} 0 0 0 {element}.ar\n"));
    }
    s.push_str("@<TRIPOS>BOND\n");
    for i in 1..=n {
        s.push_str(&format!("{i} {i} {} ar\n", i % n + 1));
    }
    s
}

#[test]
fn failed_matching_returns_an_error_without_partial_text() {
    let source = ring(5, "C");
    assert!(matches!(
        kekulize_mol2_record(&source),
        Err(Mol2KekuleError::Kekulization(_))
    ));
    assert!(source.ends_with("5 5 1 ar\n"));
}

#[test]
fn antiaromatic_and_saturated_donor_cycles_are_not_silently_written() {
    for source in [ring(4, "C"), ring(6, "O")] {
        assert!(matches!(
            kekulize_mol2_record(&source),
            Err(Mol2KekuleError::Unsupported(_))
        ));
    }
}

#[test]
fn charged_and_outside_element_profile_are_refused() {
    let charged = format!("{}@<TRIPOS>UNITY_ATOM_ATTR\n1 1\ncharge 1\n", ring(6, "N"));
    assert!(matches!(
        kekulize_mol2_record(&charged),
        Err(Mol2KekuleError::Unsupported(_))
    ));
    assert!(matches!(
        kekulize_mol2_record(&ring(6, "B")),
        Err(Mol2KekuleError::Unsupported(_))
    ));
}

#[test]
fn acyclic_aromatic_bonds_and_query_zero_bonds_are_refused() {
    let chain = "@<TRIPOS>MOLECULE\nchain\n2 1\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.ar\n2 C2 1 0 0 C.ar\n@<TRIPOS>BOND\n1 1 2 ar\n";
    assert!(matches!(
        kekulize_mol2_record(chain),
        Err(Mol2KekuleError::Unsupported(_))
    ));
    for token in ["un", "du", "nc"] {
        let s = ring(6, "C").replace("1 1 2 ar", &format!("1 1 2 {token}"));
        assert!(matches!(
            kekulize_mol2_record(&s),
            Err(Mol2KekuleError::Unsupported(_))
        ));
    }
}

#[test]
fn implicit_pyrrole_h_is_not_guessed() {
    let s = ring(5, "C").replace("1 C1 0 0 0 C.ar", "1 N1 0 0 0 N.pl3");
    assert!(kekulize_mol2_record(&s).is_err());
}

#[test]
fn graph_and_input_byte_caps_are_enforced() {
    assert!(matches!(
        kekulize_mol2_record(&ring(513, "C")),
        Err(Mol2KekuleError::Parse(_))
    ));
    let large = format!("#{}\n{}", "x".repeat(8 * 1024 * 1024), ring(6, "C"));
    assert!(matches!(
        kekulize_mol2_record(&large),
        Err(Mol2KekuleError::Parse(_))
    ));
    let mut dense =
        String::from("@<TRIPOS>MOLECULE\ndense\n64 1025\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n");
    for i in 1..=64 {
        dense.push_str(&format!("{i} C{i} 0 0 0 C\n"));
    }
    dense.push_str("@<TRIPOS>BOND\n");
    let mut id = 0;
    for a in 1..=64 {
        for b in a + 1..=64 {
            if id < 1025 {
                id += 1;
                dense.push_str(&format!("{id} {a} {b} 1\n"));
            }
        }
    }
    assert!(matches!(
        kekulize_mol2_record(&dense),
        Err(Mol2KekuleError::Parse(_))
    ));
}

#[test]
fn malformed_ids_and_multiple_records_are_rejected() {
    let s = ring(6, "C");
    assert!(kekulize_mol2_record(&s.replace("1 1 2 ar", "1 1 99 ar")).is_err());
    assert!(kekulize_mol2_record(&(s.clone() + &s)).is_err());
}

#[test]
fn lower_case_headers_tabs_crlf_and_no_final_newline_are_preserved() {
    let s = ring(6, "C")
        .replace("@<TRIPOS>", "@<tripos>")
        .replace("BOND", "bond")
        .replace(' ', "\t")
        .replace('\n', "\r\n");
    let s = s.trim_end_matches("\r\n");
    let output = kekulize_mol2_record(s).unwrap();
    assert!(output.contains("@<tripos>bond\r\n"));
    assert!(!output.ends_with('\n'));
    assert_eq!(
        output
            .lines()
            .filter(|l| l.split_whitespace().nth(3) == Some("2"))
            .count(),
        3
    );
}
