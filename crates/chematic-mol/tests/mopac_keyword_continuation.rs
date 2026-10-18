//! Open Babel #1231: '+' adds keyword lines, followed by two description lines.
use chematic_mol::{
    MopacParseLimits, parse_mopac_cartesian, parse_mopac_cartesian_with_limits,
    write_mopac_cartesian,
};

fn deck(keywords: &str) -> String {
    format!(
        "{keywords}\noriginal title\nretained comment\nC(site) 0 0 0 1 0 0\nO 1.4 1 0.2 0 -0.5 1\n\n"
    )
}

#[test]
fn two_keyword_lines_keep_description_and_geometry() {
    let text = deck("PM7 XYZ +\nCHARGE=1 SINGLET");
    let g = parse_mopac_cartesian(&text).unwrap();
    assert_eq!(g.keywords, "PM7 XYZ +\nCHARGE=1 SINGLET");
    assert_eq!(g.title, "original title");
    assert_eq!(g.comment, "retained comment");
    assert_eq!(g.source_text, text);
    assert_eq!(g.atoms.len(), 2);
    assert_eq!(g.atoms[0].label.as_deref(), Some("site"));
    assert_eq!(g.atoms[1].row.position_angstrom, [1.4, 0.2, -0.5]);
    let rewritten = write_mopac_cartesian(&g).unwrap();
    let reread = parse_mopac_cartesian(&rewritten).unwrap();
    assert_eq!(reread.keywords, g.keywords);
    assert_eq!(reread.atoms, g.atoms);
}

#[test]
fn plus_may_be_in_middle_of_each_nonfinal_line() {
    let g = parse_mopac_cartesian(&deck("PM7 + XYZ\nCHARGE=+1 + SINGLET\nPRECISE")).unwrap();
    assert_eq!(g.keywords.lines().count(), 3);
    assert_eq!(g.title, "original title");
    assert_eq!(g.atoms.len(), 2);
}

#[test]
fn charged_keyword_and_quoted_plus_are_not_continuations() {
    for keywords in [
        "PM7 XYZ CHARGE=+1",
        "PM7 TITLE=\"a + b\"",
        "PM7 TITLE=\"++ & A0\"",
    ] {
        let g = parse_mopac_cartesian(&deck(keywords)).unwrap();
        assert_eq!(g.keywords, keywords);
        assert_eq!(g.title, "original title");
    }
}

#[test]
fn crlf_source_and_indentation_are_retained_appropriately() {
    let input = deck("PM7 + \n  XYZ CHARGE=+1\t").replace('\n', "\r\n");
    let g = parse_mopac_cartesian(&input).unwrap();
    assert_eq!(g.keywords, "PM7 + \n  XYZ CHARGE=+1\t");
    assert_eq!(g.source_text, input);
    assert!(!write_mopac_cartesian(&g).unwrap().contains('\r'));
}

#[test]
fn blank_keyword_line_and_blank_descriptions_still_work() {
    let g = parse_mopac_cartesian("\n\n\nH 0 0 0\n").unwrap();
    assert_eq!((&*g.keywords, &*g.title, &*g.comment), ("", "", ""));
    let g = parse_mopac_cartesian("PM7 +\n\n\n\nH 0 0 0\n").unwrap();
    assert_eq!(g.keywords, "PM7 +\n");
    assert_eq!((&*g.title, &*g.comment), ("", ""));
}

#[test]
fn writer_rejects_unannounced_extra_keyword_lines() {
    let mut g = parse_mopac_cartesian(&deck("PM7 XYZ")).unwrap();
    for keywords in ["PM7\nXYZ", "PM7 +", "PM7 +\nXYZ +", "PM7\nXYZ\nPRECISE"] {
        g.keywords = keywords.into();
        assert!(write_mopac_cartesian(&g).is_err(), "{keywords}");
    }
    g.keywords = "PM7 +\nXYZ".into();
    assert!(write_mopac_cartesian(&g).is_ok());
}

#[test]
fn unsupported_variants_are_not_flattened() {
    for keywords in [
        "PM7 +\nA0",
        "PM7 +\n&",
        "PM7 ++\nXYZ",
        "CVB=(++\n1,2)",
        "PM7 TITLE=\"unfinished",
    ] {
        assert!(
            parse_mopac_cartesian(&deck(keywords)).is_err(),
            "{keywords}"
        );
    }
}

#[test]
fn diagnostics_use_physical_source_line_numbers() {
    let text = "PM7 +\nXYZ +\nCHARGE=1\ntitle\ncomment\nC 0 0 nan\n";
    assert_eq!(parse_mopac_cartesian(text).unwrap_err().line, 6);
    assert_eq!(parse_mopac_cartesian("PM7 +\n").unwrap_err().line, 2);
    let mut g = parse_mopac_cartesian(&deck("PM7 +\nXYZ")).unwrap();
    g.comment = "bad\ncomment".into();
    assert_eq!(write_mopac_cartesian(&g).unwrap_err().line, 4);
}

#[test]
fn keyword_line_count_is_bounded() {
    let keywords = "PM7 +\n".repeat(63) + "XYZ";
    assert!(parse_mopac_cartesian(&deck(&keywords)).is_ok());
    let keywords = "PM7 +\n".repeat(64) + "XYZ";
    assert_eq!(
        parse_mopac_cartesian(&deck(&keywords)).unwrap_err().line,
        65
    );
}

#[test]
fn keyword_bytes_and_physical_line_limits_are_bounded() {
    let keywords = format!(
        "{} +\n{} +\n{} +\n{}",
        "X".repeat(16_381),
        "X".repeat(16_381),
        "X".repeat(16_381),
        "X".repeat(16_383)
    );
    assert_eq!(keywords.len(), 65_535);
    assert!(parse_mopac_cartesian(&deck(&keywords)).is_ok());
    let exact_limit = " ".to_owned() + &keywords;
    assert_eq!(exact_limit.len(), 65_536);
    assert!(parse_mopac_cartesian(&deck(&exact_limit)).is_ok());
    let oversized = exact_limit.replacen('\n', "\n ", 1);
    assert_eq!(oversized.len(), 65_537);
    assert!(parse_mopac_cartesian(&deck(&oversized)).is_err());
    let limits = MopacParseLimits {
        max_line_len: 9,
        ..Default::default()
    };
    assert_eq!(
        parse_mopac_cartesian_with_limits("PM7 +\nLONGERTHAN9\n\n\nH 0 0 0\n", &limits)
            .unwrap_err()
            .line,
        2
    );
}

#[test]
fn emit_reference_cases() {
    for count in 1..=12 {
        for position in 0..3 {
            let mut lines = vec!["PM7 XYZ CHARGE=+1".to_owned()];
            for index in 1..count {
                lines.push(format!("PRECISE DUMMY={index}"));
            }
            for line in &mut lines[..count - 1] {
                *line = match position {
                    0 => format!("+ {line}"),
                    1 => format!("{line} +"),
                    _ => format!(
                        "{} + {}",
                        line.split_once(' ').unwrap().0,
                        line.split_once(' ').unwrap().1
                    ),
                };
            }
            let input = deck(&lines.join("\n"));
            let geometry = parse_mopac_cartesian(&input).unwrap();
            let output = write_mopac_cartesian(&geometry).unwrap();
            println!(
                "MOPAC_ORACLE {}",
                serde_json::json!({"keyword_lines": count, "plus_position": position, "input": input, "output": output, "keywords": geometry.keywords, "positions": geometry.atoms.iter().map(|a| a.row.position_angstrom).collect::<Vec<_>>() })
            );
        }
    }
}
