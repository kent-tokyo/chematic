use super::*;

#[test]
fn poscar_errors_keep_physical_lines_and_values() {
    let base = [
        "diagnostics",
        "1.0",
        "4 0 0",
        "0 4 0",
        "0 0 4",
        "C",
        "1",
        "Direct",
        "0 0 0",
    ];
    for (line, replacement, context) in [
        (1, "NaN", "non-finite"),
        (1, "0", "got 0"),
        (5, "Xx", "unknown element symbol 'Xx' at line 6"),
        (6, "2", "declared 2 atom(s), found 1 coordinate line(s)"),
        (8, "bad 0 0", "line 9"),
    ] {
        let mut lines = base;
        lines[line] = replacement;
        let error = parse_poscar(&lines.join("\n")).unwrap_err();
        let message = error.to_string();
        assert!(
            message.contains(context),
            "line {} = {replacement}: {message}",
            line + 1
        );
    }
    assert!(
        parse_poscar("")
            .unwrap_err()
            .to_string()
            .contains("empty POSCAR")
    );
    assert!(
        parse_poscar("comment")
            .unwrap_err()
            .to_string()
            .contains("missing line 2")
    );
    assert_eq!(
        parse_poscar(&base.join("\n"))
            .unwrap()
            .structure
            .sites()
            .len(),
        1
    );
}
