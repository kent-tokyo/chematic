use crate::orca::*;

#[test]
fn malformed_orca_input_reports_the_failing_construct() {
    for (text, diagnostic) in [
        ("%\n", "malformed"),
        ("%pal\nnprocs 4\n", "end"),
        ("* xyz\n", "header"),
        ("* xyz bad 1\n", "header"),
        ("* xyz 0 bad\n", "header"),
        ("* xyz 0 1 extra\n", "header"),
        ("* unknown 0 1\n", "unknown"),
        ("* xyzfile 0 1\n", "header"),
        ("* int 0 1\nC 0 0 0\n", "closing"),
        ("* xyz 0 1\nC 0 0 0\n", "closing"),
        ("* xyz 0 1\nC 0 0\n*\n", "atom line"),
        ("* xyz 0 1\nQq 0 0 0\n*\n", "element"),
        ("* xyz 0 1\nC bad 0 0\n*\n", "coordinate"),
        ("* xyz 0 1\nC NaN 0 0\n*\n", "finite"),
        ("unexpected\n", "unexpected"),
    ] {
        let error = parse_orca_input(text).unwrap_err();
        assert!(error.to_string().contains(diagnostic), "{text}: {error}");
    }
}

#[test]
fn orca_coordinate_forms_roundtrip_without_reading_external_files() {
    for text in [
        "! HF STO-3G\n* xyzfile -1 2 path with spaces.xyz\n",
        "! HF\n* gzmtfile 0 1 example.zmt\n",
        "! HF\n* int 0 1\nC 0 0 0\nH 1 0 0 1.09\n*\n",
        "! HF\n* xyz 0 1\nC 0$ 1 2 fragment=1\nH 1 2$ 3$\n*\n",
    ] {
        let parsed = parse_orca_input(text).unwrap();
        let encoded = write_orca_input(&parsed);
        let decoded = parse_orca_input(&encoded).unwrap();
        assert_eq!(parsed, decoded, "{encoded}");
    }
}
