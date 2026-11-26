use chematic_mol::{parse_orca_input, write_orca_input};

#[test]
fn blank_lines_inside_xyz_and_empty_blocks_round_trip() {
    let source = "! HF\n%pal end\n* xyz 0 1\nH 0 0 0\n\nH 0 0 0.74\n*\n";
    let parsed = parse_orca_input(source).unwrap();
    let written = write_orca_input(&parsed);
    assert!(written.contains("%pal end"));
    assert_eq!(parse_orca_input(&written).unwrap(), parsed);
}
