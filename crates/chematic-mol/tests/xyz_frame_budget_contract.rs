use chematic_mol::{XyzError, XyzParseLimits, parse_xyz_with_limits};

#[test]
fn single_frame_reader_enforces_frame_budget_over_the_whole_input() {
    let input = "1\nfirst\nH 0 0 0\n1\nsecond\nHe 1 0 0\n";
    let mut limits = XyzParseLimits {
        max_frames: 1,
        ..Default::default()
    };
    assert!(matches!(
        parse_xyz_with_limits(input, &limits),
        Err(XyzError::ResourceLimit {
            resource: "frames",
            actual: 2,
            limit: 1
        })
    ));

    limits.max_frames = 2;
    let first = parse_xyz_with_limits(input, &limits).unwrap();
    assert_eq!(first.atoms.len(), 1);
    assert_eq!(first.atoms[0].element.symbol(), "H");
}
