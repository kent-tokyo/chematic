use chematic_core::{AtomIdx, BondIdx, RenumberError, StereoGroupError};

#[test]
fn stereo_group_errors_have_stable_actionable_messages() {
    let cases = [
        (
            StereoGroupError::EmptyGroup { group: 2 },
            "stereo group 2 has no atom or bond members",
        ),
        (
            StereoGroupError::InvalidGroupNumber { group: 3 },
            "stereo group 3 requires a positive AND/OR group number",
        ),
        (
            StereoGroupError::DuplicateGroup {
                first: 1,
                second: 4,
            },
            "stereo groups 1 and 4 repeat the same kind and number",
        ),
        (
            StereoGroupError::InvalidAtom {
                group: 5,
                atom: AtomIdx(7),
            },
            "stereo group 5 references out-of-range atom 7",
        ),
        (
            StereoGroupError::InvalidBond {
                group: 6,
                bond: BondIdx(8),
            },
            "stereo group 6 references out-of-range bond 8",
        ),
        (
            StereoGroupError::AtomOverlap {
                first: 0,
                second: 2,
                atom: AtomIdx(9),
            },
            "atom 9 belongs to stereo groups 0 and 2",
        ),
        (
            StereoGroupError::BondOverlap {
                first: 1,
                second: 3,
                bond: BondIdx(10),
            },
            "bond 10 belongs to stereo groups 1 and 3",
        ),
    ];

    for (error, expected) in cases {
        assert_eq!(error.to_string(), expected);
    }
}

#[test]
fn renumber_errors_identify_the_invalid_contract() {
    let cases = [
        RenumberError::Length {
            expected: 4,
            actual: 3,
        },
        RenumberError::InvalidAtom(AtomIdx(5)),
        RenumberError::DuplicateAtom(AtomIdx(2)),
        RenumberError::InvalidMetadata {
            field: "stereo frame",
            index: 8,
        },
    ];
    for error in &cases {
        assert!(!error.to_string().is_empty());
    }
}
