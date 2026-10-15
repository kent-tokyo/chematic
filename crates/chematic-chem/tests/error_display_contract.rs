use chematic_chem::{
    LegacyMmff94ChargeError, MesoStereoGroupError, NativeDistanceMatrixError, RGroupError,
    StereoGroupAlternativesError, StereoGroupInversionError,
};
use chematic_core::{AtomIdx, BondIdx, StereoGroupError};

#[test]
fn distance_matrix_errors_have_stable_messages() {
    let cases = [
        NativeDistanceMatrixError::ResourceLimit {
            resource: "atoms",
            requested: 12,
            limit: 10,
        },
        NativeDistanceMatrixError::DimensionOverflow,
        NativeDistanceMatrixError::Allocation { buffer: "rows" },
    ];
    let messages: Vec<_> = cases.iter().map(ToString::to_string).collect();
    assert_eq!(
        messages,
        [
            "native distance matrix resource limit: 12 atoms exceeds 10",
            "native distance matrix dimension overflow",
            "native distance matrix allocation failed: rows",
        ]
    );
}

#[test]
fn rgroup_errors_keep_the_failed_stage_and_label() {
    let cases = [
        RGroupError::InvalidSmarts("[C".into()),
        RGroupError::Canonicalization("round trip".into()),
        RGroupError::InvalidLabel(0),
        RGroupError::DuplicateLabel(2),
        RGroupError::NonTerminalAttachmentLabel(3),
        RGroupError::MultipleAttachments {
            mol_idx: 4,
            label: 5,
        },
        RGroupError::MatchLimitExceeded {
            mol_idx: 6,
            limit: 100,
        },
    ];
    for error in &cases {
        assert!(!error.to_string().is_empty());
    }
}

#[test]
fn legacy_mmff_charge_error_identifies_the_atom() {
    assert_eq!(
        LegacyMmff94ChargeError {
            atom_index: 7,
            atomic_number: 26,
        }
        .to_string(),
        "legacy MMFF94 BCI charges do not support atom 7 (atomic number 26)"
    );
}

#[test]
fn stereo_group_operation_errors_have_stable_messages() {
    let group_error = StereoGroupError::EmptyGroup { group: 3 };
    let inversion_cases = [
        StereoGroupInversionError::GroupIndex {
            index: 4,
            groups: 2,
        },
        StereoGroupInversionError::AbsoluteGroup,
        StereoGroupInversionError::InvalidGroups(group_error.clone()),
        StereoGroupInversionError::UnsupportedAtom(AtomIdx(5)),
        StereoGroupInversionError::IncompleteFrame(AtomIdx(6)),
        StereoGroupInversionError::UnsupportedBond(BondIdx(7)),
        StereoGroupInversionError::UnspecifiedBond(BondIdx(8)),
        StereoGroupInversionError::SmallRingBond(BondIdx(9)),
        StereoGroupInversionError::ConflictingDirections,
        StereoGroupInversionError::UnrepresentableDirections,
    ];
    for error in &inversion_cases {
        assert!(!error.to_string().is_empty());
    }

    let alternative_cases = [
        StereoGroupAlternativesError::InvalidLimit,
        StereoGroupAlternativesError::AtomLimit {
            actual: 12,
            limit: 10,
        },
        StereoGroupAlternativesError::BondLimit {
            actual: 13,
            limit: 11,
        },
        StereoGroupAlternativesError::TooManyRelativeGroups {
            actual: 7,
            limit: 6,
        },
        StereoGroupAlternativesError::AlternativeLimit {
            required: 16,
            limit: 8,
        },
        StereoGroupAlternativesError::InvalidGroups(group_error.clone()),
        StereoGroupAlternativesError::Inversion {
            group: 2,
            cause: StereoGroupInversionError::AbsoluteGroup,
        },
        StereoGroupAlternativesError::AllocationFailed,
    ];
    for error in &alternative_cases {
        assert!(!error.to_string().is_empty());
    }

    let meso_cases = [
        MesoStereoGroupError::GroupIndex {
            index: 3,
            groups: 1,
        },
        MesoStereoGroupError::AbsoluteGroup,
        MesoStereoGroupError::GraphLimit,
        MesoStereoGroupError::InvalidVisitLimit,
        MesoStereoGroupError::BudgetExceeded {
            visits: 101,
            limit: 100,
        },
        MesoStereoGroupError::Alternatives(StereoGroupAlternativesError::InvalidLimit),
        MesoStereoGroupError::UnsupportedAtom(AtomIdx(1)),
        MesoStereoGroupError::IncompleteFrame(AtomIdx(2)),
        MesoStereoGroupError::UnspecifiedTetrahedralAtom(AtomIdx(3)),
        MesoStereoGroupError::UnsupportedBond(BondIdx(4)),
        MesoStereoGroupError::UnsupportedDirections,
    ];
    for error in &meso_cases {
        assert!(!error.to_string().is_empty());
    }
}
