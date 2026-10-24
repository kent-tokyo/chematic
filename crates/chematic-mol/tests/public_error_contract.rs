use chematic_core::{AtomIdx, BondIdx, SquarePlanarPermutation, StereoGeometry};
use chematic_mol::{
    MolFormat, MolParseError, MolStereoLoss, MolStereoWriteError, SquarePlanarRejectionReason,
    UnsupportedStereoReason,
};

#[test]
fn mol_parse_errors_keep_record_location_and_resource_context() {
    let errors = [
        MolParseError::ResourceLimit {
            resource: "records",
            actual: 4,
            limit: 3,
        },
        MolParseError::InvalidHeader {
            line: 1,
            detail: "bad title".into(),
        },
        MolParseError::InvalidCountLine {
            line: 4,
            detail: "bad count".into(),
        },
        MolParseError::InvalidAtomLine {
            line: 5,
            detail: "bad atom".into(),
        },
        MolParseError::InvalidBondLine {
            line: 6,
            detail: "bad bond".into(),
        },
        MolParseError::InvalidPropertyLine {
            line: 7,
            detail: "bad property".into(),
        },
        MolParseError::UnknownElement {
            symbol: "Xx".into(),
            line: 8,
        },
        MolParseError::UnexpectedEnd,
        MolParseError::V3000ParseError {
            line: 9,
            msg: "bad ctab".into(),
        },
        MolParseError::Io("closed".into()),
    ];
    for error in errors {
        assert!(!error.to_string().is_empty());
        assert!(std::error::Error::source(&error).is_none());
    }
}

#[test]
fn stereo_write_diagnostics_describe_every_public_refusal_reason() {
    let reasons = [
        UnsupportedStereoReason::NoConformerSupplied,
        UnsupportedStereoReason::MissingNeighborOrder,
        UnsupportedStereoReason::ImplicitHydrogenNeighbor,
        UnsupportedStereoReason::MissingConformerPosition,
        UnsupportedStereoReason::WholeMoleculeConformerFlat,
        UnsupportedStereoReason::GeometryRejected(
            SquarePlanarRejectionReason::DegenerateBondVector,
        ),
        UnsupportedStereoReason::GeometryRejected(SquarePlanarRejectionReason::NotCoplanar),
        UnsupportedStereoReason::GeometryRejected(
            SquarePlanarRejectionReason::AmbiguousTransPairing,
        ),
        UnsupportedStereoReason::GeometryTagMismatch {
            computed: SquarePlanarPermutation::SP2,
        },
    ];
    for reason in reasons {
        assert!(!reason.to_string().is_empty());
        let error = MolStereoWriteError {
            atom: AtomIdx(7),
            geometry: StereoGeometry::SquarePlanar,
            format: MolFormat::V3000,
            reason,
        };
        let message = error.to_string();
        assert!(message.contains("atom 7"));
        assert!(message.contains("V3000"));
        assert!(std::error::Error::source(&error).is_none());
    }
    assert_eq!(MolFormat::V2000.to_string(), "V2000");
    assert_eq!(MolFormat::V3000.to_string(), "V3000");
}

#[test]
fn stereo_loss_report_names_each_dropped_channel() {
    assert_eq!(MolStereoLoss::default().to_string(), "no stereo lost");
    let report = MolStereoLoss {
        centres: vec![AtomIdx(1)],
        double_bonds: vec![BondIdx(2)],
        non_tetrahedral_centres: vec![AtomIdx(3)],
        stereo_groups_dropped: true,
    };
    let message = report.to_string();
    for expected in [
        "tetrahedral centres",
        "E/Z written as either",
        "square-planar centres",
        "enhanced stereo groups",
    ] {
        assert!(message.contains(expected), "{message}");
    }
}
