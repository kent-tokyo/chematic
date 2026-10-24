use chematic_core::{AtomIdx, BondIdx, BondOrder};
use chematic_fp::{
    FpsError, PreparedFingerprintError, RdkitMorganError, RdkitRdkFingerprintError,
    RdkitTorsionDiagnosticError,
};
use chematic_perception::AromaticityError;

#[test]
fn public_fingerprint_errors_keep_context_and_sources() {
    let aromaticity = RdkitMorganError::Aromaticity(AromaticityError::KekulizationFailed {
        reason: "odd aromatic cycle".into(),
    });
    assert!(aromaticity.to_string().contains("odd aromatic cycle"));
    assert!(
        RdkitMorganError::UnsupportedBondOrder {
            bond_idx: BondIdx(7),
            order: BondOrder::QueryAny,
        }
        .to_string()
        .contains("BondIdx(7)")
    );
    assert!(
        RdkitMorganError::UnsupportedCoordinationSanitization {
            atom_idx: AtomIdx(3),
            atomic_number: 26,
            degree: 10,
        }
        .to_string()
        .contains("Z=26")
    );
    assert!(
        RdkitMorganError::InternalInvariantViolation {
            reason: "synthetic invariant".into(),
        }
        .to_string()
        .contains("synthetic invariant")
    );

    let database = PreparedFingerprintError::Database {
        index: 4,
        source: RdkitMorganError::InternalInvariantViolation {
            reason: "database".into(),
        },
    };
    assert!(database.to_string().contains("database index 4"));
    let query = PreparedFingerprintError::Query {
        source: RdkitMorganError::InternalInvariantViolation {
            reason: "query".into(),
        },
    };
    assert!(query.to_string().contains("fingerprint query"));

    for error in [
        RdkitTorsionDiagnosticError::AtomIndexOutOfRange {
            index: 5,
            atom_count: 2,
        },
        RdkitTorsionDiagnosticError::EmptyPath,
        RdkitTorsionDiagnosticError::PathTooLong {
            size: 9,
            maximum: 8,
        },
    ] {
        assert!(!error.to_string().is_empty());
    }

    assert!(
        RdkitRdkFingerprintError {
            atom: AtomIdx(8),
            atom_count: 2,
        }
        .to_string()
        .contains("index 8")
    );
}

#[test]
fn fps_error_contract_preserves_io_sources_and_parse_details() {
    let errors = [
        FpsError::MissingNumBits,
        FpsError::InvalidNumBits("zero".into()),
        FpsError::MalformedRecord("deadbeef".into()),
        FpsError::HexLengthMismatch {
            expected_hex_chars: 4,
            got_hex_chars: 2,
            num_bits: 16,
        },
        FpsError::InvalidHexDigit('x'),
    ];
    for error in errors {
        assert!(error.to_string().starts_with("fps:"));
        assert!(std::error::Error::source(&error).is_none());
    }

    let error = FpsError::from(std::io::Error::other("reader stopped"));
    assert!(error.to_string().contains("reader stopped"));
    assert!(std::error::Error::source(&error).is_some());
}
