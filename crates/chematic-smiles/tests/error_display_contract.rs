use chematic_smiles::{
    RdkitHashFunction, RdkitSmilesError, SmilesBatchCanonicalizer, SmilesBatchReaderError,
    SmilesParseLimits, StreamTerminalReason,
};

#[test]
fn rdkit_smiles_errors_distinguish_scope_from_invalid_chemistry() {
    assert_eq!(
        RdkitSmilesError::Unsupported("query bond".into()).to_string(),
        "RDKit-compatible SMILES: unsupported input: query bond"
    );
    assert_eq!(
        RdkitSmilesError::Sanitization("valence".into()).to_string(),
        "RDKit-compatible SMILES: sanitization failed: valence"
    );
}

#[test]
fn batch_terminal_reasons_limits_and_io_errors_are_stable() {
    use std::str::FromStr;

    for (reason, spelling) in [
        (StreamTerminalReason::Cancelled, "cancelled"),
        (StreamTerminalReason::TimeLimit, "time_limit"),
        (StreamTerminalReason::ResourceLimit, "resource_limit"),
        (StreamTerminalReason::ProducerError, "producer_error"),
        (StreamTerminalReason::ConsumerClosed, "consumer_closed"),
    ] {
        assert_eq!(reason.as_str(), spelling);
        assert_eq!(reason.to_string(), spelling);
        assert_eq!(StreamTerminalReason::from_str(spelling), Ok(reason));
    }
    assert!(StreamTerminalReason::from_str("other").is_err());

    let limits = SmilesParseLimits {
        max_atoms: 17,
        ..Default::default()
    };
    assert_eq!(SmilesBatchCanonicalizer::new(limits).limits().max_atoms, 17);

    let error = SmilesBatchReaderError {
        next_unprocessed_index: 4,
        source: std::io::Error::other("reader stopped"),
    };
    assert!(error.to_string().contains("input index 4"));
    assert!(std::error::Error::source(&error).is_some());
}

#[test]
fn rdkit_checked_writers_refuse_enhanced_stereo_without_dropping_it() {
    use chematic_core::{AtomIdx, StereoGroup, StereoGroupKind};

    let mut molecule = chematic_smiles::parse("C[C@H](O)F").unwrap();
    molecule.add_stereo_group(StereoGroup::new(StereoGroupKind::And(1), vec![AtomIdx(1)]));

    assert!(chematic_smiles::rdkit_stereoisomer_smiles(&molecule, 64).is_err());
    assert!(chematic_smiles::rdkit_cx_smiles(&molecule, &Default::default()).is_err());
    assert!(chematic_smiles::rdkit_cx_smarts(&molecule).is_err());
    assert!(
        chematic_smiles::rdkit_mol_hash(&molecule, RdkitHashFunction::CanonicalSmiles, true,)
            .is_err()
    );
    assert!(chematic_smiles::rdkit_stereoisomer_count(&molecule).is_err());
}
