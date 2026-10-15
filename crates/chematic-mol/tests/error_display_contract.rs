use chematic_mol::{
    CdjsonError, Dmol3Error, ExyzError, FhiaimsError, GroError, MmcifPdbError, MmcifResidueError,
    MopacError, PdbSelectionError, PqrPdbqtError, RdfileError, SdfSelectionError, TurbomoleError,
    ZMatrixError,
};

#[test]
fn interchange_errors_have_stable_actionable_messages() {
    let cdjson = [
        CdjsonError::Invalid("atom id".into()),
        CdjsonError::Unsupported("shape".into()),
        CdjsonError::ResourceLimit {
            resource: "atoms",
            actual: 9,
            limit: 8,
        },
    ];
    for error in &cdjson {
        assert!(!error.to_string().is_empty());
    }

    let gro = [
        GroError::InputTooLarge { limit: 1 },
        GroError::LineTooLong { line: 2, limit: 3 },
        GroError::TooManyAtoms { limit: 4 },
        GroError::InvalidField {
            line: 5,
            field: "x",
        },
        GroError::MissingLine {
            line: 6,
            field: "box",
        },
        GroError::TrailingData { line: 7 },
        GroError::ElementCount {
            atoms: 8,
            elements: 9,
        },
        GroError::AllocationFailed,
    ];
    for error in &gro {
        assert!(!error.to_string().is_empty());
    }

    let residues = [
        MmcifResidueError::MissingIdentity {
            atom: 1,
            field: "label_asym_id",
        },
        MmcifResidueError::DuplicateAtom {
            atom: 2,
            residue: 3,
        },
        MmcifResidueError::InvalidField {
            atom: 4,
            field: "formal charge",
        },
        MmcifResidueError::IncompleteRow { atom: 5 },
    ];
    for error in &residues {
        assert!(!error.to_string().is_empty());
    }

    let pdb = [
        MmcifPdbError::EmptyInput,
        MmcifPdbError::InvalidDocument,
        MmcifPdbError::UnsupportedModels,
        MmcifPdbError::MetadataLossRequiresConsent,
        MmcifPdbError::InvalidField {
            atom: 1,
            field: "serial",
        },
        MmcifPdbError::DuplicateSerial { atom: 2, serial: 3 },
        MmcifPdbError::ResidueCollision {
            first: 4,
            second: 5,
        },
        MmcifPdbError::AtomCollision {
            first: 6,
            second: 7,
        },
    ];
    for error in &pdb {
        assert!(!error.to_string().is_empty());
    }

    let pdbqt = [
        PqrPdbqtError::EmptyInput,
        PqrPdbqtError::AtomTypeCount { atoms: 2, types: 1 },
        PqrPdbqtError::RadiusLossRequiresConsent,
        PqrPdbqtError::InvalidField {
            atom: 3,
            field: "charge",
        },
        PqrPdbqtError::InvalidAtomType {
            atom: 4,
            atom_type: "XX".into(),
        },
        PqrPdbqtError::ElementMismatch {
            atom: 5,
            atom_type: "Cl".into(),
        },
    ];
    for error in &pdbqt {
        assert!(!error.to_string().is_empty());
    }

    assert_eq!(
        RdfileError {
            line: 12,
            detail: "bad record".into(),
        }
        .to_string(),
        "RDFile line 12: bad record"
    );
}

#[test]
fn geometry_and_selection_errors_keep_source_context() {
    assert_eq!(
        Dmol3Error {
            line: 2,
            detail: "coordinate".into(),
        }
        .to_string(),
        "DMol3 line 2: coordinate"
    );
    assert_eq!(
        ExyzError {
            line: 3,
            detail: "header".into(),
        }
        .to_string(),
        "EXYZ line 3: header"
    );
    assert_eq!(
        FhiaimsError {
            line: 4,
            detail: "lattice".into(),
        }
        .to_string(),
        "FHI-aims line 4: lattice"
    );
    assert_eq!(
        ZMatrixError {
            atom: Some(5),
            detail: "reference".into(),
        }
        .to_string(),
        "Z-matrix Some(5): reference"
    );
    assert_eq!(
        MopacError {
            line: 6,
            detail: "flag".into(),
        }
        .to_string(),
        "MOPAC line 6: flag"
    );
    assert_eq!(
        PdbSelectionError {
            line: 7,
            detail: "serial".into(),
        }
        .to_string(),
        "PDB selection line 7: serial"
    );
    assert_eq!(
        TurbomoleError {
            line: 8,
            detail: "coord".into(),
        }
        .to_string(),
        "TURBOMOLE line 8: coord"
    );

    let selection_errors = [
        SdfSelectionError::Index { index: 9, len: 2 },
        SdfSelectionError::ResourceLimit {
            resource: "records",
            actual: 10,
            limit: 3,
        },
    ];
    assert_eq!(
        selection_errors
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>(),
        [
            "SDF index 9 out of range for 2 records",
            "SDF selection records 10 exceeds 3",
        ]
    );
}
