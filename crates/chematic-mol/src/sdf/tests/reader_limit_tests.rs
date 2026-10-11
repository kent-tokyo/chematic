use super::*;

#[test]
fn single_molecule_without_delimiter_is_valid() {
    let results: Vec<_> = SdfReader::new(MOL_A).collect();
    assert_eq!(results.len(), 1);
    let (mol, _) = results[0].as_ref().expect("parse");
    assert_eq!(mol.atom_count(), 2);
}

#[test]
fn empty_input_yields_no_records() {
    let results: Vec<_> = SdfReader::new("").collect();
    assert!(results.is_empty());
}

#[test]
fn many_empty_records_do_not_recurse() {
    let mut sdf = "$$$$\n".repeat(20_000);
    sdf.push_str(MOL_A);
    let results: Vec<_> = SdfReader::new(&sdf).collect();
    assert_eq!(results.len(), 1);
    assert_eq!(results[0].as_ref().unwrap().0.atom_count(), 2);
}

#[test]
fn explicit_record_limit_is_enforced() {
    let sdf = two_mol_sdf();
    let limits = SdfParseLimits {
        max_records: 1,
        ..SdfParseLimits::default()
    };
    let mut reader = SdfReader::with_limits(&sdf, limits);
    assert!(reader.next().unwrap().is_ok());
    assert!(matches!(
        reader.next(),
        Some(Err(MolParseError::ResourceLimit {
            resource: "records",
            actual: 2,
            limit: 1,
        }))
    ));
    assert!(reader.next().is_none());
}

#[test]
fn in_memory_reader_enforces_input_line_and_record_limits() {
    let input_limits = SdfParseLimits {
        max_input_bytes: 3,
        ..SdfParseLimits::default()
    };
    let mut reader = SdfReader::with_limits(MOL_A, input_limits);
    assert!(matches!(
        reader.next(),
        Some(Err(MolParseError::ResourceLimit {
            resource: "input bytes",
            ..
        }))
    ));
    assert!(reader.next().is_none());

    for input in ["long line\n$$$$\n", "long final line"] {
        let line_limits = SdfParseLimits {
            max_line_bytes: 4,
            ..SdfParseLimits::default()
        };
        let mut reader = SdfReader::with_limits(input, line_limits);
        assert!(matches!(
            reader.next(),
            Some(Err(MolParseError::ResourceLimit {
                resource: "line bytes",
                ..
            }))
        ));
        assert!(reader.next().is_none());
    }

    let record_limits = SdfParseLimits {
        max_record_bytes: 3,
        ..SdfParseLimits::default()
    };
    let mut reader = SdfReader::with_limits("name\n$$$$\n", record_limits);
    assert!(matches!(
        reader.next(),
        Some(Err(MolParseError::ResourceLimit {
            resource: "record bytes",
            ..
        }))
    ));
    assert!(reader.next().is_none());
}

#[test]
fn delimiter_without_trailing_newline_is_skipped_iteratively() {
    let mut reader = SdfReader::new("$$$$");
    assert!(reader.next().is_none());
}
