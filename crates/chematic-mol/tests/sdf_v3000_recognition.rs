use std::io::Cursor;

use chematic_mol::{SdfFileReader, SdfReader, SdfRecordReader, read_sdf_conformer_ensembles};

const MIXED_SDF: &str =
    include_str!("../../../validation/cosmolkit_comparison/corpora/mixed-v2000-v3000-v1.sdf");
const V3000: &str = "M  END in title\n  chematic\nM  END in comment\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 2 O 1 0 0 0 CHG=-1\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n";

#[test]
fn every_sdf_reader_recognizes_mixed_v2000_v3000_and_fields() {
    for input in [MIXED_SDF.to_owned(), MIXED_SDF.replace('\n', "\r\n")] {
        let simple: Vec<_> = SdfReader::new(&input).collect::<Result<_, _>>().unwrap();
        assert_eq!(
            simple
                .iter()
                .map(|(mol, _)| mol.atom_count())
                .collect::<Vec<_>>(),
            [1, 2]
        );
        let records: Vec<_> = SdfRecordReader::new(&input)
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(records[1].mol.atom_count(), 2);
        assert_eq!(records[1].mol.bond_count(), 1);
        assert_eq!(records[1].mol.atoms().nth(1).unwrap().1.charge, -1);
        assert_eq!(records[1].meta.name, "M  END in title");
        assert_eq!(records[1].properties["id"], "new");
        assert_eq!(records[1].properties["note"], "M  END is data");

        for fast in [false, true] {
            let cursor = Cursor::new(input.as_bytes());
            let streamed: Vec<_> = if fast {
                SdfFileReader::fast(cursor).collect::<Result<_, _>>()
            } else {
                SdfFileReader::new(cursor).collect::<Result<_, _>>()
            }
            .unwrap();
            assert_eq!(
                streamed
                    .iter()
                    .map(|r| r.mol.atom_count())
                    .collect::<Vec<_>>(),
                [1, 2]
            );
            assert_eq!(streamed[1].properties["id"], "new");
            assert_eq!(streamed[1].properties["note"], "M  END is data");
        }
    }
}

#[test]
fn malformed_v3000_is_not_accepted_by_falling_back_to_v2000() {
    let malformed = V3000.replace("M  END\n", "");
    let input = format!("{malformed}$$$$\n");
    assert!(SdfReader::new(&input).next().unwrap().is_err());
    assert!(SdfRecordReader::new(&input).next().unwrap().is_err());
    assert!(
        SdfFileReader::new(Cursor::new(input.as_bytes()))
            .next()
            .unwrap()
            .is_err()
    );
    assert!(
        SdfFileReader::fast(Cursor::new(input.as_bytes()))
            .next()
            .unwrap()
            .is_err()
    );
}

#[test]
fn v3000_conformers_are_recognized_and_grouped() {
    let conformer = V3000.replace("M  V30 2 O 1 0 0 0", "M  V30 2 O 1 0 1 0");
    let input = format!("{conformer}$$$$\n{conformer}$$$$\n");
    let groups = read_sdf_conformer_ensembles(&input).unwrap();
    assert_eq!(groups.len(), 1);
    assert_eq!(groups[0].conformers.len(), 2);
    assert_eq!(groups[0].mol.atom_count(), 2);
}
