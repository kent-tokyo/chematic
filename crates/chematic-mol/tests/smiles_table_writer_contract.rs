use chematic_mol::record::MoleculeRecord;
use chematic_mol::smiles_table::{
    Delimiter, SmilesReaderOptions, SmilesRecordReader, SmilesRecordWriter, SmilesWriterOptions,
};

#[test]
fn selected_columns_roundtrip_with_literal_delimiters() {
    for delimiter in [Delimiter::Tab, Delimiter::Custom(b'|')] {
        let mut record = MoleculeRecord::new(chematic_smiles::parse("CCO").unwrap());
        record.name = "ethanol".into();
        record.properties = vec![
            ("tag".into(), "solvent".into()),
            ("score".into(), "7".into()),
        ];
        let mut writer = SmilesRecordWriter::new(
            Vec::new(),
            SmilesWriterOptions {
                delimiter,
                ..Default::default()
            },
        );
        writer.set_properties(vec!["score".into(), "tag".into()]);
        writer.write_record(&record).unwrap();
        writer.flush().unwrap();
        let bytes = writer.into_inner();
        let opts = SmilesReaderOptions {
            delimiter,
            title_line: true,
            ..Default::default()
        };
        let mut reader = SmilesRecordReader::new(bytes.as_slice(), opts);
        let restored = reader.next().unwrap().unwrap();
        assert_eq!(restored.name, "ethanol");
        assert_eq!(restored.get_property("score"), Some("7"));
        assert_eq!(restored.get_property("tag"), Some("solvent"));
        assert_eq!(
            chematic_smiles::canonical_smiles(&restored.mol),
            chematic_smiles::canonical_smiles(&record.mol)
        );
        assert!(reader.next().is_none());
    }
}

#[test]
fn reader_propagates_underlying_io_failure() {
    use std::io::{self, BufRead, Read};
    struct FailingReader;
    impl Read for FailingReader {
        fn read(&mut self, _: &mut [u8]) -> io::Result<usize> {
            Err(io::Error::other("storage failed"))
        }
    }
    impl BufRead for FailingReader {
        fn fill_buf(&mut self) -> io::Result<&[u8]> {
            Err(io::Error::other("storage failed"))
        }
        fn consume(&mut self, _: usize) {}
    }
    let mut reader = SmilesRecordReader::new(FailingReader, SmilesReaderOptions::default());
    let error = match reader.next().unwrap() {
        Ok(_) => panic!("I/O failure must not yield a record"),
        Err(error) => error,
    };
    assert!(error.to_string().contains("IO error: storage failed"));
}
