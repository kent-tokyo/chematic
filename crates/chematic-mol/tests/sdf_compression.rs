use chematic_mol::{MolParseError, SdfCompression, SdfFileReader, SdfInput, SdfParseLimits};
use flate2::{Compression, write::GzEncoder};
use std::io::{self, Cursor, Read, Write};

fn gzip(data: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(data).unwrap();
    encoder.finish().unwrap()
}

fn sdf() -> Vec<u8> {
    let mol = chematic_smiles::parse("CCO").unwrap();
    let block = chematic_mol::write_mol(&mol, &Default::default());
    format!("{block}> <Activity>\n7.2\n\n$$$$\n").into_bytes()
}

#[test]
fn auto_selection_is_independent_after_full_partial_and_failed_inputs() {
    let plain = sdf();
    let zipped = gzip(&plain);
    for _ in 0..10 {
        for data in [&zipped, &plain] {
            let input = SdfInput::new(Cursor::new(data), SdfCompression::Auto).unwrap();
            assert_eq!(
                input.compression(),
                if data == &zipped {
                    SdfCompression::Gzip
                } else {
                    SdfCompression::Plain
                }
            );
            let records = SdfFileReader::new(input)
                .collect::<Result<Vec<_>, _>>()
                .unwrap();
            assert_eq!(records.len(), 1);
            assert_eq!(records[0].mol.atom_count(), 3);
            assert_eq!(records[0].properties["Activity"], "7.2");
        }
        let mut partial = SdfInput::new(Cursor::new(&zipped), SdfCompression::Auto).unwrap();
        let mut one = [0];
        partial.read_exact(&mut one).unwrap();
        drop(partial);
        let invalid = SdfInput::new(Cursor::new(b"not gzip"), SdfCompression::Gzip).unwrap();
        assert!(SdfFileReader::new(invalid).next().unwrap().is_err());
        let fresh = SdfInput::new(Cursor::new(&plain), SdfCompression::Auto).unwrap();
        assert_eq!(
            SdfFileReader::new(fresh)
                .next()
                .unwrap()
                .unwrap()
                .mol
                .atom_count(),
            3
        );
    }
}

struct ShortReads<R> {
    inner: R,
    chunk: usize,
    interrupt: bool,
}

impl<R: Read> Read for ShortReads<R> {
    fn read(&mut self, output: &mut [u8]) -> io::Result<usize> {
        if self.interrupt {
            self.interrupt = false;
            return Err(io::ErrorKind::Interrupted.into());
        }
        let length = output.len().min(self.chunk);
        self.inner.read(&mut output[..length])
    }
}

#[test]
fn short_nonseekable_reads_and_interrupted_probe_replay_all_bytes() {
    let plain = sdf();
    for chunk in [1, 2, 7, 8192] {
        for data in [&plain, &gzip(&plain)] {
            let source = ShortReads {
                inner: Cursor::new(data),
                chunk,
                interrupt: true,
            };
            let mut input = SdfInput::new(source, SdfCompression::Auto).unwrap();
            let mut decoded = Vec::new();
            input.read_to_end(&mut decoded).unwrap();
            assert_eq!(decoded, plain);
        }
    }
}

#[test]
fn empty_and_single_byte_auto_inputs_are_not_consumed_or_reclassified() {
    for data in [b"".as_slice(), b"x", b"\x1f"] {
        let mut input = SdfInput::new(Cursor::new(data), SdfCompression::Auto).unwrap();
        assert_eq!(input.compression(), SdfCompression::Plain);
        let mut result = Vec::new();
        input.read_to_end(&mut result).unwrap();
        assert_eq!(result, data);
    }
    let input = SdfInput::new(Cursor::new(gzip(b"")), SdfCompression::Auto).unwrap();
    assert!(SdfFileReader::new(input).next().is_none());
}

#[test]
fn plain_is_an_explicit_no_decompression_policy() {
    let zipped = gzip(&sdf());
    let mut input = SdfInput::new(Cursor::new(&zipped), SdfCompression::Plain).unwrap();
    assert_eq!(input.compression(), SdfCompression::Plain);
    let mut actual = Vec::new();
    input.read_to_end(&mut actual).unwrap();
    assert_eq!(actual, zipped);
}

#[test]
fn concatenated_members_can_split_a_record_inside_a_physical_line() {
    let data = sdf();
    for split in [1, 13, data.len() / 2, data.len() - 1] {
        let mut compressed = gzip(&data[..split]);
        compressed.extend(gzip(&data[split..]));
        compressed.extend(gzip(&data));
        let input = SdfInput::new(Cursor::new(compressed), SdfCompression::Auto).unwrap();
        let records = SdfFileReader::new(input)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(records.len(), 2);
        assert!(
            records
                .iter()
                .all(|record| record.properties["Activity"] == "7.2")
        );
    }
}

#[test]
fn corrupt_truncated_or_junk_suffixed_gzip_reports_terminal_io_error() {
    let valid = gzip(&sdf());
    let mut bad_crc = valid.clone();
    bad_crc[valid.len() - 8] ^= 1;
    let mut trailing_junk = valid.clone();
    trailing_junk.extend(b"junk");
    let mut invalid_second_member = valid.clone();
    invalid_second_member.extend(b"\x1f\x8btruncated");
    for data in [
        bad_crc,
        valid[..valid.len() - 4].to_vec(),
        trailing_junk,
        invalid_second_member,
    ] {
        let input = SdfInput::new(Cursor::new(data), SdfCompression::Auto).unwrap();
        let mut reader = SdfFileReader::new(input);
        let results: Vec<_> = reader.by_ref().collect();
        assert!(matches!(results.last().unwrap(), Err(MolParseError::Io(_))));
        assert!(reader.next().is_none());
    }
}

#[test]
fn resource_limits_count_decoded_bytes_and_preserve_native_error_type() {
    let plain = sdf();
    let zipped = gzip(&plain);
    assert!(zipped.len() < plain.len() / 2);
    let limits = SdfParseLimits {
        max_input_bytes: plain.len() - 1,
        ..Default::default()
    };
    let input = SdfInput::new(Cursor::new(zipped), SdfCompression::Auto).unwrap();
    let results: Vec<_> = SdfFileReader::with_limits(input, limits).collect();
    assert!(matches!(
        results[0],
        Err(MolParseError::ResourceLimit {
            resource: "input bytes",
            ..
        })
    ));
}

#[test]
fn gzip_malformed_records_keep_positions_and_later_valid_records() {
    let mut data = sdf();
    data.extend(b"malformed\n$$$$\n");
    data.extend(sdf());
    let input = SdfInput::new(Cursor::new(gzip(&data)), SdfCompression::Auto).unwrap();
    let records: Vec<_> = SdfFileReader::new(input).collect();
    assert_eq!(records.len(), 3);
    assert!(records[0].is_ok());
    assert!(records[1].is_err());
    assert!(records[2].is_ok());
}

#[test]
fn auto_probe_failures_propagate_without_reclassifying_input() {
    struct Broken;
    impl Read for Broken {
        fn read(&mut self, _: &mut [u8]) -> io::Result<usize> {
            Err(io::ErrorKind::PermissionDenied.into())
        }
    }
    match SdfInput::new(Broken, SdfCompression::Auto) {
        Err(error) => assert_eq!(error.kind(), io::ErrorKind::PermissionDenied),
        Ok(_) => panic!("probe error must be returned"),
    }
}
