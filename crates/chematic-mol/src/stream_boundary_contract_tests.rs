//! Streaming errors, coordinate metadata, and strict grid diagnostics.
use crate::{cube::*, lammps_dump::*, opendx::*, record::MoleculeRecord, sdf::*, tdt::*};
use std::io::{self, BufReader, Cursor, Read, Write};

struct FaultReader {
    remaining: Cursor<Vec<u8>>,
}
impl Read for FaultReader {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if self.remaining.position() as usize == self.remaining.get_ref().len() {
            Err(io::Error::other("injected read failure"))
        } else {
            self.remaining.read(out)
        }
    }
}
fn fault_reader(prefix: &[u8]) -> BufReader<FaultReader> {
    BufReader::with_capacity(
        1,
        FaultReader {
            remaining: Cursor::new(prefix.to_vec()),
        },
    )
}
struct FaultWriter {
    budget: usize,
}
impl Write for FaultWriter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        if self.budget == 0 {
            return Err(io::Error::other("injected write failure"));
        }
        let count = bytes.len().min(self.budget);
        self.budget -= count;
        Ok(count)
    }
    fn flush(&mut self) -> io::Result<()> {
        Err(io::Error::other("injected flush failure"))
    }
}

const CUBE: &str = "density\ntest\n1 0 0 0\n1 1 0 0\n1 0 1 0\n1 0 0 1\n8 8 0 0 0\n0.5\n";
const DX: &str = "object 1 class gridpositions counts 1 1 1\norigin 0 0 0\ndelta 1 0 0\ndelta 0 1 0\ndelta 0 0 1\nobject 2 class gridconnections counts 1 1 1\nobject 3 class array type double rank 0 items 1 data follows\n0.5\n";
const DUMP: &str = "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\nITEM: BOX BOUNDS pp pp pp\n0 1\n0 1\n0 1\nITEM: ATOMS id x y z\n1 0 0 0\n";

#[test]
fn sdf_writer_variants_preserve_record_boundaries_and_visible_properties() {
    use crate::mol2000::*;
    let first = chematic_smiles::parse("CO").unwrap();
    let second = chematic_smiles::parse("N#N").unwrap();
    let first_meta = MolMetadata::default()
        .with_name("methanol")
        .with_comment("sample A");
    let second_meta = MolMetadata::default().with_name("nitrogen");
    let plain = write_sdf(&[(&first, &first_meta, &[]), (&second, &second_meta, &[])]);
    let records = parse_sdf(&plain).unwrap();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].1.name, "methanol");
    assert_eq!(records[1].1.name, "nitrogen");
    let charged = write_sdf_with_charges(&[
        (&first, &first_meta, &[], &[0.25, -0.25]),
        (&second, &second_meta, &[], &[]),
    ]);
    let records = SdfRecordReader::new(&charged)
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert_eq!(records[0].properties["PARTIAL_CHARGES"], "0.2500 -0.2500");
    assert!(!records[1].properties.contains_key("PARTIAL_CHARGES"));
    let props = std::collections::HashMap::from([
        ("activity".into(), "1.25\nmeasured".into()),
        ("_private".into(), "hidden".into()),
    ]);
    let mut buffer = "old record must disappear".to_string();
    write_laid_out_sdf_record_into(&mut buffer, &first, &first_meta, &props);
    let v3000 = write_sdf_record_v3000(&first, &first_meta, &[], &props);
    let record = SdfRecordReader::new(&buffer).next().unwrap().unwrap();
    assert_eq!(record.properties.len(), 1);
    assert_eq!(record.properties["activity"], "1.25\nmeasured");
    assert_eq!(
        chematic_smiles::rdkit_canonical_smiles(&record.mol).unwrap(),
        "CO"
    );
    assert_eq!(buffer.matches("$$$$").count(), 1);
    assert!(!buffer.contains("old record"));
    // The compact SDF record reader is V2000-specific; inspect the V3000
    // CTAB through its dedicated reader and check the appended data field.
    let ctab_end = v3000.find("M  END\n").unwrap() + "M  END\n".len();
    let report = crate::mol3000::read_mol_v3000_with_diagnostics(&v3000[..ctab_end]).unwrap();
    assert_eq!(report.mol.atom_count(), 2);
    assert_eq!(report.mol.bond_count(), 1);
    assert!(v3000.contains("> <activity>\n1.25\nmeasured\n\n$$$$\n"));
    assert!(!v3000.contains("_private"));
}

#[test]
fn v3000_conformer_writer_preserves_coordinates_attributes_and_enhanced_stereo_groups() {
    use crate::mol2000::MolMetadata;
    use chematic_core::{AtomIdx, Coords3D, Point3, StereoGroup, StereoGroupKind};
    let mut mol = chematic_smiles::parse("[NH3+:1][13CH2:2]C").unwrap();
    let conformer = Coords3D {
        points: vec![
            Point3::new(0.1, 0.2, 0.3),
            Point3::new(1.1, -0.2, 1.3),
            Point3::new(2.1, 1.2, -0.3),
        ],
    };
    let metadata = MolMetadata {
        name: "labeled amine".into(),
        comment: "conformer metadata".into(),
        v3000_sgroups: vec!["1 DAT 0 ATOMS=(1 1) FIELDNAME=note FIELDDATA=sample".into()],
        v3000_atom_properties: vec![(1, "INVRET=1".into())],
        v3000_bond_properties: vec![(1, "TOPO=1".into())],
    };
    for kind in [
        StereoGroupKind::Absolute,
        StereoGroupKind::Or(2),
        StereoGroupKind::And(3),
    ] {
        mol.set_stereo_groups(vec![StereoGroup::new(kind, vec![AtomIdx(1)])]);
        let text = crate::mol3000::write_mol_v3000_with_conformer(&mol, &metadata, &conformer);
        assert!(
            text.contains("INVRET=1")
                && text.contains("TOPO=1")
                && text.contains("FIELDDATA=sample")
        );
        let report = crate::mol3000::read_mol_v3000_with_diagnostics(&text).unwrap();
        assert_eq!(report.mol.stereo_groups(), mol.stereo_groups());
        assert_eq!(report.metadata.name, metadata.name);
        assert_eq!(report.metadata.comment, metadata.comment);
        assert_eq!(report.metadata.v3000_sgroups, metadata.v3000_sgroups);
        let recovered = report.conformer.unwrap();
        assert_eq!(recovered, conformer);
        assert_eq!(report.mol.atom(AtomIdx(0)).charge, 1);
        assert_eq!(report.mol.atom(AtomIdx(1)).isotope, Some(13));
    }
}

#[test]
fn v3000_count_block_and_collection_errors_fail_before_returning_a_graph() {
    use crate::{mol2000::MolMetadata, mol3000::*};
    let mol = chematic_smiles::parse("CO").unwrap();
    let base = write_mol_v3000(&mol, &MolMetadata::default(), &[]);
    for (replacement, diagnostic) in [
        ("COUNTS 2", "fewer"),
        ("COUNTS x 1 0 0 0", "atom count"),
        ("COUNTS 2 x 0 0 0", "bond count"),
        ("COUNTS 2 1 x 0 0", "SGROUP count"),
        ("COUNTS 3 1 0 0 0", "expected 3 atoms"),
        ("COUNTS 2 2 0 0 0", "bonds"),
        ("COUNTS 2 1 1 0 0", "SGROUP records"),
    ] {
        let text = base.replace("COUNTS 2 1 0 0 0", replacement);
        let message = read_mol_v3000_with_diagnostics(&text)
            .err()
            .unwrap()
            .to_string();
        assert!(message.contains(diagnostic), "{replacement}: {message}");
    }
    for (field, diagnostic) in [
        ("ATOMS=()", "empty"),
        ("ATOMS=(x 1)", "integer"),
        ("ATOMS=(2 1)", "declares"),
        ("ATOMS=(1 x)", "non-integer"),
        ("ATOMS=(1 9)", "unknown atom"),
        ("ATOMS=(0)", "no atom ids"),
        ("ATOMS=(1 1", "closed"),
        ("", "missing ATOMS"),
    ] {
        let block = format!(
            "M  V30 BEGIN COLLECTION\nM  V30 MDLV30/STEABS {field}\nM  V30 END COLLECTION\nM  V30 END CTAB"
        );
        let text = base.replace("M  V30 END CTAB", &block);
        let message = read_mol_v3000_with_diagnostics(&text)
            .err()
            .unwrap()
            .to_string();
        assert!(message.contains(diagnostic), "{field}: {message}");
    }
    for (from, to, diagnostic) in [
        (
            "M  V30 BEGIN BOND",
            "M  V30 WRONG BLOCK",
            "expected BEGIN BOND",
        ),
        ("M  V30 END BOND", "", "bond line"),
    ] {
        let text = base.replace(from, to);
        let message = read_mol_v3000_with_diagnostics(&text)
            .err()
            .unwrap()
            .to_string();
        assert!(message.contains(diagnostic), "{from}: {message}");
    }
}
fn replace_line(text: &str, index: usize, line: &str) -> String {
    let mut lines: Vec<_> = text.lines().collect();
    lines[index] = line;
    lines.join("\n") + "\n"
}

#[test]
fn cube_required_fields_report_specific_error_categories() {
    assert!(parse_cube(CUBE).is_ok());
    for (line, replacement, diagnostic) in [
        (2, "x 0 0 0", "count"),
        (2, "1 0 0 0 x", "count"),
        (2, "1 x 0 0", "origin"),
        (2, "1 NaN 0 0", "not finite"),
        (2, "-1 0 0 0", "dataset"),
        (3, "1 0 0", "axis"),
        (3, "x 1 0 0", "axis"),
        (3, "0 1 0 0", "positive"),
        (4, "0 0 1 0", "positive"),
        (5, "1 0 0 x", "axis"),
        (6, "8 8 0 0", "atom"),
        (6, "x 8 0 0 0", "atom"),
        (6, "999 8 0 0 0", "atomic number"),
        (7, "x", "voxel"),
        (7, "0.5 0.6", "extra"),
        (7, "", "values"),
    ] {
        let error = parse_cube(&replace_line(CUBE, line, replacement)).unwrap_err();
        let message = error.to_string();
        assert!(
            message.contains(diagnostic),
            "{line}: {replacement}: {message}"
        );
    }
    for count in 0..7 {
        let error =
            parse_cube(&CUBE.lines().take(count).collect::<Vec<_>>().join("\n")).unwrap_err();
        assert!(error.to_string().contains("end"));
    }
    let limits = CubeParseLimits {
        max_input_bytes: CUBE.len() - 1,
        ..Default::default()
    };
    assert!(
        parse_cube_with_limits(CUBE, &limits)
            .unwrap_err()
            .to_string()
            .contains("byte")
    );
}

#[test]
fn opendx_required_fields_and_grid_consistency_have_diagnostics() {
    assert!(parse_opendx(DX).is_ok());
    for (line, replacement, diagnostic) in [
        (
            0,
            "object 4 class gridpositions counts 1 1 1",
            "gridpositions",
        ),
        (
            0,
            "object 1 class gridpositions counts x 1 1",
            "gridpositions",
        ),
        (1, "origin 0", "origin"),
        (1, "origin x 0 0", "origin"),
        (1, "origin NaN 0 0", "not finite"),
        (2, "delta 1", "delta"),
        (
            5,
            "object 2 class gridconnections counts 2 1 1",
            "gridconnections",
        ),
        (
            6,
            "object 3 class array type float rank 0 items 1 data follows",
            "array",
        ),
        (
            6,
            "object 3 class array type double rank 0 1 data follows",
            "items",
        ),
        (
            6,
            "object 3 class array type double rank 0 items x data follows",
            "items",
        ),
        (
            6,
            "object 3 class array type double rank 0 items 2 data follows",
            "items",
        ),
        (7, "x", "value"),
        (7, "", "values"),
        (7, "0.5 0.6", "extra"),
    ] {
        let message = parse_opendx(&replace_line(DX, line, replacement))
            .unwrap_err()
            .to_string();
        assert!(
            message.contains(diagnostic),
            "{line}: {replacement}: {message}"
        );
    }
    for count in 0..7 {
        let error =
            parse_opendx(&DX.lines().take(count).collect::<Vec<_>>().join("\n")).unwrap_err();
        assert!(error.to_string().contains("end"));
    }
    let limits = OpenDxParseLimits {
        max_input_bytes: DX.len() - 1,
        ..Default::default()
    };
    assert!(
        parse_opendx_with_limits(DX, &limits)
            .unwrap_err()
            .to_string()
            .contains("byte")
    );
}

#[test]
fn lammps_dump_invalid_fields_and_stream_limits_are_reported() {
    assert!(parse_lammps_dump_frame(DUMP).is_ok());
    for (line, replacement, diagnostic) in [
        (0, "ITEM: TIME", "TIMESTEP"),
        (1, "x", "header"),
        (3, "x", "header"),
        (4, "ITEM: BOX BOUNDS pp pp", "header"),
        (4, "ITEM: BOUNDS pp pp pp", "BOX BOUNDS"),
        (5, "0 x", "header"),
        (5, "0 NaN", "non-finite"),
        (5, "0 0", "box"),
        (5, "0 1 2", "box"),
        (8, "ITEM: MOLECULES id x y z", "ATOMS"),
        (9, "x 0 0 0", "id"),
        (9, "1 NaN 0 0", "x"),
        (9, "1 0 0", "column"),
        (9, "ITEM: TIMESTEP", "ATOMS"),
    ] {
        let message = parse_lammps_dump_frame(&replace_line(DUMP, line, replacement))
            .unwrap_err()
            .to_string();
        assert!(
            message.contains(diagnostic),
            "{line}: {replacement}: {message}"
        );
    }
    for count in 1..10 {
        let error =
            parse_lammps_dump_frame(&DUMP.lines().take(count).collect::<Vec<_>>().join("\n"))
                .unwrap_err();
        assert!(error.to_string().contains("end"));
    }
    for limits in [
        LammpsDumpParseLimits {
            max_input_bytes: DUMP.len() - 1,
            ..Default::default()
        },
        LammpsDumpParseLimits {
            max_line_bytes: 10,
            ..Default::default()
        },
    ] {
        let mut reader = LammpsDumpReader::with_limits(Cursor::new(DUMP), limits);
        assert!(
            reader
                .next()
                .unwrap()
                .unwrap_err()
                .to_string()
                .contains("limit")
        );
        assert!(reader.next().is_none());
    }
}

#[test]
fn streaming_readers_propagate_io_errors_at_each_record_boundary() {
    let mol = chematic_smiles::parse("CO").unwrap();
    let sdf = crate::write_mol(&mol, &Default::default()) + "$$$$\n";
    for (input, kind) in [
        (CUBE, 0),
        (DUMP, 1),
        (sdf.as_str(), 2),
        ("$SMI<CO>\nNAME<methanol>\n|\n", 3),
    ] {
        for offset in input.match_indices('\n').map(|(i, _)| i + 1).chain([0]) {
            let message = match kind {
                0 => CubeFileReader::new(fault_reader(&input.as_bytes()[..offset]))
                    .read()
                    .unwrap_err()
                    .to_string(),
                1 => LammpsDumpReader::new(fault_reader(&input.as_bytes()[..offset]))
                    .find_map(|r| r.err())
                    .unwrap()
                    .to_string(),
                2 => SdfFileReader::new(fault_reader(&input.as_bytes()[..offset]))
                    .find_map(|r| r.err())
                    .unwrap()
                    .to_string(),
                _ => TdtRecordReader::new(
                    fault_reader(&input.as_bytes()[..offset]),
                    Default::default(),
                )
                .find_map(|r| r.err())
                .unwrap()
                .to_string(),
            };
            assert!(
                message.contains("injected read failure"),
                "{kind}: {offset}: {message}"
            );
        }
    }
}

#[test]
fn tdt_coordinates_property_selection_and_writer_failures() {
    let mut record = MoleculeRecord::new(chematic_smiles::parse("CO").unwrap());
    record.name = "methanol".into();
    record.coordinates_2d = Some(vec![[0.125, 0.25], [1.5, 2.0]]);
    record.coordinates_3d = Some(vec![[0.125, 0.25, -0.5], [1.5, 2.0, 3.0]]);
    record.properties = vec![
        ("SOURCE".into(), "first\nsecond".into()),
        ("OMIT".into(), "hidden".into()),
    ];
    let options = TdtWriterOptions {
        write_2d: true,
        ..Default::default()
    };
    let mut writer = TdtRecordWriter::new(Vec::new(), options.clone());
    writer.set_name_tag(Some("TITLE".into()));
    writer.set_precision(3);
    writer.set_properties(Some(vec!["SOURCE".into(), "MISSING".into()]));
    writer.write_record(&record).unwrap();
    writer.close().unwrap();
    writer.close().unwrap();
    let output = String::from_utf8(writer.into_inner()).unwrap();
    assert!(output.contains("TITLE<methanol>"));
    assert!(output.contains("SOURCE<first second>"));
    assert!(!output.contains("OMIT"));
    assert_eq!(output.matches("|\n").count(), 1);
    let decoded = TdtRecordReader::new(
        Cursor::new(&output),
        TdtReaderOptions {
            name_tag: Some("TITLE".into()),
            read_2d: true,
            read_3d: true,
            ..Default::default()
        },
    )
    .next()
    .unwrap()
    .unwrap();
    assert_eq!(decoded.name, record.name);
    assert_eq!(decoded.coordinates_2d, record.coordinates_2d);
    assert_eq!(decoded.coordinates_3d, record.coordinates_3d);
    let raw = TdtRecordReader::new(Cursor::new(&output), Default::default())
        .next()
        .unwrap()
        .unwrap();
    assert!(raw.get_property("2D").is_some());
    assert!(raw.get_property("3D").is_some());
    for budget in 0..output.len() + 50 {
        let mut failing = TdtRecordWriter::new(FaultWriter { budget }, options.clone());
        let error = failing
            .write_record(&record)
            .and_then(|()| failing.write_record(&record))
            .and_then(|()| failing.close())
            .unwrap_err();
        assert!(error.to_string().contains("injected"));
    }
    let mut empty = TdtRecordWriter::new(Vec::new(), Default::default());
    empty.close().unwrap();
    assert!(empty.into_inner().is_empty());
}

#[test]
fn tdt_malformed_coordinate_diagnostics_and_strict_stopping() {
    for (tag, values, expected) in [
        ("3D", "0,1", "expected 3"),
        ("3D", "x,1,2", "invalid"),
        ("2D", "0", "expected 2"),
        ("2D", "x,1", "invalid"),
    ] {
        let input = format!("$SMI<C>\n{tag}<{values};>\n|\n");
        let error = TdtRecordReader::new(
            Cursor::new(input),
            TdtReaderOptions {
                read_2d: true,
                read_3d: true,
                ..Default::default()
            },
        )
        .next()
        .unwrap()
        .err()
        .unwrap();
        let message = error.to_string();
        assert!(
            message.contains(tag) && message.contains(expected),
            "{message}"
        );
    }
    let mut strict = TdtRecordReader::new(
        Cursor::new("$SMI<invalid>\n|\n$SMI<C>\n|\n"),
        TdtReaderOptions {
            strict_parsing: true,
            ..Default::default()
        },
    );
    assert!(
        strict
            .next()
            .unwrap()
            .err()
            .unwrap()
            .to_string()
            .contains("invalid SMILES")
    );
    assert!(strict.next().is_none());
    // Recovery must leave an unconsumed $SMI boundary for the next call.
    let mut reader = TdtRecordReader::new(Cursor::new("\nBAD\n$SMI<C>\n|\n"), Default::default());
    assert!(
        reader
            .next()
            .unwrap()
            .err()
            .unwrap()
            .to_string()
            .contains("missing $SMI")
    );
    assert_eq!(reader.next().unwrap().unwrap().mol.atom_count(), 1);
}

#[test]
fn fast_sdf_batches_track_progress_and_cancellation() {
    let mol = chematic_smiles::parse("CO").unwrap();
    let record = crate::write_mol(&mol, &Default::default()) + "$$$$\n";
    let mut reader = SdfBatchReader::fast(Cursor::new(record.repeat(3)), 2);
    let first = reader.next().unwrap();
    assert_eq!(first.len(), 2);
    assert!(!first.is_empty());
    assert_eq!(reader.progress().records_emitted, 2);
    reader.cancel();
    assert!(reader.is_cancelled());
    assert!(reader.next().is_none());
    assert!(reader.manifest_json().contains("cancelled"));
    let mut all = SdfBatchReader::fast(Cursor::new(record.repeat(3)), 2);
    assert_eq!(all.next().unwrap().len(), 2);
    assert_eq!(all.next().unwrap().len(), 1);
    assert!(all.next().is_none());
    assert!(all.next().is_none());
    assert_eq!(all.progress().records_emitted, 3);
}
