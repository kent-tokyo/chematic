//! Resource limits and metadata round trips at format boundaries.
use crate::mrv::{MrvError, MrvWriteOptions, parse_mrv, write_mrv};
use crate::orca::*;
use crate::record::MoleculeRecord;

#[test]
fn orca_input_limits_accept_exact_sizes_and_reject_the_next_item() {
    let text = "! B3LYP def2-SVP Opt\n%pal nprocs 2 end\n%scf MaxIter 20 end\n* xyz 0 1\nO 0 0 0\nH 0 1 0\nH 0 -1 0\n*\n";
    let parsed = parse_orca_input(text).unwrap();
    let sizes = [
        text.len(),
        text.lines().map(str::len).max().unwrap(),
        text.lines().count(),
        parsed.keywords.len(),
        parsed.blocks.len(),
        parsed.blocks.iter().map(|b| b.raw.len()).max().unwrap(),
        3,
    ];
    let resources = [
        "input bytes",
        "line bytes",
        "lines",
        "keywords",
        "blocks",
        "block bytes",
        "atoms",
    ];
    for (which, size) in sizes.into_iter().enumerate() {
        let set = |limit| {
            let mut l = OrcaInputParseLimits::default();
            match which {
                0 => l.max_input_bytes = limit,
                1 => l.max_line_bytes = limit,
                2 => l.max_lines = limit,
                3 => l.max_keywords = limit,
                4 => l.max_blocks = limit,
                5 => l.max_block_bytes = limit,
                _ => l.max_atoms = limit,
            }
            l
        };
        assert_eq!(
            parse_orca_input_with_limits(text, &set(size)).unwrap(),
            parsed
        );
        let error = parse_orca_input_with_limits(text, &set(size - 1)).unwrap_err();
        assert!(
            matches!(&error, OrcaInputError::ResourceLimit{resource,..} if *resource == resources[which]),
            "{error}"
        );
        assert!(error.to_string().contains(resources[which]));
    }
    assert_eq!(
        parse_orca_input(&write_orca_input(&parsed)).unwrap(),
        parsed
    );
}

#[test]
fn orca_output_limits_bound_geometry_trajectories_and_frequencies() {
    let frame = "CARTESIAN COORDINATES (ANGSTROEM)\n---------------------\nO 0.0 0.0 0.0\nH 0.0 1.0 0.0\nH 0.0 -1.0 0.0\n\n";
    let text = format!(
        "{frame}{frame}VIBRATIONAL FREQUENCIES\n---------------------\n0: 123.0 cm**-1\n1: 456.0 cm**-1\n\nFINAL SINGLE POINT ENERGY -75.0\nORCA TERMINATED NORMALLY\n"
    );
    let parsed = parse_orca_output(&text).unwrap();
    assert_eq!(parsed.trajectory.len(), 2);
    assert_eq!(parsed.frequencies_cm1, [123.0, 456.0]);
    assert_eq!(parsed.final_energy_hartree, Some(-75.0));
    let sizes = [
        text.len(),
        text.lines().map(str::len).max().unwrap(),
        text.lines().count(),
        2,
        3,
        2,
    ];
    for (which, size) in sizes.into_iter().enumerate() {
        let set = |limit| {
            let mut l = OrcaOutputParseLimits::default();
            match which {
                0 => l.max_input_bytes = limit,
                1 => l.max_line_bytes = limit,
                2 => l.max_lines = limit,
                3 => l.max_geometry_frames = limit,
                4 => l.max_geometry_atoms = limit,
                _ => l.max_frequencies = limit,
            }
            l
        };
        assert!(parse_orca_output_with_limits(&text, &set(size)).is_ok());
        let error = parse_orca_output_with_limits(&text, &set(size - 1))
            .err()
            .unwrap();
        assert!(!error.to_string().is_empty());
    }
    for value in ["NaN", "inf", "-inf"] {
        for text in [
            format!("FINAL SINGLE POINT ENERGY {value}"),
            format!("CARTESIAN COORDINATES (ANGSTROEM)\n---\nC {value} 0 0\n"),
            format!("VIBRATIONAL FREQUENCIES\n---\n0: {value} cm**-1\n"),
        ] {
            let error = parse_orca_output(&text).err().unwrap();
            assert!(error.to_string().contains("finite"), "{error}");
        }
    }
}

fn mrv(bond: &str) -> String {
    format!(
        "<?xml version=\"1.0\"?>\n<!-- preamble -->\n<cml><molecule><atomArray><atom id=\"a1\" elementType=\"N\"/><atom id=\"a2\" elementType=\"C\"/></atomArray><bondArray>{bond}</bondArray></molecule></cml>"
    )
}
#[test]
fn mrv_stereo_encodings_and_coordination_conventions_have_typed_boundaries() {
    for stereo in [
        "<bondStereo dictRef=\"cml:W\"/>",
        "<bondStereo dictRef=\"cml:H\"/>",
        "<bondStereo convention=\"MDL\" conventionValue=\"1\"/>",
        "<bondStereo convention=\"MDL\" conventionValue=\"6\"/>",
        "<bondStereo convention=\"MDL\"/>",
        "<bondStereo/>",
    ] {
        let text = mrv(&format!(
            "<bond atomRefs2=\"a1 a2\" order=\"1\">{stereo}</bond>"
        ));
        let record = parse_mrv(&text).unwrap();
        assert_eq!(record.mol.bond_count(), 1);
    }
    for stereo in [
        "<bondStereo dictRef=\"unknown\"/>",
        "<bondStereo convention=\"MDL\" conventionValue=\"99\"/>",
        "<bondStereo dictRef=\"cml:W\">H</bondStereo>",
    ] {
        let error = parse_mrv(&mrv(&format!("<bond atomRefs2=\"a1 a2\">{stereo}</bond>")))
            .err()
            .unwrap();
        assert!(matches!(error, MrvError::InvalidBondStereo { .. }));
        assert!(error.to_string().contains("bondStereo"));
    }
    let record = parse_mrv(&mrv("<bond atomRefs2=\"a1 a2\" convention=\"cxn:coord\"/>")).unwrap();
    assert_eq!(
        record.mol.bonds().next().unwrap().1.order,
        chematic_core::BondOrder::Dative
    );
    let error = parse_mrv(&mrv("<bond atomRefs2=\"a1 a2\" convention=\"unknown\"/>"))
        .err()
        .unwrap();
    assert!(matches!(error, MrvError::InvalidBondOrder { .. }));
    for text in ["<?xml", "<!--", "<cml><!--", "<cml><molecule/></cml>"] {
        let error = parse_mrv(text).err().unwrap();
        assert!(!error.to_string().is_empty());
    }
}
#[test]
fn mrv_writer_preserves_three_dimensional_coordinates_and_atom_metadata() {
    let mut record = MoleculeRecord::new(chematic_smiles::parse("[15NH3+:7]C#N").unwrap());
    record.coordinates_3d = Some(vec![[0.0, 1.0, -2.0], [1.5, -0.25, 0.5], [2.7, 0.0, 3.0]]);
    let written = write_mrv(&record, &MrvWriteOptions::default()).unwrap();
    let round_trip = parse_mrv(&written).unwrap();
    assert_eq!(round_trip.coordinates_3d, record.coordinates_3d);
    assert_eq!(
        round_trip.mol.atom(chematic_core::AtomIdx(0)).isotope,
        Some(15)
    );
    assert_eq!(
        round_trip.mol.atom(chematic_core::AtomIdx(0)).atom_map,
        Some(7)
    );
    assert_eq!(round_trip.mol.atom(chematic_core::AtomIdx(0)).charge, 1);
    assert_eq!(
        round_trip.mol.bonds().last().unwrap().1.order,
        chematic_core::BondOrder::Triple
    );
}

#[test]
fn qcschema_rejections_identify_the_invalid_field_or_graph_structure() {
    use crate::qcschema::*;
    use serde_json::json;
    let base = json!({"symbols":["H","H"],"geometry":[0,0,0,0,0,1.4],"connectivity":[[0,1,1.0]]});
    assert!(parse_qcschema_molecule(&base.to_string()).is_ok());
    for (field, value, context) in [
        ("schema_name", json!("wrong-schema"), "schema_name"),
        ("symbols", json!(false), "symbols"),
        ("geometry", json!([]), "length"),
        ("connectivity", json!([[0, 2, 1.0]]), "index"),
        (
            "connectivity",
            json!([[0, 1, 1.0], [1, 0, 1.0]]),
            "duplicate",
        ),
        ("connectivity", json!([[0, 0, 1.0]]), "self"),
        ("fragments", json!([[0], [2]]), "fragment"),
        ("molecular_charge", json!("zero"), "molecular_charge"),
        ("masses", json!([1.0]), "length"),
    ] {
        let mut document = base.clone();
        document[field] = value;
        let error = parse_qcschema_molecule(&document.to_string())
            .err()
            .unwrap_or_else(|| panic!("unexpected acceptance of {field}: {document}"));
        let message = error.to_string();
        assert!(
            message.contains("QCSchema") && message.contains(context),
            "{field}: {message}"
        );
    }
    for input in ["{", "{}", "[]"] {
        let error = parse_qcschema_molecule(input).err().unwrap();
        assert!(error.to_string().contains("QCSchema"));
    }
    let mut input = json!({"molecule":base,"driver":"energy","model":{"method":"hf"}});
    assert!(parse_atomic_input(&input.to_string()).is_ok());
    input["driver"] = json!("not-a-driver");
    let error = parse_atomic_input(&input.to_string()).err().unwrap();
    assert!(error.to_string().contains("not-a-driver"));
    for (which, limit) in [(0, 1), (1, 0), (2, 1), (3, 0)] {
        let mut limits = QcSchemaParseLimits::default();
        match which {
            0 => limits.max_input_bytes = limit,
            1 => limits.max_json_depth = limit,
            2 => limits.max_array_items = limit,
            _ => limits.max_string_bytes = limit,
        }
        let error = parse_qcschema_molecule_with_limits(&base.to_string(), &limits)
            .err()
            .unwrap();
        assert!(matches!(error, QcSchemaError::ResourceLimit { .. }));
        assert!(error.to_string().contains("limit"));
    }
}
#[test]
fn volumetric_writers_reject_nonfinite_and_inconsistent_public_grids() {
    use crate::volumetric::*;
    let base = VolumetricGrid {
        origin: [0.0; 3],
        axes: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        shape: [1, 1, 1],
        values: vec![0.5],
        atoms: vec![GridAtom {
            element: chematic_core::Element::C,
            charge: 6.0,
            position: [0.0; 3],
        }],
        units: GridUnits::Bohr,
    };
    assert!(base.validate().is_ok());
    for mutation in 0..8 {
        let mut grid = base.clone();
        match mutation {
            0 => grid.shape[0] = 0,
            1 => grid.shape = [usize::MAX, 2, 2],
            2 => grid.values.clear(),
            3 => grid.origin[2] = f64::NAN,
            4 => grid.axes[1][2] = f64::INFINITY,
            5 => grid.values[0] = f64::NEG_INFINITY,
            6 => grid.atoms[0].position[1] = f64::NAN,
            _ => grid.atoms[0].charge = f64::INFINITY,
        }
        let error = grid.validate().unwrap_err();
        assert!(!error.to_string().is_empty());
        let cube = crate::cube::write_cube(&grid).unwrap_err();
        if mutation == 1 {
            assert!(matches!(
                cube,
                crate::cube::CubeError::DimensionOutOfRange { axis: 0, .. }
            ));
        } else {
            assert!(matches!(cube, crate::cube::CubeError::Grid(_)));
        }
        assert!(matches!(
            crate::opendx::write_opendx(&grid),
            Err(crate::opendx::OpenDxError::Grid(_))
        ));
    }
}

#[test]
fn json_format_budgets_report_the_resource_and_boundary() {
    use crate::{cjson::*, ket::*, moljson::*};
    let mol = chematic_smiles::parse("[13CH3+]C#N").unwrap();
    let cjson = write_cjson(&mol, &[(0., 0., 0.), (1., 0., 0.), (2., 0., 0.)]);
    let moljson = write_moljson(&mol);
    let ket = write_ket_3d(&mol, &[(0., 0., 0.), (1., 0., 0.), (2., 0., 0.)]);
    // Each independent budget must fail before returning a partial graph.
    macro_rules! budget {
        ($parse:path, $limits:ty, $error:path, $text:expr, $field:ident, $limit:expr) => {{
            let text = &$text;
            let mut limits = <$limits>::default();
            assert!($parse(text, &limits).is_ok());
            limits.$field = $limit;
            let error = $parse(text, &limits).err().unwrap();
            let $error {
                resource,
                actual,
                limit,
            } = &error
            else {
                panic!("{error}")
            };
            assert!(*actual > *limit);
            let message = error.to_string();
            assert!(
                message.contains(resource)
                    && message.contains(&actual.to_string())
                    && message.contains(&limit.to_string()),
                "{message}"
            );
        }};
    }
    macro_rules! json_budgets {
        ($parse:path, $limits:ty, $error:path, $text:expr) => {
            budget!(
                $parse,
                $limits,
                $error,
                $text,
                max_input_bytes,
                $text.len() - 1
            );
            budget!($parse, $limits, $error, $text, max_json_depth, 0);
            budget!($parse, $limits, $error, $text, max_array_items, 1);
            budget!($parse, $limits, $error, $text, max_string_bytes, 0);
            budget!($parse, $limits, $error, $text, max_atoms, 2);
            budget!($parse, $limits, $error, $text, max_bonds, 1);
        };
    }
    json_budgets!(
        parse_cjson_with_limits,
        CjsonParseLimits,
        CjsonError::ResourceLimit,
        cjson
    );
    json_budgets!(
        parse_moljson_with_limits,
        MolJsonParseLimits,
        MolJsonError::ResourceLimit,
        moljson
    );
    budget!(
        parse_ket_3d_with_limits,
        KetParseLimits,
        KetError::ResourceLimit,
        ket,
        max_input_bytes,
        ket.len() - 1
    );
    budget!(
        parse_ket_3d_with_limits,
        KetParseLimits,
        KetError::ResourceLimit,
        ket,
        max_atoms,
        2
    );
    budget!(
        parse_ket_3d_with_limits,
        KetParseLimits,
        KetError::ResourceLimit,
        ket,
        max_bonds,
        1
    );
}

#[test]
fn json_graph_errors_identify_the_missing_field_or_reference() {
    use crate::{ket::*, moljson::*};
    use serde_json::json;
    for (text, field) in [
        (json!({}), "atoms"),
        (json!({"atoms":[{}]}), "atoms[].id"),
        (json!({"atoms":[{"id":"a"}]}), "atoms[].element"),
        (
            json!({"atoms":[{"id":"a","element":"C"}],"bonds":[{}]}),
            "bonds[].source_id",
        ),
        (
            json!({"atoms":[{"id":"a","element":"C"}],"bonds":[{"source_id":"a"}]}),
            "bonds[].target_id",
        ),
    ] {
        let e = parse_moljson(&text.to_string()).err().unwrap();
        assert_eq!(e, MolJsonError::MissingField(field));
        assert!(e.to_string().contains(field));
    }
    for end in ["source_id", "target_id"] {
        let mut text = json!({"atoms":[{"id":"a","element":"C"}],"bonds":[{"id":"b17","source_id":"a","target_id":"a"}]});
        text["bonds"][0][end] = json!("unknown19");
        let e = parse_moljson(&text.to_string()).err().unwrap();
        assert!(
            matches!(&e,MolJsonError::InvalidBondRef{bond_id,ref_id} if bond_id=="b17" && ref_id=="unknown19")
        );
        assert!(e.to_string().contains("b17") && e.to_string().contains("unknown19"));
    }
    for (text, field) in [
        (json!({}), "atoms"),
        (json!({"root":{}}), "root.nodes[0].$ref"),
        (
            json!({"root":{"nodes":[{"$ref":"absent"}]}}),
            "molecule node referenced by $ref",
        ),
        (json!({"atoms":[{}]}), "atom.label"),
        (json!({"atoms":[{"label":"C"}],"bonds":[{}]}), "bond.atoms"),
        (
            json!({"atoms":[{"label":"C"}],"bonds":[{"atoms":[]}]}),
            "bond.atoms[0]",
        ),
        (
            json!({"atoms":[{"label":"C"}],"bonds":[{"atoms":[0]}]}),
            "bond.atoms[1]",
        ),
    ] {
        let e = parse_ket(&text.to_string()).err().unwrap();
        assert_eq!(e, KetError::MissingField(field));
        assert!(e.to_string().contains(field));
    }
    for references in [[0, 9], [9, 0]] {
        let text = json!({"atoms":[{"label":"C"}],"bonds":[{"atoms":references}]});
        let e = parse_ket(&text.to_string()).err().unwrap();
        assert!(matches!(
            e,
            KetError::InvalidAtomIndex {
                bond_idx: 0,
                atom_idx: 9,
                natoms: 1
            }
        ));
        assert!(e.to_string().contains('9'));
    }
    let e = parse_moljson(r#"{"atoms":[{"id":"a","element":"UnknownX"}]}"#)
        .err()
        .unwrap();
    assert!(matches!(&e,MolJsonError::UnknownElement(s) if s=="UnknownX"));
    assert!(e.to_string().contains("UnknownX"));
    let e = parse_ket(r#"{"atoms":[{"label":"UnknownX"}]}"#)
        .err()
        .unwrap();
    assert!(matches!(&e,KetError::UnknownElement(s) if s=="UnknownX"));
    assert!(e.to_string().contains("UnknownX"));
    assert!(
        parse_moljson("{")
            .err()
            .unwrap()
            .to_string()
            .contains("JSON")
    );
    assert!(parse_ket("{").err().unwrap().to_string().contains("JSON"));
}

#[test]
fn gaussian_input_diagnostics_preserve_bad_tokens_and_limits() {
    use crate::gaussian::*;
    for (input, expected) in [
        ("", GaussianError::MissingChargeMultiplicity),
        ("# test\n\ntitle\n\n0 1", GaussianError::NoAtoms),
        (
            "# test\n\ntitle\n\n0 1\nUnknownX 0 0 0",
            GaussianError::UnknownElement("UnknownX".into()),
        ),
        (
            "# test\n\ntitle\n\n0 1\n999 0 0 0",
            GaussianError::UnknownElement("999".into()),
        ),
    ] {
        let e = parse_gjf(input).err().unwrap();
        assert_eq!(e, expected);
        assert!(!e.to_string().is_empty());
    }
    for axis in 0..3 {
        for value in ["bad-token", "NaN", "inf", "-inf"] {
            let mut coords = ["0"; 3];
            coords[axis] = value;
            let text = format!(
                "# test\n\ntitle\n\n0 1\nC {} {} {}",
                coords[0], coords[1], coords[2]
            );
            let e = parse_gjf(&text).err().unwrap();
            assert!(e.to_string().contains(value));
            assert_eq!(
                e,
                if value == "bad-token" {
                    GaussianError::InvalidCoordinate(value.into())
                } else {
                    GaussianError::NonFiniteCoordinate(value.into())
                }
            );
        }
    }
    let text = "# B3LYP\n\nwater\n\n0 1\nO 0 0 0\nH 1 0 0\nH -1 0 0\n";
    for which in 0..5 {
        let mut l = GaussianParseLimits::default();
        match which {
            0 => l.max_input_bytes = text.len() - 1,
            1 => l.max_line_bytes = 1,
            2 => l.max_lines = 1,
            3 => l.max_sections = 2,
            _ => l.max_atoms = 2,
        }
        let e = parse_gjf_with_limits(text, &l).err().unwrap();
        let GaussianError::ResourceLimit {
            resource,
            actual,
            limit,
        } = &e
        else {
            panic!("{e}")
        };
        assert!(actual > limit && e.to_string().contains(resource));
    }
    let e = parse_gaussian_log("no orientation").err().unwrap();
    assert_eq!(e, GaussianError::NoStandardOrientation);
    assert!(e.to_string().contains("Standard orientation"));
}

#[test]
fn smiles_table_diagnostics_keep_physical_lines_and_record_ordinals() {
    use crate::smiles_table::*;
    use std::io::Cursor;
    for strict in [false, true] {
        for (bad, options, fragment) in [
            ("invalid", SmilesReaderOptions::default(), "invalid SMILES"),
            (
                "C",
                SmilesReaderOptions {
                    smiles_column: 2,
                    ..Default::default()
                },
                "column 2",
            ),
            (
                "\"C",
                SmilesReaderOptions {
                    delimiter: Delimiter::Comma,
                    ..Default::default()
                },
                "unterminated",
            ),
            (
                "CCCCCCCCCCCC",
                SmilesReaderOptions {
                    max_line_bytes: 10,
                    ..Default::default()
                },
                "10-byte",
            ),
            (
                "C,n,p",
                SmilesReaderOptions {
                    delimiter: Delimiter::Comma,
                    max_fields: 2,
                    ..Default::default()
                },
                "2",
            ),
        ] {
            let options = SmilesReaderOptions {
                strict_parsing: strict,
                ..options
            };
            let mut reader =
                SmilesRecordReader::new(Cursor::new(format!("# note\n{bad}\nC C C\n")), options);
            let e = reader.next().unwrap().err().unwrap();
            let message = e.to_string();
            assert!(
                message.contains("line 2")
                    && message.contains("record 0")
                    && message.contains(fragment),
                "{message}"
            );
            if strict {
                assert!(reader.next().is_none());
            }
        }
    }
    let mut reader = SmilesRecordReader::new(
        Cursor::new("C\nCC\n"),
        SmilesReaderOptions {
            max_records: 1,
            ..Default::default()
        },
    );
    assert_eq!(reader.next().unwrap().unwrap().mol.atom_count(), 1);
    let e = reader.next().unwrap().err().unwrap();
    assert_eq!(e, SmilesTableError::TooManyRecords { limit: 1 });
    assert!(e.to_string().contains('1'));
    assert!(reader.next().is_none());
}

#[test]
fn tdt_diagnostics_keep_tag_names_and_record_context() {
    use crate::tdt::*;
    use std::io::Cursor;
    for (text, fragment, options) in [
        ("NAME<x>\n|\n", "$SMI", TdtReaderOptions::default()),
        (
            "$SMI<invalid>\n|\n",
            "invalid SMILES",
            TdtReaderOptions::default(),
        ),
        (
            "$SMI<C>\nCUSTOMTAG<oops\n|\n",
            "CUSTOMTAG",
            TdtReaderOptions::default(),
        ),
        (
            "$SMI<C>\n2D<bad,0>\n|\n",
            "2D",
            TdtReaderOptions {
                read_2d: true,
                ..Default::default()
            },
        ),
        (
            "$SMI<C>\n3D<0,0>\n|\n",
            "3D",
            TdtReaderOptions {
                read_3d: true,
                ..Default::default()
            },
        ),
        (
            "$SMI<CCCCCCCC>\n|\n",
            "byte",
            TdtReaderOptions {
                max_line_bytes: 8,
                ..Default::default()
            },
        ),
        (
            "$SMI<C>\nNAME<x>\nEXTRA<y>\n|\n",
            "limit",
            TdtReaderOptions {
                max_tags_per_record: 1,
                ..Default::default()
            },
        ),
    ] {
        let mut reader = TdtRecordReader::new(Cursor::new(text), options);
        let e = reader.next().unwrap().err().unwrap();
        let message = e.to_string();
        assert!(
            message.contains(fragment) && message.contains("record 0"),
            "{message}"
        );
    }
    let mut reader = TdtRecordReader::new(
        Cursor::new("$SMI<C>\n|\n$SMI<CC>\n|\n"),
        TdtReaderOptions {
            max_records: 1,
            ..Default::default()
        },
    );
    assert!(reader.next().unwrap().is_ok());
    let e = reader.next().unwrap().err().unwrap();
    assert_eq!(e, TdtError::TooManyRecords { limit: 1 });
    assert!(e.to_string().contains('1'));
}

#[test]
fn pqr_field_errors_keep_the_line_and_original_token() {
    use crate::pqr::*;
    let fields = ["ATOM", "17", "C", "LIG", "23", "0", "0", "0", "0.1", "1.5"];
    assert!(parse_pqr(&fields.join(" ")).is_ok());
    for column in [1, 4, 5, 6, 7, 8, 9] {
        for bad in ["bad-token", "NaN", "inf", "-inf"] {
            let mut row = fields;
            row[column] = bad;
            let e = parse_pqr(&format!("REMARK test\n{}", row.join(" "))).unwrap_err();
            let message = e.to_string();
            assert!(
                message.contains("line 2") && message.contains(bad),
                "{message}"
            );
            assert!(matches!(
                e,
                PqrError::InvalidField { line: 2, .. } | PqrError::NonFiniteValue { line: 2, .. }
            ));
        }
    }
    let e = parse_pqr("ATOM 1").unwrap_err();
    assert_eq!(e, PqrError::WrongFieldCount { line: 1, found: 2 });
    assert!(e.to_string().contains("found 2"));
    let mut row = fields;
    row[2] = "123";
    let e = parse_pqr(&row.join(" ")).unwrap_err();
    assert!(matches!(&e,PqrError::UnresolvableElement{line:1,atom_name} if atom_name=="123"));
    assert!(e.to_string().contains("123"));
    for which in 0..3 {
        let mut l = PqrParseLimits::default();
        match which {
            0 => l.max_input_bytes = 1,
            1 => l.max_line_len = 1,
            _ => l.max_atoms = 0,
        }
        let e = parse_pqr_with_limits(&fields.join(" "), &l).unwrap_err();
        assert!(e.to_string().contains("limit"));
    }
    let e = parse_pqr("REMARK only").unwrap_err();
    assert_eq!(e, PqrError::NoAtomRecords);
    assert!(e.to_string().contains("ATOM/HETATM"));
}

#[test]
fn mmcif_missing_columns_and_invalid_numbers_are_actionable() {
    use crate::mmcif::*;
    let headers = [
        "type_symbol",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        "id",
        "auth_seq_id",
        "pdbx_PDB_model_num",
        "pdbx_formal_charge",
        "occupancy",
        "B_iso_or_equiv",
    ];
    let values = ["C", "0", "0", "0", "17", "23", "1", "0", "1", "0"];
    let make = |headers: &[&str], values: &[&str]| {
        format!(
            "data_test\nloop_\n{}\n{}\n",
            headers
                .iter()
                .map(|s| format!("_atom_site.{s}"))
                .collect::<Vec<_>>()
                .join("\n"),
            values.join(" ")
        )
    };
    assert!(parse_mmcif(&make(&headers, &values)).is_ok());
    for missing in 0..4 {
        let mut hs = headers.to_vec();
        hs.remove(missing);
        let mut vs = values.to_vec();
        vs.remove(missing);
        let e = parse_mmcif(&make(&hs, &vs)).unwrap_err();
        assert_eq!(
            e,
            if missing == 0 {
                MmcifError::MissingTypeSymbolColumn
            } else {
                MmcifError::MissingCoordinateColumns
            }
        );
        assert!(
            e.to_string()
                .contains(if missing == 0 { "type_symbol" } else { "Cartn" })
        );
    }
    for column in 1..headers.len() {
        let mut vs = values;
        vs[column] = "bad-token";
        let e = parse_mmcif(&make(&headers, &vs)).unwrap_err();
        assert!(e.to_string().contains("bad-token"));
        assert!(
            e.to_string()
                .to_lowercase()
                .contains(&headers[column].to_lowercase()),
            "{e}"
        );
    }
    let mut vs = values;
    vs[0] = "UnknownX";
    let e = parse_mmcif(&make(&headers, &vs)).unwrap_err();
    assert!(matches!(e, MmcifError::UnknownElement(_)));
    assert!(e.to_string().to_lowercase().contains("unknownx"));
    let text = make(&headers, &values);
    for which in 0..3 {
        let mut l = MmcifParseLimits::default();
        match which {
            0 => l.max_input_bytes = 1,
            1 => l.max_line_len = 1,
            _ => l.max_atoms = 0,
        }
        let e = parse_mmcif_with_limits(&text, &l).unwrap_err();
        assert!(e.to_string().contains("limit"));
    }
    let e = parse_mmcif("data_empty\n").unwrap_err();
    assert_eq!(e, MmcifError::NoAtomSiteLoop);
    assert!(e.to_string().contains("_atom_site"));
}

#[test]
fn xyz_and_extxyz_diagnostics_identify_bad_coordinates_and_properties() {
    use crate::xyz::*;
    for (text, fragment) in [
        ("bad-count\ncomment\n", "bad-count"),
        ("1\ncomment\nUnknownX 0 0 0\n", "UnknownX"),
        ("1\ncomment\nC 0 0\n", "line 3"),
    ] {
        assert!(parse_xyz(text).unwrap_err().to_string().contains(fragment));
    }
    for axis in 0..3 {
        for bad in ["bad-token", "NaN", "inf", "-inf"] {
            let mut xyz = ["0"; 3];
            xyz[axis] = bad;
            let text = format!("1\ncomment\nC {} {} {}\n", xyz[0], xyz[1], xyz[2]);
            let e = parse_xyz(&text).unwrap_err();
            assert!(e.to_string().contains("line 3") && e.to_string().contains(bad));
        }
    }
    for (info, row, fragment) in [
        ("=value", "C 0 0 0", "info"),
        ("name=\"unterminated", "C 0 0 0", "name"),
        ("name=a name=b", "C 0 0 0", "name"),
        ("Lattice=\"1 2\"", "C 0 0 0", "Lattice"),
        ("Properties=bad", "C 0 0 0", "Properties"),
        (
            "Properties=species:S:1:pos:R:3:charge:R:1",
            "C 0 0 0",
            "column",
        ),
        (
            "Properties=species:S:1:pos:R:3:charge:R:1",
            "C 0 0 0 bad-token",
            "charge",
        ),
        (
            "Properties=species:S:1:pos:R:3:charge:R:1",
            "C 0 0 0 NaN",
            "charge",
        ),
        (
            "Properties=species:S:1:pos:R:3:count:I:1",
            "C 0 0 0 bad-token",
            "count",
        ),
        (
            "Properties=species:S:1:pos:R:3:flag:L:1",
            "C 0 0 0 bad-token",
            "flag",
        ),
    ] {
        let e = parse_extxyz(&format!("1\n{info}\n{row}\n")).unwrap_err();
        assert!(e.to_string().contains(fragment), "{e}");
    }
    let frame = parse_xyz("1\ncomment\nC 0 0 0\n").unwrap();
    let mut bad = frame;
    bad.info.push(("invalid key".into(), "value".into()));
    let e = write_extxyz(&bad).unwrap_err();
    assert!(matches!(e, XyzError::UnwritableMetadata { .. }));
    assert!(e.to_string().contains("key"));
}

#[test]
fn conformer_mol_writer_preserves_3d_header_and_coordinate_precision() {
    use crate::mol2000::*;
    use chematic_core::{Coords3D, Point3};
    for source in ["CCO", "[13CH3+]C#N", "c1ccccc1"] {
        let mol = chematic_smiles::parse(source).unwrap();
        let conformer = Coords3D {
            points: (0..mol.atom_count())
                .map(|i| Point3::new(i as f64 * 1.25, -0.5, 0.75))
                .collect(),
        };
        let metadata = MolMetadata {
            name: "3d sample".into(),
            comment: "round trip".into(),
            ..Default::default()
        };
        let text = write_mol_with_conformer(&mol, &metadata, &conformer);
        let read = read_mol_with_diagnostics(&text).unwrap();
        assert_eq!(read.mol.atom_count(), mol.atom_count());
        assert_eq!(read.coordinate_dimension, CoordinateDimension::ThreeD);
        assert_eq!(read.conformer.unwrap(), conformer);
        assert_eq!(read.metadata.name, metadata.name);
    }
}

#[test]
fn orca_external_and_internal_coordinates_round_trip_without_interpreting_the_payload() {
    for source in [
        "! HF\n* xyzfile -1 2 external.xyz\n",
        "! HF\n* gzmtfile 1 1 external.gzmt\n",
        "! HF\n* int 0 1\nC 0 0 0 0 0 0\nH 1 0 0 1.0 0 0\n*\n",
        "! HF\n* int 0 1\n*\n",
        "! HF\n* xyz 0 1\nC 0.1$ -0.2 0.3$ fragment\n*\n",
    ] {
        let first = parse_orca_input(source).unwrap();
        let written = write_orca_input(&first);
        assert_eq!(parse_orca_input(&written).unwrap(), first);
    }
}

#[test]
fn rxn_file_failures_keep_component_counts_and_nested_mol_context() {
    use crate::rxn::*;
    let rxn = chematic_rxn::parse_reaction("CC.O>>CO").unwrap();
    let text = write_rxn_file(&rxn);
    assert!(parse_rxn_file(&text).is_ok());
    for which in 0..4 {
        let mut l = RxnFileParseLimits::default();
        match which {
            0 => l.max_input_bytes = text.len() - 1,
            1 => l.max_reactants = 1,
            2 => l.max_products = 0,
            _ => l.max_molecules = 2,
        }
        let e = parse_rxn_file_with_limits(&text, l).err().unwrap();
        let RxnParseError::ResourceLimit {
            resource,
            actual,
            limit,
        } = &e
        else {
            panic!("{e}")
        };
        assert!(actual > limit && e.to_string().contains(resource));
    }
    for (source, context) in [
        ("", "$RXN"),
        ("$RXN\nname\nprogram\ncomment\nxxx\n", "count line"),
        ("$RXN\nname\nprogram\ncomment\n  1  1\n", "declares 2"),
        (
            "$RXN\nname\nprogram\ncomment\n  1  0\n$MOL\ninvalid\n",
            "MOL parse error",
        ),
    ] {
        let e = parse_rxn_file(source).err().unwrap();
        assert!(e.to_string().contains(context), "{e}");
        let e = parse_rxn_document(source).err().unwrap();
        assert!(
            e.to_string().contains("RXN document parse error") && e.to_string().contains(context)
        );
    }
}

#[test]
fn sdf_string_readers_handle_terminal_delimiters_fields_and_each_record_budget() {
    use crate::sdf::*;
    let mol = chematic_smiles::parse("CO").unwrap();
    let block = crate::write_mol(&mol, &Default::default());
    for suffix in ["$$$$", "$$$$\r", "$$$$\n", ""] {
        let text = format!("$$$$\n{block}> <first>\none\n> <second>\ntwo\nthree\n{suffix}");
        assert_eq!(SdfReader::new(&text).count(), 1);
        let record = SdfRecordReader::new(&text).next().unwrap().unwrap();
        assert_eq!(record.mol.atom_count(), 2);
        assert_eq!(record.properties["first"], "one");
        assert_eq!(record.properties["second"], "two\nthree");
    }
    for (text, limits, resource) in [
        (
            format!("{block}$$$$\n"),
            SdfParseLimits {
                max_line_bytes: 3,
                ..Default::default()
            },
            "line bytes",
        ),
        (
            "123456789".into(),
            SdfParseLimits {
                max_line_bytes: 8,
                ..Default::default()
            },
            "line bytes",
        ),
        (
            block.clone(),
            SdfParseLimits {
                max_record_bytes: block.len() - 1,
                ..Default::default()
            },
            "record bytes",
        ),
        (
            format!("{block}$$$$\n{block}"),
            SdfParseLimits {
                max_records: 1,
                ..Default::default()
            },
            "records",
        ),
    ] {
        let e = SdfRecordReader::with_limits(&text, limits)
            .find_map(Result::err)
            .unwrap();
        assert!(matches!(e, crate::MolParseError::ResourceLimit { .. }));
        assert!(e.to_string().contains(resource), "{e}");
    }
    assert!(
        SdfRecordReader::new("broken\n$$$$")
            .next()
            .unwrap()
            .is_err()
    );
}
