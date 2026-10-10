use crate::{cdxml, cjson, lammps_data, mrv, qcschema};
use chematic_core::{BondOrder, Coords3D};

#[test]
fn qcschema_raw_conversion_errors_preserve_the_input_and_report_indices() {
    let base =
        qcschema::parse_qcschema_molecule(r#"{"symbols":["C","C"],"geometry":[0,0,0,2,0,0]}"#)
            .unwrap();
    let mut unknown = base.clone();
    unknown.symbols[1] = "Xx".into();
    let before = unknown.clone();
    let error = qcschema::qc_molecule_to_chematic(&unknown)
        .err()
        .expect("unknown symbol");
    assert_eq!(error, qcschema::QcConvertError::UnknownElement("Xx".into()));
    assert_eq!(
        error.to_string(),
        "QCSchema conversion: unknown element symbol 'Xx'"
    );
    assert_eq!(unknown, before);
    for pair in [(2, 0), (0, 3)] {
        let mut qc = base.clone();
        qc.connectivity = Some(vec![(pair.0, pair.1, 1.)]);
        let before = qc.clone();
        let error = qcschema::qc_molecule_to_chematic(&qc)
            .err()
            .expect("bad index");
        assert_eq!(
            error,
            qcschema::QcConvertError::InvalidBondIndex {
                index: pair.0.max(pair.1),
                atom_count: 2
            }
        );
        assert!(error.to_string().contains("but there are only 2 atoms"));
        assert_eq!(qc, before);
    }
    let mol = chematic_smiles::parse("CC").unwrap();
    let error =
        qcschema::chematic_to_qc_molecule(&mol, &Coords3D::new_zeroed(1), 0., 1).unwrap_err();
    assert_eq!(
        error,
        qcschema::QcConvertError::AtomCountMismatch {
            molecule: 2,
            coords: 1
        }
    );
    assert_eq!(
        error.to_string(),
        "QCSchema conversion: molecule has 2 atoms but coords has 1"
    );
}

#[test]
fn qcschema_numeric_bond_orders_round_trip_with_bohr_coordinates() {
    for (order, expected) in [
        (0., BondOrder::Zero),
        (1.5, BondOrder::Aromatic),
        (4., BondOrder::Quadruple),
    ] {
        let mut qc=qcschema::parse_qcschema_molecule(r#"{"symbols":["C","C"],"geometry":[0,0,0,2,0,0],"molecular_charge":-2,"molecular_multiplicity":3}"#).unwrap();
        qc.connectivity = Some(vec![(0, 1, order)]);
        let before = qc.clone();
        let view = qcschema::qc_molecule_to_chematic(&qc).unwrap();
        assert_eq!(
            view.molecule.bond(chematic_core::BondIdx(0)).order,
            expected
        );
        let output = qcschema::chematic_to_qc_molecule(
            &view.molecule,
            &view.coords,
            view.molecular_charge,
            view.molecular_multiplicity,
        )
        .unwrap();
        assert_eq!(output.connectivity, qc.connectivity);
        assert_eq!(output.symbols, qc.symbols);
        assert!(
            output
                .geometry
                .iter()
                .zip(&qc.geometry)
                .all(|(a, b)| (a - b).abs() < 1e-12)
        );
        assert_eq!(output.molecular_charge, -2.);
        assert_eq!(output.molecular_multiplicity, 3);
        assert_eq!(qc, before);
    }
}

#[test]
fn cjson_default_bonds_negative_indices_and_string_budgets_are_explicit() {
    let input = serde_json::json!({"atoms":{"elements":{"number":[6,8]},"coords":{"3d":[0,0,0,1.5,0,0]}},"bonds":{"connections":{"index":[0,1]}}});
    let (mol, coords) = cjson::parse_cjson(&input.to_string()).unwrap();
    assert_eq!(mol.bond_count(), 1);
    assert_eq!(mol.bond(chematic_core::BondIdx(0)).order, BondOrder::Single);
    assert_eq!(coords.len(), 2);
    for indices in [vec![-1, 1], vec![0, -1]] {
        let mut raw = input.clone();
        raw["bonds"]["connections"]["index"] = serde_json::json!(indices);
        let error = cjson::parse_cjson(&raw.to_string())
            .err()
            .expect("negative atom index");
        assert!(matches!(
            error,
            cjson::CjsonError::InvalidValue {
                field: "bonds.connections.index"
            }
        ));
        assert!(error.to_string().contains("bonds.connections.index"));
    }
    for (raw, message) in [("{", "invalid JSON"), ("{}", "missing required field")] {
        assert!(
            cjson::parse_cjson(raw)
                .err()
                .unwrap()
                .to_string()
                .contains(message)
        );
    }
    let mut raw = input;
    raw["name"] = serde_json::json!("abcdefghijkl");
    let error = cjson::parse_cjson_with_limits(
        &raw.to_string(),
        &cjson::CjsonParseLimits {
            max_string_bytes: 11,
            ..Default::default()
        },
    )
    .err()
    .unwrap();
    assert!(matches!(
        error,
        cjson::CjsonError::ResourceLimit {
            resource: "string bytes",
            actual: 12,
            limit: 11
        }
    ));
}

#[test]
fn lammps_velocities_and_image_flags_survive_writer_round_trip() {
    let input = "Example\n\n1 atoms\n0 bonds\n1 atom types\n\n0 10 xlo xhi\n0 10 ylo yhi\n0 10 zlo zhi\n\nMasses\n\n1 12.011\n\nAtoms # atomic\n\n1 1 0.5 1.5 2.5 -1 2 3\n\nVelocities\n\n1 0.125 -0.25 0.5\n";
    let data = lammps_data::parse_lammps_data(input, lammps_data::LammpsAtomStyle::Atomic).unwrap();
    assert_eq!(data.atoms[0].image, Some([-1, 2, 3]));
    assert_eq!(data.velocities.len(), 1);
    let serialized = lammps_data::write_lammps_data(&data);
    let read =
        lammps_data::parse_lammps_data(&serialized, lammps_data::LammpsAtomStyle::Atomic).unwrap();
    assert_eq!(read, data);
    let invalid = input.replace("-1 2 3", "2147483648 2 3");
    let error =
        lammps_data::parse_lammps_data(&invalid, lammps_data::LammpsAtomStyle::Atomic).unwrap_err();
    assert!(matches!(
        error,
        lammps_data::LammpsDataError::MalformedRow { .. }
    ));
    assert!(error.to_string().contains("Atoms"));
    for raw in [
        "Example\n",
        "Example\n0 10 xlo xhi\n0 10 ylo yhi\n0 10 zlo zhi\nAtoms",
    ] {
        let error =
            lammps_data::parse_lammps_data(raw, lammps_data::LammpsAtomStyle::Atomic).unwrap_err();
        assert!(matches!(
            error,
            lammps_data::LammpsDataError::TruncatedInput { .. }
        ));
        assert!(error.to_string().contains("unexpected end of input"));
    }
    let error =
        lammps_data::parse_lammps_data(input, lammps_data::LammpsAtomStyle::Other("sphere".into()))
            .unwrap_err();
    assert!(
        error
            .to_string()
            .contains("unsupported LAMMPS atom_style 'sphere'")
    );
}

#[test]
fn cdxml_empty_fragments_default_coordinates_and_missing_pseudoatom_names_are_explicit() {
    let (mol, coords) = cdxml::parse_cdxml("<CDXML><page id=\"1\"/></CDXML>").unwrap();
    assert_eq!(mol.atom_count(), 0);
    assert!(coords.is_empty());
    let input = "<CDXML><page id=\"1\"><fragment id=\"2\"><n id=\"3\" Element=\"6\"/></fragment></page></CDXML>";
    let all =
        cdxml::parse_cdxml_all_with_options(input, &cdxml::CdxmlParseOptions::default()).unwrap();
    assert_eq!(all.len(), 1);
    assert_eq!(all[0].0.atom_count(), 1);
    assert_eq!(all[0].1, vec![(0., 0.)]);
    let expanded = input.replace("><", ">\n<");
    let (expanded_mol, expanded_coords) = cdxml::parse_cdxml(&expanded).unwrap();
    assert_eq!(expanded_mol.atom_count(), all[0].0.atom_count());
    assert_eq!(expanded_coords, all[0].1);
    let invalid = input.replace("Element=\"6\"", "NodeType=\"GenericNickname\"");
    let error = cdxml::parse_cdxml(&invalid).err().unwrap();
    assert!(
        error
            .to_string()
            .contains("missing GenericNickname attribute")
    );
}

#[test]
fn mrv_missing_atom_refs_and_escaped_stereo_values_return_actionable_errors() {
    let body = "<cml><MDocument><MChemicalStruct><molecule><atomArray><atom id=\"a1\" elementType=\"C\"/><atom id=\"a2\" elementType=\"O\"/></atomArray><bondArray><bond id=\"b1\" order=\"1\"/></bondArray></molecule></MChemicalStruct></MDocument></cml>";
    let error = mrv::parse_mrv(body).err().unwrap();
    assert!(error.to_string().contains("missing atomRefs2"));
    for entity in ["&amp;", "&lt;", "&gt;", "&quot;", "&apos;"] {
        let input = body.replace(
            "order=\"1\"/>",
            &format!("order=\"1\" atomRefs2=\"a1 a2\"><bondStereo>{entity}</bondStereo></bond>"),
        );
        let error = mrv::parse_mrv(&input).err().unwrap();
        assert!(matches!(error, mrv::MrvError::InvalidBondStereo { .. }));
        assert!(!error.to_string().contains(entity));
    }
    assert!(
        mrv::parse_mrv("<wrong/>")
            .err()
            .unwrap()
            .to_string()
            .contains("no <cml> root element")
    );
}

#[test]
fn mrv_writer_reports_unkekulizable_aromatic_input_without_changing_it() {
    let mut mol = chematic_smiles::parse("C1CCCC1").unwrap();
    for i in 0..mol.atom_count() {
        mol.set_atom_aromatic(chematic_core::AtomIdx(i as u32), true);
    }
    for i in 0..mol.bond_count() {
        mol.set_bond_order(chematic_core::BondIdx(i as u32), BondOrder::Aromatic);
    }
    let record = crate::record::MoleculeRecord::new(mol);
    let before: Vec<_> = record.mol.bonds().map(|(_, b)| b.order).collect();
    let error = mrv::write_mrv(
        &record,
        &mrv::MrvWriteOptions {
            precision: 0,
            ..Default::default()
        },
    )
    .unwrap_err();
    assert!(matches!(error, mrv::MrvError::KekulizationFailed { .. }));
    assert!(error.to_string().contains("kekulization failed"));
    assert_eq!(
        record.mol.bonds().map(|(_, b)| b.order).collect::<Vec<_>>(),
        before
    );
}

#[test]
fn cdxml_minification_preserves_fragment_order_attributes_and_physical_budgets() {
    let input = "<CDXML><page id='1'><fragment id='2'><n id='3' Element='6' p='1 2' title='a > b'/><n id='4' Element='7' p='3 4'/><b B='3' E='4'/></fragment><fragment id='5'><n id='6' Element='8' p='5 6'/></fragment></page></CDXML>";
    let options = cdxml::CdxmlParseOptions::default();
    let limits = cdxml::CdxmlParseLimits {
        max_lines: 1,
        max_line_bytes: input.len(),
        ..Default::default()
    };
    let all = cdxml::parse_cdxml_all_with_options_and_limits(input, &options, &limits).unwrap();
    assert_eq!(
        all.iter().map(|(m, _)| m.atom_count()).collect::<Vec<_>>(),
        vec![2, 1]
    );
    assert_eq!(all[0].1, vec![(1., 2.), (3., 4.)]);
    assert_eq!(all[1].1, vec![(5., 6.)]);
    assert_eq!(all[0].0.bond_count(), 1);
    assert_eq!(
        all[0].0.atom(chematic_core::AtomIdx(0)).element,
        chematic_core::Element::C
    );
    let first = cdxml::parse_cdxml(input).unwrap();
    assert_eq!(first.0.atom_count(), 2);
    assert_eq!(first.1, all[0].1);
    let error = cdxml::parse_cdxml_all_with_limits(
        input,
        &cdxml::CdxmlParseLimits {
            max_line_bytes: input.len() - 1,
            ..limits
        },
    )
    .err()
    .unwrap();
    assert!(matches!(
        error,
        cdxml::CdxmlError::ResourceLimit {
            resource: "line bytes",
            ..
        }
    ));
    let pretty = input.replace("><", ">\n<");
    let pretty_all = cdxml::parse_cdxml_all(&pretty).unwrap();
    assert_eq!(
        pretty_all
            .iter()
            .map(|(m, _)| m.atom_count())
            .collect::<Vec<_>>(),
        vec![2, 1]
    );
    for ((a, ac), (b, bc)) in all.iter().zip(&pretty_all) {
        assert_eq!(
            chematic_smiles::canonical_smiles(a),
            chematic_smiles::canonical_smiles(b)
        );
        assert_eq!(ac, bc);
    }
    for source in ["C#N", "c1ccccc1", "[13CH3][NH3+]"] {
        let mol = chematic_smiles::parse(source).unwrap();
        let coords: Vec<_> = (0..mol.atom_count())
            .map(|i| (i as f64 * 1.5, (i % 2) as f64))
            .collect();
        let written = cdxml::write_cdxml(&mol, &coords);
        let (pretty_round_trip, _) = cdxml::parse_cdxml(&written).unwrap();
        let minified = written.replace('\n', "");
        let (round_trip, read_coords) = cdxml::parse_cdxml(&minified).unwrap();
        assert_eq!(round_trip.atom_count(), mol.atom_count());
        assert_eq!(round_trip.bond_count(), mol.bond_count());
        assert_eq!(read_coords, coords);
        assert_eq!(
            chematic_smiles::canonical_smiles(&round_trip),
            chematic_smiles::canonical_smiles(&pretty_round_trip),
            "{source}"
        );
    }
}

#[test]
fn cdxml_comments_cdata_and_multiline_attributes_do_not_create_phantom_atoms() {
    let input = "<?xml version='1.0'?>\n<CDXML><page id='1'><fragment id='2'><!-- x > <n id='99' Element='8'/> -->\n<![CDATA[x > <n id='98' Element='8'/>]]>\n<?ignore x > <n id='97' Element='8'/> ?>\n<n id='3'\n Element='6' p='1 2' title='a > b'/><n id='4' Element='8' p='3 4'/><b B='3' E='4'/></fragment></page></CDXML>";
    let (mol, coords) = cdxml::parse_cdxml(input).unwrap();
    assert_eq!(mol.atom_count(), 2);
    assert_eq!(mol.bond_count(), 1);
    assert_eq!(coords, vec![(1., 2.), (3., 4.)]);
    let document = crate::cdxml_document::CdxmlDocument::parse(input).unwrap();
    let atoms = document
        .pages
        .iter()
        .flat_map(|p| &p.children)
        .filter(|o| o.kind() == crate::cdxml_document::CdxmlObjectKind::Atom)
        .count();
    assert_eq!(atoms, 2);
    assert_eq!(document.write(), input);
}

#[test]
fn mol2_boundary_charge_edits_keep_status_bits_and_opaque_metadata() {
    use crate::mol2_tripos::*;
    use chematic_core::AtomIdx;
    let input = "@<TRIPOS>MOLECULE\nions\n3 1 0 0 0\nSMALL\nUSER_CHARGES\nmetadata tail\n@<TRIPOS>ATOM\n10 N 0 0 0 N.4 1 ION 0.25 BACKBONE\n20 C 1 0 0 C.3 1 ION -0.25 DICT\n30 Cl 4 0 0 Cl 2 SALT 0\n@<TRIPOS>BOND\n7 10 20 1 DICT\n@<TRIPOS>UNITY_ATOM_ATTR\n10 1\nlabel ammonium\n99 1\nlabel detached\n@<TRIPOS>COMMENT\nkeep this extension\n";
    let mut record = parse_mol2_record(input).unwrap();
    record.molecule.set_charge(AtomIdx(0), 1);
    record.molecule.set_charge(AtomIdx(2), -1);
    let written = write_mol2_record(&record);
    let restored = parse_mol2_record(&written).unwrap();
    assert_eq!(restored.molecule.atom(AtomIdx(0)).charge, 1);
    assert_eq!(restored.molecule.atom(AtomIdx(2)).charge, -1);
    assert_eq!(restored.atoms, record.atoms);
    assert_eq!(restored.bonds, record.bonds);
    assert_eq!(restored.molecule_tail, record.molecule_tail);
    assert_eq!(restored.opaque_sections, record.opaque_sections);
    assert_eq!(
        restored.unity_atom_attributes[0].attributes,
        vec![
            ("label".into(), "ammonium".into()),
            ("charge".into(), "1".into())
        ]
    );
    assert_eq!(
        restored.unity_atom_attributes[1],
        record.unity_atom_attributes[1]
    );
    assert_eq!(write_mol2_record(&restored), written);
}

#[test]
fn mol2_boundary_malformed_attribute_records_report_the_actual_field() {
    use crate::mol2_tripos::*;
    let base = "@<TRIPOS>MOLECULE\nm\n1 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C 0 0 0 C.3\n";
    for (body, message) in [
        ("1", "header must"),
        ("x 1\ncharge 1", "UNITY atom_id"),
        ("1 x\ncharge 1", "UNITY attribute count"),
        ("1 2\ncharge 1", "ends before"),
        ("1 1\ncharge invalid", "formal charge"),
    ] {
        let error = parse_mol2_record(&format!("{base}@<TRIPOS>UNITY_ATOM_ATTR\n{body}\n"))
            .err()
            .unwrap();
        assert!(matches!(error, Mol2Error::InvalidAtomLine { .. }));
        assert!(error.to_string().contains(message), "{error}");
    }
    let error = parse_mol2_record("@<TRIPOS>ATOM\n1 C 0 0 0 C.3\n")
        .err()
        .unwrap();
    assert_eq!(error, Mol2Error::MissingSection("MOLECULE".into()));
    assert!(error.to_string().contains("MOLECULE"));
    let error = parse_mol2_record(&format!("{base}{base}")).err().unwrap();
    assert_eq!(error, Mol2Error::MultipleMolecules { count: 2 });
    assert!(error.to_string().contains("2 MOLECULE"));
    let error = parse_mol2_record(&base.replace("C.3", "C.3 1 LIG wrong"))
        .err()
        .unwrap();
    assert!(error.to_string().contains("partial charge"));
    let two = base.replace("1 0\n", "2 1\n");
    let error = parse_mol2_record(&format!("{two}2 C 1 0 0 C.3\n@<TRIPOS>BOND\nx 1 2 1\n"))
        .err()
        .unwrap();
    assert!(matches!(error, Mol2Error::InvalidBondLine { .. }));
    assert!(error.to_string().contains("bond_id"));
}

#[test]
fn xyz_boundary_typed_properties_and_terminal_stream_errors() {
    use crate::xyz::*;
    let input = "2\nProperties=species:S:1:pos:R:3:label:S:1:fixed:L:1\nC 0 0 0 carbon T\nO 1 0 0 oxygen F\n";
    let frame = parse_extxyz(input).unwrap();
    assert_eq!(
        frame.properties[0].values,
        vec![
            vec![XyzValue::Str("carbon".into())],
            vec![XyzValue::Str("oxygen".into())]
        ]
    );
    assert_eq!(
        frame.properties[1].values,
        vec![
            vec![XyzValue::Logical(true)],
            vec![XyzValue::Logical(false)]
        ]
    );
    let written = write_extxyz(&frame).unwrap();
    let restored = parse_extxyz(&written).unwrap();
    assert_eq!(restored.atoms, frame.atoms);
    assert_eq!(restored.properties, frame.properties);
    for bad in ["wrong\ncomment\n", "2\ncomment\nC 0 0 0\n"] {
        let mut reader = ExtxyzReader::new(bad);
        let error = reader.next().unwrap().unwrap_err();
        assert!(!error.to_string().is_empty());
        assert!(reader.next().is_none());
    }
    let mut plain = XyzReader::new("1\ncomment\nC wrong 0 0\n1\ncomment\nC 0 0 0\n");
    assert!(plain.next().unwrap().is_err());
    assert!(plain.next().is_none());
    let error = parse_extxyz("1").unwrap_err();
    assert!(matches!(error, XyzError::MissingCommentLine));
    assert_eq!(error.to_string(), "missing comment line");
}

#[test]
fn xyz_boundary_batches_finish_and_cancel_without_reading_another_frame() {
    use crate::xyz::*;
    use std::io::Cursor;
    let input = "1\ncomment\nC 0 0 0\n";
    let mut xyz = XyzBatchReader::new(Cursor::new(input), 1);
    assert!(!xyz.is_cancelled());
    let batch = xyz.next().unwrap();
    assert!(!batch.is_empty());
    assert_eq!(batch.len(), 1);
    assert!(xyz.next().is_none());
    assert_eq!(xyz.progress().status, "complete");
    let mut ext = ExtxyzBatchReader::new(Cursor::new(input), 1);
    assert!(!ext.is_cancelled());
    assert_eq!(ext.progress().status, "running");
    let batch = ext.next().unwrap();
    assert!(!batch.is_empty());
    assert_eq!(batch.len(), 1);
    assert!(ext.next().is_none());
    assert!(ext.next().is_none());
    assert_eq!(ext.progress().status, "complete");
    let mut cancelled = ExtxyzBatchReader::new(Cursor::new("invalid"), 2);
    cancelled.cancel();
    assert!(cancelled.is_cancelled());
    assert!(cancelled.next().is_none());
    assert_eq!(cancelled.progress().frames_emitted, 0);
    assert_eq!(cancelled.progress().status, "cancelled");
}

#[test]
fn v3000_boundary_query_bonds_and_rgroup_spellings_round_trip() {
    use crate::mol3000::*;
    use chematic_core::{AtomIdx, BondIdx};
    let base = "query\n  chematic          2D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 2 C 1 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n";
    // These are query/interchange records, not force-field molecules.
    for (code, order) in [
        (0, BondOrder::Zero),
        (2, BondOrder::Double),
        (3, BondOrder::Triple),
        (5, BondOrder::QuerySingleOrDouble),
        (6, BondOrder::QuerySingleOrAromatic),
        (7, BondOrder::QueryDoubleOrAromatic),
        (8, BondOrder::QueryAny),
    ] {
        let input = base.replace("M  V30 1 1 1 2", &format!("M  V30 1 {code} 1 2"));
        let (mol, metadata) = parse_mol_v3000(&input).unwrap();
        assert_eq!(mol.bond(BondIdx(0)).order, order);
        let mut conformer = Coords3D::new_zeroed(2);
        conformer.set(AtomIdx(1), chematic_core::Point3::new(1., 0., 0.));
        for written in [
            write_mol_v3000(&mol, &metadata, &[(0., 0.), (1., 0.)]),
            write_mol_v3000_with_conformer(&mol, &metadata, &conformer),
        ] {
            let (restored, _) = parse_mol_v3000(&written).unwrap();
            assert_eq!(restored.bond(BondIdx(0)).order, order);
            assert_eq!(restored.atom_count(), 2);
        }
    }
    for (symbol, number) in [
        ("R#", None),
        ("R7", Some(7)),
        ("C 0 0 0 0 RGROUPS=(1 7)", Some(7)),
    ] {
        let line = if symbol.starts_with("C ") {
            format!("M  V30 1 {symbol}")
        } else {
            format!("M  V30 1 {symbol} 0 0 0 0")
        };
        let input = base.replace("M  V30 1 C 0 0 0 0", &line);
        let (mol, metadata) = parse_mol_v3000(&input).unwrap();
        assert!(mol.atom(AtomIdx(0)).wildcard);
        assert_eq!(mol.r_group_label(AtomIdx(0)).unwrap().number(), number);
        let written = write_mol_v3000(&mol, &metadata, &[(0., 0.), (1., 0.)]);
        let (restored, _) = parse_mol_v3000(&written).unwrap();
        assert_eq!(
            restored.r_group_label(AtomIdx(0)),
            mol.r_group_label(AtomIdx(0))
        );
    }
    for (old, new, message) in [
        ("V3000", "V2000", "version tag"),
        ("C 1 0 0 0", "C 1 nope 0 0", "y coordinate"),
        ("C 1 0 0 0", "C 1 NaN 0 0", "not finite"),
        ("M  V30 END BOND\nM  V30 END CTAB\n", "", "END BOND"),
    ] {
        let error = parse_mol_v3000(&base.replace(old, new)).err().unwrap();
        assert!(error.to_string().contains(message), "{error}");
    }
    for kind in ["EXT", "VENDOR"] {
        let group = parse_v3000_sgroup_line(&format!("1 {kind} 0 ATOMS=(1 1)")).unwrap();
        assert_eq!(group.atom_ids, vec![1]);
    }
    assert!(
        parse_v3000_sgroup_line("1 EXT")
            .unwrap_err()
            .to_string()
            .contains("needs id")
    );
}
