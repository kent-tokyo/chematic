use super::*;

#[test]
fn single_hash_centre_with_unequal_flat_bonds_reads_rdkit_parity() {
    // RDKit's layout of an artemisinin analogue (ChEMBL): the methyl
    // carbon's two O bonds are drawn 6 degrees apart with different
    // lengths. Its parity is the hash's sign and the flat bonds' angular
    // order; the triangle of their end points had turned over (#734).
    let block = r"
     RDKit          2D

 20 23  0  0  0  0  0  0  0  0999 V2000
    1.9033   -3.6288    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7339   -2.6894    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6643   -3.2324    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8337   -2.2930    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6049   -0.8105    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7743    0.1289    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1725   -0.4141    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5454    1.6113    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7148    2.5507    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1471    2.1543    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0223    1.2149    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4919    0.9148    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.4360    1.8245    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3755    2.9939    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0820    0.4707    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4264   -0.8784    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9628   -1.2069    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2066   -0.2675    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0171    0.7012    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2123    0.8558    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  6
  2  3  1  0
  3  4  1  0
  5  4  1  1
  5  6  1  0
  6  7  1  1
  6  8  1  0
  8  9  2  0
  8 10  1  0
 11 10  1  6
 11 12  1  0
 12 13  1  0
 13 14  1  6
 13 15  1  0
 15 16  1  0
 17 16  1  1
 17 18  1  0
 18 19  1  0
 19 20  1  0
 17  2  1  0
 18  5  1  0
 18 11  1  0
 20 13  1  0
M  END
";
    let (mol, _) = parse_mol(block).unwrap();
    let want =
        chematic_smiles::parse("C[C@@H]1CC[C@H]2[C@@H](C)C(=O)O[C@@H]3O[C@]4(C)CC[C@@H]1C32OO4")
            .unwrap();
    assert_eq!(
        chematic_smiles::canonical_smiles(&mol),
        chematic_smiles::canonical_smiles(&want)
    );
}

#[test]
fn push_fixed4_matches_the_float_formatter() {
    let mut values = vec![
        0.0,
        -0.0,
        0.03125,
        -0.03125,
        0.00005,
        -0.00005,
        1.23455,
        2.5e-5,
        7.5e-5,
        0.99995,
        -0.99995,
        12345.6789,
        -99999.99995,
        1e9,
        1e16,
        -1e16,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
        1.5,
        -1.5,
    ];
    // A fixed linear congruential sequence over the coordinate range.
    let mut state = 0x2545_f491_4f6c_dd1du64;
    for _ in 0..200_000 {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let unit = (state >> 11) as f64 / (1u64 << 53) as f64;
        values.push((unit - 0.5) * 2_000.0);
        values.push(((unit - 0.5) * 2e8).round() / 1e4 + 5e-5);
    }
    for x in values {
        for width in [0, 10] {
            let mut fast = String::new();
            push_fixed4(&mut fast, x, width);
            assert_eq!(fast, format!("{x:>width$.4}"), "{x:e}");
        }
    }
}

#[test]
fn aromatic_bond_atoms_read_as_aromatic() {
    for smiles in ["c1ccccc1", "Cc1ccncc1", "c1cc[nH]c1"] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let v2000 = parse_mol(&write_mol(&mol, &MolMetadata::default()))
            .unwrap()
            .0;
        let coords: Vec<_> = (0..mol.atom_count())
            .map(|i| (1.5 * i as f64, 0.0))
            .collect();
        let v3000 = crate::parse_mol_v3000(&crate::write_mol_v3000(
            &mol,
            &MolMetadata::default(),
            &coords,
        ))
        .unwrap()
        .0;
        // Pyrrole's N-H is not recoverable from type-4 bonds, so that
        // molecule is written Kekulé; compare after aromaticity
        // perception.
        let perceived = |m: &Molecule| {
            chematic_smiles::canonical_smiles(&chematic_perception::apply_aromaticity(m))
        };
        for back in [v2000, v3000] {
            assert_eq!(perceived(&back), perceived(&mol), "{smiles}");
        }
    }
}

#[test]
fn three_column_parser_matches_legacy_spelling_contract() {
    // All ASCII triples cover whitespace, signs, malformed fields and
    // every 0..999 value, including fallback-accepted "+12" and "1  ".
    for a in 0..128_u8 {
        for b in 0..128_u8 {
            for c in 0..128_u8 {
                let bytes = [a, b, c];
                let field = std::str::from_utf8(&bytes).unwrap();
                let expected = field.trim().parse::<usize>().ok();
                let actual = parse_field3(field, 0, 1, |line, detail| {
                    MolParseError::InvalidCountLine { line, detail }
                })
                .ok();
                assert_eq!(actual, expected, "field={field:?}");
            }
        }
    }
    for field in ["é1", "1é", "\u{2000}"] {
        assert!(
            parse_field3(field, 0, 1, |line, detail| {
                MolParseError::InvalidCountLine { line, detail }
            })
            .is_err()
        );
    }
}

#[test]
fn fixed_width_integer_writers_match_rust_formatting() {
    for value in 0..=1001 {
        let mut actual = String::new();
        push_right_aligned_u32(&mut actual, value, 3);
        assert_eq!(actual, format!("{value:>3}"));
    }
    for (value, width) in [
        (0_u32, 1_usize),
        (0, 3),
        (1, 3),
        (12, 3),
        (123, 3),
        (1_234, 3),
        (u32::MAX, 3),
    ] {
        let mut actual = String::new();
        push_right_aligned_u32(&mut actual, value, width);
        assert_eq!(actual, format!("{value:>width$}"));
    }

    for (value, width) in [
        (i16::MIN, 2_usize),
        (-123, 2),
        (-12, 3),
        (-1, 3),
        (0, 3),
        (12, 3),
        (123, 3),
        (i16::MAX, 2),
    ] {
        let mut actual = String::new();
        push_right_aligned_i16(&mut actual, value, width);
        assert_eq!(actual, format!("{value:>width$}"));
    }
}

#[test]
fn shared_mol_v2000_roundtrip_contract_matches() {
    let document: serde_json::Value = serde_json::from_str(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../validation/cross_binding_contract.json"
    )))
    .expect("contract JSON");
    let contract = &document["mol_v2000_contract"];
    let (mol, _) = parse_mol(contract["input"].as_str().unwrap()).unwrap();
    let serialized = write_mol(&mol, &MolMetadata::default());
    let (roundtripped, _) = parse_mol(&serialized).unwrap();
    assert_eq!(
        roundtripped.atom_count(),
        contract["expected"]["atom_count"]
    );
    assert_eq!(roundtripped.bond_count(), 1);
}

/// Minimal ethanol MOL V2000 block (CCO, 3 atoms, 2 bonds).
const ETHANOL_MOL: &str = "\
ethanol
  chematic

  3  2  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
";

#[test]
fn test_parse_ethanol_counts() {
    let (mol, meta) = parse_mol(ETHANOL_MOL).expect("parse should succeed");
    assert_eq!(mol.atom_count(), 3);
    assert_eq!(mol.bond_count(), 2);
    assert_eq!(meta.name, "ethanol");
}

#[test]
fn test_parse_elements() {
    let (mol, _) = parse_mol(ETHANOL_MOL).expect("parse should succeed");
    let atoms: Vec<_> = mol.atoms().collect();
    assert_eq!(atoms[0].1.element, Element::C);
    assert_eq!(atoms[1].1.element, Element::C);
    assert_eq!(atoms[2].1.element, Element::O);
}

#[test]
fn isotope_mass_difference_round_trips_through_v2000_and_sdf() {
    for smiles in ["[13CH4]", "[15NH3]"] {
        let mol = chematic_smiles::parse(smiles).expect("isotope SMILES");
        let meta = MolMetadata::default().with_name("isotope");

        let mol_block = write_mol(&mol, &meta);
        let (roundtrip, _) = parse_mol(&mol_block).expect("V2000 round-trip");
        assert_eq!(roundtrip.atom_count(), 1);
        assert_eq!(
            roundtrip.atom(AtomIdx(0)).isotope,
            mol.atom(AtomIdx(0)).isotope
        );

        let sdf = write_sdf_record(&mol, &meta, &[], &std::collections::HashMap::new());
        let records = parse_sdf_with_coords(&sdf).expect("SDF round-trip");
        assert_eq!(records.len(), 1);
        assert_eq!(
            records[0].0.atom(AtomIdx(0)).isotope,
            mol.atom(AtomIdx(0)).isotope
        );
    }
}

#[test]
fn atom_map_round_trips_through_v2000() {
    let mol = chematic_smiles::parse("[CH3:7]").expect("mapped SMILES");
    let written = write_mol(&mol, &MolMetadata::default());
    let (roundtrip, _) = parse_mol(&written).expect("V2000 round-trip");
    assert_eq!(roundtrip.atom(AtomIdx(0)).atom_map, Some(7));
}

#[test]
fn test_parse_bond_types() {
    // Two carbons: single, double, triple, aromatic bonds.
    let mol_str = "\
test
  chematic

  8  4  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  3  4  2  0
  5  6  3  0
  7  8  4  0
M  END
";
    let (mol, _) = parse_mol(mol_str).expect("parse should succeed");
    let bonds: Vec<_> = mol.bonds().collect();
    assert_eq!(bonds[0].1.order, BondOrder::Single);
    assert_eq!(bonds[1].1.order, BondOrder::Double);
    assert_eq!(bonds[2].1.order, BondOrder::Triple);
    assert_eq!(bonds[3].1.order, BondOrder::Aromatic);
}

/// MDL bond type 9 (dative/coordinate) must read as `BondOrder::Dative`,
/// not silently fall through to `Single` -- regression test for the
/// platinum coordination-chemistry benchmark
/// (`validation/platinum/FEASIBILITY.md`). See `mol3000.rs`'s
/// `test_dative_bond_type_9_round_trips` for the RDKit-generated V3000
/// version of the same finding (RDKit itself only ever writes bond type
/// 9 in V3000, never V2000 -- this V2000 case defensively covers any
/// other tool that does).
#[test]
fn test_parse_bond_type_9_is_dative() {
    let mol_str = "\
test
  chematic

  3  2  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 Pt  0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  9  0
  2  3  1  0
M  END
";
    let (mol, _) = parse_mol(mol_str).expect("parse should succeed");
    let bonds: Vec<_> = mol.bonds().collect();
    assert_eq!(bonds[0].1.order, BondOrder::Dative);
    assert_eq!(mol.atom(bonds[0].1.atom1).element, Element::N);
    assert_eq!(mol.atom(bonds[0].1.atom2).element, Element::PT);
    assert_eq!(bonds[1].1.order, BondOrder::Single);
}

#[test]
fn test_parse_query_bond_types_preserved() {
    let mol_str = "\
query_bonds
  chematic

  8  4  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  5  0
  3  4  6  0
  5  6  7  0
  7  8  8  0
M  END
";
    let (mol, meta) = parse_mol(mol_str).expect("parse query bonds");
    let bonds: Vec<_> = mol.bonds().collect();
    assert_eq!(bonds[0].1.order, BondOrder::QuerySingleOrDouble);
    assert_eq!(bonds[1].1.order, BondOrder::QuerySingleOrAromatic);
    assert_eq!(bonds[2].1.order, BondOrder::QueryDoubleOrAromatic);
    assert_eq!(bonds[3].1.order, BondOrder::QueryAny);

    let written = write_mol(&mol, &meta);
    assert!(written.contains("  1  2  5  0"), "{written}");
    assert!(written.contains("  7  8  8  0"), "{written}");
}

#[test]
fn test_parse_charge() {
    // Nitrogen with charge code 3 (+1 formal charge).
    let mol_str = "\
charged
  chematic

  1  0  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
M  END
";
    let (mol, _) = parse_mol(mol_str).expect("parse should succeed");
    assert_eq!(mol.atom(AtomIdx(0)).charge, 1);
}

#[test]
fn test_parse_negative_charge() {
    // Oxygen with charge code 5 (-1 formal charge).
    let mol_str = "\
negcharge
  chematic

  1  0  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
M  END
";
    let (mol, _) = parse_mol(mol_str).expect("parse should succeed");
    assert_eq!(mol.atom(AtomIdx(0)).charge, -1);
}

#[test]
fn m_chg_property_overrides_fixed_width_charge_in_all_reader_paths() {
    // RDKit uses M CHG for valid charge assignments such as S+/O- even
    // when the V2000 atom charge column is zero. The property-level
    // formal charge must be visible to both reader modes.
    let mol_str = "\
charge_property
  RDKit          2D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  CHG  2   1   1   2  -1
M  END
";
    let (ordinary, _) = parse_mol(mol_str).expect("ordinary parse");
    let (fast, _) = parse_mol_fast(mol_str).expect("fast parse");
    for mol in [&ordinary, &fast] {
        assert_eq!(mol.atom(AtomIdx(0)).charge, 1);
        assert_eq!(mol.atom(AtomIdx(1)).charge, -1);
    }
}

#[test]
fn m_iso_property_overrides_fixed_width_mass_in_all_reader_paths() {
    let mol_str = "\
isotope_property
  RDKit          2D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  ISO  1   1  13
M  END
";
    let (ordinary, _) = parse_mol(mol_str).expect("ordinary parse");
    let (fast, _) = parse_mol_fast(mol_str).expect("fast parse");
    for mol in [&ordinary, &fast] {
        assert_eq!(mol.atom(AtomIdx(0)).isotope, Some(13));
        assert_eq!(mol.atom(AtomIdx(1)).isotope, None);
    }
}

#[test]
fn malformed_m_iso_is_a_typed_property_error() {
    let mol_str = "\
bad_isotope_property
  RDKit          2D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  ISO  1   2  13
M  END
";
    assert!(matches!(
        parse_mol(mol_str),
        Err(MolParseError::InvalidPropertyLine { .. })
    ));
}

#[test]
fn rdkit_2d_mol_roundtrip_preserves_alkene_stereo_and_isotope() {
    for smiles in ["O=C(O)/C=C/C(=O)O", "O=C(O)/C=C\\C(=O)O"] {
        let original = chematic_smiles::parse(smiles).expect("SMILES parse");
        let block = chematic_smiles::rdkit_mol_block_2d(&original).expect("MOL write");
        let report = read_mol_with_diagnostics(&block).expect("MOL read");
        assert!(
            report.ez_diagnostics.is_empty(),
            "{smiles}: {:?}",
            report.ez_diagnostics
        );
        let expected = chematic_smiles::rdkit_canonical_smiles(&original);
        let actual = chematic_smiles::rdkit_canonical_smiles(&report.mol);
        assert_eq!(actual, expected, "{smiles}\n{block}");
    }

    let original = chematic_smiles::parse("[13CH3][C@H](O)Cl").expect("SMILES parse");
    let block = chematic_smiles::rdkit_mol_block_2d(&original).expect("MOL write");
    let report = read_mol_with_diagnostics(&block).expect("MOL read");
    assert_eq!(report.mol.atom(AtomIdx(0)).isotope, Some(13));
    assert_eq!(
        chematic_smiles::rdkit_canonical_smiles(&report.mol),
        chematic_smiles::rdkit_canonical_smiles(&original)
    );
}

#[test]
fn test_charged_atom_writes_implicit_h_in_smiles() {
    // MOL V2000 has no per-atom H-count field, so a charged atom read from
    // this format always has `hydrogen_count: None` (the format leaves
    // implicit H to be inferred, unlike bracket SMILES). Regression test
    // for a bug where the SMILES writer silently dropped these inferred
    // hydrogens on bracket atoms forced by charge/isotope/atom-map
    // (found via MRV oracle validation, see validation/mrv_io_parity_summary.json).
    let mol_str = "\
charged
  chematic

  1  0  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
M  END
";
    let (mol, _) = parse_mol(mol_str).expect("parse should succeed");
    assert_eq!(mol.atom(AtomIdx(0)).hydrogen_count, None);
    assert_eq!(chematic_smiles::write(&mol), "[NH4+]");
    assert_eq!(chematic_smiles::canonical_smiles(&mol), "[NH4+]");
}

#[test]
fn test_round_trip() {
    // Parse → write → parse again; atom and bond counts must match.
    let (mol1, meta1) = parse_mol(ETHANOL_MOL).expect("first parse");
    let written = write_mol(&mol1, &meta1);
    let (mol2, _meta2) = parse_mol(&written).expect("second parse");
    assert_eq!(mol1.atom_count(), mol2.atom_count());
    assert_eq!(mol1.bond_count(), mol2.bond_count());
}

#[test]
fn test_round_trip_elements_preserved() {
    let (mol1, meta1) = parse_mol(ETHANOL_MOL).expect("first parse");
    let written = write_mol(&mol1, &meta1);
    let (mol2, _) = parse_mol(&written).expect("second parse");
    for ((_, a1), (_, a2)) in mol1.atoms().zip(mol2.atoms()) {
        assert_eq!(a1.element, a2.element);
    }
}

#[test]
fn test_error_missing_v2000() {
    let bad = "\
bad
  prog

  3  2  0  0  0  0  0  0  0  0  0 V3000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
";
    assert!(matches!(
        parse_mol(bad),
        Err(MolParseError::InvalidCountLine { .. })
    ));
}

#[test]
fn test_error_truncated_input() {
    // Counts line says 3 atoms but only 1 is provided.
    let bad = "\
trunc
  prog

  3  0  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
";
    assert!(matches!(parse_mol(bad), Err(MolParseError::UnexpectedEnd)));
}

#[test]
fn test_error_invalid_counts_line() {
    // Counts line too short (no V2000 tag at all).
    let bad = "\
mol
  prog

  X  Y
M  END
";
    assert!(matches!(
        parse_mol(bad),
        Err(MolParseError::InvalidCountLine { .. })
    ));
}

#[test]
fn test_write_contains_m_end() {
    let (mol, meta) = parse_mol(ETHANOL_MOL).expect("parse");
    let written = write_mol(&mol, &meta);
    assert!(written.contains("M  END"));
}

#[test]
fn test_write_contains_v2000() {
    let (mol, meta) = parse_mol(ETHANOL_MOL).expect("parse");
    let written = write_mol(&mol, &meta);
    assert!(written.contains("V2000"));
}

#[test]
fn test_parse_stereo_up_bond() {
    // MOL V2000 with a stereo=1 (Up) bond
    let mol_str = "\n\n\n  2  1  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  1  0  0  0\nM  END\n";
    let (mol, _) = crate::parse_mol(mol_str).unwrap();
    let bond = mol.bond(chematic_core::BondIdx(0));
    assert_eq!(bond.order, chematic_core::BondOrder::Up);
}

#[test]
fn test_parse_stereo_down_bond() {
    let mol_str = "\n\n\n  2  1  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  6  0  0  0\nM  END\n";
    let (mol, _) = crate::parse_mol(mol_str).unwrap();
    let bond = mol.bond(chematic_core::BondIdx(0));
    assert_eq!(bond.order, chematic_core::BondOrder::Down);
}

#[test]
fn test_parse_stereo_either_bond_code4_not_treated_as_wedge() {
    // MDL stereo code 4 ("either"/unspecified direction) is a third,
    // distinct state from a definite wedge (1) or hash (6) -- RDKit maps
    // it to `Bond::BondDir::UNKNOWN`, never a confident tag. Collapsing
    // it into `Up` (the pre-fix behavior) would fabricate a stereocenter
    // from a bond whose direction the file explicitly declares unknown,
    // now that stereo perception actually reads `BondOrder::Up`. This is
    // a documented, accepted lossy case: the "either" information itself
    // is not preserved (round-trips as a plain, unmarked single bond),
    // but it must never masquerade as a real wedge.
    let mol_str = "\n\n\n  2  1  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  4  0  0  0\nM  END\n";
    let (mol, _) = crate::parse_mol(mol_str).unwrap();
    let bond = mol.bond(chematic_core::BondIdx(0));
    assert_eq!(bond.order, chematic_core::BondOrder::Single);

    // Round-trip: writes as stereo field 0, not a false wedge/hash.
    let (mol, meta) = crate::parse_mol(mol_str).unwrap();
    let written = write_mol(&mol, &meta);
    assert!(written.contains("  1  2  1  0"), "{written}");
}

#[test]
fn test_molmetadata_builder() {
    let meta = MolMetadata::default()
        .with_name("aspirin")
        .with_comment("test molecule");
    assert_eq!(meta.name, "aspirin");
    assert_eq!(meta.comment, "test molecule");
}

#[test]
fn valence_field_keeps_hydrogen_counts_bonds_do_not_give() {
    // Bracket atoms whose hydrogens a reader cannot infer from the bonds
    // (hypervalent iodine, a sulfonyl radical, a carbene) carry the V2000
    // valence field; 7 exposed-10k rows read back with an added H.
    for smi in [
        "Cl[I]Cl",
        "c1ccc(cc1)[I]c1ccccc1",
        "CNC[S](=O)=O",
        "C[CH]C",
        "[CH2]",
    ] {
        let mol = chematic_smiles::parse(smi).unwrap();
        let (block, _) = write_laid_out_mol(&mol, &MolMetadata::default());
        let back = read_mol_with_diagnostics(&block).unwrap().mol;
        assert_eq!(
            chematic_smiles::canonical_smiles(&back),
            chematic_smiles::canonical_smiles(&mol),
            "{smi}\n{block}"
        );
    }
    // Atoms whose hydrogens the bonds give leave the field at 0.
    let mol = chematic_smiles::parse("OCC[NH3+]").unwrap();
    let block = write_mol(&mol, &MolMetadata::default());
    for line in block.lines().skip(4).take(mol.atom_count()) {
        assert_eq!(line.get(48..51), Some("  0"), "{block}");
    }
}

#[test]
fn test_molmetadata_with_name_roundtrip() {
    // Build a two-atom molecule and write it → name appears on MOL header line 1.
    use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};
    let mut b = MoleculeBuilder::new();
    let c1 = b.add_atom(Atom::new(Element::C));
    let c2 = b.add_atom(Atom::new(Element::C));
    b.add_bond(c1, c2, BondOrder::Single).unwrap();
    let mol = b.build();

    let meta = MolMetadata::default().with_name("acetic acid");
    let molblock = crate::write_mol(&mol, &meta);
    assert!(
        molblock.starts_with("acetic acid"),
        "MOL block must start with the molecule name"
    );
}

#[test]
fn test_declared_max_atom_count_truncated_input_errors() {
    let bad = "\
max_atoms
  chematic

999  0  0  0  0  0  0  0  0  0  0 V2000
";
    assert!(matches!(parse_mol(bad), Err(MolParseError::UnexpectedEnd)));
}

#[test]
fn test_declared_large_bond_count_truncated_input_errors() {
    let bad = "\
many_bonds
  chematic

  1 999  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
";
    assert!(matches!(parse_mol(bad), Err(MolParseError::UnexpectedEnd)));
}

#[test]
fn test_write_sdf_record_props_roundtrip() {
    use crate::sdf::SdfRecordReader;
    use std::collections::HashMap;

    let (mol, meta) = parse_mol(ETHANOL_MOL).expect("parse");
    let coords = vec![(0.0, 0.0), (1.5, 0.0), (3.0, 0.0)];
    let mut props = HashMap::new();
    props.insert("Activity".to_string(), "7.2".to_string());
    props.insert("Source".to_string(), "test".to_string());
    props.insert("_Name".to_string(), "ethanol".to_string()); // internal, should be omitted

    let sdf = write_sdf_record(&mol, &meta, &coords, &props);

    // _Name must NOT appear as an SD field
    assert!(
        !sdf.contains("> <_Name>"),
        "internal prop leaked to SD block"
    );
    // Other props must appear
    assert!(sdf.contains("> <Activity>\n7.2"), "Activity missing");
    assert!(sdf.contains("> <Source>\ntest"), "Source missing");
    assert!(sdf.ends_with("$$$$\n"), "missing delimiter");

    // Parse back and verify
    let rec = SdfRecordReader::new(&sdf)
        .next()
        .expect("should have record")
        .expect("should parse");
    assert_eq!(rec.mol.atom_count(), mol.atom_count());
    assert_eq!(
        rec.properties.get("Activity").map(|s| s.as_str()),
        Some("7.2")
    );
    assert_eq!(
        rec.properties.get("Source").map(|s| s.as_str()),
        Some("test")
    );
    assert!(!rec.properties.contains_key("_Name"));
}

#[test]
fn reusable_sdf_buffer_is_byte_identical_and_cleared() {
    use std::collections::HashMap;

    let (mol, meta) = parse_mol(ETHANOL_MOL).expect("parse");
    let coords = vec![(0.0, 0.0), (1.5, 0.0), (3.0, 0.0)];
    let props = HashMap::from([
        ("Activity".to_string(), "7.2".to_string()),
        ("Multiline".to_string(), "first\nsecond".to_string()),
    ]);
    let expected = write_sdf_record(&mol, &meta, &coords, &props);

    let mut reused = String::from("stale bytes must be removed");
    write_sdf_record_into(&mut reused, &mol, &meta, &coords, &props);
    assert_eq!(reused, expected);

    write_sdf_record_into(&mut reused, &mol, &meta, &[], &HashMap::new());
    assert_eq!(reused, write_sdf_record(&mol, &meta, &[], &HashMap::new()));
    assert!(!reused.contains("stale bytes"));
    assert!(!reused.contains("Activity"));
}

// ── Issue #171: blank MOL name line must not be eaten as inter-record
// padding ──────────────────────────────────────────────────────────────
// Fixture generated by a live RDKit 2026.03.3 oracle:
//   AllChem.Compute2DCoords(Chem.MolFromSmiles("CC")); Chem.MolToMolBlock(mol)
// RDKit's own MolToMolBlock leaves the name line blank for an unnamed
// molecule -- this is the literal, spec-legal shape that previously broke
// all 3 SDF-splitting entry points (see issue #171's repro).
const RDKIT_BLANK_NAME_MOL: &str = "\
\n     RDKit          2D\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n   -0.7500    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7500   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\nM  END\n";

fn crlf(s: &str) -> String {
    s.replace('\n', "\r\n")
}

#[test]
fn test_read_sdf_with_diagnostics_blank_name_single_record() {
    let sdf = format!("{RDKIT_BLANK_NAME_MOL}$$$$\n");
    let reports = read_sdf_with_diagnostics(&sdf).expect("parse should succeed");
    assert_eq!(reports.len(), 1);
    assert_eq!(reports[0].metadata.name, "");
    assert_eq!(reports[0].mol.atom_count(), 2);
}

#[test]
fn test_read_sdf_with_diagnostics_blank_name_first_record() {
    let sdf = format!("{RDKIT_BLANK_NAME_MOL}$$$$\n{ETHANOL_MOL}$$$$\n");
    let reports = read_sdf_with_diagnostics(&sdf).expect("parse should succeed");
    assert_eq!(reports.len(), 2);
    assert_eq!(reports[0].metadata.name, "");
    assert_eq!(reports[0].mol.atom_count(), 2);
    assert_eq!(reports[1].metadata.name, "ethanol");
    assert_eq!(reports[1].mol.atom_count(), 3);
}

#[test]
fn test_read_sdf_with_diagnostics_blank_name_middle_record() {
    let sdf = format!("{ETHANOL_MOL}$$$$\n{RDKIT_BLANK_NAME_MOL}$$$$\n{ETHANOL_MOL}$$$$\n");
    let reports = read_sdf_with_diagnostics(&sdf).expect("parse should succeed");
    assert_eq!(reports.len(), 3);
    assert_eq!(reports[0].metadata.name, "ethanol");
    assert_eq!(reports[1].metadata.name, "");
    assert_eq!(reports[1].mol.atom_count(), 2);
    assert_eq!(reports[2].metadata.name, "ethanol");
}

#[test]
fn test_read_sdf_with_diagnostics_blank_name_last_record() {
    let sdf = format!("{ETHANOL_MOL}$$$$\n{RDKIT_BLANK_NAME_MOL}$$$$\n");
    let reports = read_sdf_with_diagnostics(&sdf).expect("parse should succeed");
    assert_eq!(reports.len(), 2);
    assert_eq!(reports[0].metadata.name, "ethanol");
    assert_eq!(reports[1].metadata.name, "");
    assert_eq!(reports[1].mol.atom_count(), 2);
}

#[test]
fn test_read_sdf_with_diagnostics_blank_name_crlf() {
    let sdf = crlf(&format!(
        "{ETHANOL_MOL}$$$$\n{RDKIT_BLANK_NAME_MOL}$$$$\n{ETHANOL_MOL}$$$$\n"
    ));
    let reports = read_sdf_with_diagnostics(&sdf).expect("parse should succeed (CRLF)");
    assert_eq!(reports.len(), 3);
    assert_eq!(reports[1].metadata.name, "");
    assert_eq!(reports[1].mol.atom_count(), 2);
}

#[test]
fn ctfile_default_valence_retains_nonorganic_hydrogens() {
    let atom_line = "    0.0000    0.0000    0.0000 Na  0  0  0  0  0  0  0  0  0  0  0  0";
    for (symbol, valence_code, expected_h) in [("Na", 0, 1), ("K", 0, 1), ("Na", 15, 0)] {
        let mut line = atom_line.to_string();
        line.replace_range(31..34, &format!("{symbol:<3}"));
        line.replace_range(48..51, &format!("{valence_code:>3}"));
        let block =
            format!("\n OpenBabel 0D\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n{line}\nM  END\n");
        let read = read_mol_with_diagnostics(&block).expect("isolated atom parses");
        assert_eq!(read.mol.atom(AtomIdx(0)).hydrogen_count, Some(expected_h));
    }

    let carbon = "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0";
    for (symbol, bond_order, valence_code, expected_h) in [
        ("Se", 1, 0, 1),
        ("Sn", 1, 0, 1),
        ("As", 2, 0, 1),
        ("Si", 2, 0, 2),
        ("Se", 3, 0, 1),
        ("Sn", 3, 0, 1),
        ("Se", 1, 1, 0),
    ] {
        let mut line = atom_line.to_string();
        line.replace_range(31..34, &format!("{symbol:<3}"));
        line.replace_range(48..51, &format!("{valence_code:>3}"));
        let block = format!(
            "\n OpenBabel 0D\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n{carbon}\n{line}\n  1  2  {bond_order}  0  0  0  0\nM  END\n"
        );
        let read = read_mol_with_diagnostics(&block).expect("bonded atom parses");
        assert_eq!(
            read.mol.atom(AtomIdx(1)).hydrogen_count,
            Some(expected_h),
            "{symbol} with bond order {bond_order}"
        );
    }
}

#[test]
fn test_parse_sdf_with_coords_blank_name_middle_record() {
    // parse_sdf_with_coords is a thin wrapper over read_sdf_with_diagnostics;
    // confirm the fix reaches it too.
    let sdf = format!("{ETHANOL_MOL}$$$$\n{RDKIT_BLANK_NAME_MOL}$$$$\n{ETHANOL_MOL}$$$$\n");
    let mols = parse_sdf_with_coords(&sdf).expect("parse should succeed");
    assert_eq!(mols.len(), 3);
    assert_eq!(mols[1].1.name, "");
    assert_eq!(mols[1].0.atom_count(), 2);
}

#[test]
fn test_read_sdf_with_diagnostics_malformed_input_recovery_unaffected() {
    // Existing malformed-input behavior (bad counts line) must still error,
    // confirming the blank-line-skip removal didn't loosen error recovery.
    let bad_sdf = format!("{ETHANOL_MOL}$$$$\nbad\n  prog\n\n  X  Y\nM  END\n$$$$\n");
    assert!(read_sdf_with_diagnostics(&bad_sdf).is_err());
}
