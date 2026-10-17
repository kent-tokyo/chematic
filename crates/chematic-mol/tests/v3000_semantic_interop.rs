//! A4 V3000 semantic-preservation fixtures.
//!
//! The ordering of SGROUP and COLLECTION blocks is not chemical semantics,
//! but consumers such as RDKit are sensitive to valid block placement. Keep
//! both accepted input orders in the local contract and assert that chematic
//! emits the canonical SGROUP-before-COLLECTION order without dropping the
//! typed SGROUP view.

use chematic_mol::{V3000SGroupKind, parse_mol_v3000, parse_v3000_sgroup_line, write_mol_v3000};

const HEADER: &str = "\n\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n";

#[test]
fn v3000_concrete_atoms_omit_redundant_query_hcount() {
    for (smiles, count) in [("F[C@H](Cl)Br", 1), ("F[C@](Cl)(Br)I", 0)] {
        let mol = chematic_smiles::parse(smiles).expect("SMILES parse");
        let block = write_mol_v3000(&mol, &Default::default(), &[]);
        assert!(
            !block.contains(&format!("HCOUNT={count}")),
            "valence-inferred carbon H count is not a query property: {smiles}"
        );
        parse_mol_v3000(&block).expect("concrete V3000 parse");
    }

    let radical = chematic_smiles::parse("[CH]C").expect("radical SMILES parse");
    let block = write_mol_v3000(&radical, &Default::default(), &[]);
    assert!(block.contains("HCOUNT=1"), "nonredundant H must remain");

    let aromatic = chematic_smiles::parse("c1c[nH]cn1").expect("aromatic SMILES parse");
    let block = write_mol_v3000(&aromatic, &Default::default(), &[]);
    assert!(
        !block.contains("HCOUNT=1"),
        "aromatic N H follows from the Kekule bonds"
    );
    let (reparsed, _) = parse_mol_v3000(&block).expect("aromatic V3000 parse");
    let n_h = aromatic
        .atoms()
        .find(|(_, atom)| {
            atom.element == chematic_core::Element::N && atom.hydrogen_count == Some(1)
        })
        .expect("aromatic [nH]")
        .0;
    let atom_prefix = format!("M  V30 {} N ", n_h.0 + 1);
    let with_query_h = block
        .lines()
        .map(|line| {
            if line.starts_with(&atom_prefix) {
                format!("{line} HCOUNT=1\n")
            } else {
                format!("{line}\n")
            }
        })
        .collect::<String>();
    let (baseline, _) = parse_mol_v3000(&with_query_h).expect("query HCOUNT V3000 parse");
    assert_eq!(
        chematic_smiles::canonical_smiles(&reparsed),
        chematic_smiles::canonical_smiles(&baseline)
    );
}

fn block(collection_first: bool) -> String {
    let collection =
        "M  V30 BEGIN COLLECTION\nM  V30 MDLV30/STEABS ATOMS=(1 1)\nM  V30 END COLLECTION\n";
    let sgroup =
        "M  V30 BEGIN SGROUP\nM  V30 1 SUP 0 ATOMS=(2 1 2) LABEL=core\nM  V30 END SGROUP\n";
    let blocks = if collection_first {
        format!("{collection}{sgroup}")
    } else {
        format!("{sgroup}{collection}")
    };
    format!(
        "{HEADER}M  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 1 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 2 C 1 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 END BOND\n{blocks}M  V30 END CTAB\nM  END\n"
    )
}

#[test]
fn both_v3000_block_orders_have_the_same_semantic_projection() {
    let mut projections = Vec::new();
    for collection_first in [false, true] {
        let (mol, metadata) = parse_mol_v3000(&block(collection_first)).expect("V3000 parse");
        let typed = parse_v3000_sgroup_line(&metadata.v3000_sgroups[0]).expect("typed SGROUP");
        assert_eq!(typed.kind, V3000SGroupKind::Sup);
        assert_eq!(typed.atom_ids, vec![1, 2]);
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.bond_count(), 1);

        let rewritten = write_mol_v3000(&mol, &metadata, &[]);
        let sgroup_pos = rewritten.find("BEGIN SGROUP").expect("SGROUP output");
        let collection_pos = rewritten
            .find("BEGIN COLLECTION")
            .expect("COLLECTION output");
        assert!(sgroup_pos < collection_pos, "canonical block ordering");
        assert!(rewritten.contains("M  V30 COUNTS 2 1 1 0 0"));

        let (roundtrip, roundtrip_metadata) =
            parse_mol_v3000(&rewritten).expect("rewritten V3000 parse");
        projections.push((
            roundtrip.atom_count(),
            roundtrip.bond_count(),
            roundtrip_metadata.v3000_sgroups,
            roundtrip.stereo_groups().len(),
        ));
    }
    assert_eq!(projections[0], projections[1]);
    assert_eq!(projections[0].2, vec!["1 SUP 0 ATOMS=(2 1 2) LABEL=core"]);
}

#[test]
fn malformed_known_collection_atom_count_is_rejected() {
    let input = format!(
        "{HEADER}M  V30 BEGIN CTAB\n\
M  V30 COUNTS 2 1 0 0 0\n\
M  V30 BEGIN ATOM\n\
M  V30 1 C 0 0 0 0\n\
M  V30 2 C 1 0 0 0\n\
M  V30 END ATOM\n\
M  V30 BEGIN BOND\n\
M  V30 1 1 1 2\n\
M  V30 END BOND\n\
M  V30 BEGIN COLLECTION\n\
M  V30 MDLV30/STEABS ATOMS=(2 1)\n\
M  V30 END COLLECTION\n\
M  V30 END CTAB\n\
M  END\n"
    );
    let error = match parse_mol_v3000(&input) {
        Err(error) => error,
        Ok(_) => panic!("declared count mismatch must fail"),
    };
    assert!(error.to_string().contains("declares 2 ids but contains 1"));
}

#[test]
fn known_collection_unknown_atom_is_rejected() {
    let input = format!(
        "{HEADER}M  V30 BEGIN CTAB\n\
M  V30 COUNTS 2 1 0 0 0\n\
M  V30 BEGIN ATOM\n\
M  V30 1 C 0 0 0 0\n\
M  V30 2 C 1 0 0 0\n\
M  V30 END ATOM\n\
M  V30 BEGIN BOND\n\
M  V30 1 1 1 2\n\
M  V30 END BOND\n\
M  V30 BEGIN COLLECTION\n\
M  V30 MDLV30/STEABS ATOMS=(1 99)\n\
M  V30 END COLLECTION\n\
M  V30 END CTAB\n\
M  END\n"
    );
    let error = match parse_mol_v3000(&input) {
        Err(error) => error,
        Ok(_) => panic!("unknown atom reference must fail"),
    };
    assert!(error.to_string().contains("references unknown atom id 99"));
}

#[test]
fn v3000_bond_stereo_collections_round_trip_and_merge_with_atom_members() {
    let input = format!(
        "{HEADER}M  V30 BEGIN CTAB\n\
M  V30 COUNTS 6 5 0 0 0\n\
M  V30 BEGIN ATOM\n\
M  V30 1 C 0 0 0 0\n\
M  V30 2 C 1 0 0 0\n\
M  V30 3 C 2 0 0 0\n\
M  V30 4 C 3 0 0 0\n\
M  V30 5 C 4 0 0 0\n\
M  V30 6 C 5 0 0 0\n\
M  V30 END ATOM\n\
M  V30 BEGIN BOND\n\
M  V30 4 1 1 2\n\
M  V30 7 2 2 3\n\
M  V30 9 1 3 4\n\
M  V30 11 2 4 5\n\
M  V30 15 1 5 6\n\
M  V30 END BOND\n\
M  V30 BEGIN COLLECTION\n\
M  V30 MDLV30/STERAC1 ATOMS=(1 3)\n\
M  V30 MDLV30/STEBRAC1 BONDS=(2 7 11)\n\
M  V30 END COLLECTION\n\
M  V30 END CTAB\n\
M  END\n"
    );

    let (mol, metadata) = parse_mol_v3000(&input).expect("STEB collection parse");
    assert_eq!(mol.stereo_groups().len(), 1);
    let group = &mol.stereo_groups()[0];
    assert_eq!(group.kind, chematic_core::StereoGroupKind::And(1));
    assert_eq!(group.atom_indices, vec![chematic_core::AtomIdx(2)]);
    assert_eq!(
        mol.stereo_group_bond_indices(0).unwrap(),
        vec![chematic_core::BondIdx(1), chematic_core::BondIdx(3)]
    );

    let rewritten = write_mol_v3000(&mol, &metadata, &[]);
    assert!(rewritten.contains("MDLV30/STERAC1 ATOMS=(1 3)"));
    assert!(rewritten.contains("MDLV30/STEBRAC1 BONDS=(2 2 4)"));
    let (roundtrip, _) = parse_mol_v3000(&rewritten).expect("STEB collection round trip");
    assert_eq!(roundtrip.stereo_groups(), mol.stereo_groups());
    assert_eq!(
        roundtrip.stereo_group_bond_indices(0),
        mol.stereo_group_bond_indices(0)
    );
}

#[test]
fn known_collection_unknown_bond_is_rejected() {
    let input = block(false).replace("MDLV30/STEABS ATOMS=(1 1)", "MDLV30/STEBABS BONDS=(1 99)");
    let error = match parse_mol_v3000(&input) {
        Err(error) => error,
        Ok(_) => panic!("unknown bond reference must fail"),
    };
    assert!(error.to_string().contains("references unknown bond id 99"));
}

#[test]
fn all_ctfile_2026_bond_collection_kinds_round_trip() {
    let input = format!(
        "{HEADER}M  V30 BEGIN CTAB\n\
M  V30 COUNTS 4 3 0 0 0\n\
M  V30 BEGIN ATOM\n\
M  V30 1 C 0 0 0 0\n\
M  V30 2 C 1 0 0 0\n\
M  V30 3 C 2 0 0 0\n\
M  V30 4 C 3 0 0 0\n\
M  V30 END ATOM\n\
M  V30 BEGIN BOND\n\
M  V30 1 2 1 2\n\
M  V30 2 1 2 3\n\
M  V30 3 2 3 4\n\
M  V30 END BOND\n\
M  V30 BEGIN COLLECTION\n\
M  V30 MDLV30/STEBABS BONDS=(1 1)\n\
M  V30 MDLV30/STEBREL2 BONDS=(1 2)\n\
M  V30 MDLV30/STEBRAC3 BONDS=(2 3 3)\n\
M  V30 END COLLECTION\n\
M  V30 END CTAB\n\
M  END\n"
    );
    let (mol, metadata) = parse_mol_v3000(&input).expect("all STEB kinds parse");
    assert_eq!(mol.stereo_groups().len(), 3);
    assert_eq!(mol.stereo_group_bond_indices(0).unwrap().len(), 1);
    assert_eq!(mol.stereo_group_bond_indices(1).unwrap().len(), 1);
    assert_eq!(
        mol.stereo_group_bond_indices(2).unwrap().len(),
        1,
        "deduplicated"
    );

    let rewritten = write_mol_v3000(&mol, &metadata, &[]);
    assert!(rewritten.contains("MDLV30/STEBABS BONDS=(1 1)"));
    assert!(rewritten.contains("MDLV30/STEBREL2 BONDS=(1 2)"));
    assert!(rewritten.contains("MDLV30/STEBRAC3 BONDS=(1 3)"));
    let (roundtrip, _) = parse_mol_v3000(&rewritten).expect("all STEB kinds round trip");
    assert_eq!(roundtrip.stereo_groups(), mol.stereo_groups());
    for index in 0..mol.stereo_groups().len() {
        assert_eq!(
            roundtrip.stereo_group_bond_indices(index),
            mol.stereo_group_bond_indices(index)
        );
    }
}

#[test]
fn sgroup_count_must_match_the_counts_record() {
    let input = block(false).replacen("COUNTS 2 1 1 0 0", "COUNTS 2 1 0 0 0", 1);
    let error = match parse_mol_v3000(&input) {
        Err(error) => error,
        Ok(_) => panic!("SGROUP count mismatch must fail"),
    };
    assert!(
        error
            .to_string()
            .contains("COUNTS declares 0 SGROUP records but found 1")
    );
}

#[test]
fn known_sgroup_unknown_atom_is_rejected() {
    let input = block(false).replace("ATOMS=(2 1 2)", "ATOMS=(2 1 99)");
    let error = match parse_mol_v3000(&input) {
        Err(error) => error,
        Ok(_) => panic!("known SGROUP references must resolve to CTAB atoms"),
    };
    assert!(
        error
            .to_string()
            .contains("SGROUP ATOMS references unknown atom id 99")
    );
}

#[test]
fn duplicate_known_sgroup_ids_are_rejected() {
    let second = "M  V30 1 SUP 0 ATOMS=(1 1) LABEL=duplicate\n";
    let input = block(false).replace(
        "M  V30 END SGROUP\n",
        &format!("{second}M  V30 END SGROUP\n"),
    );
    let input = input.replacen("COUNTS 2 1 1 0 0", "COUNTS 2 1 2 0 0", 1);
    let error = match parse_mol_v3000(&input) {
        Err(error) => error,
        Ok(_) => panic!("duplicate typed SGROUP ids must fail"),
    };
    assert!(error.to_string().contains("duplicate typed SGROUP id 1"));
}

#[test]
fn unknown_sgroup_kind_remains_opaque() {
    let input = block(false).replace("1 SUP 0 ATOMS=(2 1 2) LABEL=core", "1 VENDORX 0 EXTREF=99");
    let (_, metadata) = parse_mol_v3000(&input).expect("unknown SGROUP kind remains opaque");
    assert_eq!(metadata.v3000_sgroups, vec!["1 VENDORX 0 EXTREF=99"]);
}

#[test]
fn bond_count_must_match_the_counts_record() {
    let input = block(false).replacen("COUNTS 2 1 1 0 0", "COUNTS 2 2 1 0 0", 1);
    let error = match parse_mol_v3000(&input) {
        Err(error) => error,
        Ok(_) => panic!("bond count mismatch must fail"),
    };
    assert!(
        error
            .to_string()
            .contains("COUNTS declares 2 bonds but found 1")
    );
}

#[test]
fn declared_bonds_require_a_bond_block() {
    let input = format!(
        "{HEADER}M  V30 BEGIN CTAB\n\
M  V30 COUNTS 2 1 0 0 0\n\
M  V30 BEGIN ATOM\n\
M  V30 1 C 0 0 0 0\n\
M  V30 2 C 1 0 0 0\n\
M  V30 END ATOM\n\
M  V30 END CTAB\n\
M  END\n"
    );
    let error = match parse_mol_v3000(&input) {
        Err(error) => error,
        Ok(_) => panic!("declared bonds require a BOND block"),
    };
    assert!(
        error
            .to_string()
            .contains("COUNTS declares 1 bonds but the BOND block is absent")
    );
}

#[test]
fn opaque_bond_attributes_survive_v3000_round_trip() {
    let input = format!(
        "{HEADER}M  V30 BEGIN CTAB\n\
M  V30 COUNTS 2 1 0 0 0\n\
M  V30 BEGIN ATOM\n\
M  V30 1 C 0 0 0 0\n\
M  V30 2 C 1 0 0 0\n\
M  V30 END ATOM\n\
M  V30 BEGIN BOND\n\
M  V30 1 1 1 2 ENDPTS=(2 1 2) ATTACH=ANY\n\
M  V30 END BOND\n\
M  V30 END CTAB\n\
M  END\n"
    );
    let (mol, metadata) = parse_mol_v3000(&input).expect("V3000 parse");
    assert_eq!(
        metadata.v3000_bond_properties,
        vec![(1, "ENDPTS=(2 1 2) ATTACH=ANY".to_string())]
    );

    let rewritten = write_mol_v3000(&mol, &metadata, &[]);
    assert!(
        rewritten.contains("1 1 1 2 ENDPTS=(2 1 2) ATTACH=ANY"),
        "opaque bond attributes must be emitted: {rewritten}"
    );
    let (_, roundtrip_metadata) = parse_mol_v3000(&rewritten).expect("round-trip parse");
    assert_eq!(
        roundtrip_metadata.v3000_bond_properties,
        metadata.v3000_bond_properties
    );
}

#[test]
fn indigo_style_deuterium_and_atom_cfg_survive_v3000_round_trip() {
    // Indigo 1.46.0 writes deuterium as `D` and tetrahedral parity on an atom
    // line. The core representation normalizes D to H/MASS=2 but must retain
    // atom CFG so a writer that requires it can still consume the output.
    let input = format!(
        "{HEADER}M  V30 BEGIN CTAB\n\
M  V30 COUNTS 2 1 0 0 0\n\
M  V30 BEGIN ATOM\n\
M  V30 1 D 0 0 0 0\n\
M  V30 2 C 1 0 0 0 CFG=1\n\
M  V30 END ATOM\n\
M  V30 BEGIN BOND\n\
M  V30 1 1 1 2\n\
M  V30 END BOND\n\
M  V30 END CTAB\n\
M  END\n"
    );
    let (mol, metadata) = parse_mol_v3000(&input).expect("V3000 parse");
    assert_eq!(mol.atom_count(), 2);
    assert_eq!(
        metadata.v3000_atom_properties,
        vec![(2, "CFG=1".to_string())]
    );

    let rewritten = write_mol_v3000(&mol, &metadata, &[]);
    assert!(rewritten.contains("1 H 0.0000 0.0000 0.0000 0 MASS=2"));
    assert!(rewritten.contains("2 C 0.0000 0.0000 0.0000 0 CFG=1"));
    let (_, roundtrip_metadata) = parse_mol_v3000(&rewritten).expect("round-trip parse");
    assert_eq!(
        roundtrip_metadata.v3000_atom_properties,
        metadata.v3000_atom_properties
    );
}

#[test]
fn opaque_v3000_attributes_follow_writer_order_for_noncontiguous_source_ids() {
    let input = format!(
        "{HEADER}M  V30 BEGIN CTAB\n\
M  V30 COUNTS 2 1 0 0 0\n\
M  V30 BEGIN ATOM\n\
M  V30 10 C 0 0 0 0\n\
M  V30 20 C 1 0 0 0 CFG=1\n\
M  V30 END ATOM\n\
M  V30 BEGIN BOND\n\
M  V30 30 1 10 20 ENDPTS=(2 10 20) ATTACH=ANY\n\
M  V30 END BOND\n\
M  V30 END CTAB\n\
M  END\n"
    );
    let (mol, metadata) = parse_mol_v3000(&input).expect("V3000 parse");
    assert_eq!(
        metadata.v3000_atom_properties,
        vec![(2, "CFG=1".to_string())]
    );
    assert_eq!(
        metadata.v3000_bond_properties,
        vec![(1, "ENDPTS=(2 1 2) ATTACH=ANY".to_string())]
    );

    let rewritten = write_mol_v3000(&mol, &metadata, &[]);
    assert!(rewritten.contains("2 C 0.0000 0.0000 0.0000 0 CFG=1"));
    assert!(rewritten.contains("1 1 1 2 ENDPTS=(2 1 2) ATTACH=ANY"));
}
