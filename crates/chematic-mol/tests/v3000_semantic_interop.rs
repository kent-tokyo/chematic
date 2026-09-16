//! A4 V3000 semantic-preservation fixtures.
//!
//! The ordering of SGROUP and COLLECTION blocks is not chemical semantics,
//! but consumers such as RDKit are sensitive to valid block placement. Keep
//! both accepted input orders in the local contract and assert that chematic
//! emits the canonical SGROUP-before-COLLECTION order without dropping the
//! typed SGROUP view.

use chematic_mol::{V3000SGroupKind, parse_mol_v3000, parse_v3000_sgroup_line, write_mol_v3000};

const HEADER: &str = "\n\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n";

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
