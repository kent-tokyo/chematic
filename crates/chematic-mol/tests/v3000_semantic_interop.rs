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
        "{HEADER}M  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 2 C 1 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 END BOND\n{blocks}M  V30 END CTAB\nM  END\n"
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
