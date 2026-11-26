use chematic_core::BondOrder;
use chematic_mol::{MrvWriteOptions, parse_mrv, write_mrv};

#[test]
fn disabling_stereo_omits_wedges_without_changing_connectivity() {
    let input = r#"<cml><molecule><atomArray><atom id="a1" elementType="C"/><atom id="a2" elementType="N"/><atom id="a3" elementType="O"/></atomArray><bondArray><bond atomRefs2="a1 a2" order="1"><bondStereo>W</bondStereo></bond><bond atomRefs2="a1 a3" order="1"><bondStereo>H</bondStereo></bond></bondArray></molecule></cml>"#;
    let original = parse_mrv(input).unwrap();
    assert_eq!(
        original
            .mol
            .bonds()
            .map(|(_, b)| b.order)
            .collect::<Vec<_>>(),
        [BondOrder::Up, BondOrder::Down]
    );
    let options = MrvWriteOptions {
        include_stereo: false,
        ..Default::default()
    };
    let written = write_mrv(&original, &options).unwrap();
    assert!(!written.contains("bondStereo"));
    let reread = parse_mrv(&written).unwrap();
    assert_eq!(reread.mol.atom_count(), 3);
    assert_eq!(
        reread.mol.bonds().map(|(_, b)| b.order).collect::<Vec<_>>(),
        [BondOrder::Single, BondOrder::Single]
    );
}

#[test]
fn cis_trans_text_keeps_double_bonds_without_inventing_wedges() {
    for stereo in ["C", "T"] {
        let input = format!(
            r#"<cml><molecule><atomArray><atom id="a1" elementType="C"/><atom id="a2" elementType="C"/></atomArray><bondArray><bond atomRefs2="a1 a2" order="2"><bondStereo>{stereo}</bondStereo></bond></bondArray></molecule></cml>"#
        );
        let parsed = parse_mrv(&input).unwrap();
        assert_eq!(parsed.mol.bond_count(), 1);
        assert_eq!(
            parsed.mol.bond(chematic_core::BondIdx(0)).order,
            BondOrder::Double
        );
    }
}
