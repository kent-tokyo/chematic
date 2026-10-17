use super::*;

#[test]
fn a_ring_ch_cannot_be_its_own_paired_center() {
    // With explicit H, the same CH belongs to both candidate lists. There
    // is still only one substituted ring site, so it has no cis/trans pair.
    for smiles in ["[H]C1(CCCCC1)C", "CC1([H])CCCCC1", "[H][C]1(C)CCCCC1"] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        for mol2_extensions in [false, true] {
            assert!(
                additional_3d_ring_centers(&mol, mol2_extensions).is_empty(),
                "{smiles}"
            );
        }
    }
}

#[test]
fn a_distinct_ring_ch_partner_is_retained() {
    // The two exocyclic elements F/Cl make one quaternary candidate. The
    // separate methyl-substituted CH supplies the required partner.
    let mol = chematic_smiles::parse("CC1CCC(F)(Cl)CC1").unwrap();
    for mol2_extensions in [false, true] {
        let centers = additional_3d_ring_centers(&mol, mol2_extensions);
        assert_eq!(centers.len(), 2);
        assert_ne!(centers[0], centers[1]);
    }
}
