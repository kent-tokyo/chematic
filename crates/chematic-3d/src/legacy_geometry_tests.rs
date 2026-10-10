//! Keep the retained differential baseline useful: simple geometry must remain
//! deterministic, and the connectivity rewrite must solve its known fixtures.
use crate::dg::{generate_coords, generate_coords_legacy};
use chematic_smiles::parse;

#[test]
fn simple_baseline_geometry_has_correct_connectivity_and_component_separation() {
    for smiles in [
        "C",
        "CCO",
        "CC#N",
        "CC(C)C",
        "C1CCCCC1",
        "c1ccccc1C",
        "Cc1ccc(C)cc1",
        "c1ccccc1-c1ccccc1",
        "CC.O",
    ] {
        let mol = parse(smiles).unwrap();
        let first = generate_coords_legacy(&mol);
        assert_eq!(first, generate_coords_legacy(&mol), "{smiles}");
        for (_, bond) in mol.bonds() {
            let distance = first.get(bond.atom1).distance(&first.get(bond.atom2));
            assert!((0.8..1.95).contains(&distance), "{smiles}: {distance}");
        }
        if smiles == "CC.O" {
            assert!(
                first.get(chematic_core::AtomIdx(2)).x - first.get(chematic_core::AtomIdx(1)).x
                    >= 5.
            );
        }
    }
    let empty = chematic_core::MoleculeBuilder::new().build();
    assert_eq!(generate_coords_legacy(&empty).atom_count(), 0);
    assert_eq!(
        generate_coords_legacy(&parse("C").unwrap()).points,
        vec![crate::Point3::zero()]
    );
}

#[test]
fn connectivity_rewrite_resolves_legacy_fusion_and_chain_bridge_outliers() {
    for smiles in [
        "c1ccc2ccccc2c1",
        "C1CCC2(CC1)CCCC2",
        "c1ccccc1CCc1ccccc1",
        "c1ccc2c(c1)ccc1ccccc12",
    ] {
        let mol = parse(smiles).unwrap();
        let legacy = generate_coords_legacy(&mol);
        let current = generate_coords(&mol);
        let lengths = |coords: &crate::Coords3D| {
            mol.bonds()
                .map(|(_, bond)| coords.get(bond.atom1).distance(&coords.get(bond.atom2)))
                .collect::<Vec<_>>()
        };
        assert!(
            lengths(&legacy).iter().any(|&r| !(0.8..1.95).contains(&r)),
            "baseline no longer reproduces {smiles}"
        );
        assert!(
            lengths(&current).iter().all(|&r| (0.8..1.95).contains(&r)),
            "{smiles}: {:?}",
            lengths(&current)
        );
    }
}
