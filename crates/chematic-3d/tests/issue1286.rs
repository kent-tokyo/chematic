use chematic_3d::{
    add_hydrogens_with_coords, generate_and_minimize_constrained, generate_and_minimize_dreiding,
    generate_coords, generate_coords_etkdg, get_bond_angle_deg,
};
use chematic_core::AtomIdx;

#[test]
fn aromatic_substituents_and_rings() {
    for smiles in [
        "Cc1ccccc1",
        "c1ccccc1C",
        "Oc1ccccc1",
        "c1ccccc1O",
        "c1ccccc1",
    ] {
        let heavy = chematic_smiles::parse(smiles).unwrap();
        let explicit = chematic_chem::add_hydrogens(&heavy);
        let explicit_coords = generate_coords(&explicit);
        let (added, added_coords) =
            add_hydrogens_with_coords(&heavy, &generate_coords(&heavy)).unwrap();
        for (mol, coords) in [(&explicit, &explicit_coords), (&added, &added_coords)] {
            let mut checked = 0;
            for i in 0..mol.atom_count() {
                let center = AtomIdx(i as u32);
                if !mol.atom(center).aromatic {
                    continue;
                }
                let neighbors: Vec<_> = mol.neighbors(center).map(|(n, _)| n).collect();
                let ring: Vec<_> = neighbors
                    .iter()
                    .copied()
                    .filter(|n| mol.atom(*n).aromatic)
                    .collect();
                let outside: Vec<_> = neighbors
                    .iter()
                    .copied()
                    .filter(|n| !mol.atom(*n).aromatic)
                    .collect();
                if ring.len() == 2 && outside.len() == 1 {
                    checked += 1;
                    let p = coords.get(center);
                    let u = coords.get(ring[0]).sub(&p);
                    let v = coords.get(ring[1]).sub(&p);
                    let normal = u.cross(&v).normalize();
                    let height = coords.get(outside[0]).sub(&p).dot(&normal).abs();
                    assert!(height < 0.05, "{smiles}: height {height}");
                    for n in ring {
                        let angle = get_bond_angle_deg(coords, n, center, outside[0]);
                        assert!((angle - 120.0).abs() < 5.0, "{smiles}: {angle}");
                    }
                }
            }
            assert_eq!(checked, 6, "{smiles}: aromatic centers checked");
        }
    }
}

#[test]
fn etkdg_preserves_aromatic_ring_bonds() {
    for smiles in ["Cc1ccccc1", "Oc1ccccc1", "c1ccccc1"] {
        let mol = chematic_chem::add_hydrogens(&chematic_smiles::parse(smiles).unwrap());
        let coords = generate_coords_etkdg(&mol);
        for i in 0..mol.atom_count() {
            let a = AtomIdx(i as u32);
            if !mol.atom(a).aromatic {
                continue;
            }
            for (b, _) in mol.neighbors(a) {
                if mol.atom(b).aromatic {
                    let distance = coords.get(a).distance(&coords.get(b));
                    assert!((1.30..=1.50).contains(&distance), "{smiles}: {distance}");
                }
            }
        }
    }
}

#[test]
fn minimized_water_and_methane_bonds() {
    for (smiles, expected) in [("O", 0.96), ("C", 1.09)] {
        let mol = chematic_chem::add_hydrogens(&chematic_smiles::parse(smiles).unwrap());
        for coords in [
            generate_and_minimize_dreiding(&mol),
            generate_and_minimize_constrained(&mol),
        ] {
            for (h, _) in mol.neighbors(AtomIdx(0)) {
                let distance = coords.get(AtomIdx(0)).distance(&coords.get(h));
                assert!((distance - expected).abs() < 0.05, "{smiles}: {distance}");
            }
        }
    }
}

#[test]
fn biaryl_bond_is_in_plane_at_both_rings() {
    let mol = chematic_smiles::parse("c1ccccc1-c2ccccc2").unwrap();
    let coords = generate_coords(&mol);
    for (a, center, b) in [(4, 5, 6), (0, 5, 6), (5, 6, 7), (5, 6, 11)] {
        let angle = get_bond_angle_deg(&coords, AtomIdx(a), AtomIdx(center), AtomIdx(b));
        assert!((angle - 120.0).abs() < 5.0, "{a}-{center}-{b}: {angle}");
    }
}
