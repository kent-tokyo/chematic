//! Neutral three-coordinate boron is trigonal planar (120°), not tetrahedral:
//! RDKit perceives it as SP2 and types it `B_2` for UFF, even with single bonds only.

use chematic_3d::{
    Coords3D, add_hydrogens_with_coords, generate_coords_etkdg, get_bond_angle_deg, minimize_uff,
};
use chematic_core::{AtomIdx, Molecule};

/// Every X–B–X angle around the first boron atom.
fn boron_angles(mol: &Molecule, coords: &Coords3D) -> Vec<f64> {
    let (boron, _) = mol
        .atoms()
        .find(|(_, atom)| atom.element.atomic_number() == 5)
        .unwrap();
    let neighbors: Vec<AtomIdx> = mol.neighbors(boron).map(|(n, _)| n).collect();
    let mut angles = Vec::new();
    for i in 0..neighbors.len() {
        for j in i + 1..neighbors.len() {
            angles.push(get_bond_angle_deg(
                coords,
                neighbors[i],
                boron,
                neighbors[j],
            ));
        }
    }
    angles
}

fn assert_angles(smiles: &str, stage: &str, angles: &[f64], expected: f64, tolerance: f64) {
    for angle in angles {
        assert!(
            (angle - expected).abs() < tolerance,
            "{smiles}, {stage}: X-B-X {angle:.1}°, expected {expected}° ({angles:.1?})"
        );
    }
}

/// Planar (the angles add up to 360°, against 328.5° at a tetrahedral centre)
/// and each angle near 120°: the gradient-descent minimizer stops a few
/// degrees short when it starts from the placement's 97° angle.
fn assert_trigonal_planar(smiles: &str, stage: &str, angles: &[f64]) {
    let sum: f64 = angles.iter().sum();
    assert!(
        (sum - 360.0).abs() < 3.0,
        "{smiles}, {stage}: X-B-X angles add up to {sum:.1}° ({angles:.1?})"
    );
    assert_angles(smiles, stage, angles, 120.0, 6.0);
}

#[test]
fn three_coordinate_boron_is_trigonal_planar() {
    // Boron trifluoride, trimethylborane, boric acid, phenylboronic acid.
    for smiles in ["FB(F)F", "CB(C)C", "OB(O)O", "OB(O)c1ccccc1"] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        // The rule-based placement gives 120°, 120° and 97° at any trigonal
        // centre with two unplaced neighbours (acetone's carbonyl carbon too),
        // so only the minimized geometry is checked.
        let etkdg = generate_coords_etkdg(&mol);
        let minimized = minimize_uff(&mol, etkdg.clone());
        assert_trigonal_planar(smiles, "minimize_uff", &boron_angles(&mol, &minimized));
        let (explicit, coords) = add_hydrogens_with_coords(&mol, &etkdg).unwrap();
        let minimized = minimize_uff(&explicit, coords);
        let angles = boron_angles(&explicit, &minimized);
        assert_trigonal_planar(smiles, "minimize_uff with hydrogens", &angles);
    }
}

#[test]
fn four_coordinate_boron_stays_tetrahedral() {
    // Tetrafluoroborate and a trimethylamine–borane adduct.
    for smiles in ["F[B-](F)(F)F", "C[N+](C)(C)[B-](F)(F)F"] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let etkdg = generate_coords_etkdg(&mol);
        let angles = boron_angles(&mol, &etkdg);
        assert_angles(smiles, "generate_coords_etkdg", &angles, 109.5, 2.0);
        let minimized = minimize_uff(&mol, etkdg);
        let angles = boron_angles(&mol, &minimized);
        assert_angles(smiles, "minimize_uff", &angles, 109.5, 2.0);
    }
}
