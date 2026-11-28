//! Neutral three-coordinate boron is trigonal planar (120°), not tetrahedral:
//! RDKit perceives it as SP2 and types it `B_2` for UFF, even with single bonds only.

use chematic_3d::{
    Coords3D, add_hydrogens_with_coords, generate_and_minimize_dreiding, generate_coords,
    generate_coords_etkdg, get_bond_angle_deg, minimize_uff,
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

#[test]
fn acyclic_trigonal_boron_branches_start_planar() {
    // Cover a boron root, an already placed parent, and alternate atom orders.
    for smiles in ["B(F)(F)F", "FB(F)F", "CB(C)C", "OB(O)O", "ClB(F)Br"] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let coords = generate_coords(&mol);
        assert_angles(
            smiles,
            "initial placement",
            &boron_angles(&mol, &coords),
            120.0,
            0.01,
        );
        let coords = generate_coords_etkdg(&mol);
        assert_trigonal_planar(smiles, "ETKDG placement", &boron_angles(&mol, &coords));
    }
}

#[test]
fn development_boronic_acid_starts_with_planar_boron() {
    // Exposed COD development record 2877; the aromatic attachment is already
    // placed before the two acyclic B-O branches.
    let smiles = "C[N@@H+](CCNc1ccc2c3c(cccc13)C(=O)N(CCOCCO)C2=O)Cc1ccccc1B(O)O";
    let mol = chematic_smiles::parse(smiles).unwrap();
    for (stage, coords) in [
        ("initial placement", generate_coords(&mol)),
        ("ETKDG placement", generate_coords_etkdg(&mol)),
    ] {
        assert_trigonal_planar(smiles, stage, &boron_angles(&mol, &coords));
    }
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
fn dreiding_entry_point_preserves_trigonal_boron() {
    for smiles in ["FB(F)F", "CB(C)C", "OB(O)O", "OB(O)c1ccccc1"] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let coords = generate_and_minimize_dreiding(&mol);
        assert_trigonal_planar(
            smiles,
            "generate_and_minimize_dreiding",
            &boron_angles(&mol, &coords),
        );
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
