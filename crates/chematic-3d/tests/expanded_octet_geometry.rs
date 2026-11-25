//! Centres with five or six electron domains take their VSEPR shape: trigonal
//! bipyramid (PCl₅), octahedron (SF₆) and the shapes with lone pairs (SF₄,
//! ClF₃, XeF₂, XeF₄, IF₅), in rule-based placement, ETKDG and the legacy
//! minimizer.

use chematic_3d::{
    Coords3D, generate_coords, generate_coords_etkdg, get_bond_angle_deg, minimize_uff,
};
use chematic_core::{AtomIdx, Molecule};

const OCTAHEDRAL: &[(f64, usize)] = &[(90.0, 12), (180.0, 3)];
const TRIGONAL_BIPYRAMIDAL: &[(f64, usize)] = &[(90.0, 6), (120.0, 3), (180.0, 1)];

/// SMILES with the centre as atom 1, and its X–A–X angles: (angle, count).
const CASES: &[(&str, &[(f64, usize)])] = &[
    ("FS(F)(F)(F)(F)F", OCTAHEDRAL),
    ("F[Si-2](F)(F)(F)(F)F", OCTAHEDRAL),
    ("FS(F)(F)(F)(F)c1ccccc1", OCTAHEDRAL),
    ("ClP(Cl)(Cl)(Cl)Cl", TRIGONAL_BIPYRAMIDAL),
    ("FP(F)(F)(F)F", TRIGONAL_BIPYRAMIDAL),
    // Seesaw, T-shaped, linear, square planar, square pyramidal.
    ("FS(F)(F)F", &[(90.0, 4), (120.0, 1), (180.0, 1)]),
    ("FCl(F)F", &[(90.0, 2), (180.0, 1)]),
    ("F[Xe]F", &[(180.0, 1)]),
    ("F[Xe](F)(F)F", &[(90.0, 4), (180.0, 2)]),
    ("FI(F)(F)(F)F", &[(90.0, 8), (180.0, 2)]),
];

fn sorted_angles(mol: &Molecule, coords: &Coords3D) -> Vec<f64> {
    let centre = AtomIdx(1);
    let neighbors: Vec<AtomIdx> = mol.neighbors(centre).map(|(n, _)| n).collect();
    let mut angles = Vec::new();
    for i in 0..neighbors.len() {
        for j in i + 1..neighbors.len() {
            angles.push(get_bond_angle_deg(
                coords,
                neighbors[i],
                centre,
                neighbors[j],
            ));
        }
    }
    angles.sort_by(f64::total_cmp);
    angles
}

fn assert_shape(smiles: &str, stage: &str, angles: &[f64], expected: &[(f64, usize)], tol: f64) {
    let expected: Vec<f64> = expected
        .iter()
        .flat_map(|&(angle, count)| std::iter::repeat_n(angle, count))
        .collect();
    let ok = angles.len() == expected.len()
        && angles
            .iter()
            .zip(&expected)
            .all(|(a, e)| (a - e).abs() < tol);
    assert!(ok, "{smiles}, {stage}: {angles:.1?}, expected {expected:?}");
}

#[test]
fn placement_and_etkdg_follow_vsepr() {
    for &(smiles, expected) in CASES {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let placed = generate_coords(&mol);
        assert_shape(
            smiles,
            "generate_coords",
            &sorted_angles(&mol, &placed),
            expected,
            1.0,
        );
        // No two atoms on top of each other.
        for i in 0..mol.atom_count() {
            for j in i + 1..mol.atom_count() {
                let d = placed.points[i].sub(&placed.points[j]).norm();
                assert!(d > 1.0, "{smiles}: atoms {i} and {j} are {d:.2} Å apart");
            }
        }
        let etkdg = generate_coords_etkdg(&mol);
        let angles = sorted_angles(&mol, &etkdg);
        assert_shape(smiles, "generate_coords_etkdg", &angles, expected, 1.0);
    }
}

#[test]
fn legacy_minimizer_keeps_vsepr_shape() {
    for &(smiles, expected) in CASES {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let minimized = minimize_uff(&mol, generate_coords_etkdg(&mol));
        let angles = sorted_angles(&mol, &minimized);
        assert_shape(smiles, "minimize_uff", &angles, expected, 2.0);
    }
}
