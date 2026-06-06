//! `chematic-3d` — 3D coordinate generation and file formats for chematic.
//!
//! Provides:
//! - [`generate_coords`]: rule-based 3D coordinate builder.
//! - [`parse_pdb_atoms`], [`pdb_to_molecule`], [`write_pdb`]: PDB file support.
//! - [`parse_xyz`], [`write_xyz`]: XYZ file support.

#![forbid(unsafe_code)]

pub mod align;
pub mod usr;
pub mod conformer;
pub mod coords;
pub mod dg;
pub mod md;
pub mod minimize;
pub mod pdb;
pub mod shape_descriptors;
pub mod stereo3d;
pub mod xyz;

pub use align::{AlignResult, align_coords, apply_alignment, rmsd_no_align};
pub use usr::{usr_descriptors, usr_similarity};
pub use conformer::{ConformerEnsemble, ConformerError};
pub use coords::{Coords3D, Point3};
pub use dg::generate_coords;
pub use md::{MDConfig, MDFrame, MDTrajectory, Thermostat, run_md};
pub use minimize::{MinimizeConfig, minimize, minimize_uff, minimize_with_config, minimize_dreiding, minimize_dreiding_with_config};
pub use pdb::{PdbAtom, parse_pdb_atoms, pdb_to_molecule, write_pdb};
pub use shape_descriptors::{
    asphericity, eccentricity, npr1, npr2, plane_of_best_fit,
    pmi, pmi1, pmi2, pmi3, radius_of_gyration,
};
pub use stereo3d::{StereoAssignment3D, assign_stereo_from_3d};
pub use xyz::{XyzError, parse_xyz, write_xyz};

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use chematic_core::AtomIdx;
    use chematic_smiles::parse;

    use crate::{
        coords::Point3,
        dg::generate_coords,
        pdb::{parse_pdb_atoms, pdb_to_molecule, write_pdb},
        xyz::{XyzError, parse_xyz, write_xyz},
    };

    // -----------------------------------------------------------------------
    // Coords / Point3 tests
    // -----------------------------------------------------------------------

    /// Test 1: Point3 distance.
    #[test]
    fn test_point3_distance() {
        let a = Point3::new(3.0, 4.0, 0.0);
        let b = Point3::zero();
        let d = a.distance(&b);
        assert!((d - 5.0).abs() < 1e-10, "expected 5.0, got {d}");
    }

    /// Test 2: Point3 cross product — (1,0,0) × (0,1,0) = (0,0,1).
    #[test]
    fn test_point3_cross_product() {
        let x = Point3::new(1.0, 0.0, 0.0);
        let y = Point3::new(0.0, 1.0, 0.0);
        let z = x.cross(&y);
        assert!((z.x - 0.0).abs() < 1e-10);
        assert!((z.y - 0.0).abs() < 1e-10);
        assert!((z.z - 1.0).abs() < 1e-10);
    }

    // -----------------------------------------------------------------------
    // DG / generate_coords tests
    // -----------------------------------------------------------------------

    /// Test 3: Single atom placed at origin.
    #[test]
    fn test_single_atom_at_origin() {
        let mol = parse("O").expect("oxygen SMILES");
        let coords = generate_coords(&mol);
        assert_eq!(coords.atom_count(), 1);
        let p = coords.get(AtomIdx(0));
        assert!((p.x).abs() < 1e-10 && (p.y).abs() < 1e-10 && (p.z).abs() < 1e-10);
    }

    /// Test 4: Ethane — 2 distinct atoms, distance ≈ 1.54 Å (±0.1).
    #[test]
    fn test_ethane_bond_length() {
        let mol = parse("CC").expect("ethane SMILES");
        let coords = generate_coords(&mol);
        assert_eq!(coords.atom_count(), 2);
        let p0 = coords.get(AtomIdx(0));
        let p1 = coords.get(AtomIdx(1));
        let d = p0.distance(&p1);
        assert!(
            (d - 1.54).abs() < 0.1,
            "ethane C-C distance expected ~1.54, got {d}"
        );
    }

    /// Test 5: Propane — 3 distinct atoms, no two identical.
    #[test]
    fn test_propane_distinct_atoms() {
        let mol = parse("CCC").expect("propane SMILES");
        let coords = generate_coords(&mol);
        assert_eq!(coords.atom_count(), 3);
        let positions: Vec<_> = (0..3).map(|i| coords.get(AtomIdx(i))).collect();
        for i in 0..3 {
            for j in (i + 1)..3 {
                let d = positions[i].distance(&positions[j]);
                assert!(d > 0.1, "atoms {i} and {j} are too close (d={d:.4})");
            }
        }
    }

    /// Test 6: Benzene — 6 distinct atoms, all within 2.0 Å of centroid.
    #[test]
    fn test_benzene_ring() {
        let mol = parse("c1ccccc1").expect("benzene SMILES");
        let coords = generate_coords(&mol);
        assert_eq!(coords.atom_count(), 6);

        // Compute centroid.
        let cx = (0..6).map(|i| coords.get(AtomIdx(i)).x).sum::<f64>() / 6.0;
        let cy = (0..6).map(|i| coords.get(AtomIdx(i)).y).sum::<f64>() / 6.0;
        let cz = (0..6).map(|i| coords.get(AtomIdx(i)).z).sum::<f64>() / 6.0;
        let centroid = Point3::new(cx, cy, cz);

        for i in 0..6 {
            let p = coords.get(AtomIdx(i));
            let d = p.distance(&centroid);
            assert!(
                d < 2.0,
                "benzene atom {i} is {d:.3} Å from centroid, expected < 2.0"
            );
        }
    }

    /// Test 7: Water — 1 heavy atom at origin (H are implicit).
    #[test]
    fn test_water_single_atom() {
        let mol = parse("O").expect("water SMILES");
        assert_eq!(mol.atom_count(), 1, "water has 1 heavy atom");
        let coords = generate_coords(&mol);
        assert_eq!(coords.atom_count(), 1);
        let p = coords.get(AtomIdx(0));
        assert!((p.x).abs() < 1e-10 && (p.y).abs() < 1e-10 && (p.z).abs() < 1e-10);
    }

    /// Test 8: Disconnected "CC.CC" — 4 distinct atoms.
    #[test]
    fn test_disconnected_four_atoms() {
        let mol = parse("CC.CC").expect("disconnected ethane SMILES");
        assert_eq!(mol.atom_count(), 4);
        let coords = generate_coords(&mol);
        assert_eq!(coords.atom_count(), 4);

        // All four positions must be distinct.
        let positions: Vec<_> = (0..4).map(|i| coords.get(AtomIdx(i))).collect();
        for i in 0..4 {
            for j in (i + 1)..4 {
                let d = positions[i].distance(&positions[j]);
                assert!(d > 0.1, "atoms {i} and {j} overlap (d={d:.4})");
            }
        }
    }

    // -----------------------------------------------------------------------
    // XYZ tests
    // -----------------------------------------------------------------------

    /// Test 9: Write then parse roundtrip for methane — 1 atom, coord ≈ (0,0,0).
    #[test]
    fn test_xyz_roundtrip_methane() {
        let mol = parse("C").expect("methane SMILES");
        let coords = generate_coords(&mol);
        let xyz_str = write_xyz(&mol, &coords, "methane");

        let (mol2, coords2) = parse_xyz(&xyz_str).expect("roundtrip parse");
        assert_eq!(mol2.atom_count(), 1);
        let p = coords2.get(AtomIdx(0));
        assert!((p.x).abs() < 1e-6 && (p.y).abs() < 1e-6 && (p.z).abs() < 1e-6);
    }

    /// Test 10: Write ethane, parse back — 2 atoms, distance preserved (±0.01).
    #[test]
    fn test_xyz_ethane_roundtrip_distance() {
        let mol = parse("CC").expect("ethane SMILES");
        let coords = generate_coords(&mol);
        let orig_dist = coords.get(AtomIdx(0)).distance(&coords.get(AtomIdx(1)));

        let xyz_str = write_xyz(&mol, &coords, "ethane");
        let (mol2, coords2) = parse_xyz(&xyz_str).expect("roundtrip parse");
        assert_eq!(mol2.atom_count(), 2);

        let parsed_dist = coords2.get(AtomIdx(0)).distance(&coords2.get(AtomIdx(1)));
        assert!(
            (parsed_dist - orig_dist).abs() < 0.01,
            "distance changed: orig={orig_dist:.6}, parsed={parsed_dist:.6}"
        );
    }

    /// Test 11: parse_xyz returns error on invalid atom count line.
    #[test]
    fn test_xyz_invalid_atom_count() {
        let bad = "not_a_number\ncomment\n";
        let result = parse_xyz(bad);
        assert!(
            matches!(result, Err(XyzError::InvalidAtomCount)),
            "expected InvalidAtomCount error, got {:?}",
            result.err()
        );
    }

    /// Test 12: write_xyz first line is the atom count as a string.
    #[test]
    fn test_xyz_first_line_is_count() {
        let mol = parse("CCC").expect("propane SMILES");
        let coords = generate_coords(&mol);
        let xyz_str = write_xyz(&mol, &coords, "propane");
        let first_line = xyz_str.lines().next().unwrap();
        assert_eq!(first_line.trim(), "3");
    }

    // -----------------------------------------------------------------------
    // PDB tests
    // -----------------------------------------------------------------------

    /// Test 13: parse_pdb_atoms on a minimal HETATM record.
    #[test]
    fn test_pdb_parse_minimal_hetatm() {
        // Standard 80-column PDB HETATM line with known values.
        let pdb_line = "HETATM    1  C   LIG A   1       1.000   2.000   3.000  1.00  0.00           C\n";
        let atoms = parse_pdb_atoms(pdb_line);
        assert_eq!(atoms.len(), 1);
        let a = &atoms[0];
        assert_eq!(a.serial, 1);
        assert!((a.x - 1.0).abs() < 1e-3, "x={}", a.x);
        assert!((a.y - 2.0).abs() < 1e-3, "y={}", a.y);
        assert!((a.z - 3.0).abs() < 1e-3, "z={}", a.z);
        assert_eq!(a.element.trim(), "C");
    }

    /// Test 14: write_pdb then parse_pdb_atoms roundtrip preserves count and coords.
    #[test]
    fn test_pdb_write_parse_roundtrip() {
        let mol = parse("CCO").expect("ethanol SMILES");
        let coords = generate_coords(&mol);

        let pdb_str = write_pdb(&mol, &coords);
        let parsed = parse_pdb_atoms(&pdb_str);

        assert_eq!(parsed.len(), mol.atom_count());

        // Compare coordinates to within 0.001 Å.
        for i in 0..mol.atom_count() {
            let orig = coords.get(AtomIdx(i as u32));
            let p = &parsed[i];
            assert!(
                (p.x - orig.x).abs() < 0.001,
                "atom {i} x mismatch: orig={:.3} parsed={:.3}",
                orig.x,
                p.x
            );
            assert!(
                (p.y - orig.y).abs() < 0.001,
                "atom {i} y mismatch: orig={:.3} parsed={:.3}",
                orig.y,
                p.y
            );
            assert!(
                (p.z - orig.z).abs() < 0.001,
                "atom {i} z mismatch: orig={:.3} parsed={:.3}",
                orig.z,
                p.z
            );
        }
    }

    /// Test 15: pdb_to_molecule for two C atoms 1.54 Å apart — 2 atoms, 1 bond.
    #[test]
    fn test_pdb_to_molecule_bonding() {
        let pdb = "HETATM    1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n\
                   HETATM    2  C   LIG A   1       1.540   0.000   0.000  1.00  0.00           C\n\
                   END\n";
        let atoms = parse_pdb_atoms(pdb);
        let (mol, _coords) = pdb_to_molecule(&atoms);
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.bond_count(), 1);
    }
}
