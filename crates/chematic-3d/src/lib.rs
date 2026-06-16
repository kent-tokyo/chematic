//! `chematic-3d` — 3D coordinate generation and file formats for chematic.
//!
//! Provides:
//! - [`generate_coords`]: rule-based 3D coordinate builder.
//! - [`parse_pdb_atoms`], [`pdb_to_molecule`], [`write_pdb`]: PDB file support.
//! - [`parse_xyz`], [`write_xyz`]: XYZ file support.
//! - [`build_constraints`], [`satisfy_constraints`]: distance geometry constraint satisfaction.
//! - [`generate_and_minimize_constrained`]: full pipeline with constraint projection.

#![forbid(unsafe_code)]
#![allow(clippy::needless_range_loop)]

pub mod align;
pub mod conformer;
pub mod constraints;
pub mod coords;
pub mod determine_bonds;
pub mod dg;
pub mod dg_fft;
pub mod descriptors_3d;
pub mod etkdg;
pub mod etkdg_knowledge;
pub mod md;
pub mod minimize;
pub mod mol_transforms;
pub mod pharmacophore_fp_3d;
pub mod pdb;
pub mod sasa;
pub mod shape_descriptors;
pub mod stereo3d;
pub mod usr;
pub mod xyz;

pub use align::{AlignResult, align_coords, apply_alignment, rmsd_no_align};
pub use determine_bonds::{DetermineError, MAX_ATOMS as DETERMINE_BONDS_MAX_ATOMS, determine_bonds};
pub use conformer::{ConformerEnsemble, ConformerError};
pub use constraints::{
    AngleConstraint, BondConstraint, ConstraintSet, build_constraints, satisfy_constraints,
};
// Note: ConformerConfig is defined in lib.rs and exported here
pub use coords::{Coords3D, Point3};
pub use descriptors_3d::{autocorr_3d, getaway_descriptors, whim_descriptors, whim_getaway_combined};
pub use dg::generate_coords;
pub use etkdg::generate_coords_etkdg;
pub use md::{MDConfig, MDFrame, MDTrajectory, Thermostat, run_md};
pub use minimize::{
    ForceField, MinimizeConfig, minimize, minimize_dreiding, minimize_dreiding_with_config,
    minimize_mmff94, minimize_uff, minimize_with_config,
};
pub use mol_transforms::{
    center_on_origin, compute_centroid, get_bond_angle, get_bond_angle_deg, get_bond_length,
    get_dihedral, get_dihedral_deg, set_dihedral, transform_conformer,
};
pub use pharmacophore_fp_3d::{pharmacophore_fp_3d, tanimoto_pharmacophore_3d};
pub use pdb::{PdbAtom, parse_pdb_atoms, pdb_to_molecule, write_pdb};
pub use sasa::{
    sasa, sasa_per_atom, shrake_rupley_sasa, sasa_descriptor, sasa_per_element,
    sasa_from_dg, sasa_per_atom_from_dg, sasa_descriptor_from_dg, sasa_per_element_from_dg,
    calc_mol_sasa, calc_mol_sasa_with_probe, SasaDescriptor, PerElementSasa,
};
pub use shape_descriptors::{
    asphericity, eccentricity, npr1, npr2, plane_of_best_fit, pmi, pmi1, pmi2, pmi3,
    radius_of_gyration,
};
pub use stereo3d::{StereoAssignment3D, assign_stereo_from_3d};
pub use usr::{usr_descriptors, usr_similarity, usr_from_dg, shape_screen};
pub use xyz::{XyzError, parse_xyz, write_xyz};

// ---------------------------------------------------------------------------
// Configuration types
// ---------------------------------------------------------------------------

/// Configuration for conformer ensemble generation.
///
/// - `count`: number of conformers to generate (before pruning)
/// - `rmsd_threshold`: minimum RMSD (Å) between conformers. Conformers with RMSD below this
///   threshold to an already-added conformer are discarded. Set to 0.0 to disable pruning.
#[derive(Clone, Debug)]
pub struct ConformerConfig {
    pub count: usize,
    pub rmsd_threshold: f64,
}

impl Default for ConformerConfig {
    fn default() -> Self {
        Self {
            count: 1,
            rmsd_threshold: 0.5,  // Default: keep conformers at least 0.5 Å apart
        }
    }
}

// ---------------------------------------------------------------------------
// High-level 3D generation pipeline
// ---------------------------------------------------------------------------

/// Generate 3D coordinates and minimize geometry in one step.
/// Uses distance geometry for initial placement + DREIDING force field.
pub fn generate_and_minimize_dreiding(mol: &chematic_core::Molecule) -> Coords3D {
    let coords = generate_coords(mol);
    minimize_dreiding(mol, coords)
}

/// Generate 3D coordinates with constraint satisfaction and energy minimization.
///
/// Full pipeline:
/// 1. Rule-based 3D placement (`generate_coords`)
/// 2. Build bond/angle constraints from topology (`build_constraints`)
/// 3. Iterative constraint projection (`satisfy_constraints`)
/// 4. Energy minimization with DREIDING force field (`minimize_dreiding`)
pub fn generate_and_minimize_constrained(mol: &chematic_core::Molecule) -> Coords3D {
    let coords = generate_coords(mol);
    let cs = build_constraints(mol);
    let projected = satisfy_constraints(&coords, mol, &cs, 20);
    minimize_dreiding(mol, projected)
}

/// Generate 3D coordinates and minimize using UFF force field.
pub fn generate_and_minimize_uff(mol: &chematic_core::Molecule) -> Coords3D {
    let coords = generate_coords(mol);
    minimize_uff(mol, coords)
}

/// Generate multiple conformers with different initial geometries.
/// Uses distance geometry for initial placement, then minimizes with DREIDING.
/// Returns a ConformerEnsemble with all conformers.
///
/// Equivalent to `generate_conformer_ensemble_with_config(mol, ConformerConfig::default())`.
pub fn generate_conformer_ensemble(
    mol: chematic_core::Molecule,
    count: usize,
) -> Result<ConformerEnsemble, ConformerError> {
    let config = ConformerConfig {
        count,
        rmsd_threshold: 0.0,  // No pruning for backward compatibility
    };
    generate_conformer_ensemble_with_config(mol, &config)
}

/// Generate multiple conformers with RMSD-based pruning.
///
/// Generates up to `config.count` conformers via repeated distance geometry + DREIDING minimization.
/// If `rmsd_threshold > 0`, discards conformers with RMSD below the threshold to an
/// already-added conformer, effectively filtering out redundant structures.
///
/// Returns `ConformerEnsemble` with the final set of conformers (may be fewer than `count`
/// if pruning removes duplicates).
pub fn generate_conformer_ensemble_with_config(
    mol: chematic_core::Molecule,
    config: &ConformerConfig,
) -> Result<ConformerEnsemble, ConformerError> {
    if config.count == 0 {
        return Ok(ConformerEnsemble::new(mol));
    }

    let mut ensemble = ConformerEnsemble::new(mol);
    let use_pruning = config.rmsd_threshold > 0.0;

    for _ in 0..config.count {
        let coords = generate_coords(ensemble.mol());
        let minimized = minimize_dreiding(ensemble.mol(), coords);

        // Apply RMSD pruning if enabled
        if use_pruning {
            let mut is_duplicate = false;
            for i in 0..ensemble.conformer_count() {
                if let Some(existing) = ensemble.get_conformer(i) {
                    // Convert Point3 vectors to [f64; 3] arrays for rmsd_no_align
                    let min_array: Vec<[f64; 3]> = minimized.points.iter()
                        .map(|p| [p.x, p.y, p.z])
                        .collect();
                    let exist_array: Vec<[f64; 3]> = existing.points.iter()
                        .map(|p| [p.x, p.y, p.z])
                        .collect();
                    let rmsd = rmsd_no_align(&min_array, &exist_array);
                    if rmsd < config.rmsd_threshold {
                        is_duplicate = true;
                        break;
                    }
                }
            }
            if is_duplicate {
                continue;  // Skip this conformer
            }
        }

        ensemble.add_conformer(minimized)?;
    }

    Ok(ensemble)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use chematic_core::AtomIdx;
    use chematic_smiles::parse;

    use crate::{
        coords::{Coords3D, Point3},
        dg::generate_coords,
        generate_conformer_ensemble, generate_conformer_ensemble_with_config,
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
        let pdb_line =
            "HETATM    1  C   LIG A   1       1.000   2.000   3.000  1.00  0.00           C\n";
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

    // =========================================================================
    // Point3 additional tests
    // =========================================================================

    #[test]
    fn test_point3_zero() {
        let p = Point3::zero();
        assert_eq!(p.x, 0.0);
        assert_eq!(p.y, 0.0);
        assert_eq!(p.z, 0.0);
    }

    #[test]
    fn test_point3_add() {
        let p1 = Point3::new(1.0, 2.0, 3.0);
        let p2 = Point3::new(4.0, 5.0, 6.0);
        let sum = p1.add(&p2);
        assert_eq!(sum.x, 5.0);
        assert_eq!(sum.y, 7.0);
        assert_eq!(sum.z, 9.0);
    }

    #[test]
    fn test_point3_sub() {
        let p1 = Point3::new(5.0, 7.0, 9.0);
        let p2 = Point3::new(1.0, 2.0, 3.0);
        let diff = p1.sub(&p2);
        assert_eq!(diff.x, 4.0);
        assert_eq!(diff.y, 5.0);
        assert_eq!(diff.z, 6.0);
    }

    #[test]
    fn test_point3_scale() {
        let p = Point3::new(1.0, 2.0, 3.0);
        let scaled = p.scale(2.0);
        assert_eq!(scaled.x, 2.0);
        assert_eq!(scaled.y, 4.0);
        assert_eq!(scaled.z, 6.0);
    }

    #[test]
    fn test_point3_dot() {
        let p1 = Point3::new(1.0, 0.0, 0.0);
        let p2 = Point3::new(0.0, 1.0, 0.0);
        assert_eq!(p1.dot(&p2), 0.0, "perpendicular vectors have zero dot product");

        let p3 = Point3::new(1.0, 2.0, 3.0);
        let p4 = Point3::new(1.0, 2.0, 3.0);
        assert_eq!(p3.dot(&p4), 14.0); // 1 + 4 + 9
    }

    #[test]
    fn test_point3_norm() {
        let p = Point3::new(3.0, 4.0, 0.0);
        assert_eq!(p.norm(), 5.0, "3-4-5 triangle");
    }

    #[test]
    fn test_point3_normalize() {
        let p = Point3::new(3.0, 4.0, 0.0);
        let unit = p.normalize();
        assert!((unit.x - 0.6).abs() < 1e-9);
        assert!((unit.y - 0.8).abs() < 1e-9);
        assert_eq!(unit.z, 0.0);
    }

    #[test]
    #[should_panic]
    fn test_point3_normalize_zero_panics() {
        let p = Point3::zero();
        let _ = p.normalize();
    }

    // =========================================================================
    // Coords3D additional tests
    // =========================================================================

    #[test]
    fn test_coords3d_new_zeroed() {
        let coords = Coords3D::new_zeroed(5);
        assert_eq!(coords.atom_count(), 5);
        for i in 0..5 {
            let p = coords.get(AtomIdx(i as u32));
            assert_eq!(p.x, 0.0);
            assert_eq!(p.y, 0.0);
            assert_eq!(p.z, 0.0);
        }
    }

    #[test]
    fn test_coords3d_get_set_roundtrip() {
        let mut coords = Coords3D::new_zeroed(3);
        let p = Point3::new(1.5, 2.5, 3.5);
        coords.set(AtomIdx(1), p);
        let retrieved = coords.get(AtomIdx(1));
        assert_eq!(retrieved.x, 1.5);
        assert_eq!(retrieved.y, 2.5);
        assert_eq!(retrieved.z, 3.5);
    }

    #[test]
    fn test_coords3d_atom_count() {
        let coords = Coords3D::new_zeroed(10);
        assert_eq!(coords.atom_count(), 10);
    }

    // =========================================================================
    // XYZ edge cases
    // =========================================================================

    #[test]
    fn test_xyz_unknown_element() {
        let xyz = "2\n\nXx   0.0 0.0 0.0\nC    1.0 1.0 1.0\n";
        let result = parse_xyz(xyz);
        match result {
            Err(XyzError::UnknownElement(_)) => (),
            _ => panic!("expected UnknownElement error"),
        }
    }

    #[test]
    fn test_xyz_invalid_line() {
        let xyz = "2\n\nC 0.0 0.0\nC 1.0 1.0 1.0\n"; // first atom line too short
        let result = parse_xyz(xyz);
        assert!(matches!(result, Err(XyzError::InvalidLine(_))));
    }

    // =========================================================================
    // PDB edge cases
    // =========================================================================

    #[test]
    fn test_pdb_atom_record_parsed() {
        // ATOM record (not only HETATM)
        let pdb = "ATOM      1  C   ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n";
        let atoms = parse_pdb_atoms(pdb);
        assert_eq!(atoms.len(), 1);
        assert_eq!(atoms[0].element, "C");
    }

    #[test]
    fn test_pdb_remark_skipped() {
        let pdb = "REMARK This is a comment\nHETATM    1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n";
        let atoms = parse_pdb_atoms(pdb);
        assert_eq!(atoms.len(), 1, "only HETATM/ATOM records should be parsed");
    }

    #[test]
    fn test_pdb_write_ends_with_end() {
        use chematic_core::{Atom, Element, MoleculeBuilder};
        let mut builder = MoleculeBuilder::new();
        let c = Atom::new(Element::from_atomic_number(6).unwrap());
        builder.add_atom(c);
        let mol = builder.build();
        let coords = Coords3D::new_zeroed(1);
        let pdb = write_pdb(&mol, &coords);
        assert!(pdb.ends_with("END\n"), "PDB should end with 'END\\n'");
    }

    // =========================================================================
    // Conformer ensemble tests
    // =========================================================================

    #[test]
    fn test_conformer_ensemble_basic() {
        use super::ConformerConfig;
        let mol = parse("CC").expect("ethane SMILES");
        let config = ConformerConfig {
            count: 2,
            rmsd_threshold: 0.0,  // No pruning
        };
        let ensemble = generate_conformer_ensemble_with_config(mol, &config)
            .expect("should generate ensemble");
        assert_eq!(ensemble.conformer_count(), 2, "should have 2 conformers");
    }

    #[test]
    fn test_conformer_ensemble_zero_count() {
        use super::ConformerConfig;
        let mol = parse("CC").expect("ethane SMILES");
        let config = ConformerConfig {
            count: 0,
            rmsd_threshold: 0.0,
        };
        let ensemble = generate_conformer_ensemble_with_config(mol, &config)
            .expect("should create empty ensemble");
        assert_eq!(ensemble.conformer_count(), 0, "empty config should yield no conformers");
    }

    #[test]
    fn test_conformer_ensemble_rmsd_pruning() {
        use super::ConformerConfig;
        let mol = parse("C").expect("methane SMILES");
        let config = ConformerConfig {
            count: 5,
            rmsd_threshold: 1.0,  // High threshold to prune most duplicates
        };
        let ensemble = generate_conformer_ensemble_with_config(mol, &config)
            .expect("should generate ensemble with pruning");
        // With high threshold and simple molecule, should keep very few (often 1)
        assert!(
            ensemble.conformer_count() <= 3,
            "high RMSD threshold should prune duplicates; got {}",
            ensemble.conformer_count()
        );
    }

    #[test]
    fn test_conformer_backward_compatibility() {
        let mol = parse("CC").expect("ethane SMILES");
        let ensemble = generate_conformer_ensemble(mol, 2)
            .expect("should generate ensemble");
        assert_eq!(ensemble.conformer_count(), 2, "backward-compatible API should work");
    }
}
