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
pub(crate) mod clock;
pub mod conformer;
pub mod constraints;
pub mod coords;
pub mod descriptors_3d;
pub mod determine_bonds;
pub mod dg;
pub mod dg_connectivity_ordered;
pub mod dg_fft;
pub mod distance_geometry_v2;
pub mod ensemble_v2;
pub mod etkdg;
pub mod etkdg_knowledge;
#[cfg(test)]
mod issue256_255_phase2_evaluation;
pub mod md;
pub mod minimize;
pub mod mol_transforms;
pub mod o3a;
pub mod pdb;
pub mod pharmacophore_fp_3d;
pub mod pipeline_v2;
pub(crate) mod prng;
pub mod rdkit_shape_descriptors;
pub mod sasa;
pub mod shape_descriptors;
pub mod stereo3d;
pub mod stereo_constraints;
pub mod torsion_motif;
pub mod usr;
pub mod xyz;

pub use align::{AlignResult, align_coords, apply_alignment, rmsd_no_align};
pub use conformer::{ConformerEnsemble, ConformerError, rmsd_symmetric};
pub use constraints::{
    AngleConstraint, BondConstraint, ConstraintSet, build_constraints, satisfy_constraints,
};
pub use determine_bonds::{
    DetermineError, MAX_ATOMS as DETERMINE_BONDS_MAX_ATOMS, determine_bonds,
};
// Note: ConformerConfig is defined in lib.rs and exported here
pub use coords::{Coords3D, Point3};
pub use descriptors_3d::{
    autocorr_3d, getaway_descriptors, rdf_descriptors, whim_descriptors, whim_getaway_combined,
};
pub use dg::generate_coords;
pub use dg_connectivity_ordered::generate_coords_connectivity_ordered;
pub use ensemble_v2::{
    ConformerAttempt, ConformerDisposition, ConformerSuccess, EnsembleTermination,
    EnsembleV2Config, EnsembleV2ConfigError, EnsembleV2Result, embed_ensemble_v2,
};
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
pub use o3a::{O3AError, O3AResult, o3a_align};
pub use pdb::{
    PdbAtom, PdbParseLimits, PdbResourceLimitError, parse_pdb_atoms, parse_pdb_atoms_with_limits,
    pdb_to_molecule, write_pdb,
};
pub use pharmacophore_fp_3d::{pharmacophore_fp_3d, tanimoto_pharmacophore_3d};
pub use pipeline_v2::{
    PipelineV2Config, PipelineV2Failure, PipelineV2FailureCause, PipelineV2Result,
    RingTorsionApplicationPolicy, StereoPolicy, embed_pipeline_v2,
};
pub use rdkit_shape_descriptors::{
    DescriptorValue, MacrocycleStatus, MacrocycleWarning, RdkitDescriptorError,
    detect_macrocycle_status, rdkit_asphericity, rdkit_eccentricity, rdkit_inertial_shape_factor,
    rdkit_npr1, rdkit_npr2, rdkit_pbf, rdkit_pmi1, rdkit_pmi2, rdkit_pmi3,
    rdkit_radius_of_gyration, rdkit_spherocity_index,
};
pub use sasa::{
    PerElementSasa, SasaDescriptor, calc_mol_sasa, calc_mol_sasa_with_probe, sasa, sasa_descriptor,
    sasa_descriptor_from_dg, sasa_from_dg, sasa_per_atom, sasa_per_atom_from_dg, sasa_per_element,
    sasa_per_element_from_dg, shrake_rupley_sasa,
};
pub use shape_descriptors::{
    asphericity, eccentricity, npr1, npr2, plane_of_best_fit, pmi, pmi1, pmi2, pmi3,
    radius_of_gyration,
};
pub use stereo3d::{StereoAssignment3D, assign_stereo_from_3d};
pub use torsion_motif::{
    TorsionEnvironment, TorsionHistogram, TorsionMotif, TorsionProfileFit, VonMisesComponent,
    extract_torsion_motifs, fit_von_mises_mixture, motif_angles_deg, torsion_profile_distance,
    torsion_profile_to_json,
};
pub use usr::{shape_screen, usr_descriptors, usr_from_dg, usr_similarity};
pub use xyz::{XyzError, XyzParseLimits, parse_xyz, parse_xyz_with_limits, write_xyz};

// ---------------------------------------------------------------------------
// Configuration types
// ---------------------------------------------------------------------------

/// Force field used for geometry minimization during conformer ensemble generation.
#[derive(Clone, Debug, Default)]
pub enum ConformerForceField {
    /// DREIDING: fast, suitable for large-scale screening.
    #[default]
    Dreiding,
    /// MMFF94 (Halgren 1996): higher accuracy for drug-like molecules.
    /// Slower than DREIDING but produces better geometries.
    Mmff94,
}

/// Configuration for conformer ensemble generation.
///
/// - `count`: number of conformers to attempt (before RMSD pruning).
/// - `rmsd_threshold`: minimum Kabsch-aligned RMSD (Å) between retained conformers.
///   Set to 0.0 to disable pruning.
/// - `force_field`: minimization engine after ETKDG placement.
/// - `noise_sigma_deg`: standard deviation of Gaussian torsion noise (degrees).
///   Default 30°; set to 0.0 for deterministic single-conformer generation.
#[derive(Clone, Debug)]
pub struct ConformerConfig {
    pub count: usize,
    pub rmsd_threshold: f64,
    pub force_field: ConformerForceField,
    pub noise_sigma_deg: f64,
}

impl Default for ConformerConfig {
    fn default() -> Self {
        Self {
            count: 1,
            rmsd_threshold: 0.5,
            force_field: ConformerForceField::Dreiding,
            noise_sigma_deg: 30.0,
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

#[deprecated(
    note = "misnamed: this runs neither chematic-ff's real UFF nor this crate's own typed \
            DREIDING engine (`generate_and_minimize_dreiding`) — it runs chematic-3d's own \
            internal generic, element-pair-parameterized harmonic force field \
            (`minimize::minimize()`; `minimize_with_config`'s dispatch only special-cases \
            MMFF94, so `ForceField::UFF`/`ForceField::DREIDING` are indistinguishable and \
            both fall through to this same generic engine). Use \
            `generate_and_minimize_dreiding` for typed DREIDING physics, or \
            `minimize::minimize_with_policy_gated(..., minimize::ForceFieldPolicy::UffOnly, ...)` \
            / `chematic_ff::uff::{assign_uff_types, minimize_uff}` for real UFF physics. \
            See issue #204."
)]
/// Generate 3D coordinates and minimize geometry in one step.
///
/// **Despite the name, this does not run chematic-ff's UFF — and it is *not* equivalent to
/// [`generate_and_minimize_dreiding`] either.** `minimize_uff` here resolves to this crate's
/// own [`crate::minimize::minimize_uff`], which is honestly documented as an alias for
/// [`crate::minimize::minimize`] (`MinimizeConfig::default()`). That generic `minimize()`
/// dispatches on `config.force_field`, but its match arm only special-cases MMFF94 — both
/// the `UFF` and `DREIDING` enum variants fall through to the same catch-all,
/// `minimize_generic_with_config`, which uses chematic-3d's own internal, untyped,
/// element-pair bond/angle/VDW table (informally "UFF-derived", see `minimize.rs`'s
/// "UFF-derived element parameters" section) — a third, distinct engine from both
/// `chematic_ff::uff::minimize_uff` (real UFF) and [`crate::minimize::minimize_dreiding`]
/// (which assigns real per-atom DREIDING types via `assign_dreiding_types` and is what
/// [`generate_and_minimize_dreiding`] actually calls). This function is kept only so
/// existing external callers (outside this workspace) don't break; do not add new calls
/// to it.
///
/// - For typed DREIDING physics, call [`generate_and_minimize_dreiding`] directly.
/// - For real UFF physics, use
///   [`crate::minimize::minimize_with_policy_gated`] with
///   [`crate::minimize::ForceFieldPolicy::UffOnly`], or call
///   `chematic_ff::uff::{assign_uff_types, minimize_uff}` directly.
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
    generate_conformer_ensemble_with_config(
        mol,
        &ConformerConfig {
            count,
            rmsd_threshold: 0.0, // No pruning for backward compatibility
            ..ConformerConfig::default()
        },
    )
}

/// Generate multiple conformers with force-field minimization and Kabsch-RMSD pruning.
///
/// Pipeline for each attempt (up to `config.count`):
/// 1. ETKDG placement with Gaussian torsion noise (`config.noise_sigma_deg`).
/// 2. Energy minimization with the chosen force field (`config.force_field`).
/// 3. Kabsch-superposition RMSD vs all retained conformers; discard if below
///    `config.rmsd_threshold`.
///
/// Returns `ConformerEnsemble` with the retained set (may be fewer than `count`
/// after pruning).
pub fn generate_conformer_ensemble_with_config(
    mol: chematic_core::Molecule,
    config: &ConformerConfig,
) -> Result<ConformerEnsemble, ConformerError> {
    if config.count == 0 {
        return Ok(ConformerEnsemble::new(mol));
    }

    let mut ensemble = ConformerEnsemble::new(mol);
    let noise_sigma = if config.count > 1 {
        config.noise_sigma_deg
    } else {
        0.0
    };

    for _ in 0..config.count {
        let coords = etkdg::generate_coords_etkdg_with_noise(ensemble.mol(), noise_sigma);
        let minimized = match config.force_field {
            ConformerForceField::Dreiding => minimize_dreiding(ensemble.mol(), coords),
            ConformerForceField::Mmff94 => minimize_mmff94(ensemble.mol(), coords),
        };

        // Kabsch-aligned RMSD pruning: discard near-duplicates.
        if ensemble.is_duplicate(&minimized, config.rmsd_threshold) {
            continue;
        }

        ensemble.add_conformer(minimized)?;
    }

    Ok(ensemble)
}

/// Generate multiple conformers minimized with MMFF94 (Halgren 1996).
///
/// Convenience wrapper around [`generate_conformer_ensemble_with_config`] with
/// `ConformerForceField::Mmff94`.  Higher accuracy than the default DREIDING
/// pipeline, at the cost of ~3–5× longer minimization time.
///
/// ```rust,ignore
/// let ensemble = generate_conformer_ensemble_mmff94(mol, 20, 0.5)?;
/// ```
pub fn generate_conformer_ensemble_mmff94(
    mol: chematic_core::Molecule,
    count: usize,
    rmsd_threshold: f64,
) -> Result<ConformerEnsemble, ConformerError> {
    generate_conformer_ensemble_with_config(
        mol,
        &ConformerConfig {
            count,
            rmsd_threshold,
            force_field: ConformerForceField::Mmff94,
            noise_sigma_deg: 30.0,
        },
    )
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
        generate_conformer_ensemble, generate_conformer_ensemble_mmff94,
        generate_conformer_ensemble_with_config,
        pdb::{parse_pdb_atoms, pdb_to_molecule, write_pdb},
        xyz::{XyzError, XyzParseLimits, parse_xyz, parse_xyz_with_limits, write_xyz},
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
        assert_eq!(
            p1.dot(&p2),
            0.0,
            "perpendicular vectors have zero dot product"
        );

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

    #[test]
    fn test_xyz_parse_limits_reject_input_atoms_and_lines() {
        let xyz = "1\ncomment\nC 0.0 0.0 0.0\n";
        assert!(matches!(
            parse_xyz_with_limits(
                xyz,
                XyzParseLimits {
                    max_input_bytes: 4,
                    ..Default::default()
                }
            ),
            Err(XyzError::InputTooLarge { .. })
        ));
        assert!(matches!(
            parse_xyz_with_limits(
                "999\ncomment\n",
                XyzParseLimits {
                    max_atoms: 10,
                    ..Default::default()
                }
            ),
            Err(XyzError::TooManyAtoms {
                count: 999,
                limit: 10,
            })
        ));
        assert!(matches!(
            parse_xyz_with_limits(
                "1\nx\nC 0.0 0.0 0.0\n",
                XyzParseLimits {
                    max_line_bytes: 3,
                    ..Default::default()
                }
            ),
            Err(XyzError::LineTooLong { line: 3, limit: 3 })
        ));
    }

    // =========================================================================
    // PDB edge cases
    // =========================================================================

    #[test]
    fn test_pdb_atom_record_parsed() {
        // ATOM record (not only HETATM)
        let pdb =
            "ATOM      1  C   ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n";
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
            rmsd_threshold: 0.0,
            ..ConformerConfig::default()
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
            ..ConformerConfig::default()
        };
        let ensemble = generate_conformer_ensemble_with_config(mol, &config)
            .expect("should create empty ensemble");
        assert_eq!(
            ensemble.conformer_count(),
            0,
            "empty config should yield no conformers"
        );
    }

    #[test]
    fn test_conformer_ensemble_rmsd_pruning() {
        use super::ConformerConfig;
        let mol = parse("C").expect("methane SMILES");
        let config = ConformerConfig {
            count: 5,
            rmsd_threshold: 1.0,
            ..ConformerConfig::default()
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
        let ensemble = generate_conformer_ensemble(mol, 2).expect("should generate ensemble");
        assert_eq!(
            ensemble.conformer_count(),
            2,
            "backward-compatible API should work"
        );
    }

    #[test]
    fn test_conformer_ensemble_mmff94() {
        // MMFF94 ensemble must produce at least 1 conformer for a drug-like molecule.
        let mol = parse("c1ccccc1CC(=O)O").expect("phenylacetic acid");
        let ensemble = generate_conformer_ensemble_mmff94(mol, 5, 0.5)
            .expect("MMFF94 ensemble should succeed");
        assert!(
            ensemble.conformer_count() >= 1,
            "MMFF94 ensemble must produce at least 1 conformer"
        );
    }

    #[test]
    fn test_conformer_ensemble_gaussian_diversity() {
        // With Gaussian noise and n>1, we expect diverse conformers for a flexible molecule.
        let mol = parse("CCCCCC").expect("hexane");
        use super::{ConformerConfig, ConformerForceField};
        let config = ConformerConfig {
            count: 10,
            rmsd_threshold: 0.3,
            force_field: ConformerForceField::Dreiding,
            noise_sigma_deg: 30.0,
        };
        let ensemble = generate_conformer_ensemble_with_config(mol, &config).expect("ensemble ok");
        // Hexane has 3 rotatable bonds; 10 attempts with Gaussian noise should yield ≥2 unique conformers.
        assert!(
            ensemble.conformer_count() >= 2,
            "flexible molecule with Gaussian noise should produce diverse conformers, got {}",
            ensemble.conformer_count()
        );
    }

    // -----------------------------------------------------------------------
    // Regression test for issue #204
    // -----------------------------------------------------------------------

    /// `generate_and_minimize_uff`'s name has always claimed UFF, but it has never
    /// invoked `chematic_ff::uff::minimize_uff`. It also turns out NOT to be equivalent
    /// to [`generate_and_minimize_dreiding`] either (an initial version of this test
    /// asserted that and failed): `minimize::minimize()`'s dispatch only special-cases
    /// `ForceField::MMFF94`, so both `UFF` and `DREIDING` variants fall through to the
    /// same generic, untyped, element-pair harmonic engine — a third implementation,
    /// distinct from both real UFF and the typed DREIDING engine
    /// (`minimize::minimize_dreiding`, assigning real `DREIDINGType`s) that
    /// `generate_and_minimize_dreiding` actually runs.
    ///
    /// The only load-bearing, permanently-pinned contract is (1): the deprecated
    /// function's output is numerically identical (within floating-point tolerance) to
    /// calling `minimize::minimize()` directly, confirming exactly what it delegates
    /// to. (2)/(3) below (divergence from typed DREIDING / real UFF on this one
    /// molecule) are logged as diagnostic evidence, not asserted — different force
    /// fields can in principle converge to the same local minimum for a particular
    /// molecule by coincidence, so "always differs" would be a fragile regression gate
    /// unrelated to this function's actual contract.
    #[test]
    #[allow(deprecated)]
    fn generate_and_minimize_uff_delegates_to_generic_minimize() {
        use crate::minimize::{
            ForceFieldPolicy, MinimizeConfig, minimize, minimize_with_policy_gated,
        };
        use crate::{generate_and_minimize_dreiding, generate_and_minimize_uff};

        let mol = parse("c1ccccc1").expect("benzene SMILES");

        let deprecated_result = generate_and_minimize_uff(&mol);

        // (1) LOAD-BEARING: it is exactly `minimize::minimize()` on freshly-generated
        // coords — no more, no less. This is the true, honest description of its
        // current behavior, and the only property this test hard-gates.
        let generic_result = minimize(&mol, generate_coords(&mol));
        for i in 0..mol.atom_count() as u32 {
            let a = deprecated_result.get(AtomIdx(i));
            let b = generic_result.get(AtomIdx(i));
            assert!(
                a.distance(&b) < 1e-12,
                "generate_and_minimize_uff must be numerically identical (within \
                 floating-point tolerance) to minimize::minimize() — that is literally \
                 what it delegates to; atom {i} diverged"
            );
        }

        // (2)/(3) DIAGNOSTIC ONLY — evidence that today, on this one molecule, all
        // three engines happen to be distinct; not a permanent API contract, so never
        // asserted. See the doc comment above for why.
        let dreiding_result = generate_and_minimize_dreiding(&mol);
        let differs_from_dreiding = (0..mol.atom_count() as u32).any(|i| {
            deprecated_result
                .get(AtomIdx(i))
                .distance(&dreiding_result.get(AtomIdx(i)))
                > 1e-6
        });
        println!(
            "[diagnostic] generate_and_minimize_uff differs from typed DREIDING on \
             benzene: {differs_from_dreiding}"
        );

        let real_uff = minimize_with_policy_gated(
            &mol,
            generate_coords(&mol),
            ForceFieldPolicy::UffOnly,
            &MinimizeConfig::default(),
            false,
            false,
        )
        .expect("real UFF minimization of benzene should succeed");
        let differs_from_real_uff = (0..mol.atom_count() as u32).any(|i| {
            deprecated_result
                .get(AtomIdx(i))
                .distance(&real_uff.coords.get(AtomIdx(i)))
                > 1e-6
        });
        println!(
            "[diagnostic] generate_and_minimize_uff differs from real UFF on benzene: \
             {differs_from_real_uff}"
        );
    }
}
