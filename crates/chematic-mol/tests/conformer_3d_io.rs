//! Integration tests for 3D conformer support in the MOL/SDF readers and
//! writers (3D Breakthrough Program, Agent A -- see
//! `docs/rfcs/3d_breakthrough_master_plan.md` section 1a/3).
//!
//! Every MOL-block fixture below is real RDKit `2026.03.3` output (pinned
//! commit `8afba32ec539dcb2369bc84549d802aca3f7eb39`, same oracle this
//! codebase's other stereo work cites), generated via `AllChem.EmbedMolecule`
//! / `Chem.MolToMolBlock` / `Chem.MolToV3KMolBlock` in an independent venv
//! that never touched the shared repo `.venv` -- see the PR body for the
//! exact generation script. A handful of fixtures are a single, explicitly
//! called-out-in-comments digit edit of that real output (e.g. flipping a
//! wedge's stereo field from Up to Down while leaving every coordinate
//! untouched) to construct a deliberate, otherwise-unobtainable conflict
//! case -- never hand-invented geometry. The hand-authored (non-RDKit)
//! fixtures further down (typed-error inputs, the blank-z fixture) are
//! byte-exact-column-constructed, not eyeballed, matching this crate's own
//! `mol2000.rs` column layout exactly.

use chematic_core::AtomIdx;
use chematic_mol::mol2000::{
    CoordinateDimension, GeometryRank, Stereo3DDiagnostic, read_mol_with_diagnostics,
    write_mol_with_conformer,
};
use chematic_mol::mol3000::{read_mol_v3000_with_diagnostics, write_mol_v3000_with_conformer};
use chematic_mol::sdf::read_sdf_conformer_ensembles;

const ETHANOL_2D: &str = r#"
     RDKit          2D

  3  2  0  0  0  0  0  0  0  0999 V2000
   -1.2990   -0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990   -0.2500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
"#;
const ETHANOL_3D: &str = r#"
     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.8716    0.1993    0.1351 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4606   -0.5156    0.0446 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3264    0.1817   -0.8391 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
"#;
const BROMO_3D: &str = r#"
     RDKit          3D

  5  4  0  0  0  0  0  0  0  0999 V2000
   -0.0403    0.0604    0.0996 C   0  0  2  0  0  0  0  0  0  0  0  0
   -1.3475   -1.2822   -0.4144 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -0.1635    0.3239    1.4269 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3256    1.5434   -0.8199 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    1.8768   -0.6454   -0.2923 I   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  1
  1  4  1  0
  1  5  1  0
M  END
"#;
const L_ALANINE_3D: &str = r#"
     RDKit          3D

  6  5  0  0  0  0  0  0  0  0999 V2000
   -1.0297   -1.3503    0.2491 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2785   -0.2519   -0.4000 C   0  0  1  0  0  0  0  0  0  0  0  0
   -0.9835    1.0686   -0.1221 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1645   -0.2214    0.1141 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6568   -0.9800    0.9364 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.9163    0.7368   -0.4623 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  1
  2  4  1  0
  4  5  2  0
  4  6  1  0
M  END
"#;
const D_ALANINE_3D: &str = r#"
     RDKit          3D

  6  5  0  0  0  0  0  0  0  0999 V2000
   -0.9715   -1.2914   -0.4126 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2432   -0.2182    0.2939 C   0  0  2  0  0  0  0  0  0  0  0  0
   -0.8799    1.1176   -0.0660 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2483   -0.2786   -0.0635 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8615   -1.2706   -0.4330 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.9119    0.8740    0.1545 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  6
  2  4  1  0
  4  5  2  0
  4  6  1  0
M  END
"#;
const BENZENE_3D_FLAT: &str = r#"
     RDKit          3D

  6  6  0  0  0  0  0  0  0  0999 V2000
   -0.3032   -1.3614   -0.0088 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3305   -0.4183    0.0150 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0274    0.9431    0.0238 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3032    1.3614    0.0088 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3305    0.4183   -0.0150 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0274   -0.9431   -0.0238 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  1  0
  4  5  2  0
  5  6  1  0
  6  1  2  0
M  END
"#;
const CHLORO_BROMO_ETHENE_E_3D: &str = r#"
     RDKit          3D

  4  3  0  0  0  0  0  0  0  0999 V2000
   -2.1393    0.3274    0.0050 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -0.5726   -0.3650    0.0053 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5376    0.3710   -0.0052 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2494   -0.3408   -0.0053 Br  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  1  0
M  END
"#;
const CHLORO_BROMO_ETHENE_Z_3D: &str = r#"
     RDKit          3D

  4  3  0  0  0  0  0  0  0  0999 V2000
   -1.5710    1.2629    0.2063 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -0.6474   -0.1530   -0.0660 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6148   -0.1877   -0.4942 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6973    1.2603   -0.9051 Br  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  1  0
M  END
"#;
const BUTANE_ENSEMBLE_SDF: &str = r#"butane
     RDKit          3D

  4  3  0  0  0  0  0  0  0  0999 V2000
    1.8617   -0.0700    0.3396 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5248    0.4625   -0.1172 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5372   -0.4649    0.4905 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8695    0.0995    0.0112 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
M  END
$$$$
butane
     RDKit          3D

  4  3  0  0  0  0  0  0  0  0999 V2000
   -1.9282    0.0037    0.1782 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5239    0.4481   -0.0991 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5235   -0.5632    0.3245 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8763    0.0475   -0.0281 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
M  END
$$$$
butane
     RDKit          3D

  4  3  0  0  0  0  0  0  0  0999 V2000
   -1.5291   -0.5136   -0.1631 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6384    0.6532   -0.4387 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6203    0.6882    0.3875 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4982   -0.4940    0.2014 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
M  END
$$$$
"#;
const BROMO_3D_V3K: &str = r#"
     RDKit          3D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.080022 0.044032 0.115420 0 CFG=2
M  V30 2 Br -1.355208 -1.285319 -0.433136 0
M  V30 3 F -0.162340 0.349268 1.452720 0
M  V30 4 Cl -0.307369 1.530691 -0.821563 0
M  V30 5 I 1.904939 -0.638672 -0.313441 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 3 CFG=1
M  V30 3 1 1 4
M  V30 4 1 1 5
M  V30 END BOND
M  V30 END CTAB
M  END
"#;
const ENHANCED_STEREO_3D_V3K: &str = r#"enhanced_stereo_3d
  chematic

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.080022 0.044032 0.115420 0
M  V30 2 Br -1.355208 -1.285319 -0.433136 0
M  V30 3 F -0.162340 0.349268 1.452720 0
M  V30 4 Cl -0.307369 1.530691 -0.821563 0
M  V30 5 I 1.904939 -0.638672 -0.313441 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 3
M  V30 3 1 1 4
M  V30 4 1 1 5
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(1 1)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
"#;

// Hand-authored, byte-exact-column-constructed (not RDKit output, not
// eyeballed -- see the generator script referenced in the PR body).
const NAN_Z_MOL: &str = r#"bad
  chematic

  1  0  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000       NaN C   0  0  0  0  0  0  0  0  0  0  0
M  END
"#;
const INF_Z_MOL: &str = r#"bad
  chematic

  1  0  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000     1e400 C   0  0  0  0  0  0  0  0  0  0  0
M  END
"#;
const GARBLED_Z_MOL: &str = r#"bad
  chematic

  1  0  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    !!!!!! C   0  0  0  0  0  0  0  0  0  0  0
M  END
"#;
const BLANK_Z_MOL: &str = r#"blank_z
  chematic

  1  0  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000           C   0  0  0  0  0  0  0  0  0  0  0
M  END
"#;

// ---------------------------------------------------------------------------
// Root cause regression: Z coordinate is actually read now (was silently
// discarded entirely before this PR -- see PR body).
// ---------------------------------------------------------------------------

#[test]
fn v2000_3d_flagged_file_populates_conformer_with_real_z() {
    let report = read_mol_with_diagnostics(ETHANOL_3D).expect("parse");
    let conformer = report.conformer.expect("3D file must produce a conformer");
    assert_eq!(conformer.atom_count(), 3);
    // Spot-check against the literal RDKit-authored z values (not 0.0, the
    // pre-fix silent-discard behavior).
    assert!((conformer.points[0].z - 0.1351).abs() < 1e-6);
    assert!((conformer.points[1].z - 0.0446).abs() < 1e-6);
    assert!((conformer.points[2].z - (-0.8391)).abs() < 1e-6);
}

#[test]
fn v2000_2d_flagged_file_has_no_conformer_and_flatzero_rank() {
    let report = read_mol_with_diagnostics(ETHANOL_2D).expect("parse");
    assert!(report.conformer.is_none());
    assert_eq!(report.geometry_rank, GeometryRank::FlatZero);
    assert_eq!(report.coordinate_dimension, CoordinateDimension::TwoD);
    assert!(report.stereo3d_diagnostics.is_empty());
}

#[test]
fn v3000_3d_flagged_file_populates_conformer_with_real_z() {
    let report = read_mol_v3000_with_diagnostics(BROMO_3D_V3K).expect("parse");
    let conformer = report
        .conformer
        .expect("3D V3000 file must produce a conformer");
    assert_eq!(conformer.atom_count(), 5);
    assert!((conformer.points[0].z - 0.115420).abs() < 1e-6);
    assert!((conformer.points[2].z - 1.452720).abs() < 1e-6);
    assert_eq!(report.coordinate_dimension, CoordinateDimension::ThreeD);
}

// ---------------------------------------------------------------------------
// GeometryRank: distinguishing FlatZero / Coplanar / ThreeD (not one boolean)
// ---------------------------------------------------------------------------

#[test]
fn geometry_rank_threed_for_nonplanar_molecule() {
    // L-alanine's real embedded conformer is genuinely non-planar (6 atoms,
    // sp3 stereocenter) -- independently confirmed via a best-fit-plane
    // residual check in the venv (~2.08 A max deviation, see PR body).
    let report = read_mol_with_diagnostics(L_ALANINE_3D).expect("parse");
    assert_eq!(report.geometry_rank, GeometryRank::ThreeD);
    assert_eq!(report.coordinate_dimension, CoordinateDimension::ThreeD);
    // Declared and observed agree -- no dimension-mismatch diagnostic.
    assert!(!report.stereo3d_diagnostics.iter().any(|d| matches!(
        d,
        Stereo3DDiagnostic::DeclaredTwoDButNonzeroZ { .. }
            | Stereo3DDiagnostic::DeclaredThreeDButFlat { .. }
    )));
}

#[test]
fn geometry_rank_coplanar_for_a_real_but_flat_3d_embedding() {
    // Benzene's real embedded conformer is genuinely 3D-generated but the
    // molecule itself is flat -- this must NOT be conflated with "every z
    // is exactly 0" (FlatZero): independently confirmed via a best-fit-plane
    // residual check in the venv (~2.9e-6 A, see PR body) -- Coplanar, not
    // FlatZero (z values are small but nonzero).
    let report = read_mol_with_diagnostics(BENZENE_3D_FLAT).expect("parse");
    assert_eq!(report.geometry_rank, GeometryRank::Coplanar);
    assert_eq!(report.coordinate_dimension, CoordinateDimension::ThreeD);
    // A real, physically-flat 3D molecule is not a bug -- but it IS still
    // surfaced distinctly via the payload (Coplanar, not FlatZero), per the
    // Coordinator's "don't collapse to one boolean" requirement.
    assert!(report.stereo3d_diagnostics.iter().any(|d| matches!(
        d,
        Stereo3DDiagnostic::DeclaredThreeDButFlat {
            observed: GeometryRank::Coplanar
        }
    )));
}

// ---------------------------------------------------------------------------
// Header-vs-geometry dimension mismatch (cases 1/2/3 from the Coordinator's
// review, each independently identifiable via the diagnostic payload)
// ---------------------------------------------------------------------------

#[test]
fn declared_2d_but_nonzero_z_is_diagnosed() {
    // Same real RDKit 3D coordinates as ETHANOL_3D, header dimensional code
    // hand-edited from "3D" to "2D" (single substring replace, geometry
    // untouched) -- case 1.
    let hacked = ETHANOL_3D.replacen("          3D", "          2D", 1);
    let report = read_mol_with_diagnostics(&hacked).expect("parse");
    assert_eq!(report.coordinate_dimension, CoordinateDimension::TwoD);
    assert!(
        report.conformer.is_some(),
        "real z data must still be captured"
    );
    assert!(
        report
            .stereo3d_diagnostics
            .iter()
            .any(|d| matches!(d, Stereo3DDiagnostic::DeclaredTwoDButNonzeroZ { .. }))
    );
}

#[test]
fn declared_3d_but_every_z_exactly_zero_is_diagnosed_as_flatzero() {
    // Same real RDKit 2D coordinates as ETHANOL_2D (z is 0.0000 for every
    // atom), header dimensional code hand-edited from "2D" to "3D" -- case 3,
    // the "every z is literally 0" sub-case.
    let hacked = ETHANOL_2D.replacen("          2D", "          3D", 1);
    let report = read_mol_with_diagnostics(&hacked).expect("parse");
    assert_eq!(report.coordinate_dimension, CoordinateDimension::ThreeD);
    assert!(
        report.conformer.is_none(),
        "all-zero z is not a real conformer"
    );
    assert!(report.stereo3d_diagnostics.iter().any(|d| matches!(
        d,
        Stereo3DDiagnostic::DeclaredThreeDButFlat {
            observed: GeometryRank::FlatZero
        }
    )));
}

// ---------------------------------------------------------------------------
// Wedge vs. 3D geometry: tetrahedral (4 explicit neighbors)
// ---------------------------------------------------------------------------

#[test]
fn wedge_agrees_with_3d_geometry_for_real_rdkit_output_no_conflict() {
    // BROMO_3D is RDKit's OWN real, self-consistent output: the wedge on the
    // C-F bond and the real 3D geometry both encode the SAME actual
    // stereocenter. Zero conflicts expected.
    let report = read_mol_with_diagnostics(BROMO_3D).expect("parse");
    assert!(report.conformer.is_some());
    assert_ne!(
        report.mol.atom(AtomIdx(0)).chirality,
        chematic_core::Chirality::None,
        "wedge must have been perceived"
    );
    assert!(
        !report
            .stereo3d_diagnostics
            .iter()
            .any(|d| matches!(d, Stereo3DDiagnostic::WedgeVs3DParityConflict { .. })),
        "a real, self-consistent RDKit record must not conflict with itself: {:?}",
        report.stereo3d_diagnostics
    );
}

#[test]
fn wedge_disagrees_with_3d_geometry_conflict_fires() {
    // Same real geometry as BROMO_3D, wedge stereo field flipped from Up (1)
    // to Down (6) on the C-F bond line -- the ONLY edit, geometry untouched.
    // This can only disagree with the (unchanged) real 3D geometry.
    let hacked = BROMO_3D.replacen("  1  3  1  1\n", "  1  3  1  6\n", 1);
    assert_ne!(
        hacked, BROMO_3D,
        "the wedge line must actually have been edited"
    );
    let report = read_mol_with_diagnostics(&hacked).expect("parse");
    assert!(report.conformer.is_some());
    let conflicts: Vec<_> = report
        .stereo3d_diagnostics
        .iter()
        .filter(|d| matches!(d, Stereo3DDiagnostic::WedgeVs3DParityConflict { .. }))
        .collect();
    assert_eq!(
        conflicts.len(),
        1,
        "flipping only the wedge must produce exactly one conflict: {:?}",
        report.stereo3d_diagnostics
    );
}

#[test]
fn v3000_wedge_agrees_with_3d_geometry_no_conflict() {
    let report = read_mol_v3000_with_diagnostics(BROMO_3D_V3K).expect("parse");
    assert!(report.conformer.is_some());
    assert!(
        !report
            .stereo3d_diagnostics
            .iter()
            .any(|d| matches!(d, Stereo3DDiagnostic::WedgeVs3DParityConflict { .. }))
    );
}

#[test]
fn v3000_wedge_disagrees_with_3d_geometry_conflict_fires() {
    // Flip CFG=1 (Up) to CFG=3 (Down) on the C-F bond line only.
    let hacked = BROMO_3D_V3K.replacen("CFG=1", "CFG=3", 1);
    assert_ne!(hacked, BROMO_3D_V3K);
    let report = read_mol_v3000_with_diagnostics(&hacked).expect("parse");
    assert!(
        report
            .stereo3d_diagnostics
            .iter()
            .any(|d| matches!(d, Stereo3DDiagnostic::WedgeVs3DParityConflict { .. }))
    );
}

// ---------------------------------------------------------------------------
// Wedge vs. 3D geometry: implicit-H stereocenter (3 heavy neighbors + 1 H)
// ---------------------------------------------------------------------------

#[test]
fn implicit_h_stereocenter_wedge_agrees_with_3d_geometry_l_alanine() {
    // L-alanine's alpha carbon is a real 3-heavy-neighbor + implicit-H
    // stereocenter (N, CH3, COOH heavy substituents + 1 implicit H) --
    // RDKit's own real wedge (C(alpha)-CH3, stereo=1) and its own real 3D
    // embedding of the SAME molecule must agree.
    let report = read_mol_with_diagnostics(L_ALANINE_3D).expect("parse");
    let alpha_c = AtomIdx(1); // 0-indexed: N=0, C(alpha)=1, CH3=2, C(=O)=3, O=4, O=5
    assert_ne!(
        report.mol.atom(alpha_c).chirality,
        chematic_core::Chirality::None
    );
    assert_eq!(
        report.mol.stereo_neighbor_order(alpha_c).map(|o| o.len()),
        Some(4),
        "3 heavy neighbors + 1 implicit-H sentinel"
    );
    assert!(
        !report
            .stereo3d_diagnostics
            .iter()
            .any(|d| matches!(d, Stereo3DDiagnostic::WedgeVs3DParityConflict { .. })),
        "{:?}",
        report.stereo3d_diagnostics
    );
}

#[test]
fn implicit_h_stereocenter_wedge_disagrees_conflict_fires() {
    // Same real L-alanine geometry, wedge flipped Up(1) -> Down(6) on the
    // C(alpha)-CH3 bond line only.
    let hacked = L_ALANINE_3D.replacen("  2  3  1  1\n", "  2  3  1  6\n", 1);
    assert_ne!(hacked, L_ALANINE_3D);
    let report = read_mol_with_diagnostics(&hacked).expect("parse");
    let conflicts: Vec<_> = report
        .stereo3d_diagnostics
        .iter()
        .filter(|d| matches!(d, Stereo3DDiagnostic::WedgeVs3DParityConflict { .. }))
        .collect();
    assert_eq!(conflicts.len(), 1, "{:?}", report.stereo3d_diagnostics);
}

#[test]
fn implicit_h_stereocenter_d_alanine_also_self_consistent() {
    // The other enantiomer, independently embedded by RDKit -- its own
    // (different) wedge and its own (different, mirror-image) geometry must
    // still agree with each other.
    let report = read_mol_with_diagnostics(D_ALANINE_3D).expect("parse");
    assert!(
        !report
            .stereo3d_diagnostics
            .iter()
            .any(|d| matches!(d, Stereo3DDiagnostic::WedgeVs3DParityConflict { .. }))
    );
    // And the two enantiomers' wedge-perceived local parity must actually
    // differ (sanity check that this isn't vacuously true because neither
    // wedge was perceived at all).
    let l_report = read_mol_with_diagnostics(L_ALANINE_3D).expect("parse");
    assert_ne!(
        l_report.mol.atom(AtomIdx(1)).chirality,
        report.mol.atom(AtomIdx(1)).chirality
    );
}

// ---------------------------------------------------------------------------
// E/Z fixtures: parse cleanly, real conformer, no false tetrahedral
// diagnostics (E/Z-vs-3D verification itself is deferred -- see PR body).
// ---------------------------------------------------------------------------

#[test]
fn ez_fixtures_parse_with_real_conformer_and_no_tetrahedral_false_positive() {
    for block in [CHLORO_BROMO_ETHENE_E_3D, CHLORO_BROMO_ETHENE_Z_3D] {
        let report = read_mol_with_diagnostics(block).expect("parse");
        assert!(report.conformer.is_some());
        assert!(
            !report
                .stereo3d_diagnostics
                .iter()
                .any(|d| matches!(d, Stereo3DDiagnostic::WedgeVs3DParityConflict { .. }))
        );
    }
}

// ---------------------------------------------------------------------------
// Enhanced stereo group + 3D conformer coexist
// ---------------------------------------------------------------------------

#[test]
fn enhanced_stereo_group_coexists_with_3d_conformer() {
    let report = read_mol_v3000_with_diagnostics(ENHANCED_STEREO_3D_V3K).expect("parse");
    assert!(report.conformer.is_some());
    assert_eq!(report.geometry_rank, GeometryRank::ThreeD);
    let groups = report.mol.stereo_groups();
    assert_eq!(groups.len(), 1);
    assert_eq!(groups[0].kind, chematic_core::StereoGroupKind::Absolute);
    assert_eq!(groups[0].atom_indices, vec![AtomIdx(0)]);
}

// ---------------------------------------------------------------------------
// Multiple conformers: repeated SDF records for the same molecule
// ---------------------------------------------------------------------------

#[test]
fn sdf_conformer_ensemble_groups_repeated_records_of_same_molecule() {
    let ensembles = read_sdf_conformer_ensembles(BUTANE_ENSEMBLE_SDF).expect("parse");
    assert_eq!(
        ensembles.len(),
        1,
        "all 3 records are the same molecular graph"
    );
    assert_eq!(ensembles[0].mol.atom_count(), 4);
    assert_eq!(ensembles[0].conformers.len(), 3);
    // The 3 conformers are geometrically distinct (different embeddings).
    let p0 = ensembles[0].conformers[0].points[0];
    let p1 = ensembles[0].conformers[1].points[0];
    assert!(
        p0.distance(&p1) > 0.1,
        "different conformers should have different coordinates"
    );
}

#[test]
fn sdf_conformer_ensemble_omits_records_with_no_conformer() {
    // ETHANOL_2D/ETHANOL_3D both have a blank name line (RDKit's default for
    // an unnamed molecule) -- give each a real name here. This sidesteps a
    // separate, pre-existing bug in the shared SDF-block-splitting helpers
    // (`read_sdf_with_diagnostics` et al.'s leading-blank-line-skip strips a
    // legitimately-blank *first* line of a record, not just stray blank
    // lines between records), found incidentally while building this
    // fixture and reported as a known limitation in the PR body rather than
    // fixed here (unrelated to 3D coordinates; affects any blank-named 2D
    // SDF file too).
    let ethanol_2d_named = ETHANOL_2D.replacen('\n', "ethanol2d\n", 1);
    let ethanol_3d_named = ETHANOL_3D.replacen('\n', "ethanol3d\n", 1);
    let sdf = format!("{ethanol_2d_named}$$$$\n{ethanol_3d_named}$$$$\n");
    let ensembles = read_sdf_conformer_ensembles(&sdf).expect("parse");
    // ethanol_2d has no conformer (flat); ethanol_3d does -- exactly one
    // ensemble, with exactly one conformer (from the 3D record only).
    assert_eq!(ensembles.len(), 1);
    assert_eq!(ensembles[0].conformers.len(), 1);
}

// ---------------------------------------------------------------------------
// Typed errors: NaN/Inf/unparseable z, never a silent default or a panic
// ---------------------------------------------------------------------------

#[test]
fn v2000_nan_z_is_a_typed_error_not_silent_default() {
    let result = read_mol_with_diagnostics(NAN_Z_MOL);
    assert!(matches!(
        result,
        Err(chematic_mol::MolParseError::InvalidAtomLine { .. })
    ));
}

#[test]
fn v2000_infinite_z_from_overflow_is_a_typed_error() {
    // "1e400" is valid float *syntax* that overflows to +inf -- not a parse
    // failure by Rust's own rules, but not a finite coordinate either.
    let result = read_mol_with_diagnostics(INF_Z_MOL);
    assert!(matches!(
        result,
        Err(chematic_mol::MolParseError::InvalidAtomLine { .. })
    ));
}

#[test]
fn v2000_garbled_z_text_is_a_typed_error() {
    let result = read_mol_with_diagnostics(GARBLED_Z_MOL);
    assert!(matches!(
        result,
        Err(chematic_mol::MolParseError::InvalidAtomLine { .. })
    ));
}

#[test]
fn v2000_blank_z_field_still_defaults_to_zero_not_an_error() {
    // The z field's 10 columns are present but blank (all spaces) -- a
    // genuinely *missing* value, not corruption, same leniency as x/y.
    let report = read_mol_with_diagnostics(BLANK_Z_MOL).expect("blank z must not error");
    assert_eq!(report.geometry_rank, GeometryRank::FlatZero);
}

#[test]
fn v3000_nan_z_is_a_typed_error() {
    let bad = "\
bad
  chematic

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 NaN 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 END CTAB
M  END
";
    let result = read_mol_v3000_with_diagnostics(bad);
    assert!(matches!(
        result,
        Err(chematic_mol::MolParseError::InvalidAtomLine { .. })
    ));
}

#[test]
fn v3000_garbled_z_is_a_typed_error() {
    let bad = "\
bad
  chematic

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 garbage 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 END BOND
M  V30 END CTAB
M  END
";
    let result = read_mol_v3000_with_diagnostics(bad);
    assert!(matches!(
        result,
        Err(chematic_mol::MolParseError::InvalidAtomLine { .. })
    ));
}

// ---------------------------------------------------------------------------
// 3D writers: round-trip, never manufacture a fresh conflict on own output
// ---------------------------------------------------------------------------

#[test]
fn v2000_write_mol_with_conformer_round_trips_as_threed() {
    let original = read_mol_with_diagnostics(L_ALANINE_3D).expect("parse");
    let conformer = original.conformer.clone().expect("has conformer");
    let block = write_mol_with_conformer(&original.mol, &original.metadata, &conformer);
    assert!(
        block.contains("          3D"),
        "writer must stamp the 3D dimensional code"
    );

    let reparsed = read_mol_with_diagnostics(&block).expect("re-parse own output");
    assert_eq!(reparsed.coordinate_dimension, CoordinateDimension::ThreeD);
    assert_eq!(reparsed.geometry_rank, GeometryRank::ThreeD);
    assert!(
        reparsed.stereo3d_diagnostics.is_empty(),
        "writer must never manufacture a fresh conflict on its own round-trip output: {:?}",
        reparsed.stereo3d_diagnostics
    );
    // Coordinates survive round-trip within the writer's {:.4} precision.
    for (p_before, p_after) in conformer
        .points
        .iter()
        .zip(reparsed.conformer.unwrap().points.iter())
    {
        assert!((p_before.x - p_after.x).abs() < 1e-3);
        assert!((p_before.y - p_after.y).abs() < 1e-3);
        assert!((p_before.z - p_after.z).abs() < 1e-3);
    }
}

#[test]
fn v3000_write_mol_v3000_with_conformer_round_trips_as_threed() {
    let original = read_mol_v3000_with_diagnostics(BROMO_3D_V3K).expect("parse");
    let conformer = original.conformer.clone().expect("has conformer");
    let block = write_mol_v3000_with_conformer(&original.mol, &original.metadata, &conformer);
    assert!(block.contains("          3D"));

    let reparsed = read_mol_v3000_with_diagnostics(&block).expect("re-parse own output");
    assert_eq!(reparsed.coordinate_dimension, CoordinateDimension::ThreeD);
    assert!(
        reparsed.stereo3d_diagnostics.is_empty(),
        "{:?}",
        reparsed.stereo3d_diagnostics
    );
}

#[test]
fn v2000_2d_writer_is_completely_unaffected_by_3d_additions() {
    // Regression guard: the existing 2D writer's output must not change at
    // all as a side effect of this PR (additive-only constraint).
    use chematic_mol::mol2000::write_mol_with_coords;
    let (mol, meta, coords) = chematic_mol::parse_mol_with_coords(ETHANOL_2D).expect("parse");
    let written = write_mol_with_coords(&mol, &meta, &coords);
    assert!(!written.contains("3D"));
    assert!(!written.contains("2D"));
}
