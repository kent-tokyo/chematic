//! Port of RDKit 2026.09.1's distance-geometry embedding
//! (`DGeomHelpers::EmbedMolecule`, ETKDGv3). The compatibility lane targets
//! the same fixed-seed behavior while treating platform-dependent final-bit
//! coordinate drift as a versioned artifact-gate concern.

// The port keeps RDKit's index loops, mirrored if/else chains and negated
// (NaN-propagating) comparisons so it reads side by side with the C++.
#![allow(
    clippy::needless_range_loop,
    clippy::if_same_then_else,
    clippy::neg_cmp_op_on_partial_ord,
    clippy::too_many_arguments
)]

mod alkene_path;
pub mod bounds;
pub use alkene_path::{AlkenePathStereo, rdkit_alkene_path_stereo};
mod constraints;
pub mod dg;
pub mod embed;
pub mod etk;
mod molecule;
mod rng;
mod stdsort;
mod substruct;
mod torsion_prefs;

pub use bounds::{
    BondsAndAngles, BoundsInvariant, collect_bonds_and_angles, init_bounds_mat, set_topol_bounds,
};
pub use constraints::{
    BondAngleConstraint, BondLengthConstraint, InternalCoordinateConstraintError,
    InternalCoordinateConstraints, TorsionAngleConstraint,
};
pub use dg::{BoundsMatrix, ChiralSet, triangle_smooth_bounds};
pub use embed::{EmbedArgs, EmbedError, EmbedParams, embed_initial_points, embed_points};
pub use molecule::{
    RdkitEmbedError, RdkitEmbedMmffError, RdkitEmbedOptions, RdkitMmffConformer,
    RdkitMmffConformerDiagnostics, etkdg_embed_molecule_corrected,
    etkdg_embed_multiple_confs_corrected, etkdg_embed_multiple_confs_mmff_best_gradient,
    etkdg_embed_multiple_confs_mmff_corrected, rdkit_bounds_matrix,
    rdkit_bounds_matrix_with_constraints, rdkit_embed_initial_molecule, rdkit_embed_molecule,
    rdkit_embed_molecule_with_constraints, rdkit_embed_multiple_confs,
    rdkit_embed_multiple_confs_mmff, rdkit_embed_multiple_confs_mmff_adaptive,
    rdkit_embed_multiple_confs_mmff_diagnostics, rdkit_embed_multiple_confs_mmff_gradient_rescue,
    rdkit_embed_multiple_confs_mmff_rescue,
};
