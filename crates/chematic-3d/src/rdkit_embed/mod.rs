//! Port of RDKit 2026.03's distance-geometry embedding
//! (`DGeomHelpers::EmbedMolecule`, ETKDGv3) aiming at bit-identical
//! coordinates for a given `randomSeed`.

// The port keeps RDKit's index loops, mirrored if/else chains and negated
// (NaN-propagating) comparisons so it reads side by side with the C++.
#![allow(
    clippy::needless_range_loop,
    clippy::if_same_then_else,
    clippy::neg_cmp_op_on_partial_ord,
    clippy::too_many_arguments
)]

pub mod bounds;
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
pub use dg::{BoundsMatrix, ChiralSet, triangle_smooth_bounds};
pub use embed::{EmbedArgs, EmbedError, EmbedParams, embed_points};
pub use molecule::{RdkitEmbedError, RdkitEmbedOptions, rdkit_embed_molecule};
