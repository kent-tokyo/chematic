//! Port of RDKit 2026.03's distance-geometry embedding
//! (`DGeomHelpers::EmbedMolecule`, ETKDGv3) aiming at bit-identical
//! coordinates for a given `randomSeed`.

mod bfgs;
pub mod dg;
pub mod embed;
mod rng;

pub use dg::{BoundsMatrix, ChiralSet, triangle_smooth_bounds};
pub use embed::{EmbedArgs, EmbedError, EmbedParams, embed_points};
