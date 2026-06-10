//! `chematic-fp` — molecular fingerprints for chematic.
//!
//! Provides:
//! - **ECFP** (Extended Connectivity Fingerprints): Morgan-algorithm-based circular
//!   fingerprints using FNV-1a hashing for reproducibility.
//! - **MACCS 166-bit structural keys**: SMARTS-based structural key fingerprints.
//! - **Topological path fingerprints**: DFS path enumeration up to a configurable length.
//! - **AtomPair fingerprints**: atom-pair encoding with topological distances.
//! - **Topological Torsion fingerprints**: four-atom path encoding.

#![forbid(unsafe_code)]

pub mod atom_pair;
pub mod bitvec;
pub mod ecfp;
pub mod fcfp;
pub mod maccs;
pub mod path;
pub mod pharmacophore_fp;
pub mod search;
pub mod topo_path;

pub use atom_pair::{atom_pair_fp, torsion_fp};
pub use bitvec::BitVec2048;
pub use ecfp::{EcfpConfig, ecfp, ecfp4, ecfp6, morgan_fp_counts, tanimoto_ecfp4};
pub use fcfp::{fcfp, fcfp4, fcfp6, tanimoto_fcfp4};
pub use maccs::maccs;
pub use path::{RdkitPathConfig, rdkit_path_fp, rdkit_path_fp_with_config, tanimoto_rdkit_path};
pub use pharmacophore_fp::{
    pharmacophore_feature_counts, pharmacophore_fp_2d,
    tanimoto_pharmacophore_2d,
};
pub use search::{FpType, nearest_neighbors, nearest_neighbors_from_fp};
pub use topo_path::{TopoPathConfig, tanimoto_topo_path, topo_path};
