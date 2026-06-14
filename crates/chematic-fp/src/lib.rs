//! `chematic-fp` — molecular fingerprints for chematic.
//!
//! Provides:
//! - **ECFP** (Extended Connectivity Fingerprints): Morgan-algorithm-based circular
//!   fingerprints using FNV-1a hashing for reproducibility.
//! - **MACCS 166-bit structural keys**: SMARTS-based structural key fingerprints.
//! - **Topological path fingerprints**: DFS path enumeration up to a configurable length.
//! - **AtomPair fingerprints**: atom-pair encoding with topological distances.
//! - **Topological Torsion fingerprints**: four-atom path encoding.
//! - **C-Series Reaction fingerprints**: Chemical transformation encoding (Phase 1)
//! - **MHFP/SECFP**: MinHash fingerprints for fast approximate similarity searching (Phase 1)
//! - **ERG**: Extended Reduced Graph fingerprints for functional group-based similarity (Phase 1)

#![forbid(unsafe_code)]

pub mod atom_pair;
pub mod bitvec;
pub mod bulk;
pub mod map4;
pub mod ecfp;
pub mod erg;
pub mod fcfp;
pub mod layered;
pub mod lsh;
pub mod maccs;
pub mod mhfp;
pub mod path;
pub mod pattern;
pub mod pharmacophore_fp;
pub mod reaction_fp;
pub mod search;
pub mod topo_path;

pub use atom_pair::{atom_pair_fp, torsion_fp};
pub use bulk::{tanimoto_slice, tanimoto_matrix, top_k_similar};
pub use bitvec::{BitVec2048, BitVecN};
pub use ecfp::{EcfpConfig, ecfp, ecfp4, ecfp6, morgan_fp_counts, tanimoto_ecfp4};
pub use erg::{erg, erg_with_config, erg_extended, tanimoto_erg, ErgFingerprint, ErgConfig, ErgAtomType, ErgBondType,
              erg_vec, cosine_erg_vec, tanimoto_erg_vec, ERG_VEC_LEN};
pub use fcfp::{fcfp, fcfp4, fcfp6, tanimoto_fcfp4};
pub use layered::{layered_fp, layered_fp_by_layer, tanimoto_layered};
pub use maccs::maccs;
pub use mhfp::{mhfp, mhfp_with_config, mhfp_128, tanimoto_mhfp, MhfpFingerprint, MhfpConfig};
pub use path::{RdkitPathConfig, rdkit_path_fp, rdkit_path_fp_with_config, tanimoto_rdkit_path};
pub use pattern::{pattern_fp, tanimoto_pattern};
pub use pharmacophore_fp::{
    pharmacophore_feature_counts, pharmacophore_fp_2d,
    tanimoto_pharmacophore_2d,
};
pub use reaction_fp::{
    reaction_fp, reaction_fp_with_config, reaction_fp_ecfp4, tanimoto_reaction_fp,
    ReactionFingerprint, ReactionFpConfig,
};
pub use lsh::MhfpLshIndex;
pub use map4::{Map4Config, map4, map4_default, tanimoto_map4};
pub use search::{FpType, nearest_neighbors, nearest_neighbors_from_fp};
pub use topo_path::{TopoPathConfig, tanimoto_topo_path, topo_path};
