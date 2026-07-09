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
pub mod avalon;
pub mod bitvec;
pub mod bulk;
pub mod ecfp;
pub mod erg;
pub mod fcfp;
pub mod hdf;
pub mod layered;
pub mod lsh;
pub mod maccs;
pub mod map4;
pub mod mhfp;
pub mod path;
pub mod pattern;
pub mod pharmacophore_fp;
pub mod reaction_fp;
pub mod search;
pub mod topo_path;

pub use atom_pair::{atom_pair_fp, torsion_fp};
pub use avalon::{AvalonConfig, avalon_fp, avalon_fp_with_config, tanimoto_avalon};
pub use bitvec::{BitVec2048, BitVecN};
pub use bulk::{tanimoto_matrix, tanimoto_slice, top_k_similar};
pub use ecfp::{
    EcfpConfig, ecfp, ecfp_with_bitinfo, ecfp4, ecfp6, morgan_fp_counts, tanimoto_ecfp4,
};
pub use erg::{
    ERG_VEC_LEN, ErgAtomType, ErgBondType, ErgConfig, ErgFingerprint, cosine_erg_vec, erg,
    erg_extended, erg_vec, erg_with_config, tanimoto_erg, tanimoto_erg_vec,
};
pub use fcfp::{fcfp, fcfp_with_bitinfo, fcfp4, fcfp6, tanimoto_fcfp4};
pub use hdf::{HdfConfig, HdfFp, cosine_hdf, hdf, hdf_default};
pub use layered::{layered_fp, layered_fp_by_layer, tanimoto_layered};
pub use lsh::MhfpLshIndex;
pub use maccs::maccs;
pub use map4::{Map4Config, map4, map4_default, tanimoto_map4};
pub use mhfp::{MhfpConfig, MhfpFingerprint, mhfp, mhfp_128, mhfp_with_config, tanimoto_mhfp};
pub use path::{RdkitPathConfig, rdkit_path_fp, rdkit_path_fp_with_config, tanimoto_rdkit_path};
pub use pattern::{pattern_fp, tanimoto_pattern};
pub use pharmacophore_fp::{
    pharmacophore_feature_counts, pharmacophore_fp_2d, tanimoto_pharmacophore_2d,
};
pub use reaction_fp::{
    ReactionFingerprint, ReactionFpConfig, reaction_fp, reaction_fp_ecfp4, reaction_fp_with_config,
    tanimoto_reaction_fp,
};
pub use search::{FpType, nearest_neighbors, nearest_neighbors_from_fp};
pub use topo_path::{TopoPathConfig, tanimoto_topo_path, topo_path};
