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
//! - **FPS**: streaming read/write for the chemfp/OpenBabel FPS fingerprint interchange format
//!   ([`fps`])

#![forbid(unsafe_code)]

pub mod atom_pair;
pub mod avalon;
pub mod bitvec;
pub mod bulk;
pub mod ecfp;
mod ecfp_diagnostics;
pub mod erg;
pub mod fcfp;
pub mod fps;
pub mod hdf;
pub mod layered;
pub mod lsh;
pub mod maccs;
pub mod map4;
pub mod mhfp;
mod morgan_environment;
pub mod path;
pub mod pattern;
pub mod pharmacophore_fp;
pub mod rdkit_atom_pair;
mod rdkit_isotope_delta_table;
pub mod rdkit_layered;
mod rdkit_morgan_config;
mod rdkit_morgan_ecfp4;
mod rdkit_morgan_hash;
pub mod rdkit_pattern;
pub mod rdkit_rdk;
pub mod rdkit_torsion;
pub mod reaction_fp;
pub mod search;
pub mod topo_path;

pub use atom_pair::{atom_pair_fp, torsion_fp};
pub use avalon::{AvalonConfig, avalon_fp, avalon_fp_with_config, tanimoto_avalon};
pub use bitvec::{BitVec2048, BitVecN};
pub use bulk::{tanimoto_matrix, tanimoto_matrix_parallel, tanimoto_slice, top_k_similar};
pub use ecfp::{
    EcfpConfig, EcfpInvariantMode, atom_invariants, ecfp, ecfp_with_bitinfo,
    ecfp_with_bitinfo_and_mode, ecfp_with_bitinfo_rdkit_environment_experimental,
    ecfp_with_invariant_mode, ecfp4, ecfp4_rdkit_environment_experimental, ecfp4_rdkit_invariants,
    ecfp6, ecfp6_rdkit_environment_experimental, ecfp6_rdkit_invariants, morgan_fp_counts,
    tanimoto_ecfp4,
};
pub use rdkit_atom_pair::rdkit_atom_pair_fp;
pub use rdkit_layered::rdkit_layered_fp;
pub use rdkit_pattern::rdkit_pattern_fp;
pub use rdkit_rdk::rdkit_rdk_fp;
pub use rdkit_torsion::rdkit_torsion_fp;
/// Diagnostic-only APIs, not meant for production use — a per-`(atom,
/// radius)` trace of chematic's real Morgan expansion, for the RDKit
/// environment-parity oracle (see `scripts/ecfp_rdkit_environment_parity.py`),
/// and a raw-identifier suppressed-emission dump for sparse-count-shape
/// validation (see `scripts/ecfp_rdkit_suppression_parity.py`). Gated behind
/// the `diagnostics` feature.
#[cfg(feature = "diagnostics")]
#[doc(hidden)]
pub mod diagnostics {
    pub use crate::ecfp_diagnostics::{MorganTraceEntry, atom_ball, morgan_trace};
    pub use crate::morgan_environment::suppressed_environments_diagnostic;
    pub use crate::rdkit_morgan_hash::{RdkitMorganRawTraceEntry, rdkit_morgan_raw_trace};
}
pub use erg::{
    ERG_VEC_LEN, ErgAtomType, ErgBondType, ErgConfig, ErgFingerprint, cosine_erg_vec, erg,
    erg_extended, erg_vec, erg_with_config, tanimoto_erg, tanimoto_erg_vec,
};
pub use fcfp::{fcfp, fcfp_with_bitinfo, fcfp4, fcfp6, tanimoto_fcfp4};
pub use fps::{FpsError, FpsHeader, FpsReader, FpsRecord, FpsWriter};
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
pub use rdkit_morgan_config::{
    RdkitMorganConfig, RdkitMorganFingerprint, RdkitMorganFpSize, RdkitMorganRadius,
    rdkit_morgan_fingerprint,
};
pub use rdkit_morgan_ecfp4::{RdkitMorganEcfp4, RdkitMorganError, rdkit_morgan_ecfp4_experimental};
pub use reaction_fp::{
    ReactionFingerprint, ReactionFpConfig, reaction_fp, reaction_fp_ecfp4, reaction_fp_with_config,
    tanimoto_reaction_fp,
};
pub use search::{FpType, PreparedFingerprintIndex, nearest_neighbors, nearest_neighbors_from_fp};
pub use topo_path::{TopoPathConfig, tanimoto_topo_path, topo_path};
