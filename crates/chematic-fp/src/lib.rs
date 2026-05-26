//! `chematic-fp` — molecular fingerprints for chematic.
//!
//! Provides:
//! - **ECFP** (Extended Connectivity Fingerprints): Morgan-algorithm-based circular
//!   fingerprints using FNV-1a hashing for reproducibility.
//! - **MACCS 166-bit structural keys**: SMARTS-based structural key fingerprints.
//! - **Topological path fingerprints**: DFS path enumeration up to a configurable length.

#![forbid(unsafe_code)]

pub mod bitvec;
pub mod ecfp;
pub mod maccs;
pub mod topo_path;

pub use bitvec::BitVec2048;
pub use ecfp::{EcfpConfig, ecfp, ecfp4, ecfp6, tanimoto_ecfp4};
pub use maccs::maccs;
pub use topo_path::{TopoPathConfig, topo_path};
