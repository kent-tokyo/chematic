//! `chematic-smarts` — SMARTS query language and VF2 subgraph isomorphism.
//!
//! # Usage
//! ```
//! // (doc example omitted — requires chematic-smiles)
//! ```

#![forbid(unsafe_code)]

pub mod match_vf2;
pub mod mcs;
pub mod parser;
pub mod query;

pub use match_vf2::{MatchConfig, find_matches, find_matches_with_config};
pub use mcs::{AtomCompare, BondCompare, McsConfig, find_mcs, find_mcs_with_config};
pub use parser::{SmartsError, parse_smarts};
pub use query::{
    AtomPrimitive, AtomQuery, BondPrimitive, BondQuery, QueryAtom, QueryBond, QueryMolecule,
};

#[cfg(test)]
mod tests;
