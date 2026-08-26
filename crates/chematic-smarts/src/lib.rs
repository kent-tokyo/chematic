//! `chematic-smarts` — SMARTS query language and VF2 subgraph isomorphism.
//!
//! # Usage
//! ```
//! // (doc example omitted — requires chematic-smiles)
//! ```

#![forbid(unsafe_code)]

pub mod cache;
pub(crate) mod clock;
pub mod cx;
pub mod match_vf2;
pub mod mcs;
pub mod parser;
pub mod query;
pub mod rdkit_parity_match;
pub mod rdkit_ring_model;

pub use cache::{SmartsCache, named_pattern};
pub use cx::{CxQueryAtomProp, CxSmarts, parse_cxsmarts};
pub use match_vf2::{
    MatchConfig, MatchOutcome, find_matches, find_matches_with_config, find_matches_with_rings,
    find_matches_with_rings_and_config, find_matches_with_rings_and_config_checked,
    has_match_bounded,
};
pub use mcs::{
    AtomCompare, BondCompare, McsConfig, McsOutcome, find_mcs, find_mcs_with_config,
    find_mcs_with_config_checked,
};
pub use parser::{SmartsError, parse_smarts};
pub use query::{
    AtomPrimitive, AtomQuery, BondPrimitive, BondQuery, QueryAtom, QueryBond, QueryMolecule,
};
pub use rdkit_parity_match::{
    RdkitParityConfig, find_matches_rdkit_parity, has_match_rdkit_parity_bounded,
};
pub use rdkit_ring_model::{
    RdkitParityError, RdkitParityRingModel, RdkitRingModelBudget, build_rdkit_parity_ring_model,
};

#[cfg(test)]
mod tests;
