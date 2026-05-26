//! `chematic-smarts` — SMARTS query language and VF2 subgraph isomorphism.
//!
//! # Usage
//! ```
//! // (doc example omitted — requires chematic-smiles)
//! ```

#![forbid(unsafe_code)]

pub mod match_vf2;
pub mod parser;
pub mod query;

pub use match_vf2::find_matches;
pub use parser::{SmartsError, parse_smarts};
pub use query::{
    AtomPrimitive, AtomQuery, BondPrimitive, BondQuery, QueryAtom, QueryBond, QueryMolecule,
};

#[cfg(test)]
mod tests;
