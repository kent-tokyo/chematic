#![forbid(unsafe_code)]
//! `chematic-rxn` — reaction SMILES parser and writer for chematic.
//!
//! Provides:
//! - [`Reaction`]: reactants, agents, products as `Vec<Molecule>`.
//! - [`parse_reaction`]: parse a reaction SMILES string `"R>>P"` or `"R>A>P"`.
//! - [`write_reaction`]: serialize back to reaction SMILES.
//! - [`RxnError`]: parse error type.
//! - [`run_reactants`]: apply a SMIRKS template to reactant molecules.
//! - [`find_reaction_matches`]/[`apply_reaction_match`]: enumerate matches and
//!   apply one of them independently, for callers that need to accept/reject
//!   individual matches rather than an entire `run_reactants` call.
//! - [`TransformError`]: error type for SMIRKS transformation.
//! - [`enumerate_library`]: combinatorial library enumeration from SMIRKS + fragment sets.

pub mod balance;
pub mod enumerate;
pub mod green;
pub mod perf_counters;
pub mod query;
pub mod reaction;
pub mod retro;
pub mod transform;

pub use balance::{BalanceResult, balance_check};
pub use enumerate::{
    LibraryConfig, LibraryError, enumerate_library, enumerate_library_2way, enumerate_library_3way,
};
pub use green::{atom_economy, e_factor, pmi_rxn, reaction_mass_efficiency};
pub use perf_counters::PerfCounters;
pub use query::{
    BatchQueryResults, ReactionPatternLibrary, ReactionQuery, ReactionQueryError,
    batch_query_reactions, batch_query_with_library, has_reaction_substructure_match,
    parse_reaction_query, query_reaction,
};
pub use reaction::{
    Reaction, ReactionCenter, RxnError, find_reaction_center, parse_reaction, write_reaction,
};
pub use retro::{DEFAULT_TEMPLATES, RetroClass, RetroResult, RetroTemplate, retro_disconnect};
pub use transform::{
    ReactionMatch, TransformError, apply_reaction_match, find_reaction_matches, run_reactants,
    run_reactants_strict,
};
