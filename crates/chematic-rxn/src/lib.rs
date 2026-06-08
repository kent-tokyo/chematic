#![forbid(unsafe_code)]
//! `chematic-rxn` — reaction SMILES parser and writer for chematic.
//!
//! Provides:
//! - [`Reaction`]: reactants, agents, products as `Vec<Molecule>`.
//! - [`parse_reaction`]: parse a reaction SMILES string `"R>>P"` or `"R>A>P"`.
//! - [`write_reaction`]: serialize back to reaction SMILES.
//! - [`RxnError`]: parse error type.
//! - [`run_reactants`]: apply a SMIRKS template to reactant molecules.
//! - [`TransformError`]: error type for SMIRKS transformation.

pub mod balance;
pub mod green;
pub mod reaction;
pub mod transform;

pub use balance::{BalanceResult, balance_check};
pub use green::{atom_economy, e_factor, pmi_rxn, reaction_mass_efficiency};
pub use reaction::{
    Reaction, ReactionCenter, RxnError, find_reaction_center, parse_reaction, write_reaction,
};
pub use transform::{TransformError, run_reactants};
