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

pub mod reaction;
pub mod transform;

pub use reaction::{Reaction, RxnError, parse_reaction, write_reaction};
pub use transform::{TransformError, run_reactants};
