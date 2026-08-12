#![forbid(unsafe_code)]
//! `chematic-smiles` — OpenSMILES parser, writer, and canonical SMILES generator.
//!
//! # Quick start
//! ```rust
//! use chematic_smiles::{parse, write, canonical_smiles};
//!
//! let mol = parse("c1ccccc1").unwrap(); // benzene
//! assert_eq!(mol.atom_count(), 6);
//! assert_eq!(mol.bond_count(), 6);
//!
//! // Non-canonical write (DFS order).
//! let smiles = write(&mol);
//! let mol2 = parse(&smiles).unwrap();
//! assert_eq!(mol.atom_count(), mol2.atom_count());
//!
//! // Canonical SMILES (stable, unique).
//! let c1 = canonical_smiles(&mol);
//! let c2 = canonical_smiles(&parse("C1=CC=CC=C1").unwrap()); // Kekule benzene
//! // c1 and c2 differ because aromaticity differs, but both are stable.
//! assert_eq!(c1, canonical_smiles(&parse(&c1).unwrap()));
//! ```
//!
//! # Design
//! - Pure Rust: no C/C++ FFI, no unsafe.
//! - Single-pass recursive-descent parser; no separate lexer phase.
//! - WASM-compatible (no filesystem I/O, no threads).

pub mod canonical;
mod canonical_automorphism;
mod canonical_partition;
pub mod canonical_search;
pub mod cx;
pub mod error;
pub mod parser;
pub mod random_smiles;
pub mod smi_file;
pub mod writer;

pub use canonical::are_atoms_equivalent;
pub use canonical::{
    canonical_atom_order, canonical_smiles, equivalent_atom_classes, morgan_ranks,
};
pub use canonical_partition::topological_equivalence_classes;
pub use canonical_search::{
    CanonicalSearchStats, CanonicalizationError, CanonicalizationLimits,
    canonical_smiles_with_limits, reset_search_stats, search_stats_snapshot,
};
pub use cx::{CxAtomProp, CxSmiles, parse_cxsmiles, write_cxsmiles};
pub use error::SmilesError;
pub use parser::parse;
pub use random_smiles::{random_smiles, random_smiles_vect};
pub use smi_file::{parse_smi_file, write_smi_file};
pub use writer::write;
