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

pub mod batch;
pub mod canonical;
mod canonical_automorphism;
mod canonical_partition;
pub mod canonical_search;
pub mod cx;
pub mod error;
pub mod parser;
pub mod random_smiles;
pub mod rdkit;
pub mod smi_file;
pub mod writer;

pub use batch::{
    BatchCanonicalRecord, BatchCanonicalResult, BatchCanonicalStreamResult, BatchCanonicalization,
    IdentityIndexBuild, SmilesBatchCanonicalizer, SmilesBatchReader, SmilesBatchReaderError,
    SmilesBatchStream, SmilesIdentityIndex, StreamTerminalReason,
};
pub use canonical::are_atoms_equivalent;
pub use canonical::{
    canonical_atom_order, canonical_smiles, canonical_smiles_stable_key,
    canonical_smiles_with_atom_order, equivalent_atom_classes, morgan_ranks,
};
pub use canonical_partition::topological_equivalence_classes;
pub use canonical_search::{
    CanonicalSearchStats, CanonicalizationError, CanonicalizationLimits,
    canonical_smiles_with_limits, reset_search_stats, search_stats_snapshot,
};
pub use cx::{CxAtomProp, CxSmiles, attachment_point_label_number, parse_cxsmiles, write_cxsmiles};
pub use error::SmilesError;
pub use parser::{SmilesParseLimits, parse, parse_template, parse_with_limits};
pub use random_smiles::{random_smiles, random_smiles_vect};
pub use rdkit::{
    InchiOutputAtom, InchiOutputStereo0D, RdkitAlignment, RdkitHashFunction, RdkitLegacyStereo,
    RdkitMol2Molecule, RdkitMolBlock, RdkitMolBlockAtom, RdkitMolBlockBond, RdkitMolView,
    RdkitPdbMolecule, RdkitSanitizedModel, RdkitSmilesError, RdkitSmilesParams,
    RdkitTautomerEnumeration, RdkitTautomerStatus, RdkitViewAtom, RdkitViewBond, Transform3D,
    rdkit_2d_coords, rdkit_addhs_first_pattern, rdkit_align_mol, rdkit_align_points,
    rdkit_atom_stereocenter_counts, rdkit_best_rms, rdkit_calc_rms, rdkit_canonical_smiles,
    rdkit_canonical_tautomer, rdkit_chiral_centers, rdkit_crippen_contribs,
    rdkit_crippen_contribs_no_hs, rdkit_crippen_logp_mr, rdkit_cx_smarts,
    rdkit_enumerate_tautomers, rdkit_hybridizations, rdkit_hydrogen_suppressed,
    rdkit_legacy_stereo, rdkit_model_correct_view, rdkit_mol_block, rdkit_mol_block_2d,
    rdkit_mol_from_mol2_block, rdkit_mol_from_pdb_block, rdkit_mol_from_xyz_block, rdkit_mol_hash,
    rdkit_mol_view, rdkit_molecule_from_inchi_output, rdkit_murcko_scaffold, rdkit_pdb_block,
    rdkit_sanitized_model, rdkit_smarts, rdkit_smiles, rdkit_stereoisomer_count,
    rdkit_stereoisomer_smiles, rdkit_tautomer_score, register_rdkit_model_hook,
};
pub use smi_file::{
    SmiFileParseLimits, parse_smi_file, parse_smi_file_with_limits, write_smi_file,
};
pub use writer::{write, write_with_atom_order};
