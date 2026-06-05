#![forbid(unsafe_code)]
//! `chematic-mol` — SDF/MOL V2000 file format reader and writer for chematic.
//!
//! # Overview
//! - [`mol2000`]: parse and write individual MOL V2000 (Ctab) blocks.
//! - [`sdf`]: iterate over multi-molecule SDF files.
//!
//! # Quick start
//! ```rust
//! use chematic_mol::{parse_mol, write_mol, SdfReader};
//! use chematic_mol::mol2000::MolMetadata;
//!
//! let mol_str = "ethanol\n  prog\n\n  3  2  0  0  0  0  0  0  0  0  0 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\nM  END\n";
//! let (mol, meta) = parse_mol(mol_str).unwrap();
//! assert_eq!(mol.atom_count(), 3);
//! let written = write_mol(&mol, &meta);
//! ```

pub mod cdxml;
pub mod cml;
pub mod error;
pub mod mol2000;
pub mod mol3000;
pub mod sdf;

// Convenient re-exports at crate root.
pub use error::MolParseError;
// Re-export MolParseError under the alias MolError as specified by the public API.
pub use error::MolParseError as MolError;
pub use cdxml::{CdxmlError, parse_cdxml, parse_cdxml_all};
pub use cml::{CmlError, parse_cml, write_cml};
pub use mol2000::{MolMetadata, parse_mol, parse_mol_with_coords, parse_sdf_with_coords,
                  write_mol, write_mol_with_coords, write_sdf};
pub use mol3000::{parse_mol_v3000, write_mol_v3000};
pub use sdf::{SdfReader, SdfRecord, SdfRecordReader};
