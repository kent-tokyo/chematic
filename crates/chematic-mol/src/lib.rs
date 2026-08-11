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
pub mod cif;
pub mod cjson;
pub mod cml;
pub mod error;
pub mod gaussian;
pub mod ket;
pub mod mol2000;
pub mod mol2_tripos;
pub mod mol3000;
pub mod moljson;
pub mod mrv;
pub mod pdbqt;
pub mod record;
pub mod rxn;
pub mod sdf;
pub mod smiles_table;
pub mod tdt;
pub mod xyz;

// Convenient re-exports at crate root.
pub use error::MolParseError;
// Re-export MolParseError under the alias MolError as specified by the public API.
pub use cdxml::{
    CdxmlError, CdxmlParseOptions, parse_cdxml, parse_cdxml_all, parse_cdxml_all_with_options,
    parse_cdxml_with_options, write_cdxml,
};
pub use cif::{CifError, CifResult, UnitCell, parse_cif, write_cif};
pub use cjson::{CjsonError, parse_cjson, write_cjson};
pub use cml::{CmlError, parse_cml, write_cml};
pub use error::MolParseError as MolError;
pub use gaussian::{GaussianError, GaussianLogResult, parse_gaussian_log, parse_gjf, write_gjf};
pub use ket::{KetError, parse_ket, parse_ket_3d, write_ket, write_ket_3d};
pub use mol2_tripos::{Mol2Error, parse_mol2, write_mol2};
pub use mol2000::{
    CoordinateDimension, GeometryRank, MolMetadata, MolReadReport, Stereo3DDiagnostic, parse_mol,
    parse_mol_with_coords, parse_sdf_with_coords, read_mol_with_diagnostics,
    read_sdf_with_diagnostics, write_mol, write_mol_with_conformer, write_mol_with_coords,
    write_sdf, write_sdf_record, write_sdf_record_v3000, write_sdf_record_with_conformer,
    write_sdf_with_charges,
};
pub use mol3000::{
    parse_mol_v3000, parse_mol_v3000_with_coords, read_mol_v3000_with_diagnostics, write_mol_v3000,
    write_mol_v3000_with_conformer,
};
pub use moljson::{MolJsonError, parse_moljson, write_moljson};
pub use mrv::{
    MrvError, MrvParseLimits, MrvWriteOptions, parse_mrv, parse_mrv_with_limits, write_mrv,
};
pub use pdbqt::{PdbqtError, autodock_atom_type, parse_pdbqt, write_pdbqt};
pub use record::MoleculeRecord;
pub use rxn::{RxnParseError, parse_rxn_file, write_rxn_file};
pub use sdf::{
    ConformerEnsemble, SdfFileReader, SdfReader, SdfRecord, SdfRecordReader,
    read_sdf_conformer_ensembles,
};
pub use smiles_table::{
    Delimiter, SmilesReaderOptions, SmilesRecordReader, SmilesRecordWriter, SmilesTableError,
    SmilesWriterOptions,
};
pub use tdt::{TdtError, TdtReaderOptions, TdtRecordReader, TdtRecordWriter, TdtWriterOptions};
pub use xyz::{
    ExtxyzReader, ExtxyzWriter, XyzAtom, XyzError, XyzFrame, XyzProperty, XyzPropertyKind,
    XyzReader, XyzValue, XyzWriter, parse_extxyz, parse_extxyz_all, parse_xyz, parse_xyz_all,
    write_extxyz, write_xyz,
};
