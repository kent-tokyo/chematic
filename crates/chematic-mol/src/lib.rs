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
#[cfg(feature = "crystal")]
mod cif_symmetry;
pub mod cjson;
pub mod cml;
pub mod cube;
pub mod error;
pub mod gaussian;
pub mod ket;
pub mod lammps_data;
pub mod lammps_dump;
pub mod mmcif;
pub mod mol2000;
pub mod mol2_tripos;
pub mod mol3000;
pub mod moljson;
pub mod mrv;
pub mod opendx;
pub mod orca;
pub mod pdbqt;
pub mod pqr;
pub mod qcschema;
pub mod record;
pub mod rxn;
pub mod sdf;
pub mod smiles_table;
pub mod tdt;
pub mod volumetric;
pub mod xyz;

// Convenient re-exports at crate root.
pub use error::MolParseError;
// Re-export MolParseError under the alias MolError as specified by the public API.
pub use cdxml::{
    CdxmlError, CdxmlParseOptions, parse_cdxml, parse_cdxml_all, parse_cdxml_all_with_options,
    parse_cdxml_with_options, write_cdxml,
};
pub use cif::{CifError, CifResult, UnitCell, parse_cif, write_cif};
#[cfg(feature = "crystal")]
pub use cif::{
    CifPeriodicError, CifPeriodicParseOptions, CifPeriodicResult, CifSymmetryError,
    CifSymmetryStatus, parse_cif_periodic_structure, parse_cif_periodic_structure_with_options,
    write_cif_periodic_structure,
};
pub use cjson::{CjsonError, parse_cjson, write_cjson};
pub use cml::{CmlError, parse_cml, write_cml};
pub use cube::{
    CubeError, CubeFileReader, CubeParseLimits, parse_cube, parse_cube_with_limits, write_cube,
};
pub use error::MolParseError as MolError;
pub use gaussian::{GaussianError, GaussianLogResult, parse_gaussian_log, parse_gjf, write_gjf};
pub use ket::{KetError, parse_ket, parse_ket_3d, write_ket, write_ket_3d};
pub use lammps_data::{
    LammpsAtom, LammpsAtomStyle, LammpsBond, LammpsBox, LammpsData, LammpsDataError, LammpsMass,
    LammpsVelocity, parse_lammps_data, write_lammps_data,
};
pub use lammps_dump::{
    LammpsDumpError, LammpsDumpFrame, LammpsDumpReader, box_bounds_to_true,
    parse_lammps_dump_frame, true_to_box_bounds, write_lammps_dump_frame, write_lammps_trajectory,
};
pub use mmcif::{
    MmcifAtomRecord, MmcifError, MmcifParseLimits, MmcifResult, parse_mmcif,
    parse_mmcif_with_limits, write_mmcif,
};
pub use mol2_tripos::{Mol2Error, parse_mol2, write_mol2};
pub use mol2000::{
    CoordinateDimension, GeometryRank, MolFormat, MolMetadata, MolReadReport, MolStereoWriteError,
    SquarePlanarPerceptionDiagnostic, SquarePlanarRejectionReason, Stereo3DDiagnostic,
    UnsupportedStereoReason, parse_mol, parse_mol_with_coords, parse_sdf_with_coords,
    read_mol_with_diagnostics, read_sdf_with_diagnostics, validate_square_planar_for_write,
    write_mol, write_mol_with_conformer, write_mol_with_conformer_checked, write_mol_with_coords,
    write_sdf, write_sdf_record, write_sdf_record_v3000, write_sdf_record_with_conformer,
    write_sdf_record_with_conformer_checked, write_sdf_with_charges,
};
pub use mol3000::{
    parse_mol_v3000, parse_mol_v3000_with_coords, read_mol_v3000_with_diagnostics, write_mol_v3000,
    write_mol_v3000_with_conformer, write_mol_v3000_with_conformer_checked,
};
pub use moljson::{MolJsonError, parse_moljson, write_moljson};
pub use mrv::{
    MrvError, MrvParseLimits, MrvWriteOptions, parse_mrv, parse_mrv_with_limits, write_mrv,
};
pub use opendx::{
    OpenDxError, OpenDxParseLimits, parse_opendx, parse_opendx_with_limits, write_opendx,
    write_opendx_lossy,
};
pub use orca::{
    GeometryFrame, OrcaAtom, OrcaBlock, OrcaCoords, OrcaInput, OrcaInputError, OrcaOptConvergence,
    OrcaOutput, OrcaOutputError, OrcaTermination, parse_orca_input, parse_orca_output,
    write_orca_input,
};
pub use pdbqt::{PdbqtError, autodock_atom_type, parse_pdbqt, write_pdbqt};
pub use pqr::{
    PqrAtomRecord, PqrError, PqrParseLimits, PqrResult, infer_element, parse_pqr,
    parse_pqr_with_limits, write_pqr,
};
pub use qcschema::{
    AtomicInput, AtomicResult, Basis, ChematicMoleculeView, ComputeError, Connectivity, Driver,
    JsonObject, Provenance, QcConvertError, QcModel, QcMolecule, QcSchemaError, ReturnResult,
    chematic_to_qc_molecule, parse_atomic_input, parse_atomic_result, parse_qcschema_molecule,
    qc_molecule_to_chematic, write_atomic_input, write_atomic_result, write_qcschema_molecule,
};
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
pub use volumetric::{GridAtom, GridError, GridUnits, VolumetricGrid};
pub use xyz::{
    ExtxyzReader, ExtxyzWriter, XyzAtom, XyzError, XyzFrame, XyzProperty, XyzPropertyKind,
    XyzReader, XyzValue, XyzWriter, parse_extxyz, parse_extxyz_all, parse_xyz, parse_xyz_all,
    write_extxyz, write_xyz,
};
