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

pub mod ccd;
pub use ccd::{
    CcdAtom, CcdBond, CcdComponent, CcdError, CcdParseLimits, parse_ccd_component,
    parse_ccd_component_with_limits,
};

pub mod cdjson;
pub mod cdxml;
pub mod cdxml_document;
pub mod chemdraw;
pub use chemdraw::{
    ChemDrawBlock, ChemDrawError, ChemDrawFormat, ChemDrawLimits, ChemDrawMolecule,
    ChemDrawOmission, ChemDrawOptions, parse_chemdraw_block, parse_chemdraw_block_with_limits,
};
pub mod cif;
#[cfg(feature = "crystal")]
mod cif_symmetry;
pub mod cjson;
pub mod cjson_record;
pub mod cml;
pub mod cml_named;
pub mod conversion;
pub mod cube;
pub mod dmol3;
pub mod error;
pub mod exyz;
pub mod fhiaims;
pub mod gamess_zmat;
pub mod gaussian;
pub mod gro;
pub mod gro_lammps;
pub mod gzmat;
pub use gamess_zmat::{GamessZmatWriteOptions, write_gamess_zmat, xyz_to_gamess_zmat};
pub mod hin;
pub mod ket;
pub mod lammps_data;
pub mod lammps_dump;
pub mod mmcif;
pub mod mmcif_pdb;
pub mod mol2000;
pub mod mol2_kekule;
pub mod mol2_rmsd;
pub mod mol2_source;
pub mod mol2_stream;
pub mod mol2_tripos;
pub mod mol3000;
pub mod moljson;
pub mod mopac;
pub mod mrv;
pub mod nucleic_acid;
pub mod opendx;
pub mod orca;
pub mod pdb_conect;
pub mod pdb_selection;
pub mod pdbqt;
pub mod pqr;
pub mod pqr_pdbqt;
pub mod protein;
pub mod qcschema;
pub mod rdfile;
pub mod reaction_cdxml;
pub mod record;
pub mod rxn;
pub mod sdf;
pub mod semantic;
pub mod smiles_table;
pub mod stereo_depiction;
pub mod tdt;
pub mod turbomole;
pub mod volumetric;
pub mod xyz;
pub mod xyz_output;

// Convenient re-exports at crate root.
pub use cml_named::{CmlNamedError, write_cml_named_checked};
pub use conversion::{
    ConversionAtomMapping, ConversionBondMapping, ConversionError, ConversionLimits,
    ConversionObservation, ConversionPolicy, ConversionReport, PreservationStatus,
    convert_with_report, convert_with_report_json, convert_with_report_with_limits,
};
pub use error::MolParseError;
pub use reaction_cdxml::{ReactionCdxmlError, write_reaction_cdxml_checked};
// Re-export MolParseError under the alias MolError as specified by the public API.
pub use cdjson::{
    CdjsonDocument, CdjsonError, CdjsonMolecule, CdjsonParseLimits, parse_cdjson,
    parse_cdjson_with_limits, write_cdjson_graph,
};
pub use cdxml::{
    CdxmlError, CdxmlParseLimits, CdxmlParseOptions, parse_cdxml, parse_cdxml_all,
    parse_cdxml_all_with_limits, parse_cdxml_all_with_options,
    parse_cdxml_all_with_options_and_limits, parse_cdxml_with_limits, parse_cdxml_with_options,
    parse_cdxml_with_options_and_limits, write_cdxml,
};
pub use cdxml_document::{
    CdxmlDiagnostic, CdxmlDocument, CdxmlEdit, CdxmlObject, CdxmlObjectKind, CdxmlPage,
    CdxmlTextStyle, CdxmlTransform, CdxmlValue,
};
pub use cif::{
    CifBondWriteError, CifError, CifParseLimits, CifResult, UnitCell, parse_cif,
    parse_cif_with_limits, write_cif, write_cif_checked, write_cif_with_bonds_checked,
};
#[cfg(feature = "crystal")]
pub use cif::{
    CifPeriodicError, CifPeriodicParseOptions, CifPeriodicResult, CifSymmetryError,
    CifSymmetryStatus, parse_cif_periodic_structure, parse_cif_periodic_structure_with_options,
    write_cif_periodic_structure,
};
pub use cjson::{CjsonError, CjsonParseLimits, parse_cjson, parse_cjson_with_limits, write_cjson};
pub use cjson_record::{
    CjsonRecord, parse_cjson_record, parse_cjson_record_with_limits, write_cjson_record,
};
pub use cml::{
    CmlError, CmlParseLimits, parse_cml, parse_cml_strict, parse_cml_strict_with_limits,
    parse_cml_with_limits, write_cml,
};
pub use cube::{
    CubeError, CubeFileReader, CubeParseLimits, parse_cube, parse_cube_with_limits, write_cube,
};
pub use dmol3::{
    DMOL3_BOHR_ANGSTROM, Dmol3Atom, Dmol3Error, Dmol3Geometry, Dmol3ParseLimits, parse_dmol3,
    parse_dmol3_with_limits, write_dmol3,
};
pub use error::MolParseError as MolError;
pub use exyz::{
    ExyzAtom, ExyzCell, ExyzError, ExyzGeometry, ExyzParseLimits, parse_exyz,
    parse_exyz_with_limits, write_exyz,
};
pub use fhiaims::{
    FhiaimsAtom, FhiaimsCoordinates, FhiaimsError, FhiaimsGeometry, FhiaimsLatticeVector,
    FhiaimsParseLimits, parse_fhiaims, parse_fhiaims_with_limits, write_fhiaims,
};
pub use gaussian::{
    GaussianError, GaussianLogResult, GaussianParseLimits, parse_gaussian_log,
    parse_gaussian_log_with_limits, parse_gjf, parse_gjf_with_limits, write_gjf,
};
pub use gzmat::{
    GzmatWriteOptions, ZMatrixError, ZMatrixRow, cartesian_to_zmatrix, validate_zmatrix,
    write_gzmat, zmatrix_to_cartesian,
};
pub use hin::{
    HinAtom, HinBond, HinError, HinParseLimits, HinRecord, parse_hin_record,
    parse_hin_record_with_limits, write_hin_record,
};
pub use ket::{
    KetError, KetParseLimits, parse_ket, parse_ket_3d, parse_ket_3d_with_limits,
    parse_ket_with_limits, write_ket, write_ket_3d,
};
pub use lammps_data::{
    LammpsAtom, LammpsAtomStyle, LammpsBond, LammpsBox, LammpsData, LammpsDataError,
    LammpsDataParseLimits, LammpsMass, LammpsVelocity, parse_lammps_data,
    parse_lammps_data_with_limits, write_lammps_data,
};
pub use lammps_dump::{
    LammpsDumpError, LammpsDumpFrame, LammpsDumpParseLimits, LammpsDumpReader, box_bounds_to_true,
    parse_lammps_dump_frame, parse_lammps_dump_frame_with_limits, true_to_box_bounds,
    write_lammps_dump_frame, write_lammps_trajectory,
};
pub use mmcif::{
    MmcifAtomRecord, MmcifError, MmcifParseLimits, MmcifResult, parse_mmcif,
    parse_mmcif_with_limits, write_mmcif,
};
pub use mol2_tripos::{
    Mol2AtomRecord, Mol2BondRecord, Mol2Error, Mol2OpaqueSection, Mol2ParseLimits, Mol2Record,
    Mol2SectionOrder, Mol2UnityAtomAttributes, parse_mol2, parse_mol2_record,
    parse_mol2_record_with_limits, parse_mol2_with_limits, write_mol2, write_mol2_record,
    write_mol2_record_with_section_order, write_mol2_with_section_order,
};
pub use mol2000::{
    CoordinateDimension, GeometryRank, MolFormat, MolMetadata, MolReadReport, MolStereoLoss,
    MolStereoWriteError, SquarePlanarPerceptionDiagnostic, SquarePlanarRejectionReason,
    Stereo3DDiagnostic, UnsupportedStereoReason, parse_mol, parse_mol_with_coords,
    parse_sdf_with_coords, read_mol_with_diagnostics, read_sdf_with_diagnostics,
    validate_square_planar_for_write, write_laid_out_mol, write_laid_out_sdf_record_into,
    write_mol, write_mol_with_conformer, write_mol_with_conformer_checked, write_mol_with_coords,
    write_mol_with_stereo_report, write_sdf, write_sdf_record, write_sdf_record_v3000,
    write_sdf_record_with_conformer, write_sdf_record_with_conformer_checked,
    write_sdf_with_charges,
};
pub use mol3000::{
    V3000SGroup, V3000SGroupKind, parse_mol_v3000, parse_mol_v3000_with_coords,
    parse_v3000_sgroup_line, read_mol_v3000_with_diagnostics, write_mol_v3000,
    write_mol_v3000_with_conformer, write_mol_v3000_with_conformer_checked,
    write_mol_v3000_with_stereo_report,
};
pub use moljson::{
    MolJsonError, MolJsonParseLimits, parse_moljson, parse_moljson_with_limits, write_moljson,
};
pub use mopac::{
    MopacArchive, MopacArchiveAtom, MopacAtom, MopacCartesianGeometry, MopacCartesianRow,
    MopacError, MopacOptimizationFlag, MopacParseLimits, parse_mopac_archive,
    parse_mopac_archive_with_limits, parse_mopac_cartesian, parse_mopac_cartesian_with_limits,
    write_mopac_cartesian,
};
pub use mrv::{
    MrvError, MrvParseLimits, MrvWriteOptions, parse_mrv, parse_mrv_with_limits, write_mrv,
};
pub use nucleic_acid::{
    NUCLEIC_ACID_SCHEMA, NUCLEIC_ACID_VALIDATION_SCHEMA, NucleicAcidAnnotationTarget,
    NucleicAcidBase, NucleicAcidCommand, NucleicAcidDocument, NucleicAcidError, NucleicAcidLimits,
    NucleicAcidLinkage, NucleicAcidResidue, NucleicAcidStrand, NucleicAcidStrandKind,
    NucleicAcidSugar, apply_nucleic_acid_json_command, validate_nucleic_acid_json,
};
pub use opendx::{
    OpenDxError, OpenDxParseLimits, parse_opendx, parse_opendx_with_limits, write_opendx,
    write_opendx_lossy,
};
pub use orca::{
    GeometryFrame, OrcaAtom, OrcaBlock, OrcaCoords, OrcaInput, OrcaInputError,
    OrcaInputParseLimits, OrcaOptConvergence, OrcaOutput, OrcaOutputError, OrcaOutputParseLimits,
    OrcaTermination, parse_orca_input, parse_orca_input_with_limits, parse_orca_output,
    parse_orca_output_with_limits, write_orca_input,
};
pub use pdb_selection::{
    PdbSelectedAtom, PdbSelectionError, PdbSelectionLimits, PdbSelectionLoss, PdbSelectionReport,
    select_pdb_atoms, select_pdb_atoms_with_limits,
};
pub use pdbqt::{
    PdbqtError, PdbqtParseLimits, PdbqtWriteError, autodock_atom_type, parse_pdbqt,
    parse_pdbqt_with_limits, write_pdbqt, write_pdbqt_checked,
};
pub use pqr::{
    PqrAtomRecord, PqrError, PqrParseLimits, PqrResult, infer_element, parse_pqr,
    parse_pqr_with_limits, write_pqr,
};
pub use protein::{
    ProteinAtom, ProteinChain, ProteinProjection, ProteinProjectionError, ProteinResidue,
    project_protein_from_mmcif, project_protein_from_mmcif_with_limits, project_protein_from_pdb,
    project_protein_from_pdb_with_limits,
};
pub use qcschema::{
    AtomicInput, AtomicResult, Basis, ChematicMoleculeView, ComputeError, Connectivity, Driver,
    JsonObject, Provenance, QcConvertError, QcModel, QcMolecule, QcSchemaError,
    QcSchemaParseLimits, ReturnResult, chematic_to_qc_molecule, parse_atomic_input,
    parse_atomic_input_with_limits, parse_atomic_result, parse_atomic_result_with_limits,
    parse_qcschema_molecule, parse_qcschema_molecule_with_limits, qc_molecule_to_chematic,
    write_atomic_input, write_atomic_result, write_qcschema_molecule,
};
pub use record::MoleculeRecord;
pub use rxn::{
    RxnDocumentError, RxnFileParseLimits, RxnParseError, parse_rxn_document, parse_rxn_file,
    parse_rxn_file_with_limits, write_rxn_document, write_rxn_file,
};
pub use sdf::{
    ConformerEnsemble, SdfBatch, SdfBatchProgress, SdfBatchReader, SdfCompression, SdfDataset,
    SdfFileReader, SdfIndexedBatch, SdfInput, SdfParseLimits, SdfReader, SdfRecord,
    SdfRecordMetadata, SdfRecordReader, SdfSelectionError, SdfSelectionLimits,
    parse_sdf_with_limits, read_sdf_conformer_ensembles,
};
pub use semantic::{
    AtomRef, ExpandedSemantic, PolymerEndGroup, PolymerLinkage, PolymerRepeatUnit,
    RGroupDefinition, SemanticCommand, SemanticError, SemanticExpansionLimits, SemanticId,
    SemanticModel, SemanticSGroup,
};
pub use smiles_table::{
    Delimiter, SmilesReaderOptions, SmilesRecordReader, SmilesRecordWriter, SmilesTableError,
    SmilesWriterOptions,
};
pub use tdt::{TdtError, TdtReaderOptions, TdtRecordReader, TdtRecordWriter, TdtWriterOptions};
pub use volumetric::{GridAtom, GridError, GridUnits, VolumetricGrid};
pub use xyz::{
    ExtxyzBatch, ExtxyzBatchProgress, ExtxyzBatchReader, ExtxyzFileReader, ExtxyzReader,
    ExtxyzWriter, XyzAtom, XyzBatch, XyzBatchProgress, XyzBatchReader, XyzError, XyzFileReader,
    XyzFrame, XyzParseLimits, XyzProperty, XyzPropertyKind, XyzReader, XyzValue, XyzWriter,
    parse_extxyz, parse_extxyz_all, parse_extxyz_with_limits, parse_xyz, parse_xyz_all,
    parse_xyz_all_with_limits, parse_xyz_with_limits, write_extxyz, write_xyz,
};

pub use xyz_output::{XyzNotation, XyzWriteError, XyzWriteOptions, write_xyz_checked};

#[cfg(test)]
mod semantic_contract_tests;
#[cfg(test)]
mod xyz_contract_tests;

#[cfg(test)]
mod nucleic_acid_boundary_tests;
#[cfg(test)]
mod orca_boundary_tests;

#[cfg(test)]
mod xml_boundary_tests;

#[cfg(test)]
mod lammps_boundary_contract_tests;
#[cfg(test)]
mod mol_boundary_contract_tests;
#[cfg(test)]
mod stream_boundary_contract_tests;

#[cfg(test)]
mod docking_boundary_contract_tests;

#[cfg(test)]
mod cdxml_boundary_contract_tests;

#[cfg(test)]
mod format_limit_contract_tests;

#[cfg(test)]
mod format_bridge_boundary_contract_tests;

pub use mol2_kekule::{Mol2KekuleError, kekulize_mol2_record};
pub use mol2_rmsd::{Mol2ReferenceRmsd, Mol2RmsdError, mol2_reference_rmsds};
pub use mol2_source::{
    Mol2SourceAtomRecord, Mol2SourceBondRecord, Mol2SourceRecord, parse_mol2_source_record,
    parse_mol2_source_record_with_limits, write_mol2_source_record,
};
pub use mol2_stream::{Mol2SourceReader, Mol2StreamError, Mol2StreamLimits};

pub use pqr_pdbqt::{PqrPdbqtError, PqrPdbqtReport, pqr_to_pdbqt};

pub use mmcif::residues::{
    MmcifAtomIdentity, MmcifResidue, MmcifResidueDocument, MmcifResidueError, MmcifResidueIdentity,
    parse_mmcif_residues, parse_mmcif_residues_with_limits,
};

pub use mmcif_pdb::{MmcifPdbError, MmcifPdbReport, PdbNamespace, write_mmcif_pdb_checked};

pub use gro::{
    GroAtomRecord, GroError, GroFrame, GroParseLimits, parse_gro, parse_gro_with_limits,
};

pub use turbomole::{
    TURBOMOLE_BOHR_ANGSTROM, TurbomoleAtom, TurbomoleError, TurbomoleGeometry,
    TurbomoleParseLimits, parse_turbomole, parse_turbomole_with_limits, write_turbomole,
};

pub use rdfile::{
    RdfileData, RdfileDocument, RdfileError, RdfileKind, RdfileParseLimits, RdfileRecord,
    RdfileWriteOptions, parse_rdfile, parse_rdfile_with_limits, write_rdfile,
    write_rdfile_with_options,
};

pub use pdb_conect::{
    PdbConectError, PdbConectFrame, PdbConectLimits, PdbConectOptions, parse_pdb_conect,
    parse_pdb_conect_with_options, parse_pdb_conect_with_options_and_limits,
};
