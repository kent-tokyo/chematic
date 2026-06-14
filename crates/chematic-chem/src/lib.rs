//! `chematic-chem` — molecular descriptors for chematic.
//!
//! All descriptor functions take a `&Molecule` reference.
//! Values are approximate; calibrated against RDKit defaults.

#![forbid(unsafe_code)]

pub mod abbreviations;
pub mod alerts;
pub mod atropisomer;
pub mod brics;
pub mod cache;
pub mod cip;
pub mod condensed;
pub mod descriptors;
pub mod diversity;
pub mod esol;
pub mod estate;
pub mod gasteiger;
pub mod hash;
pub mod hydrogen;
pub mod ifg;
pub mod isotope_distribution;
pub mod logd;
pub mod mmp;
pub mod mmff94_bci;
pub mod named_groups;
pub mod qed;
pub mod recap;
pub mod sa_score;
pub mod scaffold;
pub mod standardize;
pub mod stereo;
pub mod tautomer;
pub mod topo_descriptors;
pub mod vsa;
pub mod workflow;
pub mod iupac_stereo;
pub mod xlogp3;

pub use cip::{CipAssignment, assign_cip};
pub use iupac_stereo::iupac_name_stereo;
pub use descriptors::{
    aromatic_ring_count, autocorr_2d, balaban_j, calc_mol_formula, egan_passes, exact_mass,
    formal_charge_sum, fsp3, ghose_passes, hba_count, hbd_count, hall_kier_alpha,
    heavy_atom_count, ipc, lipinski_passes, logp_crippen, logp_crippen_per_atom, mmff94_charges,
    molar_refractivity, molecular_weight, mqn, mr_per_atom, num_aliphatic_heterocycles,
    num_aliphatic_rings, num_amide_bonds, num_aromatic_heterocycles, num_bridgehead_atoms,
    num_bromines, num_carbons, num_chlorines, num_ester_bonds, num_fluorines, num_heteroatoms,
    num_hydrogens, num_iodines, num_nitrogens, num_oxygens, num_phosphorus, num_saturated_heterocycles,
    num_saturated_rings, num_spiro_atoms, num_stereocenters, num_sulfurs,
    num_unspecified_stereocenters, reos_passes, ring_count, rotatable_bond_count, tpsa, usrcat,
    veber_passes,
};

pub use abbreviations::{abbreviations, expand_abbreviation};
pub use alerts::{brenk_matches, brenk_passes, pains_matches, pains_passes};
pub use atropisomer::{
    AtropisomerType, assign_atropisomer_chirality, detect_atropisomers,
};
pub use brics::{BricsConfig, brics_bonds, brics_fragments, brics_fragments_with_config};
pub use cache::{DescriptorCache, DescriptorEntry};
pub use condensed::{CondensedError, parse_condensed};
pub use diversity::{butina_cluster, maxmin_picks};
pub use esol::esol_solubility;
pub use estate::{estate_indices, max_estate, min_estate, sum_estate};
pub use gasteiger::gasteiger_charges;
pub use hash::{are_identical, mol_hash};
pub use hydrogen::{add_hydrogens, remove_hydrogens};
pub use ifg::{FunctionalGroup, identify_functional_groups};
pub use isotope_distribution::isotope_distribution;
pub use logd::{logd_profile, logd_simple};
pub use mmp::{MmpPair, find_mmp};
pub use named_groups::{NamedGroup, detect_named_functional_groups};
pub use qed::qed;
pub use recap::{recap_breakable_bond_count, recap_fragment};
pub use sa_score::sa_score;
pub use scaffold::{
    generic_murcko_scaffold, murcko_scaffold, scaffold_network, scaffold_network_with_counts,
    schuffenhauer_parents, ScaffoldNetwork,
};
pub use standardize::{
    MoleculeSnapshot, PipelineStatus, StandardizationPipeline, StandardizationReport,
    StandardizationStep, StandardizationStepReport, StandardizationWarning, StandardizeOptions,
    ZwitterionHandling, has_zwitterion, largest_fragment, neutralize_charges, normalize_zwitterion,
    remove_isotopes, remove_stereo, standardize,
};
pub use stereo::{assign_complete_stereochemistry, enumerate_stereoisomers, invert_stereocenter};
pub use tautomer::{
    TautomerConfig, canonical_tautomer, canonical_tautomer_with_config, enumerate_tautomers,
    enumerate_tautomers_with_config,
};
pub use mmff94_bci::{mmff94_charges_bci, mmff94_charges_typed, assign_mmff94_type, MmffType};
pub use topo_descriptors::{
    bertz_ct, chi0, chi0v, chi1, chi1v, chi2, chi2v, chi3, chi3v, chi4, chi4v, kappa1, kappa2,
    kappa3, labute_asa, labute_asa_per_atom, randic_index, topological_distance_matrix,
    wiener_index, zagreb_index_m1,
};
pub use vsa::{estate_vsa, peoe_vsa, slogp_vsa, smr_vsa};
pub use workflow::{
    CompareOptions, DescriptorDelta, DescriptorSummary, FilterSummary, FunctionalGroupSummary,
    MoleculeComparison, MoleculeReport, NamedGroupSummary, PairwiseComparison, ReportOptions,
    ScreenOptions, ScreeningRecord, ScreeningReport, SimilaritySummary, WorkflowError,
    WorkflowLimits, compare_molecules, compare_molecules_with_options, molecule_report,
    molecule_report_with_options, screen_smiles, screen_smiles_with_options,
};
pub use xlogp3::{xlogp3, xlogp3_per_atom};
