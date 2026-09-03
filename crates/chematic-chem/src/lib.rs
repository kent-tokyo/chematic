//! `chematic-chem` — molecular descriptors for chematic.
//!
//! All descriptor functions take a `&Molecule` reference.
//! Values are approximate; calibrated against RDKit defaults.

#![forbid(unsafe_code)]

pub mod abbreviations;
pub mod activity_cliff;
pub mod admet;
pub mod alerts;
pub mod atropisomer;
pub mod brics;
pub mod cache;
pub mod canonical;
pub mod cip;
pub mod condensed;
pub mod descriptors;
pub mod diversity;
pub mod drug_score;
pub mod esol;
pub mod estate;
pub mod formula;
pub mod gasteiger;
pub mod hash;
pub mod hydrogen;
pub mod ifg;
pub mod isotope_distribution;
pub mod iupac_stereo;
pub mod logd;
pub mod mlp;
pub mod mmff94_bci;
pub mod mmp;
pub mod named_groups;
pub mod parent;
pub mod pka;
pub mod qed;
pub mod recap;
pub mod rgroup;
pub mod sa_score;
pub mod scaffold;
pub mod standardize;
pub mod stereo;
pub mod tautomer;
pub mod topo_descriptors;
pub mod vsa;
pub mod workflow;
pub mod xlogp3;

pub use cip::{
    CipAssignment, CipMode, CipModeAssignment, CipModeError, CipUnresolvedReason, EzCompleteness,
    assign_cip, assign_cip_with_mode, ez_completeness, tetrahedral_stereo_neighbors,
};
pub use descriptors::{
    Bcut2D, CarbonTypes, InformationContent, RingBundle, aromatic_ring_count, autocorr_2d,
    balaban_j, bcut2d, calc_mol_formula, carbon_types, cns_mpo_from_parts, cns_mpo_score,
    egan_passes, exact_mass, formal_charge_per_atom, formal_charge_sum, fraction_rotatable_bonds,
    fsp3, geary_autocorr, ghose_passes, hall_kier_alpha, hba_count, hba_count_lipinski, hbd_count,
    heavy_atom_count, hybridization_per_atom, implicit_hcount_per_atom, information_content, ipc,
    lead_like_passes, lipinski_passes, logp_and_mr, logp_crippen, logp_crippen_per_atom,
    mcf_passes, mde_carbon, mmff94_charges, molar_refractivity, molecular_weight, moran_autocorr,
    mqn, mr_per_atom, num_aliphatic_heterocycles, num_aliphatic_rings, num_amide_bonds,
    num_aromatic_heterocycles, num_bridgehead_atoms, num_bromines, num_carbons, num_chlorines,
    num_ester_bonds, num_fluorines, num_heteroatoms, num_hydrogens, num_iodines, num_nitrogens,
    num_oxygens, num_phosphorus, num_saturated_heterocycles, num_saturated_rings, num_spiro_atoms,
    num_stereocenters, num_sulfurs, num_unspecified_stereocenters, pfizer_3_75_passes, reos_passes,
    ring_bundle, ring_count, ring_system_count, ro3_passes, rotatable_bond_atom_pairs,
    rotatable_bond_count, tpsa, tpsa_per_atom, usrcat, veber_passes,
};
pub use iupac_stereo::iupac_name_stereo;

pub use abbreviations::{abbreviations, expand_abbreviation};
pub use activity_cliff::{ActivityCliff, activity_cliffs};
pub use admet::{
    AdmetProfile, BoiledEggProfile, ClearanceClass, admet_profile, ames_alerts, ames_passes,
    ames_risk_score, bbb_passes, bbb_score, bbb_score_from_parts, boiled_egg, boiled_egg_from,
    caco2_permeability, caco2_precomputed, clearance_class, clearance_score,
    cyp3a4_inhibition_risk, cyp3a4_precomputed, herg_risk_precomputed, herg_risk_score,
    ppb_percent,
};
pub use alerts::{
    brenk_matches, brenk_matches_checked, brenk_matches_detailed, brenk_passes,
    brenk_passes_and_matches, pains_matches, pains_matches_checked, pains_matches_detailed,
    pains_passes, pains_passes_and_matches,
};
pub use atropisomer::{AtropisomerType, assign_atropisomer_chirality, detect_atropisomers};
pub use brics::{BricsConfig, brics_bonds, brics_fragments, brics_fragments_with_config};
pub use cache::{DescriptorCache, DescriptorEntry};
pub use canonical::{CanonicalMode, canonical_smiles_mode};
pub use condensed::{CondensedError, parse_condensed};
pub use diversity::{butina_cluster, maxmin_picks};
pub use drug_score::drug_score;
pub use esol::esol_solubility;
pub use estate::{estate_all, estate_indices, max_estate, min_estate, sum_estate};
pub use formula::{FormulaParseError, parse_formula};
pub use gasteiger::gasteiger_charges;
pub use hash::{are_identical, mol_hash, stable_are_identical};
pub use hydrogen::{add_hydrogens, remove_hydrogens};
pub use ifg::{FunctionalGroup, identify_functional_groups};
pub use isotope_distribution::isotope_distribution;
pub use logd::{logd_from_logp, logd_profile, logd_simple};
pub use mlp::{MLP_SOLUBILITY_TRAINED, mlp_solubility};
pub use mmff94_bci::{MmffType, assign_mmff94_type, mmff94_charges_bci, mmff94_charges_typed};
pub use mmp::{MmpPair, MmsMember, MmsSeries, find_mmp, find_mms};
pub use named_groups::{NamedGroup, detect_named_functional_groups};
pub use parent::{
    AbstainReason, InvalidInputReason, ParentAudit, ParentComputationStatus, ParentResult,
    super_parent,
};
pub use pka::{PkaSite, PkaSiteType, pka_acid, pka_base, pka_both, predict_pka};
pub use qed::{qed, qed_with_bundle};
pub use recap::{recap_breakable_bond_count, recap_fragment};
pub use rgroup::{RGroupError, RGroupResult, rgroup_decompose};
pub use sa_score::{sa_score, sa_score_with_bundle};
pub use scaffold::{
    ScaffoldNetwork, generic_murcko_scaffold, murcko_scaffold, scaffold_network,
    scaffold_network_with_counts, schuffenhauer_parents,
};
pub use standardize::{
    FragmentDecision, FragmentPolicy, FragmentRecord, FragmentSnapshot, MoleculeSnapshot,
    PipelineStatus, StandardizationPipeline, StandardizationReport, StandardizationStep,
    StandardizationStepReport, StandardizationWarning, StandardizeOptions, TransformationRecord,
    ZwitterionHandling, charge_parent, fragment_parent, has_zwitterion, isotope_parent,
    largest_fragment, neutralize_charges, normalize_groups, normalize_zwitterion, prefer_organic,
    reionize, remove_isotopes, remove_stereo, select_fragment, standardize, stereo_parent,
    uncharge,
};
pub use stereo::{assign_complete_stereochemistry, enumerate_stereoisomers, invert_stereocenter};
pub use tautomer::{
    AppliedTransform, ScoreContribution, TautomerAuditRecord, TautomerConfig, TautomerLimits,
    TautomerRuleId, TautomerScoreTerm, canonical_tautomer, canonical_tautomer_with_config,
    enumerate_tautomers, enumerate_tautomers_with_config, tautomer_parent,
};
pub use topo_descriptors::{
    bertz_ct, chi_all, chi0, chi0v, chi1, chi1v, chi2, chi2v, chi3, chi3v, chi4, chi4v,
    eccentric_connectivity_index, graph_diameter, graph_eccentricities, graph_radius,
    gravitational_index, gutman_mti, hosoya_index, kappa_all, kappa1, kappa2, kappa3, labute_asa,
    labute_asa_per_atom, num_valence_electrons, padmakar_ivan_index, petitjean_index, randic_index,
    schultz_mti, topological_distance_matrix, vabc, wiener_index, zagreb_index_m1, zagreb_index_m2,
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
