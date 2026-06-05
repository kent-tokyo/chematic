//! `chematic-chem` — molecular descriptors for chematic.
//!
//! All descriptor functions take a `&Molecule` reference.
//! Values are approximate; calibrated against RDKit defaults.

#![forbid(unsafe_code)]

pub mod alerts;
pub mod brics;
pub mod cip;
pub mod named_groups;
pub mod descriptors;
pub mod diversity;
pub mod estate;
pub mod gasteiger;
pub mod hydrogen;
pub mod ifg;
pub mod mmp;
pub mod qed;
pub mod sa_score;
pub mod scaffold;
pub mod standardize;
pub mod tautomer;
pub mod topo_descriptors;
pub mod vsa;

pub use cip::{assign_cip, CipAssignment};
pub use descriptors::{
    aromatic_ring_count, egan_passes, exact_mass, formal_charge_sum, fsp3,
    ghose_passes, heavy_atom_count, hba_count, hbd_count, lipinski_passes,
    logp_crippen, logp_crippen_per_atom, molar_refractivity, mr_per_atom,
    molecular_weight, num_aliphatic_heterocycles,
    num_aliphatic_rings, num_aromatic_heterocycles, num_bridgehead_atoms,
    num_heteroatoms, num_saturated_heterocycles, num_saturated_rings,
    num_spiro_atoms, num_stereocenters, num_unspecified_stereocenters,
    reos_passes, ring_count, rotatable_bond_count, tpsa, veber_passes,
};

pub use alerts::{pains_matches, pains_passes};
pub use brics::{brics_bonds, brics_fragments};
pub use mmp::{MmpPair, find_mmp};
pub use diversity::{butina_cluster, maxmin_picks};
pub use gasteiger::gasteiger_charges;
pub use hydrogen::{add_hydrogens, remove_hydrogens};
pub use ifg::{identify_functional_groups, FunctionalGroup};
pub use named_groups::{detect_named_functional_groups, NamedGroup};
pub use qed::qed;
pub use sa_score::sa_score;
pub use scaffold::{generic_murcko_scaffold, murcko_scaffold};
pub use standardize::{largest_fragment, neutralize_charges};
pub use tautomer::{canonical_tautomer, enumerate_tautomers};
pub use topo_descriptors::{
    bertz_ct,
    chi0, chi1, chi2, chi3, chi4,
    chi0v, chi1v, chi2v, chi3v, chi4v,
    kappa1, kappa2, kappa3,
    labute_asa, labute_asa_per_atom,
    wiener_index,
};
pub use estate::{estate_indices, max_estate, min_estate, sum_estate};
pub use vsa::{peoe_vsa, slogp_vsa, smr_vsa};
