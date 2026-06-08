//! VSA (van der Waals surface area) binned descriptors.
//!
//! Bin boundaries from RDKit `rdkit/Chem/MolSurf.py`.
//!
//! - SlogP_VSA1–12: bins per-atom Labute ASA by Crippen LogP contribution.
//! - SMR_VSA1–10:   bins per-atom Labute ASA by Crippen MR contribution.
//! - PEOE_VSA1–14:  bins per-atom Labute ASA by Gasteiger partial charge.

use crate::descriptors::{logp_crippen_per_atom, mr_per_atom};
use crate::gasteiger::gasteiger_charges;
use crate::topo_descriptors::labute_asa_per_atom;
use chematic_core::Molecule;

/// SlogP_VSA bin boundaries (11 cuts → 12 bins).
/// Source: RDKit MolSurf.py _slogpVSA_bins
const SLOGP_CUTS: &[f64] = &[-0.4, -0.2, 0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6];

/// SMR_VSA bin boundaries (9 cuts → 10 bins).
/// Source: RDKit MolSurf.py _smrVSA_bins
const SMR_CUTS: &[f64] = &[1.29, 1.82, 2.24, 2.45, 2.75, 3.05, 3.63, 3.8, 4.0];

/// PEOE_VSA bin boundaries (13 cuts → 14 bins).
/// Source: RDKit MolSurf.py _peoe_VSA_bins
const PEOE_CUTS: &[f64] = &[
    -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3,
];

fn bin_idx(value: f64, cuts: &[f64]) -> usize {
    cuts.partition_point(|&c| value >= c)
}

fn vsa_bins(mol: &Molecule, contrib: &[f64], cuts: &[f64]) -> Vec<f64> {
    let asa = labute_asa_per_atom(mol);
    let nbins = cuts.len() + 1;
    let mut bins = vec![0.0f64; nbins];
    for (i, (&c, &a)) in contrib.iter().zip(asa.iter()).enumerate() {
        let _ = i;
        bins[bin_idx(c, cuts)] += a;
    }
    bins
}

/// SlogP_VSA descriptors (12 values).
///
/// Returns sum of Labute ASA contributions for atoms whose Crippen LogP
/// contribution falls in each bin. Bin boundaries (RDKit):
/// -0.4, -0.2, 0.0, 0.1, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60
pub fn slogp_vsa(mol: &Molecule) -> Vec<f64> {
    vsa_bins(mol, &logp_crippen_per_atom(mol), SLOGP_CUTS)
}

/// SMR_VSA descriptors (10 values).
///
/// Returns sum of Labute ASA contributions for atoms whose Crippen MR
/// contribution falls in each bin. Bin boundaries (RDKit):
/// 1.29, 1.82, 2.24, 2.45, 2.75, 3.05, 3.63, 3.80, 4.00
pub fn smr_vsa(mol: &Molecule) -> Vec<f64> {
    vsa_bins(mol, &mr_per_atom(mol), SMR_CUTS)
}

/// PEOE_VSA descriptors (14 values).
///
/// Returns sum of Labute ASA contributions for atoms whose Gasteiger partial
/// charge falls in each bin. Bin boundaries (RDKit):
/// -0.30, -0.25, -0.20, -0.15, -0.10, -0.05, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30
pub fn peoe_vsa(mol: &Molecule) -> Vec<f64> {
    vsa_bins(mol, &gasteiger_charges(mol), PEOE_CUTS)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> Molecule {
        parse(s).unwrap()
    }

    #[test]
    fn slogp_vsa_length_is_12() {
        assert_eq!(slogp_vsa(&mol("CC(=O)O")).len(), 12);
    }

    #[test]
    fn smr_vsa_length_is_10() {
        assert_eq!(smr_vsa(&mol("CC(=O)O")).len(), 10);
    }

    #[test]
    fn peoe_vsa_length_is_14() {
        assert_eq!(peoe_vsa(&mol("CC(=O)O")).len(), 14);
    }

    #[test]
    fn slogp_vsa_sum_equals_labute_asa() {
        use crate::topo_descriptors::labute_asa;
        let mol = mol("c1ccccc1");
        let total: f64 = slogp_vsa(&mol).iter().sum();
        let expected = labute_asa(&mol);
        assert!(
            (total - expected).abs() < 1e-6,
            "SlogP_VSA sum {total:.4} != LabuteASA {expected:.4}"
        );
    }

    #[test]
    fn smr_vsa_sum_equals_labute_asa() {
        use crate::topo_descriptors::labute_asa;
        let mol = mol("CC(=O)Oc1ccccc1C(=O)O");
        let total: f64 = smr_vsa(&mol).iter().sum();
        let expected = labute_asa(&mol);
        assert!(
            (total - expected).abs() < 1e-6,
            "SMR_VSA sum {total:.4} != LabuteASA {expected:.4}"
        );
    }

    #[test]
    fn peoe_vsa_sum_equals_labute_asa() {
        use crate::topo_descriptors::labute_asa;
        let mol = mol("c1ccccc1");
        let total: f64 = peoe_vsa(&mol).iter().sum();
        let expected = labute_asa(&mol);
        assert!(
            (total - expected).abs() < 1e-6,
            "PEOE_VSA sum {total:.4} != LabuteASA {expected:.4}"
        );
    }

    #[test]
    fn empty_molecule_zero_bins() {
        let mol = mol("C");
        assert!(
            slogp_vsa(&mol).iter().any(|&v| v > 0.0),
            "methane should have nonzero VSA"
        );
    }
}
