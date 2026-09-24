//! QED — Quantitative Estimate of Druglikeness.
//!
//! Implements the score from Bickerton et al. 2012 (Nature Chemistry 4, 90–98).
//!
//! The score is the weighted geometric mean of 8 desirability functions using
//! the WEIGHT_MEAN array from the paper.  Each desirability uses the
//! 7-parameter ADS (Asymmetric Double Sigmoidal) function:
//!
//! ```text
//! d(x) = (A + B / (1 + exp(-(x - C + D/2) / E))
//!              * (1 - 1 / (1 + exp(-(x - C - D/2) / F)))) / DMAX
//! ```
//!
//! Parameters and structural alerts are taken directly from RDKit's QED.py,
//! which implements Bickerton 2012 Supplementary Table S3.

use std::sync::OnceLock;

use chematic_core::{AtomIdx, BondOrder, Molecule};
use chematic_smarts::{MatchConfig, QueryMolecule, find_matches, parse_smarts};

use crate::descriptors::{
    RingBundle, hbd_count, logp_crippen, rdkit_molecular_weight, rdkit_tpsa, rotatable_bond_count,
};

/// ADS (Asymmetric Double Sigmoidal) desirability function from Bickerton 2012.
#[inline]
#[allow(clippy::too_many_arguments)]
fn ads(x: f64, a: f64, b: f64, c: f64, d: f64, e: f64, f: f64, dmax: f64) -> f64 {
    let exp1 = 1.0 + (-(x - c + d / 2.0) / e).exp();
    let exp2 = 1.0 + (-(x - c - d / 2.0) / f).exp();
    let dx = a + b / exp1 * (1.0 - 1.0 / exp2);
    dx / dmax
}

// ADS parameters — Bickerton 2012 Table S3 (via RDKit QED.py).
// Order: MW, ALOGP, HBA, HBD, PSA, ROTB, AROM, ALERTS.
// Each entry: (A, B, C, D, E, F, DMAX).

const ADS_PARAMS: [(f64, f64, f64, f64, f64, f64, f64); 8] = [
    (
        2.817_065_973,
        392.575_495_3,
        290.748_976_4,
        2.419_764_353,
        49.223_256_77,
        65.370_517_07,
        104.980_556_1,
    ), // MW
    (
        3.172_690_585,
        137.862_475_1,
        2.534_937_431,
        4.581_497_897,
        0.822_739_154,
        0.576_295_591,
        131.318_660_4,
    ), // ALOGP
    (
        2.948_620_388,
        160.460_597_2,
        3.615_294_657,
        4.435_986_202,
        0.290_141_953,
        1.300_669_958,
        148.776_304_6,
    ), // HBA
    (
        1.618_662_227,
        1_010.051_101,
        0.985_094_388,
        0.000_000_001,
        0.713_820_843,
        0.920_922_555,
        258.163_261_6,
    ), // HBD
    (
        1.876_861_559,
        125.223_265_7,
        62.907_735_54,
        87.833_666_14,
        12.019_998_24,
        28.513_247_32,
        104.568_616_7,
    ), // PSA
    (
        0.010_000_000,
        272.412_142_7,
        2.558_379_970,
        1.565_547_684,
        1.271_567_166,
        2.758_063_707,
        105.442_040_3,
    ), // ROTB
    (
        3.217_788_970,
        957.737_410_8,
        2.274_627_939,
        0.000_000_001,
        1.317_690_384,
        0.375_760_881,
        312.337_261_0,
    ), // AROM
    (
        0.010_000_000,
        1_199.094_025,
        -0.090_028_83,
        0.000_000_001,
        0.185_904_477,
        0.875_193_782,
        417.725_314_0,
    ), // ALERTS
];

// Weighted geometric mean weights (WEIGHT_MEAN from Bickerton 2012).
// Order: MW, ALOGP, HBA, HBD, PSA, ROTB, AROM, ALERTS.
const WEIGHTS_MEAN: [f64; 8] = [0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95];

// Structural alerts — Brenk 2008 / Bickerton 2012 (via RDKit QED.py).
// Patterns that fail to parse (isotopes, disconnected, unsupported syntax)
// are silently skipped at init time via filter_map.

static STRUCTURAL_ALERT_SMARTS: &[&str] = &[
    // Verbatim `rdkit.Chem.QED.StructuralAlertSmarts` (RDKit 2026.03.6), in order.
    "*1[O,S,N]*1",
    "[S,C](=[O,S])[F,Br,Cl,I]",
    "[CX4][Cl,Br,I]",
    "[#6]S(=O)(=O)O[#6]",
    "[$([CH]),$(CC)]#CC(=O)[#6]",
    "[$([CH]),$(CC)]#CC(=O)O[#6]",
    "n[OH]",
    "[$([CH]),$(CC)]#CS(=O)(=O)[#6]",
    "C=C(C=O)C=O",
    "n1c([F,Cl,Br,I])cccc1",
    "[CH1](=O)",
    "[#8][#8]",
    "[C;!R]=[N;!R]",
    "[N!R]=[N!R]",
    "[#6](=O)[#6](=O)",
    "[#16][#16]",
    "[#7][NH2]",
    "C(=O)N[NH2]",
    "[#6]=S",
    "[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]=[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]",
    "C1(=[O,N])C=CC(=[O,N])C=C1",
    "C1(=[O,N])C(=[O,N])C=CC=C1",
    "a21aa3a(aa1aaaa2)aaaa3",
    "a31a(a2a(aa1)aaaa2)aaaa3",
    "a1aa2a3a(a1)A=AA=A3=AA=A2",
    "c1cc([NH2])ccc1",
    "[Hg,Fe,As,Sb,Zn,Se,se,Te,B,Si,Na,Ca,Ge,Ag,Mg,K,Ba,Sr,Be,Ti,Mo,Mn,Ru,Pd,Ni,Cu,Au,Cd,Al,Ga,Sn,Rh,Tl,Bi,Nb,Li,Pb,Hf,Ho]",
    "I",
    "OS(=O)(=O)[O-]",
    "[N+](=O)[O-]",
    "C(=O)N[OH]",
    "C1NC(=O)NC(=O)1",
    "[SH]",
    "[S-]",
    "c1ccc([Cl,Br,I,F])c([Cl,Br,I,F])c1[Cl,Br,I,F]",
    "c1cc([Cl,Br,I,F])cc([Cl,Br,I,F])c1[Cl,Br,I,F]",
    "[CR1]1[CR1][CR1][CR1][CR1][CR1][CR1]1",
    "[CR1]1[CR1][CR1]cc[CR1][CR1]1",
    "[CR2]1[CR2][CR2][CR2][CR2][CR2][CR2][CR2]1",
    "[CR2]1[CR2][CR2]cc[CR2][CR2][CR2]1",
    "[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1",
    "[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1",
    "C#C",
    "[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]",
    "[$([N+R]),$([n+R]),$([N+]=C)][O-]",
    "[#6]=N[OH]",
    "[#6]=NOC=O",
    "[#6](=O)[CX4,CR0X3,O][#6](=O)",
    "c1ccc2c(c1)ccc(=O)o2",
    "[O+,o+,S+,s+]",
    "N=C=O",
    "[NX3,NX4][F,Cl,Br,I]",
    "c1ccccc1OC(=O)[#6]",
    "[CR0]=[CR0][CR0]=[CR0]",
    "[C+,c+,C-,c-]",
    "N=[N+]=[N-]",
    "C12C(NC(N1)=O)CSC2",
    "c1c([OH])c([OH,NH2,NH])ccc1",
    "P",
    "[N,O,S]C#N",
    "C=C=O",
    "[Si][F,Cl,Br,I]",
    "[SX2]O",
    "[SiR0,CR0](c1ccccc1)(c2ccccc2)(c3ccccc3)",
    "O1CCCCC1OC2CCC3CCCCC3C2",
    "N=[CR0][N,n,O,S]",
    "[cR2]1[cR2][cR2]([Nv3X3,Nv4X4])[cR2][cR2][cR2]1[cR2]2[cR2][cR2][cR2]([Nv3X3,Nv4X4])[cR2][cR2]2",
    "C=[C!r]C#N",
    "[cR2]1[cR2]c([N+0X3R0,nX3R0])c([N+0X3R0,nX3R0])[cR2][cR2]1",
    "[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2]c([N+0X3R0,nX3R0])[cR2]1",
    "[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2][cR2]c1([N+0X3R0,nX3R0])",
    "[OH]c1ccc([OH,NH2,NH])cc1",
    "c1ccccc1OC(=O)O",
    "[SX2H0][N]",
    "c12ccccc1(SC(S)=N2)",
    "c12ccccc1(SC(=S)N2)",
    "c1nnnn1C=O",
    "s1c(S)nnc1NC=O",
    "S1C=CSC1=S",
    "C(=O)Onnn",
    "OS(=O)(=O)C(F)(F)F",
    "N#CC[OH]",
    "N#CC(=O)",
    "S(=O)(=O)C#N",
    "N[CH2]C#N",
    "C1(=O)NCC1",
    "S(=O)(=O)[O-,OH]",
    "NC[F,Cl,Br,I]",
    "C=[C!r]O",
    "[NX2+0]=[O+0]",
    "[OR0,NR0][OR0,NR0]",
    "C(=O)O[C,H1].C(=O)O[C,H1].C(=O)O[C,H1]",
    "[CX2R0][NX3R0]",
    "c1ccccc1[C;!R]=[C;!R]c2ccccc2",
    "[NX3R0,NX4R0,OR0,SX2R0][CX4][NX3R0,NX4R0,OR0,SX2R0]",
    "[s,S,c,C,n,N,o,O]~[n+,N+](~[s,S,c,C,n,N,o,O])(~[s,S,c,C,n,N,o,O])~[s,S,c,C,n,N,o,O]",
    "[s,S,c,C,n,N,o,O]~[nX3+,NX3+](~[s,S,c,C,n,N])~[s,S,c,C,n,N]",
    "[*]=[N+]=[*]",
    "[SX3](=O)[O-,OH]",
    "N#N",
    "F.F.F.F",
    "[R0;D2][R0;D2][R0;D2][R0;D2]",
    "[cR,CR]~C(=O)NC(=O)~[cR,CR]",
    "C=!@CC=[O,S]",
    "[#6,#8,#16][#6](=O)O[#6]",
    "c[C;R0](=[O,S])[#6]",
    "c[SX2][C;!R]",
    "C=C=C",
    "c1nc([F,Cl,Br,I,S])ncc1",
    "c1ncnc([F,Cl,Br,I,S])c1",
    "c1nc(c2c(n1)nc(n2)[F,Cl,Br,I])",
    "[#6]S(=O)(=O)c1ccc(cc1)F",
    "[15N]",
    "[13C]",
    "[18O]",
    "[34S]",
];

/// Verbatim `rdkit.Chem.QED.AcceptorSmarts` (RDKit 2026.03.6), in order.
static ACCEPTOR_SMARTS: &[&str] = &[
    "[oH0;X2]",
    "[OH1;X2;v2]",
    "[OH0;X2;v2]",
    "[OH0;X1;v2]",
    "[O-;X1]",
    "[SH0;X2;v2]",
    "[SH0;X1;v2]",
    "[S-;X1]",
    "[nH0;X2]",
    "[NH0;X1;v3]",
    "[$([N;+0;X3;v3]);!$(N[C,S]=O)]",
];

static ALERT_QUERIES: OnceLock<Vec<QueryMolecule>> = OnceLock::new();
static ACCEPTOR_QUERIES: OnceLock<Vec<QueryMolecule>> = OnceLock::new();

fn alert_queries() -> &'static [QueryMolecule] {
    ALERT_QUERIES.get_or_init(|| {
        STRUCTURAL_ALERT_SMARTS
            .iter()
            .filter_map(|s| parse_smarts(s).ok())
            .collect()
    })
}

fn acceptor_queries() -> &'static [QueryMolecule] {
    ACCEPTOR_QUERIES.get_or_init(|| {
        ACCEPTOR_SMARTS
            .iter()
            .filter_map(|s| parse_smarts(s).ok())
            .collect()
    })
}

/// `sum(1 for alert in StructuralAlerts if mol.HasSubstructMatch(alert))`.
/// RDKit SMARTS always honour isotope labels (`[15N]` only matches ¹⁵N), so
/// the isotope alerts are matched with isotopes enforced.
fn structural_alert_count(view: &Molecule) -> usize {
    let config = MatchConfig {
        max_matches: Some(1),
        uniquify: false,
        use_isotopes: true,
        ..Default::default()
    };
    alert_queries()
        .iter()
        .filter(|q| chematic_smarts::has_match_with_config(q, view, &config))
        .count()
}

/// `sum(len(mol.GetSubstructMatches(p)) for p in Acceptors if mol.HasSubstructMatch(p))`
/// (uniquified matches, as `GetSubstructMatches` defaults).
fn qed_acceptor_count(view: &Molecule) -> usize {
    acceptor_queries()
        .iter()
        .map(|q| find_matches(q, view).len())
        .sum()
}

/// `len(Chem.GetSSSR(Chem.DeleteSubstructs(mol, AliphaticRings)))` with
/// `AliphaticRings = '[$([A;R][!a])]'`: delete every non-aromatic ring atom
/// that has a non-aromatic neighbour, then count SSSR rings of what remains
/// (the cycle rank, E - V + C over ring-eligible bonds).
fn qed_aromatic_ring_count(view: &Molecule) -> usize {
    let n = view.atom_count();
    let ring_atom = chematic_perception::ring_atom_flags(view);
    let deleted: Vec<bool> = (0..n)
        .map(|i| {
            let idx = AtomIdx(i as u32);
            ring_atom[i]
                && !view.atom(idx).aromatic
                && view.neighbors(idx).any(|(nb, _)| !view.atom(nb).aromatic)
        })
        .collect();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    let mut edges = 0usize;
    for (_, bond) in view.bonds() {
        if matches!(bond.order, BondOrder::Zero | BondOrder::Dative) {
            continue;
        }
        let (a, b) = (bond.atom1.0 as usize, bond.atom2.0 as usize);
        if deleted[a] || deleted[b] {
            continue;
        }
        edges += 1;
        let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
        if ra != rb {
            parent[ra] = rb;
        }
    }
    let kept: Vec<usize> = (0..n).filter(|&i| !deleted[i]).collect();
    let components = kept.iter().filter(|&&i| find(&mut parent, i) == i).count();
    (edges + components).saturating_sub(kept.len())
}

/// The eight QED properties exactly as `rdkit.Chem.QED.properties` defines
/// them, computed on the RDKit-perceived aromatic view of `mol`.
fn qed_properties(mol: &Molecule) -> [f64; 8] {
    let perceived = chematic_perception::apply_aromaticity_rdkit_parity_experimental(mol).ok();
    let view = perceived.as_ref().unwrap_or(mol);
    [
        rdkit_molecular_weight(mol),
        logp_crippen(mol),
        qed_acceptor_count(view) as f64,
        hbd_count(mol) as f64,
        rdkit_tpsa(mol),
        rotatable_bond_count(mol) as f64,
        qed_aromatic_ring_count(view) as f64,
        structural_alert_count(view) as f64,
    ]
}

fn qed_compute(props: [f64; 8]) -> f64 {
    let mut t = 0.0_f64;
    for (i, &pi) in props.iter().enumerate() {
        let (a, b, c, d, e, f, dmax) = ADS_PARAMS[i];
        let d_i = ads(pi, a, b, c, d, e, f, dmax).max(1e-10);
        t += WEIGHTS_MEAN[i] * d_i.ln();
    }
    let w_sum: f64 = WEIGHTS_MEAN.iter().sum();
    (t / w_sum).exp()
}

/// Compute the QED (Quantitative Estimate of Druglikeness) score.
///
/// Returns a value in `[0, 1]` where higher is more drug-like.
/// Uses the 7-parameter ADS desirability function and WEIGHT_MEAN from
/// Bickerton et al. 2012 (Nature Chemistry 4, 90–98).
///
/// # Example
/// ```
/// # use chematic_smiles::parse;
/// # use chematic_chem::qed;
/// let aspirin = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
/// let score = qed(&aspirin);
/// assert!(score > 0.4 && score < 0.7, "aspirin QED={score:.3}");
/// ```
pub fn qed(mol: &Molecule) -> f64 {
    qed_compute(qed_properties(mol))
}

/// Same value as [`qed`]. The ring bundle is no longer needed: every QED
/// property now follows `rdkit.Chem.QED.properties` directly. Kept for API
/// compatibility with callers that already hold a [`RingBundle`].
pub fn qed_with_bundle(mol: &Molecule, _rb: &RingBundle) -> f64 {
    qed(mol)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> chematic_core::Molecule {
        parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"))
    }

    #[test]
    fn test_qed_aspirin_rdkit_range() {
        // RDKit reference: ~0.553
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        let score = qed(&m);
        assert!(
            score > 0.40 && score < 0.72,
            "aspirin QED={score:.3} expected ~0.55"
        );
    }

    #[test]
    fn test_qed_caffeine_rdkit_range() {
        // RDKit reference: ~0.540
        let m = mol("Cn1cnc2c1c(=O)n(c(=O)n2C)C");
        let score = qed(&m);
        assert!(
            score > 0.40 && score < 0.72,
            "caffeine QED={score:.3} expected ~0.54"
        );
    }

    #[test]
    fn test_qed_ibuprofen_rdkit_range() {
        // RDKit reference: ~0.677
        let m = mol("CC(C)Cc1ccc(cc1)C(C)C(=O)O");
        let score = qed(&m);
        assert!(
            score > 0.50 && score < 0.85,
            "ibuprofen QED={score:.3} expected ~0.68"
        );
    }

    #[test]
    fn test_qed_paracetamol_rdkit_range() {
        // RDKit reference: ~0.551
        let m = mol("CC(=O)Nc1ccc(O)cc1");
        let score = qed(&m);
        assert!(
            score > 0.40 && score < 0.72,
            "paracetamol QED={score:.3} expected ~0.55"
        );
    }

    #[test]
    fn test_qed_benzene_range() {
        let m = mol("c1ccccc1");
        let score = qed(&m);
        assert!(score > 0.0 && score <= 1.0, "benzene QED={score:.3}");
    }

    #[test]
    fn test_qed_valid_range_for_common_molecules() {
        for smiles in &[
            "C",
            "CC",
            "c1ccccc1",
            "CCO",
            "CC(=O)O",
            "CC(=O)Oc1ccccc1C(=O)O",
            "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
        ] {
            let m = mol(smiles);
            let score = qed(&m);
            assert!(
                score > 0.0 && score <= 1.0,
                "QED for '{smiles}' = {score:.4} out of range (0, 1]"
            );
        }
    }

    #[test]
    fn test_structural_alert_count_aspirin() {
        // Aspirin has structural alerts: phenyl ester (c1ccccc1OC(=O)[#6])
        // and generic ester ([#6,#8,#16][#6](=O)O[#6]) — consistent with RDKit.
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        // RDKit 2026.03.6: QED.properties(aspirin) -> ALERTS=2, AROM=1, HBA=4.
        assert_eq!(structural_alert_count(&m), 2);
        assert_eq!(qed_aromatic_ring_count(&m), 1);
        assert_eq!(qed_acceptor_count(&m), 4);
    }

    #[test]
    fn test_structural_alert_nitro_compound() {
        // Nitrobenzene: [N+](=O)[O-] matches pattern 30
        let m = mol("c1ccc([N+](=O)[O-])cc1");
        // RDKit 2026.03.6: QED.properties(nitrobenzene).ALERTS == 2.
        assert_eq!(structural_alert_count(&m), 2);
    }

    #[test]
    fn test_structural_alert_queries_loaded() {
        // At least 80 patterns should parse successfully
        let queries = alert_queries();
        assert!(
            queries.len() >= 80,
            "expected >= 80 alert queries, got {}",
            queries.len()
        );
    }

    #[test]
    fn test_ads_function_smoke() {
        // ADS for MW=290 (near optimum) should be close to 1.0
        let (a, b, c, d, e, f, dmax) = ADS_PARAMS[0];
        let d_val = ads(290.0, a, b, c, d, e, f, dmax);
        assert!(d_val > 0.9 && d_val <= 1.1, "ADS(MW=290)={d_val:.4}");
    }
}
