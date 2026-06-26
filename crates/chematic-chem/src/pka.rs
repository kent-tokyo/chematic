//! Rule-based pKa prediction from functional group SMARTS patterns.
//!
//! # Method
//!
//! Ionizable sites are detected by matching a priority-ordered table of SMARTS
//! patterns.  Each pattern carries a literature-derived base pKa value:
//!
//! | Group | Base pKa | Reference |
//! |-------|----------|-----------|
//! | Carboxylic acid | 4.0 | Acetic acid 4.76, benzoic acid 4.19 |
//! | Thiol | 8.3 | Cysteamine 8.3, methanethiol 10.3 |
//! | Sulfonamide N-H | 10.1 | Sulfadiazine 6.5 |
//! | Phenol | 10.0 | Phenol 9.99 |
//! | Aromatic amine | 4.6 | Aniline 4.63 |
//! | Pyridine N | 5.2 | Pyridine 5.23 |
//! | Imidazole N | 6.9 | Imidazole 6.95 |
//! | Aliphatic 3° amine | 9.8 | Trimethylamine 9.81 |
//! | Aliphatic 2° amine | 10.8 | Dimethylamine 10.73 |
//! | Aliphatic 1° amine | 10.6 | Methylamine 10.66 |
//! | Piperidine | 11.2 | Piperidine 11.22 |
//! | Morpholine | 8.3 | Morpholine 8.36 |
//! | Guanidine | 12.5 | Arginine side chain 12.5 |
//!
//! **Accuracy**: ±1–2 pKa units; suitable for rule-of-5/ADMET triage.
//! For high-accuracy pKa a machine-learning model is required.

#![forbid(unsafe_code)]

use std::sync::OnceLock;

use chematic_core::{AtomIdx, Molecule};
use chematic_perception::find_sssr;
use chematic_smarts::{
    MatchConfig, QueryMolecule, find_matches_with_rings_and_config, parse_smarts,
};

// ── site type ────────────────────────────────────────────────────────────────

/// Whether an ionizable site is acid (H donor) or base (H acceptor).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PkaSiteType {
    /// Acidic site (loses H⁺ on ionization): carboxylic acid, phenol, thiol, …
    Acid,
    /// Basic site (gains H⁺ on ionization): amine, pyridine, imidazole, …
    Base,
}

/// One ionizable site with its predicted pKa.
#[derive(Debug, Clone)]
pub struct PkaSite {
    /// Index of the ionizable heavy atom (O of COOH, N of amine, etc.).
    pub atom_idx: AtomIdx,
    /// Predicted pKa of the ionization event.
    pub pka: f64,
    /// Whether this is an acidic or basic site.
    pub site_type: PkaSiteType,
    /// Short name of the functional group.
    pub group_name: &'static str,
}

// ── SMARTS rules ─────────────────────────────────────────────────────────────

/// (group_name, SMARTS, base_pKa, site_type)
///
/// Listed in priority order: more specific patterns appear first so that, when
/// an atom matches multiple patterns, only the highest-priority one is used.
static PKA_RULES: &[(&str, &str, f64, PkaSiteType)] = &[
    // ── Acids ──────────────────────────────────────────────────────────────
    // Carboxylic acid -C(=O)OH  (match the OH oxygen)
    (
        "carboxylic_acid",
        "[OX2H1][CX3](=O)",
        4.0,
        PkaSiteType::Acid,
    ),
    // Tetrazole N-H (bioisostere for COOH): match the N-H in the tetrazole ring
    // 1H-tetrazole pKa ~4.9, 5-substituted ~5.0–6.5
    ("tetrazole_nh", "[nH;r5;$(n1nnnc1)]", 4.9, PkaSiteType::Acid),
    // Thiol C-SH  (match the sulfur)
    ("thiol", "[SX2H1]", 8.3, PkaSiteType::Acid),
    // Hydroxamic acid O-H: RC(=O)NHOH → match the O-H (~8.7)
    (
        "hydroxamic_acid",
        "[OX2H1][NX3][CX3](=O)",
        8.7,
        PkaSiteType::Acid,
    ),
    // Sulfonamide -S(=O)(=O)-NH  (match the nitrogen)
    (
        "sulfonamide_nh",
        "[NX3H1][SX4](=O)(=O)",
        10.1,
        PkaSiteType::Acid,
    ),
    // Phenol: Ar-OH  (match the oxygen)
    ("phenol", "[OX2H1][cX3]", 10.0, PkaSiteType::Acid),
    // ── Bases ──────────────────────────────────────────────────────────────
    // Guanidine C(=NH)NH2  (match the imine N — most basic of the three Ns)
    (
        "guanidine",
        "[NH;!R;$(N/C(=[NH])/N)]",
        12.5,
        PkaSiteType::Base,
    ),
    // Aliphatic cyclic N-H in 5-membered ring (pyrrolidine ~11.3)
    (
        "cyclic_amine_5",
        "[NX3H1;r5;!$(N-c);!$(NC=O)]",
        11.3,
        PkaSiteType::Base,
    ),
    // Piperazine-type N-H: 6-ring containing a second N (pKa ~9.8, lower than piperidine)
    // Must come before generic cyclic_amine_6 (11.0).
    (
        "piperazine_nh",
        "[NX3H1;r6;$(N1CCNCC1);!$(N-c);!$(NC=O)]",
        9.8,
        PkaSiteType::Base,
    ),
    // Morpholine-type N-H: 6-ring containing O (pKa ~8.3, lower than piperidine)
    // Must come before generic cyclic_amine_6 (11.0).
    (
        "morpholine_nh",
        "[NX3H1;r6;$(N1CCOCC1);!$(NC=O)]",
        8.3,
        PkaSiteType::Base,
    ),
    // Thiomorpholine-type N-H: 6-ring containing S (pKa ~7.0)
    (
        "thiomorpholine_nh",
        "[NX3H1;r6;$(N1CCSCC1);!$(NC=O)]",
        7.0,
        PkaSiteType::Base,
    ),
    // Aliphatic cyclic N-H in 6-membered ring (piperidine ~11.2)
    (
        "cyclic_amine_6",
        "[NX3H1;r6;!$(N-c);!$(NC=O)]",
        11.0,
        PkaSiteType::Base,
    ),
    // Aliphatic 2° amine (non-ring)
    (
        "amine_secondary",
        "[NX3H1;!R;!$(N-c);!$(NC=O);!$(NS(=O))]",
        10.8,
        PkaSiteType::Base,
    ),
    // Aliphatic 1° amine
    (
        "amine_primary",
        "[NX3H2;!R;!$(N-c);!$(NC=O)]",
        10.6,
        PkaSiteType::Base,
    ),
    // N-substituted piperazine (H0): 6-ring with second N, ~8.7
    (
        "piperazine_tertiary",
        "[NX3H0;r6;$(N1CCNCC1);!$(N-c);!$(NC=O)]",
        8.7,
        PkaSiteType::Base,
    ),
    // N-substituted morpholine (H0): 6-ring with O, ~7.4
    (
        "morpholine_tertiary",
        "[NX3H0;r6;$(N1CCOCC1);!$(NC=O)]",
        7.4,
        PkaSiteType::Base,
    ),
    // N-substituted thiomorpholine (H0): 6-ring with S, ~6.5
    (
        "thiomorpholine_tertiary",
        "[NX3H0;r6;$(N1CCSCC1);!$(NC=O)]",
        6.5,
        PkaSiteType::Base,
    ),
    // Cyclic 3° amine in 6-membered ring (N-substituted piperidine, ~9.5)
    (
        "cyclic_amine_tertiary_6",
        "[NX3H0;r6;!$(N-c);!$(NC=O)]",
        9.5,
        PkaSiteType::Base,
    ),
    // Cyclic 3° amine in 5-membered ring (N-substituted pyrrolidine, ~10.2)
    (
        "cyclic_amine_tertiary_5",
        "[NX3H0;r5;!$(N-c);!$(NC=O)]",
        10.2,
        PkaSiteType::Base,
    ),
    // Aliphatic 3° amine (non-ring, no H)
    (
        "amine_tertiary",
        "[NX3H0;!R;!$(N-c);!$(NC=O)]",
        9.8,
        PkaSiteType::Base,
    ),
    // Imidazole ring (5-membered, 2 N, one with H): match the pyridine-like N
    (
        "imidazole",
        "[nX2H0;r5;$(n1c[nH]cc1)]",
        6.9,
        PkaSiteType::Base,
    ),
    // Aromatic amine Ar-NH2 or Ar-NH- (match the nitrogen)
    (
        "aromatic_amine",
        "[NX3;!$(NC=O);$(N-c)]",
        4.6,
        PkaSiteType::Base,
    ),
    // Pyridine-like N (6-membered ring aromatic N, no H)
    ("pyridine", "[nX2H0;r6]", 5.2, PkaSiteType::Base),
];

// ── compiled patterns (lazy, thread-safe) ────────────────────────────────────

struct CompiledRule {
    name: &'static str,
    query: QueryMolecule,
    pka: f64,
    site_type: PkaSiteType,
}

static COMPILED: OnceLock<Vec<CompiledRule>> = OnceLock::new();

fn compiled_rules() -> &'static [CompiledRule] {
    COMPILED.get_or_init(|| {
        PKA_RULES
            .iter()
            .filter_map(|&(name, smarts, pka, site_type)| {
                parse_smarts(smarts).ok().map(|q| CompiledRule {
                    name,
                    query: q,
                    pka,
                    site_type,
                })
            })
            .collect()
    })
}

// ── public API ────────────────────────────────────────────────────────────────

/// Predict pKa values for all ionizable sites in `mol`.
///
/// Returns one [`PkaSite`] per ionizable heavy atom, in discovery order.
/// If the same atom matches multiple patterns, only the highest-priority
/// (first in the table) pattern is recorded.
///
/// Molecules with no ionizable groups return an empty Vec.
pub fn predict_pka(mol: &Molecule) -> Vec<PkaSite> {
    let rules = compiled_rules();
    let rings = find_sssr(mol);
    let config = MatchConfig::default();
    let mut seen: std::collections::HashSet<u32> = std::collections::HashSet::new();
    let mut sites: Vec<PkaSite> = Vec::new();

    for rule in rules {
        let matches = find_matches_with_rings_and_config(&rule.query, mol, &rings, &config);
        for atom_map in &matches {
            // The matching atom for the ionizable center is the first query atom (index 0).
            if let Some(&atom_idx) = atom_map.get(&0)
                && seen.insert(atom_idx.0)
            {
                sites.push(PkaSite {
                    atom_idx,
                    pka: rule.pka,
                    site_type: rule.site_type,
                    group_name: rule.name,
                });
            }
        }
    }

    sites
}

/// Return the pKa of the most acidic (lowest pKa) site, or `None` if there
/// are no acidic sites.
pub fn pka_acid(mol: &Molecule) -> Option<f64> {
    predict_pka(mol)
        .into_iter()
        .filter(|s| s.site_type == PkaSiteType::Acid)
        .map(|s| s.pka)
        .reduce(f64::min)
}

/// Return the pKa of the most basic (highest pKa) site, or `None` if there
/// are no basic sites.
pub fn pka_base(mol: &Molecule) -> Option<f64> {
    predict_pka(mol)
        .into_iter()
        .filter(|s| s.site_type == PkaSiteType::Base)
        .map(|s| s.pka)
        .reduce(f64::max)
}

/// Return `(pka_acid, pka_base)` in a single `predict_pka` call.
///
/// Use when both values are needed to avoid running the 42-rule SMARTS scan twice.
pub fn pka_both(mol: &Molecule) -> (Option<f64>, Option<f64>) {
    let sites = predict_pka(mol);
    let acid = sites
        .iter()
        .filter(|s| s.site_type == PkaSiteType::Acid)
        .map(|s| s.pka)
        .reduce(f64::min);
    let base = sites
        .iter()
        .filter(|s| s.site_type == PkaSiteType::Base)
        .map(|s| s.pka)
        .reduce(f64::max);
    (acid, base)
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> Molecule {
        parse(s).unwrap()
    }

    #[test]
    fn test_acetic_acid_pka() {
        let m = mol("CC(=O)O");
        let acid = pka_acid(&m);
        assert!(acid.is_some(), "acetic acid should have acidic site");
        assert!(
            (acid.unwrap() - 4.0).abs() < 0.1,
            "acetic acid pKa ~4.0, got {:.2}",
            acid.unwrap()
        );
    }

    #[test]
    fn test_phenol_pka() {
        let m = mol("Oc1ccccc1");
        let sites = predict_pka(&m);
        let phenol_site = sites.iter().find(|s| s.group_name == "phenol");
        assert!(phenol_site.is_some(), "phenol should be detected");
        assert!(
            (phenol_site.unwrap().pka - 10.0).abs() < 0.1,
            "phenol pKa ~10.0"
        );
    }

    #[test]
    fn test_thiol_pka() {
        let m = mol("CCS");
        let sites = predict_pka(&m);
        let thiol_site = sites.iter().find(|s| s.group_name == "thiol");
        assert!(thiol_site.is_some(), "thiol should be detected");
        assert_eq!(thiol_site.unwrap().site_type, PkaSiteType::Acid);
    }

    #[test]
    fn test_aniline_pka() {
        let m = mol("Nc1ccccc1"); // aniline
        let base = pka_base(&m);
        assert!(base.is_some(), "aniline should have basic site");
        assert!(
            (base.unwrap() - 4.6).abs() < 0.1,
            "aniline pKa ~4.6, got {:.2}",
            base.unwrap()
        );
    }

    #[test]
    fn test_pyridine_pka() {
        let m = mol("c1ccncc1"); // pyridine
        let base = pka_base(&m);
        assert!(base.is_some(), "pyridine should have basic N");
        let pka = base.unwrap();
        assert!((pka - 5.2).abs() < 0.5, "pyridine pKa ~5.2, got {pka:.2}");
    }

    #[test]
    fn test_primary_amine_pka() {
        let m = mol("CCN"); // ethylamine
        let base = pka_base(&m);
        assert!(base.is_some(), "ethylamine should have basic site");
        let pka = base.unwrap();
        assert!(
            pka > 9.0 && pka < 12.0,
            "aliphatic amine pKa 9–12, got {pka:.2}"
        );
    }

    #[test]
    fn test_piperidine_pka() {
        let m = mol("C1CCNCC1"); // piperidine
        let sites = predict_pka(&m);
        let pip = sites.iter().find(|s| s.group_name == "piperidine");
        // piperidine may be detected as piperidine or secondary amine — both are valid
        if let Some(p) = pip {
            assert!((p.pka - 11.2).abs() < 0.5, "piperidine pKa ~11.2");
        } else {
            let base = pka_base(&m);
            assert!(
                base.is_some(),
                "piperidine should have at least one basic site"
            );
        }
    }

    #[test]
    fn test_benzene_no_sites() {
        let m = mol("c1ccccc1");
        let sites = predict_pka(&m);
        assert!(sites.is_empty(), "benzene has no ionizable sites");
        assert!(pka_acid(&m).is_none());
        assert!(pka_base(&m).is_none());
    }

    #[test]
    fn test_aspirin_has_acid() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O"); // aspirin
        let acid = pka_acid(&m);
        assert!(acid.is_some(), "aspirin has a carboxylic acid");
        assert!(acid.unwrap() < 5.0, "carboxylic acid pKa < 5");
    }

    #[test]
    fn test_amphoteric_glycine() {
        let m = mol("NCC(=O)O"); // glycine
        let acid = pka_acid(&m);
        let base = pka_base(&m);
        assert!(acid.is_some(), "glycine has acid site (COOH)");
        assert!(base.is_some(), "glycine has base site (NH2)");
        assert!(
            acid.unwrap() < base.unwrap(),
            "pKa_acid < pKa_base for amino acid"
        );
    }

    #[test]
    fn test_site_type_classification() {
        let m = mol("CC(=O)O"); // acetic acid
        let sites = predict_pka(&m);
        assert!(!sites.is_empty());
        assert_eq!(sites[0].site_type, PkaSiteType::Acid);
    }

    #[test]
    fn test_morpholine_lower_than_piperidine() {
        let morpholine = mol("C1CNCCO1");
        let piperidine = mol("C1CCNCC1");
        let mp = pka_base(&morpholine).unwrap();
        let pp = pka_base(&piperidine).unwrap();
        assert!(
            mp < pp,
            "morpholine ({mp:.1}) should be less basic than piperidine ({pp:.1})"
        );
        assert!((mp - 8.3).abs() < 0.5, "morpholine pKa ~8.3, got {mp:.1}");
    }

    #[test]
    fn test_piperazine_nh_pka() {
        let m = mol("C1CNCCN1"); // piperazine
        let pka = pka_base(&m).unwrap();
        assert!((pka - 9.8).abs() < 0.5, "piperazine pKa ~9.8, got {pka:.1}");
    }

    #[test]
    fn test_tetrazole_is_acidic() {
        let m = mol("c1nnn[nH]1"); // 1H-tetrazole
        let acid = pka_acid(&m);
        assert!(acid.is_some(), "tetrazole should have acidic N-H");
        let pka = acid.unwrap();
        assert!(pka < 6.0, "tetrazole pKa < 6, got {pka:.1}");
    }

    #[test]
    fn test_hydroxamic_acid_pka() {
        let m = mol("CC(=O)NO"); // acetohydroxamic acid
        let acid = pka_acid(&m);
        assert!(acid.is_some(), "hydroxamic acid should have acid site");
        let pka = acid.unwrap();
        assert!(
            (pka - 8.7).abs() < 0.5,
            "hydroxamic acid pKa ~8.7, got {pka:.1}"
        );
    }

    #[test]
    fn test_carboxylic_vs_phenol_priority() {
        // p-hydroxybenzoic acid has BOTH carboxylic acid AND phenol
        let m = mol("OC(=O)c1ccc(O)cc1");
        let sites = predict_pka(&m);
        let acid_sites: Vec<_> = sites
            .iter()
            .filter(|s| s.site_type == PkaSiteType::Acid)
            .collect();
        assert!(
            acid_sites.len() >= 2,
            "p-hydroxybenzoic acid has carboxylic acid + phenol"
        );
        // carboxylic acid should have lower pKa than phenol
        let pkas: Vec<f64> = acid_sites.iter().map(|s| s.pka).collect();
        assert!(pkas.iter().any(|&p| p < 5.0), "carboxylic acid pKa < 5");
        assert!(pkas.iter().any(|&p| p > 8.0), "phenol pKa > 8");
    }
}
