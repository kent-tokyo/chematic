//! Explainable structural reactivity findings related to genotoxicity triage.
//!
//! This module deliberately reports motifs and evidence only.  It is not a
//! genotoxicity predictor, score, classifier, or biological risk assessment.

use std::sync::OnceLock;

use chematic_core::{AtomIdx, Molecule};
use chematic_smarts::{QueryMolecule, find_matches, parse_smarts};

/// Confidence attached to a finding produced by a deterministic SMARTS rule.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum GenotoxConfidence {
    /// The finding is supported by a structural pattern match only.
    PatternOnly,
}

/// Applicability boundary for the current rule set.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum GenotoxApplicability {
    /// A small organic molecule represented as a molecular graph.
    OrganicSmallMolecule,
}

/// Structural motif associated with electrophilic reactivity triage.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum GenotoxMotifKind {
    Epoxide,
    Aziridine,
    MichaelAcceptor,
}

/// One explainable structural finding.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenotoxFinding {
    pub motif: GenotoxMotifKind,
    pub pattern_name: &'static str,
    /// Target atom indices in deterministic ascending order.
    pub matched_atoms: Vec<AtomIdx>,
    pub evidence: &'static str,
    pub interpretation: &'static str,
    pub confidence: GenotoxConfidence,
    pub applicability: GenotoxApplicability,
}

/// Results of structural genotoxicity-reactivity triage.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct GenotoxReactivityReport {
    pub findings: Vec<GenotoxFinding>,
}

const MOTIFS: &[(GenotoxMotifKind, &str, &str, &str, &str)] = &[
    (
        GenotoxMotifKind::Epoxide,
        "epoxide",
        "[C;!a]1O[C;!a]1",
        "three-membered oxygen heterocycle match",
        "compatible with electrophilic ring-opening reactivity; not a biological prediction",
    ),
    (
        GenotoxMotifKind::Aziridine,
        "aziridine",
        "[C;!a]1N[C;!a]1",
        "three-membered nitrogen heterocycle match",
        "compatible with electrophilic ring-opening reactivity; not a biological prediction",
    ),
    (
        GenotoxMotifKind::MichaelAcceptor,
        "michael_acceptor",
        "C=CC(=O)",
        "alpha,beta-unsaturated carbonyl-like match",
        "compatible with conjugate-addition reactivity; not a biological prediction",
    ),
];

type CompiledPattern = (
    GenotoxMotifKind,
    &'static str,
    QueryMolecule,
    &'static str,
    &'static str,
);

fn patterns() -> &'static [CompiledPattern] {
    static CACHE: OnceLock<Vec<CompiledPattern>> = OnceLock::new();
    CACHE.get_or_init(|| {
        MOTIFS
            .iter()
            .filter_map(|(motif, name, smarts, evidence, interpretation)| {
                parse_smarts(smarts)
                    .ok()
                    .map(|query| (*motif, *name, query, *evidence, *interpretation))
            })
            .collect()
    })
}

/// Find explainable electrophilic-reactivity motifs without assigning risk.
pub fn genotox_reactivity(mol: &Molecule) -> GenotoxReactivityReport {
    let mut findings = patterns()
        .iter()
        .flat_map(|(motif, name, query, evidence, interpretation)| {
            find_matches(query, mol).into_iter().map(|mapping| {
                let mut matched_atoms: Vec<_> = mapping.values().copied().collect();
                matched_atoms.sort_unstable();
                matched_atoms.dedup();
                GenotoxFinding {
                    motif: *motif,
                    pattern_name: name,
                    matched_atoms,
                    evidence,
                    interpretation,
                    confidence: GenotoxConfidence::PatternOnly,
                    applicability: GenotoxApplicability::OrganicSmallMolecule,
                }
            })
        })
        .collect::<Vec<_>>();
    findings.sort_by(|a, b| (&a.motif, &a.matched_atoms).cmp(&(&b.motif, &b.matched_atoms)));
    GenotoxReactivityReport { findings }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn reports_epoxide_with_explainable_boundary() {
        let report = genotox_reactivity(&parse("C1CO1").unwrap());
        assert_eq!(report.findings.len(), 1);
        assert_eq!(report.findings[0].motif, GenotoxMotifKind::Epoxide);
        assert_eq!(
            report.findings[0].matched_atoms,
            vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)]
        );
        assert_eq!(
            report.findings[0].confidence,
            GenotoxConfidence::PatternOnly
        );
        assert!(
            report.findings[0]
                .interpretation
                .contains("not a biological prediction")
        );
    }

    #[test]
    fn reports_aziridine_and_michael_acceptor() {
        assert_eq!(
            genotox_reactivity(&parse("C1CN1").unwrap()).findings[0].motif,
            GenotoxMotifKind::Aziridine
        );
        assert_eq!(
            genotox_reactivity(&parse("C=CC(=O)C").unwrap()).findings[0].motif,
            GenotoxMotifKind::MichaelAcceptor
        );
    }

    #[test]
    fn does_not_flag_unrelated_aromatic_structure() {
        assert!(
            genotox_reactivity(&parse("c1ccccc1").unwrap())
                .findings
                .is_empty()
        );
    }
}
