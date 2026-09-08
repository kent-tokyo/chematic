//! Explainable structural reactivity findings related to genotoxicity triage.
//!
//! This module deliberately reports motifs and evidence only.  It is not a
//! genotoxicity predictor, score, classifier, or biological risk assessment.

use std::collections::VecDeque;
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
    BifunctionalElectrophile,
}

/// One explainable structural finding.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenotoxFinding {
    pub motif: GenotoxMotifKind,
    pub pattern_name: &'static str,
    /// Target atom indices in deterministic ascending order.
    pub matched_atoms: Vec<AtomIdx>,
    /// Independent reactive-site groups used by a bifunctional finding.
    /// Empty for single-motif findings.
    pub site_groups: Vec<Vec<AtomIdx>>,
    /// Minimum heavy-atom bond distance between the two site groups. This is
    /// a topological spacer estimate, not a 3D geometry or biological claim.
    pub through_bond_distance: Option<u32>,
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

const MAX_BIFUNCTIONAL_SITE_PAIRS: usize = 64;

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
                    site_groups: Vec::new(),
                    through_bond_distance: None,
                    evidence,
                    interpretation,
                    confidence: GenotoxConfidence::PatternOnly,
                    applicability: GenotoxApplicability::OrganicSmallMolecule,
                }
            })
        })
        .collect::<Vec<_>>();
    let mut sites: Vec<(GenotoxMotifKind, Vec<AtomIdx>)> = findings
        .iter()
        .filter(|finding| finding.motif != GenotoxMotifKind::BifunctionalElectrophile)
        .map(|finding| (finding.motif, finding.matched_atoms.clone()))
        .collect();
    sites.sort_by(|(left_motif, left_atoms), (right_motif, right_atoms)| {
        (left_motif, left_atoms).cmp(&(right_motif, right_atoms))
    });
    sites.dedup();
    let mut pair_count = 0;
    'pairs: for left in 0..sites.len() {
        for right in (left + 1)..sites.len() {
            if sites[left]
                .1
                .iter()
                .any(|atom| sites[right].1.contains(atom))
            {
                continue;
            }
            let Some(distance) = site_distance(mol, &sites[left].1, &sites[right].1) else {
                continue;
            };
            let mut matched_atoms = sites[left].1.clone();
            matched_atoms.extend_from_slice(&sites[right].1);
            matched_atoms.sort_unstable();
            matched_atoms.dedup();
            findings.push(GenotoxFinding {
                motif: GenotoxMotifKind::BifunctionalElectrophile,
                pattern_name: "bifunctional_electrophile",
                matched_atoms,
                site_groups: vec![sites[left].1.clone(), sites[right].1.clone()],
                through_bond_distance: Some(distance),
                evidence: "two independent electrophilic motif matches with a bounded through-bond spacer estimate",
                interpretation: "structurally compatible with two-site electrophilic reactivity; not a biological or 3D-geometry prediction",
                confidence: GenotoxConfidence::PatternOnly,
                applicability: GenotoxApplicability::OrganicSmallMolecule,
            });
            pair_count += 1;
            if pair_count == MAX_BIFUNCTIONAL_SITE_PAIRS {
                break 'pairs;
            }
        }
    }
    findings.sort_by(|a, b| (&a.motif, &a.matched_atoms).cmp(&(&b.motif, &b.matched_atoms)));
    GenotoxReactivityReport { findings }
}

fn site_distance(mol: &Molecule, left: &[AtomIdx], right: &[AtomIdx]) -> Option<u32> {
    let mut distances = vec![u32::MAX; mol.atom_count()];
    let mut queue = VecDeque::new();
    for &atom in left {
        distances[atom.0 as usize] = 0;
        queue.push_back(atom);
    }
    while let Some(atom) = queue.pop_front() {
        let next_distance = distances[atom.0 as usize].saturating_add(1);
        if right.contains(&atom) {
            return Some(distances[atom.0 as usize]);
        }
        for (neighbor, _) in mol.neighbors(atom) {
            let slot = &mut distances[neighbor.0 as usize];
            if *slot == u32::MAX {
                *slot = next_distance;
                queue.push_back(neighbor);
            }
        }
    }
    None
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
        assert!(report.findings[0].site_groups.is_empty());
        assert_eq!(report.findings[0].through_bond_distance, None);
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

    #[test]
    fn reports_bifunctional_sites_with_bounded_topological_spacer() {
        let report = genotox_reactivity(&parse("C1CO1CC2CO2").unwrap());
        let finding = report
            .findings
            .iter()
            .find(|finding| finding.motif == GenotoxMotifKind::BifunctionalElectrophile)
            .expect("two epoxides should produce a pair finding");
        assert_eq!(finding.site_groups.len(), 2);
        assert_eq!(finding.through_bond_distance, Some(2));
        assert!(
            finding
                .interpretation
                .contains("not a biological or 3D-geometry prediction")
        );
    }

    #[test]
    fn single_site_does_not_produce_bifunctional_finding() {
        assert!(
            !genotox_reactivity(&parse("C1CO1").unwrap())
                .findings
                .iter()
                .any(|finding| finding.motif == GenotoxMotifKind::BifunctionalElectrophile)
        );
    }

    #[test]
    fn source_referenced_structures_cover_external_smoke_fixtures() {
        let manifest: serde_json::Value = serde_json::from_str(include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../validation/genotox_structural_fixtures.json"
        )))
        .expect("genotoxicity fixture manifest should be valid JSON");
        for fixture in manifest["fixtures"].as_array().expect("fixture list") {
            let molecule = parse(fixture["isomeric_smiles"].as_str().expect("SMILES"))
                .expect("source fixture SMILES should parse");
            let report = genotox_reactivity(&molecule);
            for expected in fixture["expected_motifs"]
                .as_array()
                .expect("expected motifs")
            {
                let name = expected.as_str().expect("motif name");
                assert!(
                    report
                        .findings
                        .iter()
                        .any(|finding| finding.pattern_name == name),
                    "{} should report {name}",
                    fixture["id"]
                );
            }
        }
    }
}
