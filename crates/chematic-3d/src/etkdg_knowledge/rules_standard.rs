//! Standard (acyclic) experimental torsion rules --
//! [`TorsionKnowledgeSource::StandardExperimental`].
//!
//! A **curated subset**, not a wholesale reproduction, of RDKit's
//! `torsionPreferences_v1.in`/`torsionPreferences_v2.in` (592 and 378 lines
//! respectively at the pinned commit -- see the sources manifest). Every
//! rule below cites its exact source file and 1-based line number; SMARTS
//! text and Fourier `(periodicity, sign, V)` coefficients are copied
//! verbatim from the fetched file, only the `!@` (non-ring) bond qualifier
//! is dropped from the SMARTS text itself (chematic-smarts has no trouble
//! parsing `!@`, but ring-exclusion for this tier is instead enforced
//! structurally by `matcher.rs`, which only applies `StandardExperimental`
//! to bonds `classify_bond` reports as `acyclic_rotatable` -- functionally
//! equivalent, more testable in one place).
//!
//! Copyright note: RDKit is BSD-3-Clause
//! ("Copyright (C) 2017-2023 Sereina Riniker and other RDKit contributors",
//! per `TorsionPreferences.cpp`'s own header, fetched and hashed in the
//! sources manifest). This file preserves that attribution here and in the
//! per-rule `source_line` field rather than redistributing RDKit source
//! files themselves.

use chematic_smarts::{QueryMolecule, parse_smarts};

use super::types::FourierTorsionTerm;

/// One standard-experimental torsion rule: a 4-atom-mapped SMARTS (`:1:2:3:4`
/// mark the A-B-C-D dihedral atoms, B-C is the central bond) plus its
/// Fourier terms.
pub struct StandardTorsionRule {
    pub rule_id: &'static str,
    pub smarts: &'static str,
    /// `(periodicity, sign, amplitude)` triples, RDKit's own encoding (see
    /// `FourierTorsionTerm::from_rdkit`). Zero-amplitude terms are omitted
    /// (RDKit's data always lists all 6, most `0.0` -- keeping only the
    /// nonzero ones is a lossless simplification, not a data change).
    pub terms: &'static [(u8, i8, f64)],
    /// Exact source citation, e.g. `"torsionPreferences_v2.in:142"`.
    pub source_line: &'static str,
}

/// Curated standard-experimental rules. See module docs for the
/// non-wholesale-reproduction rationale. Each `source_line` is independently
/// checkable against the fetched file recorded in
/// `validation/manifests/etkdg_torsion_knowledge_sources.json`.
pub static STANDARD_TORSION_RULES: &[StandardTorsionRule] = &[
    // Secondary amide O=C-N(H)-R': the real, dominant amide-bond rotational
    // preference (this is why N-methylacetamide is ~99% trans). A single
    // n=1 term with V=100 (very stiff) reproduces that near-total
    // trans/cis population skew from real data, not an invented number.
    StandardTorsionRule {
        rule_id: "standard:secondary_amide",
        smarts: "[O:1]=[CX3:2][NX3H1:3][!#1:4]",
        terms: &[(1, -1, 100.0)],
        source_line: "torsionPreferences_v2.in:142 (identical in torsionPreferences_v1.in:215)",
    },
    // Tertiary (N,N-disubstituted) amide O=C-N(R)(R'), generic non-aromatic
    // acyl case.
    StandardTorsionRule {
        rule_id: "standard:tertiary_amide",
        smarts: "[O:1]=[CX3:2][NX3H0:3][!#1:4]",
        terms: &[(2, -1, 8.0)],
        source_line: "torsionPreferences_v1.in:214",
    },
    // Tertiary amide on an aromatic acyl carbon (Ar-C(=O)-NR2), stronger
    // restriction than the generic case above.
    StandardTorsionRule {
        rule_id: "standard:aromatic_tertiary_amide",
        smarts: "[O:1]=[CX3:2](a)[NX3H0:3][!#1:4]",
        terms: &[(2, -1, 13.9)],
        source_line: "torsionPreferences_v1.in:212",
    },
    // Aryl ester Ar-O-C(=O): the O-Caryl bond.
    StandardTorsionRule {
        rule_id: "standard:aryl_ester",
        smarts: "[$(C=O):1][O:2][c:3][*:4]",
        terms: &[(2, 1, 0.8)],
        source_line: "torsionPreferences_v2.in:15",
    },
    // Unsubstituted biphenyl (all-cH1 ortho positions): the REAL multi-term,
    // two-minima potential the legacy heuristic's single-45°-angle rule
    // could not represent (see docs/rfcs/3d_torsion_knowledge_audit.md's L202
    // entry for the concrete contrast).
    StandardTorsionRule {
        rule_id: "standard:biphenyl_unsubstituted",
        smarts: "[cH1:1][c:2]([cH1])[c:3]([cH1:4])[cH1]",
        terms: &[(1, -1, 0.7), (2, 1, 8.0), (4, 1, 4.4), (6, 1, 1.5)],
        source_line: "torsionPreferences_v2.in:270",
    },
    // ortho,ortho'-disubstituted ("hindered") biphenyl: both flanking
    // positions on both rings are cH0 (fully substituted). A single n=2
    // term whose minimum sits at 90/270 degrees -- the real RDKit-sourced
    // analogue of the legacy `SMARTS_TORSION_RULES` "hindered biaryl -> 90°"
    // guess, now with an actual citable source.
    StandardTorsionRule {
        rule_id: "standard:biphenyl_hindered_22prime",
        smarts: "[cH0:1][c:2]([cH0])[c:3]([cH0:4])[cH0]",
        terms: &[(2, 1, 3.6)],
        source_line: "torsionPreferences_v2.in:265",
    },
    // Secondary-amide alpha carbon, methyl-flanked: N(H)-CH2-CH3 context.
    StandardTorsionRule {
        rule_id: "standard:secondary_amide_alpha_methyl",
        smarts: "[$([CX3]=O):1][NX3H1:2][CX4H2:3][C:4]",
        terms: &[(3, 1, 4.0)],
        source_line: "torsionPreferences_v2.in:64",
    },
    // Sulfonamide-flanked alpha carbon: a genuinely multi-term (real, not
    // invented) potential from real data.
    StandardTorsionRule {
        rule_id: "standard:sulfonamide_alpha_carbon",
        smarts: "[$(S(=O)(=O)):1][NX3H1:2][CX4H2:3][!#1:4]",
        terms: &[
            (1, 1, 17.9),
            (2, -1, 13.3),
            (3, 1, 9.2),
            (4, 1, 4.7),
            (5, -1, 2.3),
            (6, -1, 0.9),
        ],
        source_line: "torsionPreferences_v2.in:66",
    },
    // Secondary-amide alpha carbon, generic (any non-H substituent).
    StandardTorsionRule {
        rule_id: "standard:secondary_amide_alpha_generic",
        smarts: "[$(C=O):1][NX3H1:2][CX4H2:3][!#1:4]",
        terms: &[(3, 1, 2.0)],
        source_line: "torsionPreferences_v2.in:108",
    },
];

/// Parse `rule.smarts` and convert `rule.terms` into
/// [`FourierTorsionTerm`]s. Returns `None` (never panics) on a SMARTS parse
/// failure so the caller (`matcher.rs`) can record a typed diagnostic
/// instead of the legacy code's silent `continue` (audit doc §3.7).
pub fn parse_rule(rule: &StandardTorsionRule) -> Option<QueryMolecule> {
    parse_smarts(rule.smarts).ok()
}

pub fn rule_fourier_terms(rule: &StandardTorsionRule) -> Vec<FourierTorsionTerm> {
    rule.terms
        .iter()
        .map(|&(n, s, v)| FourierTorsionTerm::from_rdkit(n, s, v))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_rule_smarts_parses() {
        for rule in STANDARD_TORSION_RULES {
            assert!(
                parse_rule(rule).is_some(),
                "rule {} SMARTS failed to parse: {}",
                rule.rule_id,
                rule.smarts
            );
        }
    }

    #[test]
    fn every_rule_has_atom_map_1_through_4() {
        for rule in STANDARD_TORSION_RULES {
            let q = parse_rule(rule).unwrap();
            for want in 1u16..=4 {
                assert!(
                    q.atoms.iter().any(|a| a.atom_map == Some(want)),
                    "rule {} is missing atom map :{want}",
                    rule.rule_id
                );
            }
        }
    }

    #[test]
    fn biphenyl_unsubstituted_is_genuinely_multi_term() {
        let rule = STANDARD_TORSION_RULES
            .iter()
            .find(|r| r.rule_id == "standard:biphenyl_unsubstituted")
            .unwrap();
        let terms = rule_fourier_terms(rule);
        assert!(
            terms.len() > 1,
            "biphenyl potential must be multi-modal, not a single term"
        );
    }

    #[test]
    fn rule_ids_are_unique() {
        let mut seen = std::collections::HashSet::new();
        for rule in STANDARD_TORSION_RULES {
            assert!(
                seen.insert(rule.rule_id),
                "duplicate rule_id {}",
                rule.rule_id
            );
        }
    }

    #[test]
    fn every_rule_cites_a_real_source_file_and_line() {
        // Every rule_id must be traceable to a specific line in one of the
        // RDKit files hashed in the sources manifest -- never an uncited
        // number (spec §1/§2). Checked here (not just eyeballed) so a
        // future edit that drops a citation fails CI, not just review.
        for rule in STANDARD_TORSION_RULES {
            assert!(
                rule.source_line.starts_with("torsionPreferences_v1.in:")
                    || rule.source_line.starts_with("torsionPreferences_v2.in:"),
                "rule {} has an uncitable source_line: {}",
                rule.rule_id,
                rule.source_line
            );
        }
    }
}
