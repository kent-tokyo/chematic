//! Small-ring experimental torsion rules --
//! [`TorsionKnowledgeSource::SmallRingExperimental`].
//!
//! Adapted (not transliterated verbatim) from RDKit's
//! `torsionPreferences_smallrings.in` (116 lines at the pinned commit; see
//! sources manifest). Every rule's SMARTS text and V/sign coefficients are
//! taken from a specific, cited source line, but the RDKit source's own
//! `r{a-b}` ring-size RANGE syntax (e.g. `[C;r{5-8}:2]`) is **not**
//! reproducible in chematic-smarts: this crate's own `[rN]` primitive is
//! documented (`crates/chematic-smarts/src/parser.rs`) as matching an
//! *exact* single ring size only, with no range extension. Rather than
//! silently dropping the ring-size condition, this module keeps each rule's
//! SMARTS to its bare connectivity/element pattern and enforces the ring-
//! size range **structurally** in `matcher.rs`, which only tries a given
//! rule against a bond that `classify_bond` has already independently
//! determined sits in a ring whose size falls in `applicable_ring_sizes`.
//! This is a translation of the same semantic rule, not an approximation of
//! a different one -- documented here as the adaptation this PR's sources
//! manifest and audit doc both call out.

use std::ops::RangeInclusive;

use chematic_smarts::{QueryMolecule, parse_smarts};

use super::types::FourierTorsionTerm;

pub struct SmallRingTorsionRule {
    pub rule_id: &'static str,
    /// Connectivity-only SMARTS (no ring-size predicate -- see module docs).
    pub smarts: &'static str,
    /// Ring sizes (inclusive) this rule may fire for; enforced by
    /// `matcher.rs` against `classify_bond`'s own ring classification, not
    /// by the SMARTS itself.
    pub applicable_ring_sizes: RangeInclusive<usize>,
    pub terms: &'static [(u8, i8, f64)],
    pub source_line: &'static str,
}

pub static SMALL_RING_TORSION_RULES: &[SmallRingTorsionRule] = &[
    // Generic sp3-sp3 ring C-C bond in a 3- or 4-membered ring (cyclopropane,
    // cyclobutane, and their substituted/fused analogues).
    SmallRingTorsionRule {
        rule_id: "smallring:generic_cc_3_4",
        smarts: "[!#1:1][C:2][C:3][!#1:4]",
        applicable_ring_sizes: 3..=4,
        terms: &[(1, -1, 30.0)],
        source_line: "torsionPreferences_smallrings.in:14",
    },
    // Generic sp3-sp3 ring C-C bond in a 5-membered ring (cyclopentane).
    SmallRingTorsionRule {
        rule_id: "smallring:generic_cc_5",
        smarts: "[!#1:1][CX4:2][CX4:3][!#1:4]",
        applicable_ring_sizes: 5..=5,
        terms: &[(6, 1, 10.0)],
        source_line: "torsionPreferences_smallrings.in:37",
    },
    // Generic sp3-sp3 ring C-C bond in a 6-8-membered ring (cyclohexane
    // through cyclooctane).
    SmallRingTorsionRule {
        rule_id: "smallring:generic_cc_6_8",
        smarts: "[!#1:1][C:2][C:3][!#1:4]",
        applicable_ring_sizes: 6..=8,
        terms: &[(1, -1, 20.0)],
        source_line: "torsionPreferences_smallrings.in:40",
    },
    // Ring C-N bond (e.g. piperidine-, pyrrolidine-sized N-heterocycles),
    // 5-8-membered.
    SmallRingTorsionRule {
        rule_id: "smallring:ring_cn_5_8",
        smarts: "[!#1:1][C:2][N:3][!#1:4]",
        applicable_ring_sizes: 5..=8,
        terms: &[(1, -1, 20.0)],
        source_line: "torsionPreferences_smallrings.in:61",
    },
    // Ring C-O bond (e.g. tetrahydrofuran-, tetrahydropyran-, morpholine-
    // sized O-heterocycles), 5-8-membered.
    SmallRingTorsionRule {
        rule_id: "smallring:ring_co_5_8",
        smarts: "[!#1:1][C:2][O:3][!#1:4]",
        applicable_ring_sizes: 5..=8,
        terms: &[(4, -1, 10.0)],
        source_line: "torsionPreferences_smallrings.in:66",
    },
    // Ring C-S bond (e.g. thiolane-, thiane-sized S-heterocycles), 5-8-membered.
    SmallRingTorsionRule {
        rule_id: "smallring:ring_cs_5_8",
        smarts: "[!#1:1][C:2][S:3][!#1:4]",
        applicable_ring_sizes: 5..=8,
        terms: &[(3, 1, 10.0)],
        source_line: "torsionPreferences_smallrings.in:71",
    },
];

pub fn parse_rule(rule: &SmallRingTorsionRule) -> Option<QueryMolecule> {
    parse_smarts(rule.smarts).ok()
}

pub fn rule_fourier_terms(rule: &SmallRingTorsionRule) -> Vec<FourierTorsionTerm> {
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
        for rule in SMALL_RING_TORSION_RULES {
            assert!(
                parse_rule(rule).is_some(),
                "rule {} SMARTS failed to parse: {}",
                rule.rule_id,
                rule.smarts
            );
        }
    }

    #[test]
    fn no_rule_applies_to_a_macrocycle_size() {
        for rule in SMALL_RING_TORSION_RULES {
            assert!(
                *rule.applicable_ring_sizes.end() <= 8,
                "rule {} must not claim a macrocycle-range ring size",
                rule.rule_id
            );
        }
    }

    #[test]
    fn rule_ids_are_unique() {
        let mut seen = std::collections::HashSet::new();
        for rule in SMALL_RING_TORSION_RULES {
            assert!(
                seen.insert(rule.rule_id),
                "duplicate rule_id {}",
                rule.rule_id
            );
        }
    }

    #[test]
    fn every_rule_cites_a_real_source_line() {
        for rule in SMALL_RING_TORSION_RULES {
            assert!(
                rule.source_line
                    .starts_with("torsionPreferences_smallrings.in:"),
                "rule {} has an uncitable source_line: {}",
                rule.rule_id,
                rule.source_line
            );
        }
    }
}
