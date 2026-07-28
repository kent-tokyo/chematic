//! Macrocycle-specific experimental torsion rules --
//! [`TorsionKnowledgeSource::MacrocycleAdaptation`].
//!
//! Adapted from RDKit's `torsionPreferences_macrocycles.in` (379 lines at
//! the pinned commit; see sources manifest), same `r{9-}` range-syntax
//! adaptation caveat as `rules_smallring.rs` (chematic-smarts has no ring-
//! size-range primitive; the >=9 condition is enforced in `matcher.rs`
//! against `classify_bond`'s own classification, not in the SMARTS text).

use std::ops::RangeInclusive;

use chematic_smarts::{QueryMolecule, parse_smarts};

use super::classify::MACROCYCLE_MIN;
use super::types::FourierTorsionTerm;

pub struct MacrocycleTorsionRule {
    pub rule_id: &'static str,
    pub smarts: &'static str,
    pub applicable_ring_sizes: RangeInclusive<usize>,
    pub terms: &'static [(u8, i8, f64)],
    pub source_line: &'static str,
}

pub static MACROCYCLE_TORSION_RULES: &[MacrocycleTorsionRule] = &[
    // Lactam (ring amide) bond: the C(=O)-N central bond of an amide whose
    // carbonyl carbon sits in a >=9-membered ring (a macrolactam).
    MacrocycleTorsionRule {
        rule_id: "macrocycle:lactam_amide",
        smarts: "[C:1][C:2](=O)[NX3H0:3][CX4H1:4]",
        applicable_ring_sizes: MACROCYCLE_MIN..=usize::MAX,
        terms: &[(2, -1, 8.0)],
        source_line: "torsionPreferences_macrocycles.in:12",
    },
    // Macrolactone ester bond: O=C-O-C in a >=9-membered ring.
    MacrocycleTorsionRule {
        rule_id: "macrocycle:lactone_ester",
        smarts: "[O:1]=[C:2][O:3][CH0:4]",
        applicable_ring_sizes: MACROCYCLE_MIN..=usize::MAX,
        terms: &[(1, -1, 78.2)],
        source_line: "torsionPreferences_macrocycles.in:16",
    },
    // Acetal-like O-C-O-C chain within a macrocycle ring (e.g. crown-ether-
    // adjacent contexts).
    MacrocycleTorsionRule {
        rule_id: "macrocycle:ring_o_c_o_chain",
        smarts: "[O:1][CX4:2][O:3][CX4:4]",
        applicable_ring_sizes: MACROCYCLE_MIN..=usize::MAX,
        terms: &[(3, 1, 8.8)],
        source_line: "torsionPreferences_macrocycles.in:27",
    },
    // Aromatic-N-adjacent ring CH2: real, genuinely multi-term data is not
    // present at this exact line (single term here); kept as its own rule
    // since the macrocyclic CH2-n(aromatic) environment is chemically
    // distinct from the macrocyclic sp3-amine chain rule below.
    MacrocycleTorsionRule {
        rule_id: "macrocycle:ring_ch2_aromatic_n",
        smarts: "[!#1:1][CH2:2][n:3][cH0:4]",
        applicable_ring_sizes: MACROCYCLE_MIN..=usize::MAX,
        terms: &[(2, 1, 6.0)],
        source_line: "torsionPreferences_macrocycles.in:173",
    },
    // Macrocyclic amine chain C-C-N-C: a genuinely multi-term (real, not
    // invented) potential -- the macrocycle analogue of standard.rs's
    // biphenyl multi-term evidence for why a single-angle model is wrong.
    MacrocycleTorsionRule {
        rule_id: "macrocycle:ring_amine_chain",
        smarts: "[CX4:1][CX4H2:2][NX3:3][CX4:4]",
        applicable_ring_sizes: MACROCYCLE_MIN..=usize::MAX,
        terms: &[
            (1, 1, 4.0),
            (2, 1, 3.1),
            (3, 1, 3.9),
            (4, -1, 0.8),
            (6, 1, 0.7),
        ],
        source_line: "torsionPreferences_macrocycles.in:176",
    },
];

pub fn parse_rule(rule: &MacrocycleTorsionRule) -> Option<QueryMolecule> {
    parse_smarts(rule.smarts).ok()
}

pub fn rule_fourier_terms(rule: &MacrocycleTorsionRule) -> Vec<FourierTorsionTerm> {
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
        for rule in MACROCYCLE_TORSION_RULES {
            assert!(
                parse_rule(rule).is_some(),
                "rule {} SMARTS failed to parse: {}",
                rule.rule_id,
                rule.smarts
            );
        }
    }

    #[test]
    fn no_rule_applies_below_macrocycle_min() {
        for rule in MACROCYCLE_TORSION_RULES {
            assert!(
                *rule.applicable_ring_sizes.start() >= MACROCYCLE_MIN,
                "rule {} must not claim a small-ring-range ring size",
                rule.rule_id
            );
        }
    }

    #[test]
    fn amine_chain_rule_is_genuinely_multi_term() {
        let rule = MACROCYCLE_TORSION_RULES
            .iter()
            .find(|r| r.rule_id == "macrocycle:ring_amine_chain")
            .unwrap();
        assert!(rule_fourier_terms(rule).len() > 1);
    }

    #[test]
    fn rule_ids_are_unique() {
        let mut seen = std::collections::HashSet::new();
        for rule in MACROCYCLE_TORSION_RULES {
            assert!(
                seen.insert(rule.rule_id),
                "duplicate rule_id {}",
                rule.rule_id
            );
        }
    }

    #[test]
    fn every_rule_cites_a_real_source_line() {
        for rule in MACROCYCLE_TORSION_RULES {
            assert!(
                rule.source_line
                    .starts_with("torsionPreferences_macrocycles.in:"),
                "rule {} has an uncitable source_line: {}",
                rule.rule_id,
                rule.source_line
            );
        }
    }
}
