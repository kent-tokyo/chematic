//! Priority-tiered rule application (Wave 2 spec §4).
//!
//! Explicit priority order, highest first:
//! 1. specific validated SMARTS (currently empty -- see PR body's judgment-
//!    call section: no rule in this PR clears the bar of "individually,
//!    specifically validated beyond membership in the bulk statistical
//!    library" that would distinguish it from tier 4; the tier is wired and
//!    tested, just unpopulated, rather than populated with something
//!    dishonestly labeled)
//! 2. small-ring torsion
//! 3. macrocycle-specific torsion
//! 4. standard experimental torsion
//! 5. basic chemical knowledge
//! 6. legacy heuristic (opt-in only)
//!
//! This ordering is a **chematic-specific design decision**, not a
//! reproduction of RDKit's own iteration order: RDKit's real
//! `getExperimentalTorsions` loop (`TorsionPreferences.cpp`, fetched and
//! hashed in the sources manifest) is a single first-match-not-yet-done pass
//! over the concatenated standard+smallring+macrocycle parameter string, in
//! that concatenation order -- meaning **standard** rules actually take
//! priority over small-ring/macrocycle ones in real RDKit, the reverse of
//! this file's tier order. The Wave 2 spec explicitly mandates small-ring
//! and macrocycle above standard, so this file follows the spec, not the
//! RDKit source -- documented here so nobody mistakes this for a bug
//! relative to the oracle (spec §12 only asks this PR's rule-family/ring-
//! classification/minima *differential* to be checked against RDKit, not
//! the internal iteration order).
//!
//! For bonds skipped as non-torsional (terminal / double-triple / aromatic-
//! ring), no tier is tried at all -- recorded once in `skipped_bonds`, never
//! silently absent from the report.
//!
//! Within one tier, multiple rules matching the same bond are: deduplicated
//! if their term lists are identical; composed (merged) if their term lists
//! occupy disjoint Fourier periodicities; and reported as a typed
//! [`TorsionKnowledgeDiagnostic`] (never silently resolved) if they disagree
//! on the same periodicity.

use std::collections::HashMap;

use chematic_core::{AtomIdx, Molecule};
use chematic_smarts::{QueryMolecule, find_matches};

use super::classify::{
    BondClassification, RingMembershipIndex, candidate_central_bonds, classify_bond,
};
use super::rules_basic;
use super::rules_macrocycle;
use super::rules_smallring;
use super::rules_standard;
use super::types::{
    FourierTorsionTerm, TorsionKnowledgeConfig, TorsionKnowledgeDiagnostic,
    TorsionKnowledgeDiagnosticKind, TorsionKnowledgeReport, TorsionKnowledgeSource,
    TorsionPotential,
};
use crate::etkdg_knowledge::{build_smarts_torsion_map, get_torsion_preference};
use chematic_perception::find_sssr;

type BondKey = (u32, u32);

fn bond_key(a: AtomIdx, b: AtomIdx) -> BondKey {
    (a.0.min(b.0), a.0.max(b.0))
}

/// One candidate match within a tier, before same-tier conflict resolution.
struct Candidate {
    rule_id: String,
    /// The concrete A-B-C-D atoms this specific match found (B-C is the
    /// central bond). Needed so the resulting [`TorsionPotential`] can
    /// actually compute a real dihedral angle from coordinates later --
    /// carrying only the central bond would lose the substituent choice.
    atoms: [AtomIdx; 4],
    terms: Vec<FourierTorsionTerm>,
    source: TorsionKnowledgeSource,
    ring_size: Option<usize>,
}

/// Extract the atoms mapped `:1`, `:2`, `:3`, `:4` from one VF2 match,
/// generic over the match map's hasher (chematic-smarts's `find_matches`
/// returns an `FxHashMap`, but this function only needs `.get`, so it stays
/// decoupled from that specific hasher type).
fn atoms_by_map4<S: std::hash::BuildHasher>(
    query: &QueryMolecule,
    m: &HashMap<usize, AtomIdx, S>,
) -> Option<[AtomIdx; 4]> {
    let mut out: [Option<AtomIdx>; 4] = [None; 4];
    for (qi, qatom) in query.atoms.iter().enumerate() {
        if let Some(map) = qatom.atom_map
            && (1..=4).contains(&map)
            && let Some(&aidx) = m.get(&qi)
        {
            out[(map - 1) as usize] = Some(aidx);
        }
    }
    Some([out[0]?, out[1]?, out[2]?, out[3]?])
}

fn effective_ring_size(bc: &BondClassification) -> Option<usize> {
    use super::classify::RingMembership;
    match bc.ring {
        RingMembership::NotInRing => None,
        RingMembership::SmallRing(n) | RingMembership::Macrocycle(n) => Some(n),
        RingMembership::FusedOrBridged { chosen_size } => Some(chosen_size),
    }
}

/// Build the full torsion-knowledge report for `mol` under `config`. Pure
/// function: reads `mol`, never mutates it; never touches coordinates.
///
/// All-flags-false is a genuine no-op: returns
/// [`TorsionKnowledgeReport::default()`] (empty in every field) without
/// even enumerating candidate bonds -- spec §11/§15's "true no-op" for the
/// disabled-flags negative control.
pub fn build_torsion_knowledge(
    mol: &Molecule,
    config: &TorsionKnowledgeConfig,
) -> TorsionKnowledgeReport {
    let mut report = TorsionKnowledgeReport::default();

    if !config.use_exp_torsions
        && !config.use_small_ring_torsions
        && !config.use_macrocycle_torsions
        && !config.include_legacy_heuristic
    {
        // `use_macrocycle_14_bounds` is deliberately excluded from this
        // check: it gates a wholly separate API (`bounds14.rs`), not this
        // torsion-potential report.
        return report;
    }

    let ring_index = RingMembershipIndex::build(mol);
    let central_bonds = candidate_central_bonds(mol);

    let mut classifications: HashMap<BondKey, BondClassification> = HashMap::new();
    let mut rotatable_bonds: Vec<(AtomIdx, AtomIdx)> = Vec::new();

    for (a, b) in &central_bonds {
        let bc = classify_bond(mol, &ring_index, *a, *b);
        if bc.terminal || bc.double_or_triple || bc.aromatic_ring_bond {
            report.skipped_bonds.push(TorsionKnowledgeDiagnostic {
                central_bond: (*a, *b),
                kind: TorsionKnowledgeDiagnosticKind::NonRotatableBondSkipped,
                message: bc.reasoning.clone(),
                candidate_rule_ids: vec![],
            });
            continue;
        }
        if bc.ring.is_fused_or_bridged() {
            report.ambiguous_matches.push(TorsionKnowledgeDiagnostic {
                central_bond: (*a, *b),
                kind: TorsionKnowledgeDiagnosticKind::FusedOrBridgedRingBoundary,
                message: bc.reasoning.clone(),
                candidate_rule_ids: vec![],
            });
            // Not a fatal skip -- fused/bridged bonds are still eligible for
            // rule matching using the chosen (smallest-ring) bucket; this
            // diagnostic exists purely so the boundary decision is visible,
            // per spec §5.
        }
        rotatable_bonds.push((*a, *b));
        classifications.insert(bond_key(*a, *b), bc);
    }

    // Tier 2: small-ring.
    let mut tier2: HashMap<BondKey, Vec<Candidate>> = HashMap::new();
    if config.use_small_ring_torsions {
        for rule in rules_smallring::SMALL_RING_TORSION_RULES {
            let Some(query) = rules_smallring::parse_rule(rule) else {
                report.skipped_bonds.push(parse_failure_diagnostic(
                    rule.rule_id,
                    rule.smarts,
                    rule.source_line,
                ));
                continue;
            };
            for m in find_matches(&query, mol) {
                let Some(atoms) = atoms_by_map4(&query, &m) else {
                    continue;
                };
                let (b, c) = (atoms[1], atoms[2]);
                if mol.bond_between(b, c).is_none() {
                    continue;
                }
                let key = bond_key(b, c);
                let Some(bc) = classifications.get(&key) else {
                    continue;
                };
                let Some(size) = effective_ring_size(bc) else {
                    continue;
                };
                if !rule.applicable_ring_sizes.contains(&size) {
                    continue;
                }
                tier2.entry(key).or_default().push(Candidate {
                    rule_id: rule.rule_id.to_string(),
                    atoms,
                    terms: rules_smallring::rule_fourier_terms(rule),
                    source: TorsionKnowledgeSource::SmallRingExperimental,
                    ring_size: Some(size),
                });
            }
        }
    }

    // Tier 3: macrocycle.
    let mut tier3: HashMap<BondKey, Vec<Candidate>> = HashMap::new();
    if config.use_macrocycle_torsions {
        for rule in rules_macrocycle::MACROCYCLE_TORSION_RULES {
            let Some(query) = rules_macrocycle::parse_rule(rule) else {
                report.skipped_bonds.push(parse_failure_diagnostic(
                    rule.rule_id,
                    rule.smarts,
                    rule.source_line,
                ));
                continue;
            };
            for m in find_matches(&query, mol) {
                let Some(atoms) = atoms_by_map4(&query, &m) else {
                    continue;
                };
                let (b, c) = (atoms[1], atoms[2]);
                if mol.bond_between(b, c).is_none() {
                    continue;
                }
                let key = bond_key(b, c);
                let Some(bc) = classifications.get(&key) else {
                    continue;
                };
                let Some(size) = effective_ring_size(bc) else {
                    continue;
                };
                if !rule.applicable_ring_sizes.contains(&size) {
                    continue;
                }
                tier3.entry(key).or_default().push(Candidate {
                    rule_id: rule.rule_id.to_string(),
                    atoms,
                    terms: rules_macrocycle::rule_fourier_terms(rule),
                    source: TorsionKnowledgeSource::MacrocycleAdaptation,
                    ring_size: Some(size),
                });
            }
        }
    }

    // Tier 4: standard experimental (acyclic bonds only).
    let mut tier4: HashMap<BondKey, Vec<Candidate>> = HashMap::new();
    if config.use_exp_torsions {
        for rule in rules_standard::STANDARD_TORSION_RULES {
            let Some(query) = rules_standard::parse_rule(rule) else {
                report.skipped_bonds.push(parse_failure_diagnostic(
                    rule.rule_id,
                    rule.smarts,
                    rule.source_line,
                ));
                continue;
            };
            for m in find_matches(&query, mol) {
                let Some(atoms) = atoms_by_map4(&query, &m) else {
                    continue;
                };
                let (b, c) = (atoms[1], atoms[2]);
                if mol.bond_between(b, c).is_none() {
                    continue;
                }
                let key = bond_key(b, c);
                let Some(bc) = classifications.get(&key) else {
                    continue;
                };
                if bc.ring.is_in_any_ring() {
                    continue; // standard tier is acyclic-only, mirrors RDKit's `!@`
                }
                tier4.entry(key).or_default().push(Candidate {
                    rule_id: rule.rule_id.to_string(),
                    atoms,
                    terms: rules_standard::rule_fourier_terms(rule),
                    source: TorsionKnowledgeSource::StandardExperimental,
                    ring_size: None,
                });
            }
        }
    }

    // Tier 5: basic chemical knowledge (flat 4-6-membered all-sp2 rings;
    // linear sp centers). Built directly against the rotatable bonds rather
    // than a SMARTS table, since both sub-rules are structural.
    let mut tier5: HashMap<BondKey, Vec<Candidate>> = HashMap::new();
    if config.use_exp_torsions {
        for &(a, b) in &rotatable_bonds {
            let key = bond_key(a, b);
            if let Some((rule_id, atoms, term)) =
                basic_knowledge_term_for_bond(mol, &ring_index, a, b)
            {
                tier5.entry(key).or_default().push(Candidate {
                    rule_id,
                    atoms,
                    terms: vec![term],
                    source: TorsionKnowledgeSource::BasicChemicalKnowledge,
                    ring_size: None,
                });
            }
        }
    }

    // Tier 6: legacy heuristic (opt-in only). Wraps the pre-existing,
    // behaviorally-unchanged `get_torsion_preference`/`build_smarts_torsion_map`
    // output as a single-term potential -- an APPROXIMATION, not an exact
    // translation: the legacy model is a single angle + LINEAR penalty,
    // which has no exact periodic-cosine equivalent. This maps the
    // preferred angle to a cosine minimum at the same angle (n=1) and scales
    // `amplitude = penalty_per_degree * 50.0` as a documented, arbitrary-but-
    // stated proportionality constant -- never claimed to reproduce the
    // legacy scoring exactly, only to make its *preferred angle* visible
    // through the same report type when a caller explicitly opts in. This
    // is a one-directional wrapper (legacy -> v2 report); it does not change
    // `get_torsion_preference`'s own behavior or make it delegate to v2 (see
    // spec §8), so no new regression fixtures are needed for the legacy API
    // itself.
    let mut tier6: HashMap<BondKey, Vec<Candidate>> = HashMap::new();
    if config.include_legacy_heuristic {
        let ring_bond_set: std::collections::HashSet<(u32, u32)> = find_sssr(mol)
            .rings()
            .iter()
            .flat_map(|r| {
                let n = r.len();
                (0..n).flat_map(move |i| {
                    let x = r[i].0;
                    let y = r[(i + 1) % n].0;
                    [(x, y), (y, x)]
                })
            })
            .collect();
        let smarts_map = build_smarts_torsion_map(mol, &ring_bond_set);
        for &(a, b) in &rotatable_bonds {
            let key = bond_key(a, b);
            // A real A-D pair is needed regardless of whether the SMARTS map
            // or the atom-type cascade ultimately supplies the preference,
            // since the resulting `TorsionPotential.atoms` must be a genuine
            // 4-distinct-atom quadruple (unlike some of the legacy test
            // suite's own duplicate-index calls -- see audit doc §4).
            let d_candidate = mol.neighbors(b).map(|(n, _)| n).find(|&n| n != a);
            let a_candidate = mol.neighbors(a).map(|(n, _)| n).find(|&n| n != b);
            let (Some(a_atom), Some(d_atom)) = (a_candidate, d_candidate) else {
                continue;
            };
            let legacy_pref = smarts_map
                .get(&(a.0, b.0))
                .cloned()
                .or_else(|| get_torsion_preference(mol, a_atom, a, b, d_atom));
            if let Some(pref) = legacy_pref {
                let phase_deg = (pref.angle_deg - 180.0).rem_euclid(360.0);
                let amplitude = pref.penalty_per_degree * 50.0;
                tier6.entry(key).or_default().push(Candidate {
                    rule_id: "legacy:get_torsion_preference".to_string(),
                    atoms: [a_atom, a, b, d_atom],
                    terms: vec![FourierTorsionTerm::new(1, phase_deg, amplitude)],
                    source: TorsionKnowledgeSource::LegacyHeuristic,
                    ring_size: None,
                });
            }
        }
    }

    // Resolve, per bond, the first tier (in priority order) with any
    // candidate, applying same-tier dedup/compose/conflict resolution.
    for &(a, b) in &rotatable_bonds {
        let key = bond_key(a, b);
        let tiers: [&HashMap<BondKey, Vec<Candidate>>; 5] =
            [&tier2, &tier3, &tier4, &tier5, &tier6];
        let mut resolved = false;
        for tier_map in tiers {
            let Some(candidates) = tier_map.get(&key) else {
                continue;
            };
            if candidates.is_empty() {
                continue;
            }
            match resolve_same_tier(candidates) {
                ResolvedTier::Single {
                    rule_ids,
                    atoms,
                    terms,
                    source,
                    ring_size,
                } => {
                    report.matched_rule_ids.extend(rule_ids.iter().cloned());
                    report.potentials.push(TorsionPotential {
                        atoms,
                        central_bond: (a, b),
                        source,
                        rule_id: rule_ids.join("+"),
                        terms,
                        ring_size,
                    });
                    resolved = true;
                }
                ResolvedTier::Conflict { rule_ids } => {
                    report.ambiguous_matches.push(TorsionKnowledgeDiagnostic {
                        central_bond: (a, b),
                        kind: TorsionKnowledgeDiagnosticKind::AmbiguousSameTierConflict,
                        message: format!(
                            "{} rules matched bond (atom{},atom{}) at the same priority tier with conflicting terms",
                            rule_ids.len(), a.0, b.0
                        ),
                        candidate_rule_ids: rule_ids,
                    });
                    resolved = true; // conflict is a resolution (a diagnosed one), do not fall through to a lower tier
                }
            }
            break;
        }
        if !resolved {
            report.unmatched_rotatable_bonds.push((a, b));
        }
    }

    report
}

fn parse_failure_diagnostic(
    rule_id: &str,
    smarts: &str,
    source_line: &str,
) -> TorsionKnowledgeDiagnostic {
    TorsionKnowledgeDiagnostic {
        // Rule-level (not bond-specific) failure: sentinel atom pair, noted
        // in the message so it reads as global, not a real bond.
        central_bond: (AtomIdx(u32::MAX), AtomIdx(u32::MAX)),
        kind: TorsionKnowledgeDiagnosticKind::SmartsParseFailure,
        message: format!(
            "rule {rule_id} (source: {source_line}): SMARTS failed to parse: {smarts}"
        ),
        candidate_rule_ids: vec![rule_id.to_string()],
    }
}

enum ResolvedTier {
    Single {
        rule_ids: Vec<String>,
        atoms: [AtomIdx; 4],
        terms: Vec<FourierTorsionTerm>,
        source: TorsionKnowledgeSource,
        ring_size: Option<usize>,
    },
    Conflict {
        rule_ids: Vec<String>,
    },
}

/// Dedup/compose/conflict resolution among multiple same-tier candidates
/// for one bond (spec §4). Two candidates are:
/// - equivalent if their term lists are identical (same periodicity, phase,
///   amplitude for every term, any order) -- deduplicated to one.
/// - composable if their term sets touch disjoint periodicities -- merged
///   into one potential carrying every term.
/// - conflicting if they share a periodicity with a different phase or
///   amplitude -- reported as [`ResolvedTier::Conflict`], never silently
///   resolved by picking one side.
fn resolve_same_tier(candidates: &[Candidate]) -> ResolvedTier {
    if candidates.len() == 1 {
        return ResolvedTier::Single {
            rule_ids: vec![candidates[0].rule_id.clone()],
            atoms: candidates[0].atoms,
            terms: candidates[0].terms.clone(),
            source: candidates[0].source,
            ring_size: candidates[0].ring_size,
        };
    }

    // Merge one at a time, checking for conflicts against the accumulated set.
    let mut merged_terms: Vec<FourierTorsionTerm> = Vec::new();
    let mut rule_ids: Vec<String> = Vec::new();
    let source = candidates[0].source;
    let mut ring_size = candidates[0].ring_size;

    for cand in candidates {
        rule_ids.push(cand.rule_id.clone());
        for &term in &cand.terms {
            if let Some(existing) = merged_terms
                .iter()
                .find(|t| t.periodicity == term.periodicity)
            {
                let same = (existing.phase_deg - term.phase_deg).abs() < 1e-9
                    && (existing.amplitude - term.amplitude).abs() < 1e-9;
                if !same {
                    return ResolvedTier::Conflict {
                        rule_ids: candidates.iter().map(|c| c.rule_id.clone()).collect(),
                    };
                }
                // Equivalent term already present -- skip (dedup).
            } else {
                merged_terms.push(term);
            }
        }
        if cand.source != source {
            // Different sources should never share a tier in this
            // implementation (each tier is homogeneous by construction),
            // but guard defensively rather than silently picking one.
            return ResolvedTier::Conflict {
                rule_ids: candidates.iter().map(|c| c.rule_id.clone()).collect(),
            };
        }
        if cand.ring_size != ring_size {
            ring_size = ring_size.or(cand.ring_size);
        }
    }

    ResolvedTier::Single {
        rule_ids,
        atoms: candidates[0].atoms,
        terms: merged_terms,
        source,
        ring_size,
    }
}

/// Tier 5 (basic chemical knowledge) for one bond: flat-ring planarity term
/// if the bond sits in a 4-6-membered all-sp2 ring, else a linear-sp term if
/// either endpoint is a linear sp center, else `None`.
fn basic_knowledge_term_for_bond(
    mol: &Molecule,
    ring_index: &RingMembershipIndex,
    a: AtomIdx,
    b: AtomIdx,
) -> Option<(String, [AtomIdx; 4], FourierTorsionTerm)> {
    for ring in ring_index.rings() {
        let n = ring.len();
        if !rules_basic::FLAT_RING_SIZES.contains(&n) {
            continue;
        }
        if let Some(i) = ring.iter().position(|&x| x == a) {
            let next = ring[(i + 1) % n];
            let prev = ring[(i + n - 1) % n];
            if next == b {
                let d = ring[(i + 2) % n];
                if rules_basic::flat_ring_applies(mol, n, [prev, a, b, d]) {
                    return Some((
                        "basic:flat_ring".to_string(),
                        [prev, a, b, d],
                        rules_basic::flat_ring_term(),
                    ));
                }
            } else if prev == b {
                let d = ring[(i + n - 2) % n];
                if rules_basic::flat_ring_applies(mol, n, [d, b, a, next]) {
                    return Some((
                        "basic:flat_ring".to_string(),
                        [d, b, a, next],
                        rules_basic::flat_ring_term(),
                    ));
                }
            }
        }
    }
    if rules_basic::is_linear_sp_center(mol, a) || rules_basic::is_linear_sp_center(mol, b) {
        // A genuine A-D pair for the linear-sp diagnostic case: any other
        // neighbor of a/b distinct from the bond itself.
        let a_other = mol.neighbors(a).map(|(n, _)| n).find(|&n| n != b)?;
        let d_other = mol.neighbors(b).map(|(n, _)| n).find(|&n| n != a)?;
        return Some((
            "basic:linear_sp_center".to_string(),
            [a_other, a, b, d_other],
            rules_basic::linear_sp_term(),
        ));
    }
    None
}

// ---------------------------------------------------------------------------
// Tests, including spec §13's required negative controls that need access to
// private internals (`resolve_same_tier`, `Candidate`) not reachable from an
// integration test. See `tests/torsion_knowledge_negative_controls.rs` for
// the negative controls expressible through the public API instead.
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn fake_candidate(rule_id: &str, terms: Vec<FourierTorsionTerm>) -> Candidate {
        Candidate {
            rule_id: rule_id.to_string(),
            atoms: [AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3)],
            terms,
            source: TorsionKnowledgeSource::StandardExperimental,
            ring_size: None,
        }
    }

    /// Negative control (spec §13, "swapping rule application order changes
    /// the result" must FAIL to occur): two candidates with genuinely
    /// different (conflicting) terms must resolve to `Conflict` regardless
    /// of the order they're fed in -- order must never change the verdict.
    #[test]
    fn resolve_same_tier_is_order_independent_for_a_genuine_conflict() {
        let cand_a = fake_candidate("rule:a", vec![FourierTorsionTerm::from_rdkit(2, 1, 10.0)]);
        let cand_b = fake_candidate("rule:b", vec![FourierTorsionTerm::from_rdkit(2, -1, 10.0)]);

        let forward = resolve_same_tier(&[
            fake_candidate("rule:a", cand_a.terms.clone()),
            fake_candidate("rule:b", cand_b.terms.clone()),
        ]);
        let backward = resolve_same_tier(&[
            fake_candidate("rule:b", cand_b.terms.clone()),
            fake_candidate("rule:a", cand_a.terms.clone()),
        ]);

        assert!(
            matches!(forward, ResolvedTier::Conflict { .. }),
            "forward order must detect the conflict"
        );
        assert!(
            matches!(backward, ResolvedTier::Conflict { .. }),
            "backward order must detect the SAME conflict -- order must not change the verdict"
        );
    }

    /// Negative control (spec §13, "silently accepting same-tier
    /// conflicting rules" must FAIL to occur): two candidates disagreeing on
    /// the same periodicity's phase/amplitude must never resolve to
    /// `Single` (a silent pick of one side).
    #[test]
    fn resolve_same_tier_never_silently_picks_a_side_of_a_real_conflict() {
        let candidates = vec![
            fake_candidate("rule:x", vec![FourierTorsionTerm::from_rdkit(1, 1, 50.0)]),
            fake_candidate("rule:y", vec![FourierTorsionTerm::from_rdkit(1, -1, 50.0)]),
        ];
        match resolve_same_tier(&candidates) {
            ResolvedTier::Conflict { rule_ids } => {
                assert_eq!(rule_ids.len(), 2);
            }
            ResolvedTier::Single { .. } => {
                panic!("a genuine same-periodicity phase disagreement must never resolve to Single")
            }
        }
    }

    /// Positive control paired with the above: candidates that are truly
    /// equivalent (identical terms) DO dedupe to `Single` -- confirms the
    /// conflict detection isn't simply "always Conflict for 2+ candidates."
    #[test]
    fn resolve_same_tier_dedupes_genuinely_equivalent_candidates() {
        let terms = vec![FourierTorsionTerm::from_rdkit(2, 1, 5.0)];
        let candidates = vec![
            fake_candidate("rule:x", terms.clone()),
            fake_candidate("rule:y", terms),
        ];
        match resolve_same_tier(&candidates) {
            ResolvedTier::Single {
                rule_ids, terms, ..
            } => {
                assert_eq!(rule_ids.len(), 2);
                assert_eq!(terms.len(), 1);
            }
            ResolvedTier::Conflict { .. } => panic!("identical terms must dedupe, not conflict"),
        }
    }

    /// Composable (disjoint-periodicity) candidates merge into one
    /// multi-term potential rather than conflicting.
    #[test]
    fn resolve_same_tier_composes_disjoint_periodicities() {
        let candidates = vec![
            fake_candidate("rule:x", vec![FourierTorsionTerm::from_rdkit(1, 1, 5.0)]),
            fake_candidate("rule:y", vec![FourierTorsionTerm::from_rdkit(3, -1, 2.0)]),
        ];
        match resolve_same_tier(&candidates) {
            ResolvedTier::Single { terms, .. } => assert_eq!(terms.len(), 2),
            ResolvedTier::Conflict { .. } => {
                panic!("disjoint periodicities must compose, not conflict")
            }
        }
    }

    /// Negative control (spec §13, "applying an acyclic torsion rule to an
    /// aromatic ring bond" must FAIL to occur): benzene must produce zero
    /// `StandardExperimental` potentials -- every one of its bonds is an
    /// aromatic ring bond, skipped before any tier is even tried.
    #[test]
    fn standard_tier_never_fires_on_an_aromatic_ring_bond() {
        let mol = parse("c1ccccc1").unwrap();
        let config = TorsionKnowledgeConfig {
            use_exp_torsions: true,
            ..TorsionKnowledgeConfig::default()
        };
        let report = build_torsion_knowledge(&mol, &config);
        assert!(
            report
                .potentials
                .iter()
                .all(|p| p.source != TorsionKnowledgeSource::StandardExperimental),
            "no StandardExperimental potential should ever target an aromatic ring bond: {:?}",
            report
                .potentials
                .iter()
                .map(|p| &p.rule_id)
                .collect::<Vec<_>>()
        );
    }

    /// Negative control (spec §13, "applying a macrocycle rule to a
    /// 3-8-membered ring" must FAIL to occur): cyclohexane (6-membered) must
    /// produce zero `MacrocycleAdaptation` potentials even with the
    /// macrocycle flag enabled.
    #[test]
    fn macrocycle_tier_never_fires_on_a_small_ring() {
        let mol = parse("C1CCCCC1").unwrap();
        let config = TorsionKnowledgeConfig {
            use_macrocycle_torsions: true,
            ..TorsionKnowledgeConfig::default()
        };
        let report = build_torsion_knowledge(&mol, &config);
        assert!(
            report
                .potentials
                .iter()
                .all(|p| p.source != TorsionKnowledgeSource::MacrocycleAdaptation),
            "no MacrocycleAdaptation potential should ever target a 6-membered ring"
        );
    }

    /// Negative control (spec §13, "treating a 9+-membered ring as
    /// small-ring" must FAIL to occur): cyclononane (9-membered, the exact
    /// macrocycle boundary) must produce zero `SmallRingExperimental`
    /// potentials even with the small-ring flag enabled.
    #[test]
    fn small_ring_tier_never_fires_on_a_macrocycle() {
        let mol = parse("C1CCCCCCCC1").unwrap(); // cyclononane
        let config = TorsionKnowledgeConfig {
            use_small_ring_torsions: true,
            ..TorsionKnowledgeConfig::default()
        };
        let report = build_torsion_knowledge(&mol, &config);
        assert!(
            report
                .potentials
                .iter()
                .all(|p| p.source != TorsionKnowledgeSource::SmallRingExperimental),
            "no SmallRingExperimental potential should ever target a 9-membered ring"
        );
    }

    /// Negative control (spec §13, "coordinates changing when all flags are
    /// false" must FAIL to occur): all-flags-false must be a true report-
    /// level no-op (empty in every field), the precondition for the
    /// coordinate-level no-op the gap-check example measures.
    #[test]
    fn all_flags_false_report_is_genuinely_empty() {
        let mol = parse("CC(=O)Nc1ccc(O)cc1").unwrap(); // paracetamol, exercises many tiers
        let report = build_torsion_knowledge(&mol, &TorsionKnowledgeConfig::default());
        assert!(report.potentials.is_empty());
        assert!(report.matched_rule_ids.is_empty());
        assert!(report.unmatched_rotatable_bonds.is_empty());
        assert!(report.ambiguous_matches.is_empty());
        assert!(report.skipped_bonds.is_empty());
    }
}
