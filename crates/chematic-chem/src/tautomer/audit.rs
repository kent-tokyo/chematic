/// Budget limits for [`tautomer_parent`].
///
/// `max_transforms`/`max_tautomers` are deterministic budgets, counted
/// exactly as [`TautomerConfig::max_iter`]/[`TautomerConfig::max_tautomers`]
/// always have been: `max_transforms` counts outer-loop iterations in which
/// a transform was actually applied (a rule that matches but produces an
/// already-seen fingerprint does not consume budget); `max_tautomers`
/// counts distinct direct-aromatic-shift candidates found while choosing
/// among them at the end (duplicates rejected by fingerprint do not consume
/// budget either). Both are reproducible: same input, same limits, same
/// result, always.
///
/// `timeout_ms` is a different kind of bound: a wall-clock timeout depends
/// on machine speed and load, so a `TimedOut` result is explicitly outside
/// the "same canonical tautomer regardless of input" determinism guarantee
/// the other two limits carry. `None` (the default) performs no wall-clock
/// check.
///
/// `max_tautomers: 0` is treated as "skip the direct-aromatic-shift
/// comparison step entirely" (status stays `Completed`, using whatever the
/// rule-based loop converged to) rather than reporting `MaxTautomersReached`
/// on every input regardless of whether any alternate form actually
/// exists -- the latter would be a boundary artifact of the search
/// algorithm, not a meaningful signal.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TautomerLimits {
    pub max_transforms: usize,
    pub max_tautomers: usize,
    pub timeout_ms: Option<u64>,
}

impl Default for TautomerLimits {
    fn default() -> Self {
        Self {
            max_transforms: 16,
            max_tautomers: 32,
            timeout_ms: None,
        }
    }
}

/// Stable identifier for one of the built-in [`TautomerRule`]s, used in
/// [`AppliedTransform::rule_id`]. `#[non_exhaustive]` so new rules can be
/// added without a breaking change; `Other` is a defensive fallback for a
/// rule name not yet mapped to a variant here -- should not occur for any
/// rule currently in `RULES` (see `all_rules_map_to_a_named_tautomer_rule_id`
/// in this module's tests).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum TautomerRuleId {
    KetoEnol,
    AmideIminol,
    IminolAmide,
    ImineEnamine,
    OneThreeNToO,
    OneThreeNToN,
    Thioamide,
    ThioIminolAmide,
    ThioKetoEnol,
    ThioEnolKetone,
    OneThreeNToS,
    OneThreeSToO,
    OneThreeSToN,
    OneThreeOToS,
    OneThreeSToS,
    OneThreeOToNAnyBridge,
    OneThreeOToOAnyBridge,
    OneThreeNToCAnyBridge,
    OneThreeCToOAnyBridge,
    NitrosoOxime,
    OneThreeCToNAnyBridge,
    OneFiveOToOBetaDiketone,
    OneFiveOToN,
    OneFiveNToN,
    OneFiveNToO,
    OneFiveCToO,
    OneFiveOToONBridge,
    OneFiveOToOSBridge,
    OneFiveNToONBridge,
    OneFiveNToOSBridge,
    OneFiveNToNNBridge,
    OneFiveNToNSBridge,
    OneFiveOToNNBridge,
    OneFiveOToNSBridge,
    OneFiveSToONBridge,
    OneFiveSToOSBridge,
    OneFiveSToNCBridge,
    OneFiveSToNNBridge,
    OneFiveCToNCBridge,
    OneFiveCToNNBridge,
    OneFiveCToNSBridge,
    OneFiveCToSNBridge,
    OneFiveCToSSBridge,
    /// The round 2C-2 exocyclic lactam/lactim shift -- not a `RULES` table
    /// entry (it's a directional step applied alongside them, see
    /// `apply_exocyclic_lactam_shift_tracked`), so it is never produced by
    /// `rule_id_for` and is exempt from
    /// `all_rules_map_to_a_named_tautomer_rule_id`'s drift check.
    AromaticExocyclicLactamLactim,
    /// H transfer between the two aromatic nitrogens flanking a carbonyl.
    AromaticDualFlank,
    /// A rule name not yet mapped to a named variant above.
    Other(&'static str),
}

fn rule_id_for(name: &'static str) -> TautomerRuleId {
    match name {
        "keto-enol" => TautomerRuleId::KetoEnol,
        "amide-iminol" => TautomerRuleId::AmideIminol,
        "iminol-amide" => TautomerRuleId::IminolAmide,
        "imine-enamine" => TautomerRuleId::ImineEnamine,
        "1,3-N-to-O" => TautomerRuleId::OneThreeNToO,
        "1,3-N-to-N" => TautomerRuleId::OneThreeNToN,
        "thioamide" => TautomerRuleId::Thioamide,
        "thio-iminol-amide" => TautomerRuleId::ThioIminolAmide,
        "thio-keto-enol" => TautomerRuleId::ThioKetoEnol,
        "thio-enol-ketone" => TautomerRuleId::ThioEnolKetone,
        "1,3-N-to-S" => TautomerRuleId::OneThreeNToS,
        "1,3-S-to-O" => TautomerRuleId::OneThreeSToO,
        "1,3-S-to-N" => TautomerRuleId::OneThreeSToN,
        "1,3-O-to-S" => TautomerRuleId::OneThreeOToS,
        "1,3-S-to-S" => TautomerRuleId::OneThreeSToS,
        "1,3-O-to-N-any-bridge" => TautomerRuleId::OneThreeOToNAnyBridge,
        "1,3-O-to-O-any-bridge" => TautomerRuleId::OneThreeOToOAnyBridge,
        "1,3-N-to-C-any-bridge" => TautomerRuleId::OneThreeNToCAnyBridge,
        "1,3-C-to-O-any-bridge" => TautomerRuleId::OneThreeCToOAnyBridge,
        "nitroso-oxime" => TautomerRuleId::NitrosoOxime,
        "1,3-C-to-N-any-bridge" => TautomerRuleId::OneThreeCToNAnyBridge,
        "1,5-O-to-O-beta-diketone" => TautomerRuleId::OneFiveOToOBetaDiketone,
        "1,5-O-to-N" => TautomerRuleId::OneFiveOToN,
        "1,5-N-to-N" => TautomerRuleId::OneFiveNToN,
        "1,5-N-to-O" => TautomerRuleId::OneFiveNToO,
        "1,5-C-to-O" => TautomerRuleId::OneFiveCToO,
        "1,5-O-to-O-N-bridge" => TautomerRuleId::OneFiveOToONBridge,
        "1,5-O-to-O-S-bridge" => TautomerRuleId::OneFiveOToOSBridge,
        "1,5-N-to-O-N-bridge" => TautomerRuleId::OneFiveNToONBridge,
        "1,5-N-to-O-S-bridge" => TautomerRuleId::OneFiveNToOSBridge,
        "1,5-N-to-N-N-bridge" => TautomerRuleId::OneFiveNToNNBridge,
        "1,5-N-to-N-S-bridge" => TautomerRuleId::OneFiveNToNSBridge,
        "1,5-O-to-N-N-bridge" => TautomerRuleId::OneFiveOToNNBridge,
        "1,5-O-to-N-S-bridge" => TautomerRuleId::OneFiveOToNSBridge,
        "1,5-S-to-O-N-bridge" => TautomerRuleId::OneFiveSToONBridge,
        "1,5-S-to-O-S-bridge" => TautomerRuleId::OneFiveSToOSBridge,
        "1,5-S-to-N-C-bridge" => TautomerRuleId::OneFiveSToNCBridge,
        "1,5-S-to-N-N-bridge" => TautomerRuleId::OneFiveSToNNBridge,
        "1,5-C-to-N-C-bridge" => TautomerRuleId::OneFiveCToNCBridge,
        "1,5-C-to-N-N-bridge" => TautomerRuleId::OneFiveCToNNBridge,
        "1,5-C-to-N-S-bridge" => TautomerRuleId::OneFiveCToNSBridge,
        "1,5-C-to-S-N-bridge" => TautomerRuleId::OneFiveCToSNBridge,
        "1,5-C-to-S-S-bridge" => TautomerRuleId::OneFiveCToSSBridge,
        other => TautomerRuleId::Other(other),
    }
}

/// One named contribution to a candidate tautomer's [`tautomer_score`],
/// decomposed for [`TautomerAuditRecord::score_breakdown`].
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum TautomerScoreTerm {
    AromaticRing,
    /// `element` is an atomic number (today: 8=O, 7=N, 16=S).
    HeteroatomHydrogen {
        element: u8,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ScoreContribution {
    pub term: TautomerScoreTerm,
    pub value: i32,
}

/// One transformation applied by the rule-based 1,3-/1,5-shift loop while
/// computing a [`TautomerAuditRecord`].
///
/// Does not implement `Serialize`/`Deserialize` even under the `serde`
/// feature: `AtomIdx`/`BondIdx` (from `chematic-core`) do not implement
/// them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AppliedTransform {
    pub rule_id: TautomerRuleId,
    pub affected_atoms: Vec<AtomIdx>,
    pub affected_bonds: Vec<BondIdx>,
}

/// Explainable audit record for one [`tautomer_parent`] call.
///
/// `applied_transforms` records only the rule-based 1,3-/1,5-shift loop (the
/// same loop [`canonical_tautomer_with_config`] runs); the final
/// direct-aromatic-shift tie-break contributes to `score_breakdown` and
/// `candidate_count` instead, since that search does not track the specific
/// shift sequence taken to reach each candidate -- a disclosed scope
/// limitation, not an oversight.
///
/// `lost_stereo`/`affected_isotopes` are atom-index lists (empty = none
/// affected), not booleans: every transform this module applies preserves
/// `chirality`/`isotope` (see the RFC's non-bug finding, section 1.2), so
/// both are expected to always be empty today; the check is real (not
/// hardcoded) so it stays correct if that ever changes.
///
/// Does not implement `Serialize`/`Deserialize` even under the `serde`
/// feature: it embeds [`AppliedTransform`] (see its own doc comment).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TautomerAuditRecord {
    pub selected: MoleculeSnapshot,
    pub candidate_count: usize,
    pub score_breakdown: Vec<ScoreContribution>,
    pub applied_transforms: Vec<AppliedTransform>,
    pub lost_stereo: Vec<AtomIdx>,
    pub affected_isotopes: Vec<AtomIdx>,
}

type ParentDeadline = Option<(Instant, std::time::Duration)>;

struct ParentSearchState {
    current: Molecule,
    seen: HashSet<u64>,
    applied_transforms: Vec<AppliedTransform>,
    transforms_applied: usize,
    status: ParentComputationStatus,
}

impl ParentSearchState {
    fn new(mol: &Molecule) -> Self {
        let current = mol.clone();
        let mut seen = HashSet::new();
        seen.insert(mol_fingerprint(&current));
        Self {
            current,
            seen,
            applied_transforms: Vec::new(),
            transforms_applied: 0,
            status: ParentComputationStatus::Completed,
        }
    }

    fn has_more_transforms(&self, config: &TautomerConfig, rank: &[u32]) -> bool {
        active_rules(config)
            .into_iter()
            .filter(|rule| rule.prefer_forward)
            .any(|rule| {
                apply_first_match_tracked(&self.current, rule, config, rank)
                    .is_some_and(|(next, ..)| !self.seen.contains(&mol_fingerprint(&next)))
            })
            || apply_exocyclic_lactam_shift_tracked(&self.current, config, rank)
                .is_some_and(|(next, ..)| !self.seen.contains(&mol_fingerprint(&next)))
            || apply_dual_flank_tracked(&self.current, config, rank)
                .is_some_and(|(next, ..)| !self.seen.contains(&mol_fingerprint(&next)))
    }

    fn apply_rule_round(
        &mut self,
        config: &TautomerConfig,
        rank: &[u32],
        max_transforms: usize,
    ) -> bool {
        for rule in active_rules(config)
            .into_iter()
            .filter(|rule| rule.prefer_forward)
        {
            let remaining = max_transforms - self.transforms_applied;
            let (next, applied) = apply_nonconflicting_matches_tracked(
                &self.current,
                rule,
                config,
                rank,
                remaining,
            );
            if applied.is_empty() || !self.seen.insert(mol_fingerprint(&next)) {
                continue;
            }
            self.record_rule_matches(rule_id_for(rule.name), &applied);
            self.current = next;
            self.transforms_applied += applied.len();
            return true;
        }
        false
    }

    fn record_rule_matches(
        &mut self,
        rule_id: TautomerRuleId,
        applied: &[(AtomIdx, AtomIdx, AtomIdx)],
    ) {
        for &(donor, bridge, acceptor) in applied {
            let mut affected_bonds = Vec::new();
            if let Some((bond, _)) = self.current.bond_between(donor, bridge) {
                affected_bonds.push(bond);
            }
            if let Some((bond, _)) = self.current.bond_between(bridge, acceptor) {
                affected_bonds.push(bond);
            }
            self.applied_transforms.push(AppliedTransform {
                rule_id,
                affected_atoms: vec![donor, bridge, acceptor],
                affected_bonds,
            });
        }
    }

    fn apply_fallback_round(&mut self, config: &TautomerConfig, rank: &[u32]) -> bool {
        if let Some((next, donor, bridge, acceptor)) =
            apply_exocyclic_lactam_shift_tracked(&self.current, config, rank)
            && self.seen.insert(mol_fingerprint(&next))
        {
            self.record_rule_matches(
                TautomerRuleId::AromaticExocyclicLactamLactim,
                &[(donor, bridge, acceptor)],
            );
            self.current = next;
            self.transforms_applied += 1;
            return true;
        }
        if let Some((next, donor, bridge, acceptor)) =
            apply_dual_flank_tracked(&self.current, config, rank)
            && self.seen.insert(mol_fingerprint(&next))
        {
            self.applied_transforms.push(AppliedTransform {
                rule_id: TautomerRuleId::AromaticDualFlank,
                affected_atoms: vec![donor, bridge, acceptor],
                affected_bonds: Vec::new(),
            });
            self.current = next;
            self.transforms_applied += 1;
            return true;
        }
        false
    }
}

fn parent_deadline_expired(deadline: &ParentDeadline) -> bool {
    deadline
        .as_ref()
        .is_some_and(|(start, budget)| start.elapsed() >= *budget)
}

fn empty_parent_result(mol: &Molecule) -> ParentResult {
    ParentResult {
        molecule: mol.clone(),
        status: ParentComputationStatus::InvalidInput(InvalidInputReason::EmptyMolecule),
        audit: ParentAudit::Tautomer(TautomerAuditRecord {
            selected: MoleculeSnapshot::from_mol(mol),
            candidate_count: 0,
            score_breakdown: Vec::new(),
            applied_transforms: Vec::new(),
            lost_stereo: Vec::new(),
            affected_isotopes: Vec::new(),
        }),
    }
}

fn select_parent_candidate(
    state: &mut ParentSearchState,
    limits: &TautomerLimits,
    rank: &[u32],
    deadline: &ParentDeadline,
) -> usize {
    if !matches!(state.status, ParentComputationStatus::Completed) || limits.max_tautomers == 0 {
        return 1;
    }
    if parent_deadline_expired(deadline) {
        state.status = ParentComputationStatus::TimedOut;
        return 1;
    }
    let (extra, truncated) = enumerate_direct_aromatic_forms_tracked(
        &state.current,
        &HashSet::new(),
        limits.max_tautomers,
        rank,
    );
    let mut candidates = vec![state.current.clone()];
    candidates.extend(extra);
    let candidate_count = candidates.len();
    if candidates.len() > 1 {
        candidates.sort_by(|a, b| {
            tautomer_score(b).cmp(&tautomer_score(a)).then_with(|| {
                chematic_smiles::canonical_smiles(a).cmp(&chematic_smiles::canonical_smiles(b))
            })
        });
        state.current = candidates.remove(0);
    }
    if truncated {
        state.status = ParentComputationStatus::MaxTautomersReached;
    }
    candidate_count
}

fn finish_parent_result(
    original: &Molecule,
    state: ParentSearchState,
    candidate_count: usize,
) -> ParentResult {
    let common_atoms = original.atom_count().min(state.current.atom_count());
    let lost_stereo = (0..common_atoms)
        .map(|i| AtomIdx(i as u32))
        .filter(|&i| {
            original.atom(i).chirality != Chirality::None
                && state.current.atom(i).chirality == Chirality::None
        })
        .collect();
    let affected_isotopes = (0..common_atoms)
        .map(|i| AtomIdx(i as u32))
        .filter(|&i| original.atom(i).isotope != state.current.atom(i).isotope)
        .collect();
    ParentResult {
        status: state.status,
        audit: ParentAudit::Tautomer(TautomerAuditRecord {
            selected: MoleculeSnapshot::from_mol(&state.current),
            candidate_count,
            score_breakdown: score_breakdown(&state.current),
            applied_transforms: state.applied_transforms,
            lost_stereo,
            affected_isotopes,
        }),
        molecule: state.current,
    }
}

/// Like [`canonical_tautomer`], but budget-limited via [`TautomerLimits`]
/// and returns a [`ParentResult`] (status + explainable
/// [`TautomerAuditRecord`]) instead of a bare `Molecule` -- see
/// `docs/rfcs/tautomer_parent_identity_phase2_rfc.md` sections 4.2-4.3 and
/// 4.4a.
///
/// Includes the section 4.4a directional exocyclic lactam/lactim shift
/// (round 2C-2, `apply_exocyclic_lactam_shift_tracked`), wired into this
/// function's own loop the same way it's wired into
/// [`canonical_tautomer_with_config`]'s -- `MaxTransformsReached` accounts
/// for it. The former Section 1.7 cytosine/guanine dual-flank residual is
/// handled by rooted N-H normalization and is covered by the
/// `tp2_07_09_dual_flank_cytosine_and_guanine_converge` regression. The
/// former Section 1.6 nitroso/oxime residual (`CCN=O`/`CC=NO`) is handled by
/// the dedicated forward `nitroso-oxime` rule and is covered by the
/// `tp2_04_nitroso_oxime_converges_without_enolizing_ketones` regression.
/// The asymmetric `tp2-39` and N9-methylhypoxanthine holdout remain explicit
/// residuals because their broader ring-bond tautomer components are not yet
/// fully enumerated by the canonical search.
///
/// Independent matches of one rule are applied together in canonical-rank
/// order, subject to `limits.max_transforms`; each applied match is retained
/// in `TautomerAuditRecord::applied_transforms`. This keeps bounded parent
/// selection atom-order invariant without weakening the transform limit.
pub fn tautomer_parent(mol: &Molecule, limits: &TautomerLimits) -> ParentResult {
    if mol.atom_count() == 0 {
        return empty_parent_result(mol);
    }

    // Every candidate/ring tie broken via `rank` (issue #415 residual --
    // see the "Atom-order canonicalization" section above), so the result
    // is independent of `mol`'s caller-supplied atom order.
    let rank = atom_rank(mol);

    let config = TautomerConfig {
        max_iter: limits.max_transforms,
        max_tautomers: limits.max_tautomers,
        ..TautomerConfig::default()
    };
    let deadline = limits
        .timeout_ms
        .map(|ms| (Instant::now(), std::time::Duration::from_millis(ms)));
    let mut state = ParentSearchState::new(mol);

    loop {
        if parent_deadline_expired(&deadline) {
            state.status = ParentComputationStatus::TimedOut;
            break;
        }
        if state.transforms_applied >= limits.max_transforms {
            if state.has_more_transforms(&config, &rank) {
                state.status = ParentComputationStatus::MaxTransformsReached;
            }
            break;
        }
        let changed = state.apply_rule_round(&config, &rank, limits.max_transforms)
            || state.apply_fallback_round(&config, &rank);
        if !changed {
            break; // converged: no forward rule produces an unseen form
        }
    }
    let candidate_count = select_parent_candidate(&mut state, limits, &rank, &deadline);
    finish_parent_result(mol, state, candidate_count)
}
