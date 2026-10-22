//! Molecular standardization routines.
//!
//! Provides utilities for cleaning up molecular representations:
//! - Selecting the largest connected fragment.
//! - Neutralizing simple formal charges.

#![forbid(unsafe_code)]

use std::collections::{HashMap, VecDeque};

use chematic_core::{AtomIdx, BondIdx, Element, Molecule, MoleculeBuilder};

use chematic_smarts::{MatchConfig, parse_smarts};

mod cleanup;
mod groups;
mod reionize;
pub use cleanup::cleanup;
pub use groups::normalize_groups;

/// Salt removal catalog: common salt patterns (counterions and solvates).
///
/// Each pattern is a (name, SMARTS) tuple for organic and inorganic salts.
/// Used by [`remove_salts`] to filter out counterions and keep drug-like fragments.
#[derive(Clone, Debug)]
pub struct SaltCatalog {
    /// (name, SMARTS) pairs for salt patterns
    patterns: Vec<(&'static str, &'static str)>,
}

impl Default for SaltCatalog {
    fn default() -> Self {
        Self::new()
    }
}

impl SaltCatalog {
    /// Create a default salt catalog with common counterions and solvates.
    pub fn new() -> Self {
        Self {
            patterns: vec![
                // Organic salts (carboxylates, sulfonates, etc.)
                ("acetate", "[#6](-[#1])(-[#1])-[#6](=[#8])[O-]"),
                ("formate", "[#6](=[#8])[O-]"),
                (
                    "propionate",
                    "[#6](-[#1])(-[#1])-[#6](-[#1])-[#6](=[#8])[O-]",
                ),
                ("benzoate", "c1ccccc1-[#6](=[#8])[O-]"),
                (
                    "trifluoroacetate",
                    "[#9]-[#6](-[#9])(-[#9])-[#6](=[#8])[O-]",
                ),
                (
                    "mesylate",
                    "[#16](=[#8])(=[#8])-[#8]-[#6](-[#1])(-[#1])-[#1]",
                ),
                ("tosylate", "c1ccc(cc1)-[#16](=[#8])(=[#8])-[#8]"),
                (
                    "nosylate",
                    "[#8]-[#6](-[#1])(-[#1])-[#8]-[#16](=[#8])(=[#8])-c1ccc([N+](=O)[O-])cc1",
                ),
                ("sulfate", "[#16](=[#8])(=[#8])(-[#8])-[#8]"),
                ("phosphate", "[#15](=[#8])(-[#8])(-[#8])-[#8]"),
                (
                    "citrate",
                    "[#6](-[#6](=[#8])[O-])(-[#6](=[#8])[O-])-[#6](-[#8])-[#6](=[#8])[O-]",
                ),
                (
                    "tartrate",
                    "[#6](-[#8])(-[#6](-[#8])-[#6](=[#8])[O-])-[#6](=[#8])[O-]",
                ),
                // Inorganic salts (single atoms/small molecules)
                ("sodium_cation", "[Na+]"),
                ("potassium_cation", "[K+]"),
                ("lithium_cation", "[Li+]"),
                ("calcium_cation", "[Ca+2]"),
                ("magnesium_cation", "[Mg+2]"),
                ("chloride_anion", "[Cl-]"),
                ("bromide_anion", "[Br-]"),
                ("iodide_anion", "[I-]"),
                ("fluoride_anion", "[F-]"),
                ("oxide_anion", "[O-2]"),
                ("sulfate_anion", "[#16](=[#8])(=[#8])(-[#8])-[#8-]"),
                ("phosphate_anion", "[#15](=[#8])(-[#8])(-[#8-])-[#8]"),
                // Solvates and additives
                ("water", "[#8](-[#1])-[#1]"),
                (
                    "dmso",
                    "[#16](=[#8])(-[#6](-[#1])(-[#1])-[#1])-[#6](-[#1])(-[#1])-[#1]",
                ),
                ("methanol", "[#6](-[#1])(-[#1])-[#8]-[#1]"),
                ("ethanol", "[#6](-[#1])(-[#1])-[#6](-[#1])(-[#1])-[#8]-[#1]"),
                (
                    "isopropanol",
                    "[#6](-[#1])(-[#1])-[#6](-[#8]-[#1])(-[#1])-[#6](-[#1])(-[#1])-[#1]",
                ),
                // Rare but important salts
                ("borate", "[#5](-[#8])(-[#8])-[#8]"),
                ("ammonium", "[#7+;H0,H1,H2,H3]"),
            ],
        }
    }

    /// Add a custom salt pattern to this catalog.
    pub fn add(&mut self, name: &'static str, smarts: &'static str) {
        self.patterns.push((name, smarts));
    }

    /// Check if a molecule fragment matches any salt pattern.
    pub fn is_salt(&self, frag: &Molecule) -> bool {
        // max_matches: Some(1) stops each pattern's VF2 search at the first
        // embedding instead of enumerating every match — this loop only
        // needs to know whether a match exists, not how many.
        let config = MatchConfig {
            max_visit_budget: Some(1_000_000),
            max_matches: Some(1),
            uniquify: false,
            ..Default::default()
        };
        for (_, smarts_str) in &self.patterns {
            if let Ok(query) = parse_smarts(smarts_str)
                && chematic_smarts::has_match_with_config(&query, frag, &config)
            {
                return true;
            }
        }
        false
    }
}

/// Find all connected components of `mol` via BFS, sorted descending by size.
fn connected_components(mol: &Molecule) -> Vec<Vec<AtomIdx>> {
    let n = mol.atom_count();
    let mut visited = vec![false; n];
    let mut components: Vec<Vec<AtomIdx>> = Vec::new();

    for start in 0..n {
        if visited[start] {
            continue;
        }
        visited[start] = true;
        let mut component = Vec::new();
        let mut queue: VecDeque<AtomIdx> = VecDeque::new();
        queue.push_back(AtomIdx(start as u32));

        while let Some(current) = queue.pop_front() {
            component.push(current);
            for (neighbor, _) in mol.neighbors(current) {
                let ni = neighbor.0 as usize;
                if !visited[ni] {
                    visited[ni] = true;
                    queue.push_back(neighbor);
                }
            }
        }
        components.push(component);
    }

    components.sort_by_key(|b| std::cmp::Reverse(b.len()));
    components
}

/// Copy bonds from `mol` into `builder` when both endpoints are remapped.
fn copy_bonds(mol: &Molecule, builder: &mut MoleculeBuilder, remap: &HashMap<AtomIdx, AtomIdx>) {
    for i in 0..mol.bond_count() {
        let bond = mol.bond(BondIdx(i as u32));
        if let (Some(&new_a), Some(&new_b)) = (remap.get(&bond.atom1), remap.get(&bond.atom2)) {
            let _ = builder.add_bond(new_a, new_b, bond.order);
        }
    }
}

/// Check if a fragment is a common inorganic salt or counterion.
///
/// Returns true if fragment matches patterns like: NaCl, KCl, Na+, K+, Cl-, Br-, I-, etc.
fn is_salt_fragment(frag: &Molecule) -> bool {
    let n = frag.atom_count();

    // Single atom: check if it's a common counterion
    if n == 1 {
        let atom = frag.atom(AtomIdx(0));
        return matches!(
            atom.element.atomic_number(),
            11 | 19 | 37 | 55 |  // Na, K, Rb, Cs (alkali metals)
            17 | 35 | 53 |       // Cl, Br, I (halogens)
            8 // O (oxide)
        );
    }

    // Two atoms: check for common binary salts (NaCl, KBr, etc.)
    if n == 2 {
        let a0 = frag.atom(AtomIdx(0)).element.atomic_number();
        let a1 = frag.atom(AtomIdx(1)).element.atomic_number();
        let bond_count = frag.bond_count();

        // Ionic pair (no bond between them) — cation + anion
        if bond_count == 0 {
            let metals = [11, 19, 37, 55]; // Na, K, Rb, Cs
            let nonmetals = [17, 35, 53, 8]; // Cl, Br, I, O
            return (metals.contains(&a0) && nonmetals.contains(&a1))
                || (metals.contains(&a1) && nonmetals.contains(&a0));
        }
    }

    // Small molecules with only metal/nonmetal atoms (common solvate salts)
    if n <= 4 {
        let has_organic = frag.atoms().any(|(_, a)| a.element.atomic_number() == 6);
        if !has_organic {
            // Pure inorganic salt (no carbons)
            return true;
        }
    }

    false
}

/// Return a new `Molecule` with salts/solvents removed, keeping the fragment
/// selected by [`FragmentPolicy::default()`] (see [`select_fragment`]).
///
/// Uses a small, structural classification (monatomic counterions, water, a
/// no-carbon/small-fragment heuristic) rather than [`SaltCatalog`]'s named
/// pattern list — see `docs/rfcs/explainable_standardization_phase1_rfc.md`
/// for why. Callers who specifically want the legacy named-catalog behavior
/// can still call [`remove_salts_with_catalog`] directly.
///
/// If no non-salt fragment exists or molecule is empty, returns the largest
/// fragment by the same policy's ranking (never panics, never drops the
/// input to nothing).
pub fn remove_salts(mol: &Molecule) -> Molecule {
    select_fragment_molecule(mol, &FragmentPolicy::default())
}

/// The molecule [`select_fragment`] would return, without building the
/// audit record.
///
/// `select_fragment` canonicalizes every fragment and snapshots the input and
/// output for its [`TransformationRecord`]; callers that only want the kept
/// fragment paid for all of that. Here canonical SMILES are computed only for
/// fragments that tie on `(rank_size, has_carbon)`, which is the only place
/// the canonical-string tie-break can matter. Selection (including
/// `max_by_key`'s last-of-equals rule) is identical.
fn select_fragment_molecule(mol: &Molecule, policy: &FragmentPolicy) -> Molecule {
    if mol.atom_count() == 0 {
        return MoleculeBuilder::new().build();
    }
    let components = connected_components(mol);
    if components.len() == 1 {
        return extract_fragment(mol, &components[0]);
    }
    let frags: Vec<Molecule> = components
        .iter()
        .map(|c| extract_fragment(mol, c))
        .collect();
    let primary: Vec<(usize, bool)> = frags
        .iter()
        .map(|f| {
            let rank_size = if policy.count_heavy_atoms_only {
                heavy_atom_count(f)
            } else {
                f.atom_count()
            };
            (rank_size, policy.prefer_organic && fragment_has_carbon(f))
        })
        .collect();
    let kept: Vec<usize> = (0..frags.len())
        .filter(|&i| classify_fragment(&frags[i], policy) == FragmentClass::Kept)
        .collect();
    let candidates: Vec<usize> = if kept.is_empty() {
        (0..frags.len()).collect()
    } else {
        kept
    };
    let best = candidates.iter().map(|&i| primary[i]).max().unwrap();
    let tied: Vec<usize> = candidates
        .into_iter()
        .filter(|&i| primary[i] == best)
        .collect();
    let winner = if tied.len() == 1 {
        tied[0]
    } else {
        let canon: Vec<(usize, String)> = tied
            .iter()
            .map(|&i| (i, chematic_smiles::canonical_smiles(&frags[i])))
            .collect();
        canon
            .iter()
            .max_by_key(|(_, c)| std::cmp::Reverse(c.clone()))
            .unwrap()
            .0
    };
    frags.into_iter().nth(winner).unwrap()
}

/// Remove salts using a custom catalog.
///
/// # Arguments
/// - `mol`: the molecule to process
/// - `catalog`: custom salt catalog (use `SaltCatalog::new()` for default)
pub fn remove_salts_with_catalog(mol: &Molecule, catalog: &SaltCatalog) -> Molecule {
    if mol.atom_count() == 0 {
        return MoleculeBuilder::new().build();
    }

    let components = connected_components(mol);

    // Find largest non-salt fragment
    let mut largest_non_salt: Option<&Vec<AtomIdx>> = None;
    let mut largest_non_salt_size = 0;

    for component in &components {
        // Extract fragment molecule temporarily to check if it's a salt.
        // `extract_fragment` preserves stereo_neighbor_order/bond_directions/
        // stereo_groups via `Molecule::remove_atom` -- an earlier version of
        // this function built fragments with a bare MoleculeBuilder remap
        // that silently dropped stereo_neighbor_order, flipping @/@@ on any
        // stereocenter-bearing fragment (found while adding Phase 1 tests).
        let frag = extract_fragment(mol, component);

        // Check with catalog first, fall back to heuristic
        let is_salt = catalog.is_salt(&frag) || is_salt_fragment(&frag);

        // If not a salt and larger than current best, use it
        if !is_salt && component.len() > largest_non_salt_size {
            largest_non_salt = Some(component);
            largest_non_salt_size = component.len();
        }
    }

    // Fall back to largest fragment if no non-salt found
    let component = largest_non_salt.unwrap_or(&components[0]);
    extract_fragment(mol, component)
}

/// Return a new `Molecule` containing only the largest connected fragment.
///
/// If the molecule is empty, an empty `Molecule` is returned.
/// This is an alias for `remove_salts()` for backward compatibility.
pub fn largest_fragment(mol: &Molecule) -> Molecule {
    remove_salts(mol)
}

/// Policy controlling fragment ranking and salt/solvent classification for
/// [`select_fragment`].
///
/// Classification is purely structural (monatomic-ion identity, water,
/// no-carbon/small-fragment shape) — it never depends on matching a named
/// substance pattern, unlike [`SaltCatalog`]. See
/// `docs/rfcs/explainable_standardization_phase1_rfc.md` section 3.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FragmentPolicy {
    /// Rank fragments by heavy-atom count rather than total atom count, so
    /// an explicit-hydrogen spelling of a small fragment can't outrank a
    /// larger heavy-atom-bearing fragment.
    pub count_heavy_atoms_only: bool,
    /// Prefer a carbon-containing fragment when ranks would otherwise tie.
    pub prefer_organic: bool,
    /// Never classify an isotopically-labeled fragment as salt/solvent, even
    /// if it would otherwise match the structural heuristic.
    pub preserve_isotopes: bool,
    /// Reserved for a future policy refinement (prefer keeping a counterion
    /// when the pipeline can't confidently identify an organic parent at
    /// all). Not consulted by [`select_fragment`] yet.
    pub preserve_counterion_if_required: bool,
}

impl Default for FragmentPolicy {
    fn default() -> Self {
        Self {
            count_heavy_atoms_only: true,
            prefer_organic: true,
            preserve_isotopes: true,
            preserve_counterion_if_required: false,
        }
    }
}

/// Per-fragment atom/formula/canonical-SMILES snapshot, used in
/// [`FragmentRecord`].
#[derive(Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FragmentSnapshot {
    /// Total atom count (including explicit hydrogens, if any).
    pub atom_count: usize,
    /// Heavy-atom (non-hydrogen) count.
    pub heavy_atom_count: usize,
    /// Hill-order molecular formula.
    pub formula: String,
    /// Canonical SMILES of this fragment alone.
    pub canonical_smiles: String,
}

/// Why a fragment was kept or removed by [`select_fragment`].
#[derive(Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum FragmentDecision {
    /// This fragment is the pipeline's output. `rank_key` is a
    /// human-readable rendering of the ranking fields that won.
    Kept {
        /// Human-readable ranking-field summary, e.g.
        /// `"rank_size=7, has_carbon=true, canonical_smiles=..."`.
        rank_key: String,
    },
    /// This fragment was excluded from the output.
    Removed {
        /// Stable, machine-readable rule identifier, e.g.
        /// `"monatomic_always_strip_ion"`, `"water_always_strip"`,
        /// `"structural_no_carbon_small_fragment"`, `"ranked_out_by_size"`,
        /// `"duplicate_of_kept_fragment"`.
        rule_id: String,
        /// Human-readable explanation.
        reason: String,
    },
}

/// One fragment's classification/ranking outcome within a
/// [`TransformationRecord`].
#[derive(Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FragmentRecord {
    /// Structure summary of this fragment as it appeared in the input.
    pub snapshot: FragmentSnapshot,
    /// Why this fragment was kept or removed.
    pub decision: FragmentDecision,
}

/// Explainable audit record for one [`select_fragment`] call: every input
/// fragment, its classification/ranking decision, and before/after
/// structure — see
/// `docs/rfcs/explainable_standardization_phase1_rfc.md` section 4.
#[derive(Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TransformationRecord {
    /// Stable rule identifier for this transformation as a whole.
    pub rule_id: String,
    /// Rule version, bumped whenever the classification/ranking policy's
    /// observable behavior changes.
    pub rule_version: u32,
    /// Every input fragment (kept and removed), in a deterministic,
    /// canonical-SMILES-ordered sequence — never input/discovery order.
    pub fragments: Vec<FragmentRecord>,
    /// Set when no fragment classified as a confident non-salt candidate
    /// (e.g. a purely inorganic input like `NaCl`), explaining that the
    /// kept fragment was chosen by falling back to ranking among
    /// salt-classified fragments rather than from a real organic parent.
    pub abstained: Option<String>,
    /// Molecule summary before this transformation.
    pub before: MoleculeSnapshot,
    /// Molecule summary after this transformation.
    pub after: MoleculeSnapshot,
    /// Warnings raised while processing.
    pub warnings: Vec<StandardizationWarning>,
}

/// Monatomic ions always classified as salt/counterion, regardless of what
/// else is present in the input. Deliberately tiny (fixed-identity,
/// commonly-encountered single-atom counterions) rather than a general
/// "known salts" list — see
/// `docs/rfcs/explainable_standardization_phase1_rfc.md` section 3.2.
/// Atomic numbers: Li, Na, K, F, Cl, Br, I.
const ALWAYS_STRIP_MONATOMIC_IONS: &[u8] = &[3, 11, 19, 9, 17, 35, 53];

/// Upper bound on heavy-atom count for the structural "no carbon → likely
/// salt/solvent" fallback. Inherited from the pre-existing `is_salt_fragment`
/// threshold (this module); flagged in the RFC as open to revisiting via
/// fixture evidence, not re-derived from first principles here.
const MAX_SALT_HEAVY_ATOMS_NO_CARBON: usize = 4;

fn heavy_atom_count(frag: &Molecule) -> usize {
    frag.atoms()
        .filter(|(_, a)| a.element.atomic_number() != 1)
        .count()
}

fn fragment_has_carbon(frag: &Molecule) -> bool {
    frag.atoms().any(|(_, a)| a.element.atomic_number() == 6)
}

fn fragment_has_isotope(frag: &Molecule) -> bool {
    frag.atoms().any(|(_, a)| a.isotope.is_some())
}

/// A lone, neutral oxygen fragment — water, however many (implicit or
/// explicit) hydrogens accompany it in the graph.
fn is_water_fragment(frag: &Molecule) -> bool {
    heavy_atom_count(frag) == 1
        && frag
            .atoms()
            .any(|(_, a)| a.element.atomic_number() == 8 && a.charge == 0)
}

fn is_always_strip_monatomic_ion(frag: &Molecule) -> bool {
    if heavy_atom_count(frag) != 1 {
        return false;
    }
    frag.atoms().any(|(_, a)| {
        a.charge != 0 && ALWAYS_STRIP_MONATOMIC_IONS.contains(&a.element.atomic_number())
    })
}

fn is_structural_salt_fallback(frag: &Molecule) -> bool {
    heavy_atom_count(frag) <= MAX_SALT_HEAVY_ATOMS_NO_CARBON && !fragment_has_carbon(frag)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum FragmentClass {
    Kept,
    AlwaysStripWater,
    AlwaysStripIon,
    StructuralSalt,
}

struct ScoredFragment {
    heavy: usize,
    rank_size: usize,
    has_carbon: bool,
    canon: String,
    class: FragmentClass,
}

fn score_fragments(frags: &[Molecule], policy: &FragmentPolicy) -> Vec<ScoredFragment> {
    frags
        .iter()
        .map(|fragment| {
            let heavy = heavy_atom_count(fragment);
            ScoredFragment {
                heavy,
                rank_size: if policy.count_heavy_atoms_only {
                    heavy
                } else {
                    fragment.atom_count()
                },
                has_carbon: fragment_has_carbon(fragment),
                canon: chematic_smiles::canonical_smiles(fragment),
                class: classify_fragment(fragment, policy),
            }
        })
        .collect()
}

fn fragment_rank_key<'a>(
    score: &'a ScoredFragment,
    policy: &FragmentPolicy,
) -> (usize, bool, std::cmp::Reverse<&'a str>) {
    (
        score.rank_size,
        policy.prefer_organic && score.has_carbon,
        std::cmp::Reverse(score.canon.as_str()),
    )
}

fn select_fragment_winner(
    scored: &[ScoredFragment],
    policy: &FragmentPolicy,
) -> (usize, Option<String>) {
    let kept = (0..scored.len()).filter(|&i| scored[i].class == FragmentClass::Kept);
    if let Some(winner) = kept.max_by_key(|&i| fragment_rank_key(&scored[i], policy)) {
        return (winner, None);
    }

    let winner = (0..scored.len())
        .max_by_key(|&i| fragment_rank_key(&scored[i], policy))
        .expect("non-empty molecule has at least one fragment");
    let reason = format!(
        "no_fragment_classified_kept: all {} fragment(s) matched a salt/ion/water \
         classification rule; falling back to the size-ranked choice among them",
        scored.len()
    );
    (winner, Some(reason))
}

fn removed_fragment_decision(score: &ScoredFragment) -> FragmentDecision {
    let (rule_id, reason) = match score.class {
        FragmentClass::AlwaysStripWater => (
            "water_always_strip",
            "single neutral oxygen fragment (water)".to_string(),
        ),
        FragmentClass::AlwaysStripIon => (
            "monatomic_always_strip_ion",
            "single charged atom on the fixed always-strip element list".to_string(),
        ),
        FragmentClass::StructuralSalt => (
            "structural_no_carbon_small_fragment",
            format!("no carbon and <= {MAX_SALT_HEAVY_ATOMS_NO_CARBON} heavy atoms"),
        ),
        FragmentClass::Kept => (
            "ranked_out_by_size",
            "outranked by the kept fragment under the fragment policy".to_string(),
        ),
    };
    FragmentDecision::Removed {
        rule_id: rule_id.to_string(),
        reason,
    }
}

fn fragment_records(
    frags: &[Molecule],
    scored: &[ScoredFragment],
    winner: usize,
) -> Vec<FragmentRecord> {
    let winner_canon = &scored[winner].canon;
    frags
        .iter()
        .enumerate()
        .map(|(i, fragment)| {
            let snapshot = FragmentSnapshot {
                atom_count: fragment.atom_count(),
                heavy_atom_count: scored[i].heavy,
                formula: fragment.formula(),
                canonical_smiles: scored[i].canon.clone(),
            };
            let decision = if i == winner {
                FragmentDecision::Kept {
                    rank_key: format!(
                        "rank_size={}, has_carbon={}, canonical_smiles={}",
                        scored[i].rank_size, scored[i].has_carbon, scored[i].canon
                    ),
                }
            } else if scored[i].canon == *winner_canon {
                FragmentDecision::Removed {
                    rule_id: "duplicate_of_kept_fragment".to_string(),
                    reason: "identical structure to the kept fragment (same canonical SMILES)"
                        .to_string(),
                }
            } else {
                removed_fragment_decision(&scored[i])
            };
            FragmentRecord { snapshot, decision }
        })
        .collect()
}

fn classify_fragment(frag: &Molecule, policy: &FragmentPolicy) -> FragmentClass {
    if policy.preserve_isotopes && fragment_has_isotope(frag) {
        return FragmentClass::Kept;
    }
    if is_water_fragment(frag) {
        FragmentClass::AlwaysStripWater
    } else if is_always_strip_monatomic_ion(frag) {
        FragmentClass::AlwaysStripIon
    } else if is_structural_salt_fallback(frag) {
        FragmentClass::StructuralSalt
    } else {
        FragmentClass::Kept
    }
}

/// Select which connected fragment of `mol` to keep, per `policy`, returning
/// both the kept fragment and a full per-fragment audit record.
///
/// Ranking is a pure function of each fragment's own structure — heavy-atom
/// count, carbon presence, and (as an intrinsic tie-break) canonical SMILES
/// — never of atom index or input spelling order, so two differently-spelled
/// inputs describing the same set of fragments always select the same
/// output. See `docs/rfcs/explainable_standardization_phase1_rfc.md`.
pub fn select_fragment(
    mol: &Molecule,
    policy: &FragmentPolicy,
) -> (Molecule, TransformationRecord) {
    let before = MoleculeSnapshot::from_mol(mol);

    if mol.atom_count() == 0 {
        let empty = MoleculeBuilder::new().build();
        let after = MoleculeSnapshot::from_mol(&empty);
        return (
            empty,
            TransformationRecord {
                rule_id: "fragment_policy_v1".to_string(),
                rule_version: 1,
                fragments: Vec::new(),
                abstained: None,
                before,
                after,
                warnings: Vec::new(),
            },
        );
    }

    let components = connected_components(mol);
    let frags: Vec<Molecule> = components
        .iter()
        .map(|c| extract_fragment(mol, c))
        .collect();

    let scored = score_fragments(&frags, policy);
    let (winner, abstained) = select_fragment_winner(&scored, policy);
    let mut warnings = Vec::new();
    if let Some(reason) = &abstained {
        warnings.push(StandardizationWarning::new(
            "fragment_selection_abstained",
            reason.clone(),
        ));
    }

    let fragments = fragment_records(&frags, &scored, winner);

    let output = frags[winner].clone();
    let after = MoleculeSnapshot::from_mol(&output);

    (
        output,
        TransformationRecord {
            rule_id: "fragment_policy_v1".to_string(),
            rule_version: 1,
            fragments,
            abstained,
            before,
            after,
            warnings,
        },
    )
}

// ---------------------------------------------------------------------------
// Parent identity (ROADMAP.md Phase 2 round 2B) -- see
// docs/rfcs/tautomer_parent_identity_phase2_rfc.md section 4.3.
// ---------------------------------------------------------------------------

/// [`select_fragment`] with the default [`FragmentPolicy`] -- the "Parent"
/// framing of Phase 1's fragment-selection logic.
pub fn fragment_parent(mol: &Molecule) -> (Molecule, TransformationRecord) {
    select_fragment(mol, &FragmentPolicy::default())
}

/// Select the fragment parent, then neutralize all formal charges on it.
///
/// **Not** the same as [`neutralize_charges`], which neutralizes every
/// fragment in a possibly-multi-fragment molecule and leaves them all in
/// the output. `charge_parent` first resolves fragment-selection ambiguity
/// down to one representative structure, matching this module's definition
/// of a Parent: an idempotent reduction to *one* structure, not a
/// multi-fragment mechanical transform.
///
/// The returned record's `fragments` (inherited from `fragment_parent`)
/// describe each candidate fragment **as selected**, before neutralization
/// -- a record of what was chosen and why. Only the record's own `after`
/// snapshot is updated to reflect the neutralized result.
pub fn charge_parent(mol: &Molecule) -> (Molecule, TransformationRecord) {
    let (selected, mut record) = fragment_parent(mol);
    let neutralized = neutralize_charges(&selected);
    record.rule_id = "charge_parent_v1".to_string();
    record.after = MoleculeSnapshot::from_mol(&neutralized);
    (neutralized, record)
}

/// [`remove_isotopes`] with an explainable audit record.
///
/// `fragments` is always empty: no fragment-selection decision is made
/// here (unlike `fragment_parent`/`charge_parent`), only isotope removal
/// on the whole input molecule.
pub fn isotope_parent(mol: &Molecule) -> (Molecule, TransformationRecord) {
    let before = MoleculeSnapshot::from_mol(mol);
    let output = remove_isotopes(mol);
    let after = MoleculeSnapshot::from_mol(&output);
    (
        output,
        TransformationRecord {
            rule_id: "isotope_parent_v1".to_string(),
            rule_version: 1,
            fragments: Vec::new(),
            abstained: None,
            before,
            after,
            warnings: Vec::new(),
        },
    )
}

/// [`remove_stereo`] with an explainable audit record.
///
/// `fragments` is always empty -- see [`isotope_parent`]'s doc comment.
pub fn stereo_parent(mol: &Molecule) -> (Molecule, TransformationRecord) {
    let before = MoleculeSnapshot::from_mol(mol);
    let output = remove_stereo(mol);
    let after = MoleculeSnapshot::from_mol(&output);
    (
        output,
        TransformationRecord {
            rule_id: "stereo_parent_v1".to_string(),
            rule_version: 1,
            fragments: Vec::new(),
            abstained: None,
            before,
            after,
            warnings: Vec::new(),
        },
    )
}

/// Extract a connected fragment (given its atom set) as a standalone
/// `Molecule`, preserving stereo data.
///
/// Built via repeated [`Molecule::remove_atom`] rather than a fresh
/// `MoleculeBuilder` + manual atom/bond remap: `remove_atom` already
/// correctly remaps `stereo_neighbor_order`/`bond_directions`/`stereo_groups`
/// (see its own doc comment), which a from-scratch builder copy silently
/// drops, flipping `@`/`@@` on any stereocenter-bearing fragment. Atoms not
/// in `component` are removed in descending original-index order, which
/// keeps every not-yet-removed index — including every atom still to be
/// kept or removed — stable until it is itself removed (each `remove_atom`
/// call only shifts indices *above* the one removed).
fn extract_fragment(mol: &Molecule, component: &[AtomIdx]) -> Molecule {
    let keep: std::collections::HashSet<u32> = component.iter().map(|a| a.0).collect();
    let mut result = mol.clone();
    let mut to_remove: Vec<u32> = (0..mol.atom_count() as u32)
        .filter(|i| !keep.contains(i))
        .collect();
    to_remove.sort_unstable_by(|a, b| b.cmp(a));
    for idx in to_remove {
        result.remove_atom(AtomIdx(idx));
    }
    result
}

/// Detect if a molecule contains a zwitterion (internal salt).
///
/// A zwitterion is defined as having both positive and negative formal charges.
/// Examples: amino acids in zwitterionic form ([NH3+][COO-]).
pub fn has_zwitterion(mol: &Molecule) -> bool {
    let mut has_positive = false;
    let mut has_negative = false;

    for (_, atom) in mol.atoms() {
        if atom.charge > 0 {
            has_positive = true;
        } else if atom.charge < 0 {
            has_negative = true;
        }
        if has_positive && has_negative {
            return true;
        }
    }
    false
}

/// Normalize a zwitterion to neutral form by proton transfer.
///
/// For each negatively-charged atom, find the nearest positively-charged atom.
/// A proton is transferred (donor's H count -1, acceptor's H count +1, both
/// charges move one step toward neutral) only if the chosen positive atom
/// actually has a hydrogen to give -- `has_zwitterion`'s "some + and some -
/// charge exists somewhere" check is necessary but not sufficient for a real
/// protonation-state zwitterion (amino-acid-style `[NH3+]`.../`[COO-]`...);
/// permanently-charge-separated groups with no available proton on either
/// side (e.g. a diazo-`N,N'`-dioxide, `[N+]([O-])=[N+]...[O-]`) are left
/// untouched atom-by-atom, never invented from nowhere. Every atom not
/// involved in an actual transfer keeps its original charge and H count
/// exactly -- atom count, per-element counts, and isotopes are always
/// preserved (see issue #407).
///
/// Returns a new molecule with normalized charges.
pub fn normalize_zwitterion(mol: &Molecule) -> Molecule {
    if !has_zwitterion(mol) {
        return clone_molecule(mol);
    }

    let mut modifications: HashMap<AtomIdx, (i8, Option<u8>)> = HashMap::new();

    // Collect positive and negative charge atoms
    let mut positive_atoms: Vec<AtomIdx> = Vec::new();
    let mut negative_atoms: Vec<AtomIdx> = Vec::new();

    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let atom = mol.atom(idx);
        if atom.charge > 0 {
            positive_atoms.push(idx);
        } else if atom.charge < 0 {
            negative_atoms.push(idx);
        }
    }

    // For each negative atom, transfer proton from nearest positive atom
    for &neg_idx in &negative_atoms {
        if positive_atoms.is_empty() {
            continue;
        }

        // Find nearest positive charge (by BFS distance)
        let mut closest_pos_idx = positive_atoms[0];
        let mut closest_distance = i32::MAX;

        for &pos_idx in &positive_atoms {
            if let Some(dist) = bfs_distance(mol, neg_idx, pos_idx)
                && dist < closest_distance
            {
                closest_distance = dist;
                closest_pos_idx = pos_idx;
            }
        }

        // Transfer proton: N+ loses H, O- gains H. Both sides of the pair must
        // move together or not at all -- a positive atom with no available H
        // (e.g. the fully-substituted N of a diazo-N,N'-dioxide, which is not
        // a real protonation-state zwitterion) has no proton to give, so the
        // pair is left untouched entirely. Previously the negative atom was
        // neutralized unconditionally regardless of whether the positive atom
        // actually donated anything, inventing a hydrogen and an unbalanced
        // charge/atom-count change from nowhere (issue #407).
        let pos_atom = mol.atom(closest_pos_idx);
        let pos_h = pos_atom.hydrogen_count.unwrap_or(0);
        if pos_h == 0 {
            continue;
        }

        let neg_atom = mol.atom(neg_idx);
        let new_neg_charge = neg_atom.charge + 1;
        let neg_h = neg_atom.hydrogen_count.unwrap_or(0);
        modifications.insert(neg_idx, (new_neg_charge, Some(neg_h + 1)));

        let new_pos_charge = pos_atom.charge - 1;
        modifications.insert(closest_pos_idx, (new_pos_charge, Some(pos_h - 1)));
    }

    // Reconstruct molecule with modified charges
    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();

    for i in 0..mol.atom_count() {
        let old_idx = AtomIdx(i as u32);
        let mut atom = mol.atom(old_idx).clone();
        if let Some(&(new_charge, new_h)) = modifications.get(&old_idx) {
            atom.charge = new_charge;
            atom.hydrogen_count = new_h;
        }
        let new_idx = builder.add_atom(atom);
        remap.insert(old_idx, new_idx);
    }
    copy_bonds(mol, &mut builder, &remap);
    // Identity-preserving rebuild (every atom carried forward at the same
    // relative order, remap is the identity) -- see issue #399.
    builder.copy_stereo_from(mol);
    builder.copy_bond_directions_from(mol);
    builder.copy_stereo_groups_from(mol);
    builder.build()
}

/// BFS distance between two atoms in a molecule.
fn bfs_distance(mol: &Molecule, start: AtomIdx, end: AtomIdx) -> Option<i32> {
    if start == end {
        return Some(0);
    }

    let n = mol.atom_count();
    let mut visited = vec![false; n];
    let mut queue = std::collections::VecDeque::new();
    queue.push_back((start, 0));
    visited[start.0 as usize] = true;

    while let Some((current, dist)) = queue.pop_front() {
        for (neighbor, _) in mol.neighbors(current) {
            if neighbor == end {
                return Some(dist + 1);
            }
            let ni = neighbor.0 as usize;
            if !visited[ni] {
                visited[ni] = true;
                queue.push_back((neighbor, dist + 1));
            }
        }
    }
    None
}

/// Neutralize simple formal charges in a molecule.
///
/// Rules applied:
/// - `[O-]` with a carbon neighbor → charge 0, +1 H (carboxylate → carboxylic acid).
/// - `[N+]` with at least one explicit H → charge 0, −1 H (ammonium → amine).
/// - `[O+]` with at least one explicit H → charge 0, −1 H (protonated ether → ether).
pub fn neutralize_charges(mol: &Molecule) -> Molecule {
    let mut modifications: HashMap<AtomIdx, (i8, Option<u8>)> = HashMap::new();

    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let atom = mol.atom(idx);
        let h = atom.hydrogen_count.unwrap_or(0);

        match (atom.element, atom.charge) {
            (Element::O, -1) => {
                let has_c_neighbor = mol
                    .neighbors(idx)
                    .any(|(nb, _)| mol.atom(nb).element == Element::C);
                if has_c_neighbor {
                    modifications.insert(idx, (0, Some(h + 1)));
                }
            }
            (Element::N, 1) | (Element::O, 1) if h > 0 => {
                modifications.insert(idx, (0, Some(h - 1)));
            }
            _ => {}
        }
    }

    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();
    for i in 0..mol.atom_count() {
        let old_idx = AtomIdx(i as u32);
        let mut atom = mol.atom(old_idx).clone();
        if let Some(&(new_charge, new_h)) = modifications.get(&old_idx) {
            atom.charge = new_charge;
            atom.hydrogen_count = new_h;
        }
        let new_idx = builder.add_atom(atom);
        remap.insert(old_idx, new_idx);
    }
    copy_bonds(mol, &mut builder, &remap);
    // Identity-preserving rebuild -- see issue #399.
    builder.copy_stereo_from(mol);
    builder.copy_bond_directions_from(mol);
    builder.copy_stereo_groups_from(mol);
    builder.build()
}

/// Remove all isotope labels from atoms.
///
/// Returns a new molecule with `atom.isotope = None` for all atoms.
pub fn remove_isotopes(mol: &Molecule) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();

    for i in 0..mol.atom_count() {
        let old_idx = AtomIdx(i as u32);
        let mut atom = mol.atom(old_idx).clone();
        atom.isotope = None;
        let new_idx = builder.add_atom(atom);
        remap.insert(old_idx, new_idx);
    }
    copy_bonds(mol, &mut builder, &remap);
    // Identity-preserving rebuild -- see issue #399.
    builder.copy_stereo_from(mol);
    builder.copy_bond_directions_from(mol);
    builder.copy_stereo_groups_from(mol);
    builder.build()
}

/// Remove all stereochemistry from a molecule.
///
/// Sets `atom.chirality = Chirality::None` for all atoms and converts
/// wedge/wedge-hash bonds to single bonds.
pub fn remove_stereo(mol: &Molecule) -> Molecule {
    use chematic_core::{BondOrder, Chirality};

    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();

    for i in 0..mol.atom_count() {
        let old_idx = AtomIdx(i as u32);
        let mut atom = mol.atom(old_idx).clone();
        atom.chirality = Chirality::None;
        let new_idx = builder.add_atom(atom);
        remap.insert(old_idx, new_idx);
    }

    for i in 0..mol.bond_count() {
        let bond = mol.bond(BondIdx(i as u32));
        if let (Some(&new_a), Some(&new_b)) = (remap.get(&bond.atom1), remap.get(&bond.atom2)) {
            let order = match bond.order {
                BondOrder::Up | BondOrder::Down => BondOrder::Single,
                other => other,
            };
            let _ = builder.add_bond(new_a, new_b, order);
        }
    }

    builder.build()
}

/// Remove atoms from stereo groups that are no longer chiral, and drop empty groups.
///
/// Analog of RDKit PR #9051 ("cleanup of stereogroups and wedges for non-chiral sites").
/// When molecular operations (fragment removal, chirality clearing, …) leave atoms
/// without chirality flags, their stereo group membership becomes invalid.  This
/// function filters each [`StereoGroup`]'s atom list to only those atoms where
/// `atom.chirality != Chirality::None`, and discards any group that becomes empty.
pub fn clean_stereo_groups(mol: &Molecule) -> Molecule {
    use chematic_core::{Chirality, StereoGroup};

    let cleaned = mol
        .stereo_groups()
        .iter()
        .enumerate()
        .filter_map(|(group_index, g)| {
            let chiral_atoms: Vec<AtomIdx> = g
                .atom_indices
                .iter()
                .copied()
                .filter(|&idx| mol.atom(idx).chirality != Chirality::None)
                .collect();
            // Bond members can be double-bond geometry or single-bond
            // atropisomer axes. This atom-centre cleanup has no basis for
            // removing either, so preserve every existing bond member.
            let stereo_bonds = mol
                .stereo_group_bond_indices(group_index)
                .unwrap_or_default()
                .to_vec();
            if chiral_atoms.is_empty() && stereo_bonds.is_empty() {
                None
            } else {
                Some((StereoGroup::new(g.kind.clone(), chiral_atoms), stereo_bonds))
            }
        })
        .collect();

    let mut out = MoleculeBuilder::from_molecule(mol).build();
    out.set_stereo_groups_with_bonds(cleaned);
    out
}

/// Keep only the largest organic (carbon-containing) fragment.
///
/// Removes all inorganic fragments (those without carbon atoms).
/// Useful for removing metal ions, salts, and other counterions.
/// Falls back to largest fragment if no organic fragment exists.
pub fn prefer_organic(mol: &Molecule) -> Molecule {
    if mol.atom_count() == 0 {
        return MoleculeBuilder::new().build();
    }

    let components = connected_components(mol);

    // Find largest organic fragment
    let mut largest_organic: Option<&Vec<AtomIdx>> = None;
    let mut largest_organic_size = 0;

    for component in &components {
        // Check if fragment contains carbon (organic)
        let has_carbon = component
            .iter()
            .any(|&idx| mol.atom(idx).element.atomic_number() == 6);

        if has_carbon && component.len() > largest_organic_size {
            largest_organic = Some(component);
            largest_organic_size = component.len();
        }
    }

    // Fall back to largest fragment if no organic found
    let target_component = largest_organic.or_else(|| components.first());

    if let Some(component) = target_component {
        // extract_fragment remaps stereo_neighbor_order/bond_directions/
        // stereo_groups correctly for a genuine atom-subset extraction; a
        // fresh MoleculeBuilder copy here would silently drop them (#399).
        extract_fragment(mol, component)
    } else {
        MoleculeBuilder::new().build()
    }
}

/// Move an existing negative charge to a stronger acid site when the acid/base
/// classes are supported. Neutral acids and bases are not independently ionized.
///
/// This bounded RDKit-compatible operation currently handles a single simple
/// carbanion or oxygen anion and unique oxygen-acid donors. Other charge and
/// resonance patterns are returned unchanged, with their stereo metadata intact.
pub fn reionize(mol: &Molecule) -> Molecule {
    reionize::reionize(mol)
}

/// Remove formal charges by protonation/deprotonation, as RDKit's
/// `Uncharger` does.
///
/// A cation loses its charge only when it has hydrogens to give up
/// (`[NH4+]` → `N`); one without (`[N+](C)(C)(C)C`, `[Na+]`, the N of a nitro
/// group) keeps it. An anion bonded to such a cation stays charge-separated
/// (nitro, N-oxide); other anions gain hydrogens, except as many as are
/// needed to balance the cations that keep their charge (`[Na+].[O-]C(=O)C`
/// stays a salt), and except anions whose neutral form would exceed a common
/// valence (`[B-](F)(F)(F)F`). Which balancing anions keep their charge
/// follows canonical atom order, so the result does not depend on input
/// atom order.
///
/// # Note
/// This differs from [`neutralize_charges`] which uses specific rules.
pub fn uncharge(mol: &Molecule) -> Molecule {
    use chematic_core::implicit_hcount;
    let n = mol.atom_count();
    let bond_sum = |idx: AtomIdx| -> i32 {
        mol.neighbors(idx)
            .map(|(_, b)| match mol.bond(b).order {
                chematic_core::BondOrder::Aromatic => 1,
                o => i32::from(o.order_int()),
            })
            .sum()
    };
    let max_neutral_valence = |z: u8| -> Option<i32> {
        match z {
            6 => Some(4),
            7 => Some(3),
            8 => Some(2),
            15 => Some(5),
            16 => Some(6),
            9 | 17 | 35 | 53 => Some(1),
            _ => None,
        }
    };
    let mut charge = vec![0i8; n];
    let mut hcount: Vec<Option<u8>> = vec![None; n];
    let mut keep = vec![false; n];
    for (idx, atom) in mol.atoms() {
        let i = idx.0 as usize;
        charge[i] = atom.charge;
        hcount[i] = atom.hydrogen_count;
    }
    // Cations: give up H where there is H to give.
    let mut kept_positive = 0i32;
    for (idx, atom) in mol.atoms() {
        let i = idx.0 as usize;
        if atom.charge <= 0 {
            continue;
        }
        let h = implicit_hcount(mol, idx);
        if i32::from(h) >= i32::from(atom.charge) {
            charge[i] = 0;
            if atom.hydrogen_count.is_some() || atom.aromatic {
                hcount[i] = Some(h - atom.charge as u8);
            }
        } else {
            keep[i] = true;
            kept_positive += i32::from(atom.charge);
        }
    }
    // Anions bonded to a kept cation stay charge-separated.
    let mut anions = Vec::new();
    for (idx, atom) in mol.atoms() {
        if atom.charge >= 0 {
            continue;
        }
        if mol.neighbors(idx).any(|(nb, _)| keep[nb.0 as usize]) {
            keep[idx.0 as usize] = true;
            kept_positive += i32::from(atom.charge);
        } else {
            anions.push(idx);
        }
    }
    // The remaining cationic charge is balanced by anions in canonical order.
    if kept_positive > 0 && !anions.is_empty() {
        let (_, order) = chematic_smiles::canonical_smiles_with_atom_order(mol);
        let mut rank = vec![usize::MAX; n];
        for (r, a) in order.iter().enumerate() {
            rank[a.0 as usize] = r;
        }
        anions.sort_by_key(|a| rank[a.0 as usize]);
    }
    for idx in anions {
        let i = idx.0 as usize;
        let atom = mol.atom(idx);
        let q = -i32::from(atom.charge);
        if kept_positive > 0 {
            kept_positive -= q;
            continue;
        }
        let h = i32::from(implicit_hcount(mol, idx));
        let fits = max_neutral_valence(atom.element.atomic_number())
            .is_some_and(|max| bond_sum(idx) + h + q <= max);
        if !fits {
            continue;
        }
        charge[i] = 0;
        if atom.hydrogen_count.is_some() || atom.aromatic {
            hcount[i] = Some((h + q) as u8);
        }
    }
    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();
    for (idx, atom) in mol.atoms() {
        let mut atom = atom.clone();
        atom.charge = charge[idx.0 as usize];
        atom.hydrogen_count = hcount[idx.0 as usize];
        let new_idx = builder.add_atom(atom);
        remap.insert(idx, new_idx);
    }
    copy_bonds(mol, &mut builder, &remap);
    // Identity-preserving rebuild -- see issue #399.
    builder.copy_stereo_from(mol);
    builder.copy_bond_directions_from(mol);
    builder.copy_stereo_groups_from(mol);
    builder.build()
}

/// Set every formal charge to zero, keeping all other atom data: a topology
/// key, not a chemical transformation (`CanonicalMode::Backbone`).
pub(crate) fn strip_charges(mol: &Molecule) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();
    for (idx, atom) in mol.atoms() {
        let mut atom = atom.clone();
        atom.charge = 0;
        remap.insert(idx, builder.add_atom(atom));
    }
    copy_bonds(mol, &mut builder, &remap);
    builder.copy_stereo_from(mol);
    builder.copy_bond_directions_from(mol);
    builder.copy_stereo_groups_from(mol);
    builder.build()
}

fn clone_molecule(mol: &Molecule) -> Molecule {
    MoleculeBuilder::from_molecule(mol).build()
}

fn is_metal(element: Element) -> bool {
    matches!(
        element.atomic_number(),
        3 | 4
            | 11 | 12 | 13
            | 19 | 20 | 21 | 22 | 23 | 24 | 25 | 26 | 27 | 28 | 29 | 30 | 31
            | 37 | 38 | 39 | 40 | 41 | 42 | 43 | 44 | 45 | 46 | 47 | 48 | 49 | 50
            | 55 | 56 | 57..=71
            | 72 | 73 | 74 | 75 | 76 | 77 | 78 | 79 | 80 | 81 | 82 | 83
            | 87 | 88 | 89..=103
            | 104 | 105 | 106 | 107 | 108 | 109 | 110 | 111 | 112 | 113 | 114 | 115 | 116
    )
}

mod pipeline;
pub use pipeline::{
    MoleculeSnapshot, PipelineStatus, StandardizationPipeline, StandardizationReport,
    StandardizationStep, StandardizationStepReport, StandardizationWarning, StandardizeOptions,
    ZwitterionHandling, standardize,
};

#[cfg(test)]
mod tests;

#[cfg(test)]
mod validation_boundary_tests;
