//! Molecular standardization routines.
//!
//! Provides utilities for cleaning up molecular representations:
//! - Selecting the largest connected fragment.
//! - Neutralizing simple formal charges.

#![forbid(unsafe_code)]

use std::collections::{HashMap, HashSet, VecDeque};

use chematic_core::{
    AtomIdx, BondIdx, Element, Molecule, MoleculeBuilder, valence_inferred_hcount, validate_valence,
};

use crate::{hash::mol_hash, hydrogen::remove_hydrogens, tautomer::canonical_tautomer};
use chematic_smarts::{MatchConfig, find_matches_with_config, parse_smarts};

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
                && !find_matches_with_config(&query, frag, &config).is_empty()
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
    select_fragment(mol, &FragmentPolicy::default()).0
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

    struct Scored {
        heavy: usize,
        rank_size: usize,
        has_carbon: bool,
        canon: String,
        class: FragmentClass,
    }

    let scored: Vec<Scored> = frags
        .iter()
        .map(|f| {
            let heavy = heavy_atom_count(f);
            let rank_size = if policy.count_heavy_atoms_only {
                heavy
            } else {
                f.atom_count()
            };
            Scored {
                heavy,
                rank_size,
                has_carbon: fragment_has_carbon(f),
                canon: chematic_smiles::canonical_smiles(f),
                class: classify_fragment(f, policy),
            }
        })
        .collect();

    let rank_key = |i: usize| -> (usize, bool, std::cmp::Reverse<String>) {
        let s = &scored[i];
        let carbon_key = policy.prefer_organic && s.has_carbon;
        (s.rank_size, carbon_key, std::cmp::Reverse(s.canon.clone()))
    };

    let kept_candidates: Vec<usize> = (0..frags.len())
        .filter(|&i| scored[i].class == FragmentClass::Kept)
        .collect();

    let (winner, abstained) = if kept_candidates.is_empty() {
        let w = (0..frags.len()).max_by_key(|&i| rank_key(i)).unwrap();
        (
            w,
            Some(format!(
                "no_fragment_classified_kept: all {} fragment(s) matched a salt/ion/water \
                 classification rule; falling back to the size-ranked choice among them",
                frags.len()
            )),
        )
    } else {
        let w = *kept_candidates
            .iter()
            .max_by_key(|&&i| rank_key(i))
            .unwrap();
        (w, None)
    };

    let winner_canon = scored[winner].canon.clone();
    let mut warnings = Vec::new();
    if let Some(reason) = &abstained {
        warnings.push(StandardizationWarning::new(
            "fragment_selection_abstained",
            reason.clone(),
        ));
    }

    let fragments: Vec<FragmentRecord> = frags
        .iter()
        .enumerate()
        .map(|(i, f)| {
            let snapshot = FragmentSnapshot {
                atom_count: f.atom_count(),
                heavy_atom_count: scored[i].heavy,
                formula: f.formula(),
                canonical_smiles: scored[i].canon.clone(),
            };
            let decision = if i == winner {
                FragmentDecision::Kept {
                    rank_key: format!(
                        "rank_size={}, has_carbon={}, canonical_smiles={}",
                        scored[i].rank_size, scored[i].has_carbon, scored[i].canon
                    ),
                }
            } else if scored[i].canon == winner_canon {
                FragmentDecision::Removed {
                    rule_id: "duplicate_of_kept_fragment".to_string(),
                    reason: "identical structure to the kept fragment (same canonical SMILES)"
                        .to_string(),
                }
            } else {
                let (rule_id, reason) = match scored[i].class {
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
            };
            FragmentRecord { snapshot, decision }
        })
        .collect();

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

/// Normalize chemical groups (nitro groups, etc.).
///
/// Transforms:
/// - `[N+](=O)[O-]` → `N(=O)=O` (nitro normalization: N charge 0, O- → double bond)
///
/// Returns a new molecule with normalized groups.
pub fn normalize_groups(mol: &Molecule) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();
    let mut nitro_atoms = std::collections::HashSet::new();
    let mut oxide_atoms = std::collections::HashSet::new();
    let mut azide_atoms = std::collections::HashSet::new();
    let mut sulfoxide_atoms = std::collections::HashSet::new();

    // First pass: identify functional groups via per-group detectors.
    for (idx, atom) in mol.atoms() {
        if atom.element.atomic_number() == 7 && atom.charge == 1 {
            detect_nitro(mol, idx, atom, &mut nitro_atoms, &mut oxide_atoms);
            detect_azide(mol, idx, &mut azide_atoms);
        }
        if atom.element.atomic_number() == 16 {
            detect_sulfoxide(mol, idx, &mut sulfoxide_atoms);
        }
    }

    // Second pass: copy atoms with normalized charges
    for (idx, atom) in mol.atoms() {
        let mut new_atom = atom.clone();

        if nitro_atoms.contains(&idx) {
            // Neutralize the N and O in nitro group
            if atom.element.atomic_number() == 7 || atom.element.atomic_number() == 8 {
                new_atom.charge = 0;
            }
        }

        if azide_atoms.contains(&idx) {
            // Neutralize all N in azide group [N-][N+]#N -> N=N=N
            if atom.element.atomic_number() == 7 {
                new_atom.charge = 0;
            }
        }

        // Sulfoxide: keep as is (S=O is already correct form)

        let new_idx = builder.add_atom(new_atom);
        remap.insert(idx, new_idx);
    }

    // Third pass: copy bonds, normalizing functional groups
    for i in 0..mol.bond_count() {
        let bond = mol.bond(chematic_core::BondIdx(i as u32));
        let mut new_order = bond.order;

        // Nitro groups: convert single N-O (where O is negative) to double
        if nitro_atoms.contains(&bond.atom1) && nitro_atoms.contains(&bond.atom2) {
            let a1_is_n = mol.atom(bond.atom1).element.atomic_number() == 7;
            let a2_is_o = mol.atom(bond.atom2).element.atomic_number() == 8;
            let a1_is_o = mol.atom(bond.atom1).element.atomic_number() == 8;
            let a2_is_n = mol.atom(bond.atom2).element.atomic_number() == 7;

            if (a1_is_n
                && a2_is_o
                && bond.order == chematic_core::BondOrder::Single
                && mol.atom(bond.atom2).charge == -1)
                || (a1_is_o
                    && a2_is_n
                    && bond.order == chematic_core::BondOrder::Single
                    && mol.atom(bond.atom1).charge == -1)
            {
                new_order = chematic_core::BondOrder::Double;
            }
        }

        // AZIDE normalization: [N-][N+]#N -> N=N=N (convert single to double)
        if azide_atoms.contains(&bond.atom1) && azide_atoms.contains(&bond.atom2) {
            let a1_is_n = mol.atom(bond.atom1).element.atomic_number() == 7;
            let a2_is_n = mol.atom(bond.atom2).element.atomic_number() == 7;

            if a1_is_n && a2_is_n && bond.order == chematic_core::BondOrder::Single {
                // Convert single bonds in azide to double
                new_order = chematic_core::BondOrder::Double;
            }
        }

        // N-oxide: keep as single bond (already correct after atom charge normalization)
        if oxide_atoms.contains(&bond.atom1) || oxide_atoms.contains(&bond.atom2) {
            // No bond order change needed — already single
        }

        // Sulfoxide: keep as is (S=O is already correct form)

        if let (Some(&new_a1), Some(&new_a2)) = (remap.get(&bond.atom1), remap.get(&bond.atom2)) {
            let _ = builder.add_bond(new_a1, new_a2, new_order);
        }
    }

    // Every atom above is copied 1:1 in original index order (never
    // filtered), so `remap` is the identity and bond indices are assigned in
    // the same relative order too -- `copy_stereo_from`/
    // `copy_bond_directions_from`'s verbatim-clone semantics are exactly
    // valid here. Without this, any stereocenter/declared-E-Z-bond
    // surviving nitro/azide/sulfoxide normalization silently lost its
    // `stereo_neighbor_order`/`bond_directions` entry -- issue #399's root
    // cause, the same defect class #392 fixed in `remove_hydrogens`.
    builder.copy_stereo_from(mol);
    builder.copy_bond_directions_from(mol);
    builder.copy_stereo_groups_from(mol);
    builder.build()
}

// ---------------------------------------------------------------------------
// Functional group detectors for normalize_groups
// ---------------------------------------------------------------------------

/// Mark nitro [N+](=O)[O-] and aromatic N-oxide atoms.
fn detect_nitro(
    mol: &Molecule,
    idx: AtomIdx,
    atom: &chematic_core::Atom,
    nitro_atoms: &mut std::collections::HashSet<AtomIdx>,
    oxide_atoms: &mut std::collections::HashSet<AtomIdx>,
) {
    let o_nbrs: Vec<_> = mol
        .neighbors(idx)
        .filter(|(n, _)| mol.atom(*n).element.atomic_number() == 8)
        .collect();

    if o_nbrs.len() == 2 {
        let mut has_double_o = false;
        let mut has_single_neg_o = false;
        for (o_idx, bid) in &o_nbrs {
            let o = mol.atom(*o_idx);
            let b = mol.bond(*bid);
            if b.order == chematic_core::BondOrder::Double && o.charge == 0 {
                has_double_o = true;
            }
            if b.order == chematic_core::BondOrder::Single && o.charge == -1 {
                has_single_neg_o = true;
                nitro_atoms.insert(*o_idx);
            }
        }
        if has_double_o && has_single_neg_o {
            nitro_atoms.insert(idx);
        }
    } else if let Some((o_idx, bid)) = o_nbrs.first() {
        let o = mol.atom(*o_idx);
        let b = mol.bond(*bid);
        if atom.aromatic && b.order == chematic_core::BondOrder::Single && o.charge == -1 {
            nitro_atoms.insert(idx);
            oxide_atoms.insert(*o_idx);
        }
    }
}

/// Mark azide [N-][N+]#N atoms.
fn detect_azide(
    mol: &Molecule,
    idx: AtomIdx,
    azide_atoms: &mut std::collections::HashSet<AtomIdx>,
) {
    let n_nbrs: Vec<_> = mol
        .neighbors(idx)
        .filter(|(n, _)| mol.atom(*n).element.atomic_number() == 7)
        .collect();

    for (n_idx, bid) in &n_nbrs {
        let n = mol.atom(*n_idx);
        let b = mol.bond(*bid);
        if b.order == chematic_core::BondOrder::Triple && n.charge == 0 {
            for (other_idx, other_bid) in n_nbrs.iter() {
                if other_idx == n_idx {
                    continue;
                }
                let other = mol.atom(*other_idx);
                let other_b = mol.bond(*other_bid);
                if other_b.order == chematic_core::BondOrder::Single && other.charge == -1 {
                    azide_atoms.insert(idx);
                    azide_atoms.insert(*n_idx);
                    azide_atoms.insert(*other_idx);
                }
            }
        }
    }
}

/// Mark sulfoxide S=O atoms.
fn detect_sulfoxide(
    mol: &Molecule,
    idx: AtomIdx,
    sulfoxide_atoms: &mut std::collections::HashSet<AtomIdx>,
) {
    for (o_idx, bid) in mol.neighbors(idx) {
        if mol.atom(o_idx).element.atomic_number() == 8
            && mol.bond(bid).order == chematic_core::BondOrder::Double
        {
            sulfoxide_atoms.insert(idx);
            sulfoxide_atoms.insert(o_idx);
        }
    }
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

    let cleaned: Vec<StereoGroup> = mol
        .stereo_groups()
        .iter()
        .filter_map(|g| {
            let chiral_atoms: Vec<AtomIdx> = g
                .atom_indices
                .iter()
                .copied()
                .filter(|&idx| mol.atom(idx).chirality != Chirality::None)
                .collect();
            if chiral_atoms.is_empty() {
                None
            } else {
                Some(StereoGroup::new(g.kind.clone(), chiral_atoms))
            }
        })
        .collect();

    let mut out = MoleculeBuilder::from_molecule(mol).build();
    out.set_stereo_groups(cleaned);
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

/// Reionize a molecule by adjusting protonation to favored forms.
///
/// Simple heuristic approach that adjusts charges on common acidic and basic groups:
/// - Carboxylic acids with OH: deprotonate to COO- (carboxylate)
/// - Phenols: deprotonate to phenoxide (O-)
/// - Primary amines: protonate to ammonium (NH3+)
/// - Imidazoles: protonate to imidazolium (N+)
///
/// This is a simplified version suitable for typical organic molecules.
pub fn reionize(mol: &Molecule) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();

    // Copy all atoms, adjusting charges
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let mut atom = mol.atom(idx).clone();
        let an = atom.element.atomic_number();

        // Check for carboxylic acid or phenol: C=O with O-H or Ar-O-H
        if an == 8 {
            // Oxygen: check if it's OH bonded to C
            if let Some((c_idx, _)) = mol.neighbors(idx).find(|(neighbor, bond_idx)| {
                mol.bond(*bond_idx).order == chematic_core::BondOrder::Single
                    && mol.atom(*neighbor).element.atomic_number() == 6
            }) {
                // Check if C is aromatic (phenol) or has a double-bonded O (carboxylic acid)
                let is_aromatic = mol.atom(c_idx).aromatic;
                let has_double_bonded_o = mol.neighbors(c_idx).any(|(other, bond_idx)| {
                    mol.bond(bond_idx).order == chematic_core::BondOrder::Double
                        && mol.atom(other).element.atomic_number() == 8
                        && other != idx
                });

                // Only deprotonate if it's a phenol or carboxylic acid, not aliphatic OH
                if (is_aromatic || has_double_bonded_o) && atom.charge >= 0 {
                    atom.charge -= 1; // Deprotonate: OH → O-
                }
            }
        }

        // Check for primary/secondary amines (but NOT amides)
        if an == 7 {
            // Check if this N is NOT part of an amide (C(=O)-N)
            let is_amide = mol.neighbors(idx).any(|(neighbor, bond_idx)| {
                mol.bond(bond_idx).order == chematic_core::BondOrder::Single
                    && mol.atom(neighbor).element.atomic_number() == 6
                    && mol.neighbors(neighbor).any(|(o_neighbor, o_bond)| {
                        mol.bond(o_bond).order == chematic_core::BondOrder::Double
                            && (mol.atom(o_neighbor).element.atomic_number() == 8
                                || mol.atom(o_neighbor).element.atomic_number() == 16)
                    })
            });

            if !is_amide {
                let h_count = chematic_core::implicit_hcount(mol, idx);
                // Protonate free amines only (not amides)
                if (h_count == 2 || h_count == 1) && atom.charge <= 0 {
                    atom.charge += 1; // Protonate: NH2 → NH3+
                }
            }
        }

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

/// Remove all charges from a molecule by protonation/deprotonation.
///
/// Neutralizes positively charged atoms by removing protons and
/// negatively charged atoms by adding protons. This is an aggressive
/// neutralization that may create chemically unrealistic structures.
///
/// # Note
/// This differs from [`neutralize_charges`] which uses specific rules.
/// `uncharge` is a brute-force approach suitable for structure cleanup.
pub fn uncharge(mol: &Molecule) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();

    // Copy all atoms, removing charges
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let mut atom = mol.atom(idx).clone();
        atom.charge = 0; // Force neutral
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

/// One transformation stage in the standardization pipeline.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum StandardizationStep {
    /// Remove metal coordination bonds before charge and fragment processing.
    DisconnectMetals,
    /// Select the largest connected component.
    LargestFragment,
    /// Apply simple neutralization rules for common formal charges.
    NeutralizeCharges,
    /// Normalize chemical groups (nitro groups, etc.).
    NormalizeGroups,
    /// Normalize zwitterionic forms to neutral molecules.
    ZwitterionNormalization,
    /// Remove explicit hydrogen atoms.
    RemoveExplicitHydrogens,
    /// Canonicalize supported tautomer systems.
    CanonicalTautomer,
    /// Keep only the largest non-salt fragment.
    FragmentParent,
    /// Neutralize all formal charges.
    ChargeParent,
    /// Remove all isotope labels.
    IsotopeParent,
    /// Remove all stereochemistry (wedge/wedge-hash bonds, chirality).
    StereoParent,
}

impl StandardizationStep {
    /// Stable machine-readable stage name.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::DisconnectMetals => "disconnect_metals",
            Self::LargestFragment => "largest_fragment",
            Self::NeutralizeCharges => "neutralize_charges",
            Self::NormalizeGroups => "normalize_groups",
            Self::ZwitterionNormalization => "zwitterion_normalization",
            Self::RemoveExplicitHydrogens => "remove_explicit_hydrogens",
            Self::CanonicalTautomer => "canonical_tautomer",
            Self::FragmentParent => "fragment_parent",
            Self::ChargeParent => "charge_parent",
            Self::IsotopeParent => "isotope_parent",
            Self::StereoParent => "stereo_parent",
        }
    }
}

/// High-level status for a standardization run.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum PipelineStatus {
    /// Pipeline completed and the output is structurally identical to the input.
    Unchanged,
    /// Pipeline completed and at least one enabled stage changed the molecule.
    Modified,
    /// Pipeline completed, but warnings indicate unsupported or suspicious input features.
    CompletedWithWarnings,
}

/// Warning emitted during standardization.
#[derive(Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct StandardizationWarning {
    /// Stable machine-readable warning code.
    pub code: String,
    /// Human-readable detail.
    pub message: String,
}

impl StandardizationWarning {
    fn new(code: &str, message: String) -> Self {
        Self {
            code: code.to_string(),
            message,
        }
    }
}

/// Atom/bond/hash summary before or after a pipeline stage.
#[derive(Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MoleculeSnapshot {
    /// Number of atoms in the molecule.
    pub atoms: usize,
    /// Number of bonds in the molecule.
    pub bonds: usize,
    /// Deterministic structure hash based on canonical SMILES.
    pub hash: u64,
}

impl MoleculeSnapshot {
    pub(crate) fn from_mol(mol: &Molecule) -> Self {
        Self {
            atoms: mol.atom_count(),
            bonds: mol.bond_count(),
            hash: mol_hash(mol),
        }
    }
}

/// Per-stage audit entry for a standardization run.
#[derive(Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct StandardizationStepReport {
    /// Pipeline stage.
    pub step: StandardizationStep,
    /// Whether the stage was enabled in the config.
    pub enabled: bool,
    /// Whether the stage changed the molecule hash.
    pub changed: bool,
    /// Molecule summary before the stage.
    pub before: MoleculeSnapshot,
    /// Molecule summary after the stage.
    pub after: MoleculeSnapshot,
}

/// Full standardization audit result.
#[derive(Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct StandardizationReport {
    /// Overall status.
    pub status: PipelineStatus,
    /// Summary of the input molecule.
    pub input: MoleculeSnapshot,
    /// Summary of the output molecule.
    pub output: MoleculeSnapshot,
    /// Ordered per-stage results.
    pub steps: Vec<StandardizationStepReport>,
    /// Validation and unsupported-feature warnings.
    pub warnings: Vec<StandardizationWarning>,
}

impl StandardizationReport {
    /// Returns `true` if the molecule changed at any enabled stage.
    pub fn changed(&self) -> bool {
        self.input.hash != self.output.hash
    }
}

/// Handling strategy for zwitterions (internal salts).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum ZwitterionHandling {
    /// Keep zwitterionic form as-is.
    Keep,
    /// Normalize to neutral form via proton transfer.
    #[default]
    Normalize,
}

/// Options for molecular standardization.
///
/// Controls which cleaning transformations are applied in a standardization pipeline.
#[derive(Clone, Debug)]
pub struct StandardizeOptions {
    /// Convert to canonical tautomer. Default: `true`.
    pub canonical_tautomer: bool,
    /// Neutralize simple formal charges. Default: `true`.
    pub neutralize_charges: bool,
    /// Remove explicit hydrogen atoms. Default: `true`.
    pub remove_explicit_h: bool,
    /// Keep only the largest connected fragment. Default: `false`.
    pub largest_fragment_only: bool,
    /// Handle zwitterions (internal salts). Default: `Normalize`.
    pub zwitterion_handling: ZwitterionHandling,
}

impl Default for StandardizeOptions {
    fn default() -> Self {
        Self {
            canonical_tautomer: true,
            neutralize_charges: true,
            remove_explicit_h: true,
            largest_fragment_only: false,
            zwitterion_handling: ZwitterionHandling::Normalize,
        }
    }
}

/// RDKit-style standardization pipeline with an auditable report.
#[derive(Clone, Debug, Default)]
pub struct StandardizationPipeline {
    options: StandardizeOptions,
}

impl StandardizationPipeline {
    /// Create a pipeline from explicit options.
    pub fn new(options: StandardizeOptions) -> Self {
        Self { options }
    }

    /// Borrow the pipeline options.
    pub fn options(&self) -> &StandardizeOptions {
        &self.options
    }

    /// Standardize a molecule and return both the output molecule and an audit report.
    pub fn run(&self, mol: &Molecule) -> (Molecule, StandardizationReport) {
        let input = MoleculeSnapshot::from_mol(mol);
        let mut current = clone_molecule(mol);
        let mut steps = Vec::new();
        let mut warnings = detect_initial_warnings(mol);

        // Disconnect metals early (remove dative/coordinate bonds). Keep this
        // implicit chemistry step in the report: otherwise the input snapshot
        // and the first reported stage describe different molecules.
        let has_metals = current.atoms().any(|(_, a)| is_metal(a.element));
        current = self.apply_stage(
            current,
            StandardizationStep::DisconnectMetals,
            has_metals,
            disconnect_metals,
            &mut steps,
            &mut warnings,
        );

        // Apply NeutralizeCharges BEFORE LargestFragment to ensure predictable fragment selection.
        // Example: [NH3+].[Cl-] should be neutralized first to [NH3].[Cl-], then largest fragment.
        current = self.apply_stage(
            current,
            StandardizationStep::NeutralizeCharges,
            self.options.neutralize_charges,
            neutralize_charges,
            &mut steps,
            &mut warnings,
        );
        current = self.apply_stage(
            current,
            StandardizationStep::LargestFragment,
            self.options.largest_fragment_only,
            largest_fragment,
            &mut steps,
            &mut warnings,
        );
        let zwitterion_enabled = self.options.zwitterion_handling == ZwitterionHandling::Normalize;
        current = self.apply_stage(
            current,
            StandardizationStep::ZwitterionNormalization,
            zwitterion_enabled,
            normalize_zwitterion,
            &mut steps,
            &mut warnings,
        );
        current = self.apply_stage(
            current,
            StandardizationStep::RemoveExplicitHydrogens,
            self.options.remove_explicit_h,
            remove_hydrogens,
            &mut steps,
            &mut warnings,
        );
        current = self.apply_stage(
            current,
            StandardizationStep::CanonicalTautomer,
            self.options.canonical_tautomer,
            canonical_tautomer,
            &mut steps,
            &mut warnings,
        );

        let output = MoleculeSnapshot::from_mol(&current);
        // Status depends only on structure change (hash), not on warnings.
        // Warnings are reported separately for the user to inspect.
        let status = if input.hash == output.hash {
            PipelineStatus::Unchanged
        } else if !warnings.is_empty() {
            PipelineStatus::CompletedWithWarnings
        } else {
            PipelineStatus::Modified
        };

        (
            current,
            StandardizationReport {
                status,
                input,
                output,
                steps,
                warnings,
            },
        )
    }

    fn apply_stage(
        &self,
        current: Molecule,
        step: StandardizationStep,
        enabled: bool,
        f: fn(&Molecule) -> Molecule,
        steps: &mut Vec<StandardizationStepReport>,
        warnings: &mut Vec<StandardizationWarning>,
    ) -> Molecule {
        let before = MoleculeSnapshot::from_mol(&current);
        let next = if enabled {
            f(&current)
        } else {
            clone_molecule(&current)
        };
        let after = MoleculeSnapshot::from_mol(&next);
        steps.push(StandardizationStepReport {
            step,
            enabled,
            changed: before.hash != after.hash,
            before,
            after,
        });
        if enabled {
            append_valence_warnings(step, &next, warnings);
        }
        next
    }
}

fn clone_molecule(mol: &Molecule) -> Molecule {
    MoleculeBuilder::from_molecule(mol).build()
}

fn detect_initial_warnings(mol: &Molecule) -> Vec<StandardizationWarning> {
    let mut warnings = Vec::new();
    // Metal disconnection is now handled in the pipeline, so we don't warn about it
    let valence_errors = validate_valence(mol);
    if !valence_errors.is_empty() {
        warnings.push(StandardizationWarning::new(
            "input_valence_validation_failed",
            format!(
                "input molecule has {} valence validation issue(s)",
                valence_errors.len()
            ),
        ));
    }
    warnings
}

fn append_valence_warnings(
    step: StandardizationStep,
    mol: &Molecule,
    warnings: &mut Vec<StandardizationWarning>,
) {
    let errors = validate_valence(mol);
    if errors.is_empty() {
        return;
    }
    warnings.push(StandardizationWarning::new(
        "valence_validation_failed",
        format!(
            "{} produced {} valence validation issue(s)",
            step.as_str(),
            errors.len()
        ),
    ));
}

/// Disconnect metal-nonmetal bonds by removing dative/coordinate bonds to metals.
///
/// Iterates through all atoms; if a metal is found, removes all bonds between
/// that metal and organic/inorganic atoms. Returns the molecule with metal
/// coordination bonds severed.
///
/// For every non-metal atom losing a bond this way, an explicitly-stored
/// `hydrogen_count` is recomputed by valence inference from the *surviving*
/// bonds -- a dative M-O/M-N bond is commonly written with a formal charge on
/// the non-metal atom that exactly balances the bond (e.g. `[O+]`
/// single-bonded to the metal, satisfying O+'s valence-3), and leaving the
/// stale, now-too-low H count in place would let the very next pipeline
/// stage, `neutralize_charges`, see `h == 0` (its guard is `h > 0`), skip
/// neutralizing it, and leave a dangling formal charge with no bond left to
/// justify it -- a real, confirmed idempotency bug (issue #403): the charge
/// only got neutralized on a *second* standardize pass, once a fresh parse of
/// the first pass's (incorrectly charged) output SMILES stored the H count
/// explicitly instead of inferring it. Recomputing here (not merely resetting
/// to `None` -- `neutralize_charges` reads the raw stored field, not
/// `implicit_hcount`, so an unresolved `None` would look identical to a
/// stored `Some(0)` to its guard) means `neutralize_charges` sees the true
/// post-disconnection valence on the very first pass and can correctly
/// neutralize (or correctly decide not to) right away, converging to the
/// same result every time.
fn disconnect_metals(mol: &Molecule) -> Molecule {
    // Atoms adjacent to a metal, on the *non-metal* side of the bond -- their
    // hydrogen_count needs recomputing once that bond is dropped.
    let mut affected: HashSet<AtomIdx> = HashSet::new();
    for i in 0..mol.bond_count() {
        let bond = mol.bond(BondIdx(i as u32));
        let atom1_is_metal = is_metal(mol.atom(bond.atom1).element);
        let atom2_is_metal = is_metal(mol.atom(bond.atom2).element);
        if atom1_is_metal && !atom2_is_metal {
            affected.insert(bond.atom2);
        } else if atom2_is_metal && !atom1_is_metal {
            affected.insert(bond.atom1);
        }
    }

    // Pass 1: copy atoms unchanged, drop metal bonds -- builds the
    // post-disconnection *topology* so valence inference below has the
    // surviving bonds (not the original, metal-including ones) to work from.
    let mut builder = MoleculeBuilder::new();
    for i in 0..mol.atom_count() {
        builder.add_atom(mol.atom(AtomIdx(i as u32)).clone());
    }
    for i in 0..mol.bond_count() {
        let old_bidx = BondIdx(i as u32);
        let bond = mol.bond(old_bidx);
        let atom1_is_metal = is_metal(mol.atom(bond.atom1).element);
        let atom2_is_metal = is_metal(mol.atom(bond.atom2).element);
        if !atom1_is_metal
            && !atom2_is_metal
            && let Ok(new_bidx) = builder.add_bond(bond.atom1, bond.atom2, bond.order)
            && let Some(direction) = mol.bond_direction(old_bidx)
        {
            builder.set_bond_direction(new_bidx, direction);
        }
    }
    let disconnected = builder.build();

    // Pass 2: recompute the affected atoms' hydrogen_count against the
    // now-disconnected topology (charge unchanged; `neutralize_charges` runs
    // next and decides whether the charge itself should move).
    let mut builder = MoleculeBuilder::new();
    for i in 0..disconnected.atom_count() {
        let idx = AtomIdx(i as u32);
        let mut atom = disconnected.atom(idx).clone();
        if affected.contains(&idx) && atom.hydrogen_count.is_some() {
            atom.hydrogen_count = Some(valence_inferred_hcount(&disconnected, idx));
        }
        builder.add_atom(atom);
    }
    for i in 0..disconnected.bond_count() {
        let bond = disconnected.bond(BondIdx(i as u32));
        if let Ok(new_bidx) = builder.add_bond(bond.atom1, bond.atom2, bond.order)
            && let Some(direction) = disconnected.bond_direction(BondIdx(i as u32))
        {
            builder.set_bond_direction(new_bidx, direction);
        }
    }
    builder.copy_stereo_from(&disconnected);
    builder.copy_stereo_groups_from(&disconnected);
    builder.build()
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

/// Apply a series of standardization steps to a molecule.
///
/// Transformations are applied in this order:
/// 1. If `largest_fragment_only`, select the largest connected component.
/// 2. If `neutralize_charges`, neutralize simple charges.
/// 3. If `remove_explicit_h`, remove explicit H atoms.
/// 4. If `canonical_tautomer`, convert to the canonical tautomer.
///
/// Useful for cleaning pasted structures or database entries.
///
/// With `canonical_tautomer` enabled (the default), a stereocenter's CIP
/// (R/S) label can legitimately change even though its real spatial
/// configuration never moves -- see [`crate::tautomer::canonical_tautomer`]'s
/// own doc comment (issue #402) before treating a pre- vs. post-standardize
/// CIP-label difference as a bug on its own.
pub fn standardize(mol: &Molecule, opts: &StandardizeOptions) -> Molecule {
    StandardizationPipeline::new(opts.clone()).run(mol).0
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn largest_fragment_two_fragments_picks_larger() {
        // "CC.CCC" — ethane (2 C) and propane (3 C)
        let mol = parse("CC.CCC").unwrap();
        let result = largest_fragment(&mol);
        assert_eq!(result.atom_count(), 3, "should keep propane (3 C)");
    }

    #[test]
    fn largest_fragment_single_fragment_unchanged() {
        // "CC" — ethane, only one fragment
        let mol = parse("CC").unwrap();
        let result = largest_fragment(&mol);
        assert_eq!(result.atom_count(), 2);
    }

    #[test]
    fn largest_fragment_keeps_benzene_over_ethane() {
        // "CC.c1ccccc1" — ethane (2 C) vs benzene (6 C)
        let mol = parse("CC.c1ccccc1").unwrap();
        let result = largest_fragment(&mol);
        assert_eq!(result.atom_count(), 6, "should keep benzene (6 atoms)");
    }

    #[test]
    fn largest_fragment_ionic_pair_keeps_one_atom() {
        // "[Na+].[Cl-]" — both fragments are single atoms; either is fine
        let mol = parse("[Na+].[Cl-]").unwrap();
        let result = largest_fragment(&mol);
        assert_eq!(result.atom_count(), 1);
    }

    #[test]
    fn neutralize_neutral_molecule_unchanged() {
        // "CC" is already neutral; no atom should gain/lose charge
        let mol = parse("CC").unwrap();
        let result = neutralize_charges(&mol);
        for i in 0..result.atom_count() {
            let atom = result.atom(AtomIdx(i as u32));
            assert_eq!(atom.charge, 0, "all atoms should remain neutral");
        }
    }

    #[test]
    fn neutralize_acetate_oxygen() {
        // "CC(=O)[O-]" — acetate; the [O-] should become neutral with H added
        let mol = parse("CC(=O)[O-]").unwrap();
        let result = neutralize_charges(&mol);

        // Find the oxygen that was originally [O-]: it should now have charge 0
        // and hydrogen_count == Some(1).
        let neutralized_o = (0..result.atom_count())
            .map(|i| result.atom(AtomIdx(i as u32)))
            .find(|a| a.element == Element::O && a.hydrogen_count == Some(1));

        assert!(
            neutralized_o.is_some(),
            "neutralized [O-] should have hydrogen_count == Some(1)"
        );
        assert_eq!(
            neutralized_o.unwrap().charge,
            0,
            "neutralized [O-] should have charge == 0"
        );
    }

    #[test]
    fn standardize_with_defaults() {
        // "CC(=O)[O-]" — acetate ion
        let mol = parse("CC(=O)[O-]").unwrap();
        let opts = StandardizeOptions::default();
        let result = standardize(&mol, &opts);

        // With default options, remove_explicit_h will be applied.
        // Check that [O-] was neutralized (charge should be 0).
        let has_neutral_o = (0..result.atom_count())
            .map(|i| result.atom(AtomIdx(i as u32)))
            .any(|a| a.element == Element::O && a.charge == 0);
        assert!(has_neutral_o, "acetate oxygen should be neutralized");

        // Should have at least 3 atoms (C, C, O) with no explicit H
        assert!(
            result.atom_count() >= 3,
            "should have at least 3 atoms after standardization"
        );
    }

    #[test]
    fn standardize_skip_largest_fragment() {
        // "CC.CCC" — ethane and propane
        let mol = parse("CC.CCC").unwrap();
        let opts = StandardizeOptions {
            zwitterion_handling: ZwitterionHandling::Normalize,
            largest_fragment_only: false,
            ..Default::default()
        };
        let result = standardize(&mol, &opts);

        // Should keep both fragments
        assert_eq!(
            result.atom_count(),
            5,
            "should keep both fragments when largest_fragment_only=false"
        );
    }

    #[test]
    fn pipeline_report_tracks_enabled_stage_changes() {
        let mol = parse("CC.CCC").unwrap();
        let pipeline = StandardizationPipeline::new(StandardizeOptions {
            largest_fragment_only: true,
            neutralize_charges: false,
            remove_explicit_h: false,
            canonical_tautomer: false,
            zwitterion_handling: ZwitterionHandling::Keep,
        });

        let (result, report) = pipeline.run(&mol);

        assert_eq!(result.atom_count(), 3);
        assert_eq!(report.status, PipelineStatus::Modified);
        assert!(report.changed());
        assert_eq!(report.steps.len(), 6);
        // Metal disconnection is the first (not enabled, so no change)
        assert_eq!(report.steps[0].step, StandardizationStep::DisconnectMetals);
        assert!(!report.steps[0].enabled);
        // NeutralizeCharges is applied second (not enabled, so no change)
        assert_eq!(report.steps[1].step, StandardizationStep::NeutralizeCharges);
        assert!(!report.steps[1].enabled);
        // LargestFragment is applied third and is enabled
        assert_eq!(report.steps[2].step, StandardizationStep::LargestFragment);
        assert!(report.steps[2].enabled);
        assert!(report.steps[2].changed);
    }

    #[test]
    fn pipeline_report_marks_unchanged_clean_molecule() {
        let mol = parse("CC").unwrap();
        let pipeline = StandardizationPipeline::new(StandardizeOptions {
            canonical_tautomer: false,
            neutralize_charges: false,
            remove_explicit_h: false,
            largest_fragment_only: false,
            zwitterion_handling: ZwitterionHandling::Keep,
        });

        let (_result, report) = pipeline.run(&mol);

        assert_eq!(report.status, PipelineStatus::Unchanged);
        assert!(!report.changed());
        assert!(report.warnings.is_empty());
        assert!(report.steps.iter().all(|s| !s.enabled && !s.changed));
    }

    #[test]
    fn pipeline_report_disconnects_metal_bonds() {
        let mol = parse("[Na]OC").unwrap();
        assert_eq!(mol.bond_count(), 2, "input has Na-O and O-C bonds");

        let pipeline = StandardizationPipeline::new(StandardizeOptions {
            canonical_tautomer: false,
            neutralize_charges: false,
            remove_explicit_h: false,
            largest_fragment_only: false,
            zwitterion_handling: ZwitterionHandling::Keep,
        });

        let (result, report) = pipeline.run(&mol);

        // Metal disconnection should run automatically, removing the Na-O bond
        assert_eq!(result.bond_count(), 1, "Na-O bond should be disconnected");
        // The remaining bond should be O-C
        assert!(
            result.bond(BondIdx(0)).atom1.0 < 3 && result.bond(BondIdx(0)).atom2.0 < 3,
            "remaining bond should connect organic atoms"
        );
        let metal_step = &report.steps[0];
        assert_eq!(metal_step.step, StandardizationStep::DisconnectMetals);
        assert!(metal_step.enabled);
        assert!(metal_step.changed);
        assert_eq!(metal_step.before.bonds, 2);
        assert_eq!(metal_step.after.bonds, 1);
    }

    #[test]
    fn bug3_ionic_pair_neutralize_before_largest_fragment() {
        // BUG #3 Fix Verification: NeutralizeCharges must run BEFORE LargestFragment
        // Example: [NH4+].[OH-] (ammonium hydroxide)
        // NH4+ will be neutralized (reduce H by 1)
        // After neutralization: [NH3].[OH-] - still two fragments
        // LargestFragment then picks the larger one
        // Key: stage order affects which fragment is selected
        let mol = parse("[NH4+].[OH-]").unwrap();
        let pipeline = StandardizationPipeline::new(StandardizeOptions {
            largest_fragment_only: true,
            neutralize_charges: true,
            remove_explicit_h: false,
            canonical_tautomer: false,
            zwitterion_handling: ZwitterionHandling::Normalize,
        });

        let (_result, report) = pipeline.run(&mol);

        // Verify step order: NeutralizeCharges MUST come before LargestFragment
        assert_eq!(report.steps.len(), 6, "Should have 6 steps in pipeline");
        assert_eq!(
            report.steps[0].step,
            StandardizationStep::DisconnectMetals,
            "DisconnectMetals must be step 0"
        );
        assert_eq!(
            report.steps[1].step,
            StandardizationStep::NeutralizeCharges,
            "NeutralizeCharges must be step 1"
        );
        assert_eq!(
            report.steps[2].step,
            StandardizationStep::LargestFragment,
            "LargestFragment must be step 2"
        );

        // The test passes if step order is correct and the pipeline runs without error
        assert!(report.changed(), "Pipeline should report changes");
        assert_eq!(
            report.status,
            PipelineStatus::Modified,
            "Should be marked as Modified"
        );
        assert!(
            report.warnings.is_empty(),
            "valid ammonium hydroxide should not warn"
        );
    }

    // ── Parent structure extraction tests ────────────────────────────────────

    #[test]
    fn remove_isotopes_strips_isotope_labels() {
        // "[13C]CC" has one 13C isotope label
        let mol = parse("[13C]CC").unwrap();
        let result = remove_isotopes(&mol);
        for i in 0..result.atom_count() {
            assert_eq!(
                result.atom(chematic_core::AtomIdx(i as u32)).isotope,
                None,
                "atom {} should have no isotope",
                i
            );
        }
    }

    #[test]
    fn remove_isotopes_preserves_structure() {
        // "[13C]CC" and "CC" should be structurally identical after isotope removal
        let mol = parse("[13C]CC").unwrap();
        let result = remove_isotopes(&mol);
        assert_eq!(result.atom_count(), 3, "atom count preserved");
        assert_eq!(result.bond_count(), 2, "bond count preserved");
    }

    #[test]
    fn remove_stereo_strips_chirality() {
        // "N[C@@H](C)C(=O)O" — alanine with (S) stereochemistry
        let mol = parse("N[C@@H](C)C(=O)O").unwrap();
        let result = remove_stereo(&mol);
        for i in 0..result.atom_count() {
            use chematic_core::Chirality;
            assert_eq!(
                result.atom(chematic_core::AtomIdx(i as u32)).chirality,
                Chirality::None,
                "atom {} should have no chirality",
                i
            );
        }
    }

    #[test]
    fn remove_stereo_converts_wedge_bonds_to_single() {
        // "C[C@H](O)C" — methylcarbinol with wedge stereochemistry
        // After remove_stereo, Up/Down bonds should become Single
        let mol = parse("C[C@H](O)C").unwrap();
        let result = remove_stereo(&mol);
        for i in 0..result.bond_count() {
            use chematic_core::BondOrder;
            let bond = result.bond(chematic_core::BondIdx(i as u32));
            assert_ne!(
                bond.order,
                BondOrder::Up,
                "bond {} should not be Up after stereo removal",
                i
            );
            assert_ne!(
                bond.order,
                BondOrder::Down,
                "bond {} should not be Down after stereo removal",
                i
            );
        }
    }

    #[test]
    fn remove_stereo_preserves_structure() {
        // "N[C@@H](C)C(=O)O" (alanine) should keep 6 atoms, 5 bonds after stereo removal
        let mol = parse("N[C@@H](C)C(=O)O").unwrap();
        let result = remove_stereo(&mol);
        assert_eq!(result.atom_count(), 6, "atom count preserved");
        assert_eq!(result.bond_count(), 5, "bond count preserved");
    }

    #[test]
    fn parent_variant_step_names_distinct() {
        // Verify that the 4 new parent variants have distinct step names
        let frag_parent = StandardizationStep::FragmentParent;
        let charge_parent = StandardizationStep::ChargeParent;
        let isotope_parent = StandardizationStep::IsotopeParent;
        let stereo_parent = StandardizationStep::StereoParent;

        assert_eq!(frag_parent.as_str(), "fragment_parent");
        assert_eq!(charge_parent.as_str(), "charge_parent");
        assert_eq!(isotope_parent.as_str(), "isotope_parent");
        assert_eq!(stereo_parent.as_str(), "stereo_parent");
    }

    // ── Cleanup transform tests (B3) ────────────────────────────────────────

    #[test]
    fn prefer_organic_removes_inorganic_salts() {
        // "CCO.[Na+].[Cl-]" — ethanol + sodium chloride
        let mol = parse("CCO.[Na+].[Cl-]").unwrap();
        assert_eq!(mol.atom_count(), 5, "input has CCO + Na + Cl");

        let result = prefer_organic(&mol);

        // Should keep only the organic ethanol fragment
        assert_eq!(result.atom_count(), 3, "should keep only ethanol (C, C, O)");
    }

    #[test]
    fn prefer_organic_keeps_organic_if_no_inorganic() {
        // "CC" — ethane only
        let mol = parse("CC").unwrap();
        let result = prefer_organic(&mol);
        assert_eq!(result.atom_count(), 2, "ethane unchanged");
    }

    #[test]
    fn prefer_organic_falls_back_to_largest() {
        // "C.C.C" — three separate carbons (all organic)
        let mol = parse("C.C.C").unwrap();
        let result = prefer_organic(&mol);
        // Should keep the largest organic fragment (any one of them, but all are size 1)
        assert_eq!(
            result.atom_count(),
            1,
            "falls back to largest fragment (one C)"
        );
    }

    #[test]
    fn uncharge_neutralizes_all_charges() {
        // "[NH4+].[OH-]" — ammonium hydroxide
        let mol = parse("[NH4+].[OH-]").unwrap();
        assert!(
            mol.atoms().any(|(_, a)| a.charge != 0),
            "input has charged atoms"
        );

        let result = uncharge(&mol);

        // All atoms should be neutral
        for (_, atom) in result.atoms() {
            assert_eq!(atom.charge, 0, "all atoms should be neutral");
        }
    }

    #[test]
    fn reionize_deprotonates_carboxylic_acids() {
        // "CC(=O)O" — acetic acid
        let mol = parse("CC(=O)O").unwrap();

        let result = reionize(&mol);

        // Should have a negatively charged oxygen (carboxylate anion)
        let has_negative_oxygen = result
            .atoms()
            .any(|(_, a)| a.element.atomic_number() == 8 && a.charge < 0);

        assert!(
            has_negative_oxygen,
            "reionize should deprotonate carboxylic acids"
        );
    }

    #[test]
    fn reionize_protonates_amines() {
        // "CC(N)C" — secondary amine
        let mol = parse("CC(N)C").unwrap();

        let result = reionize(&mol);

        // Should have a positively charged nitrogen (ammonium)
        let has_positive_nitrogen = result
            .atoms()
            .any(|(_, a)| a.element.atomic_number() == 7 && a.charge > 0);

        assert!(has_positive_nitrogen, "reionize should protonate amines");
    }

    #[test]
    fn reionize_protects_amide_nitrogen() {
        // BUG FIX #2: Amide nitrogen should NOT be protonated
        // "CC(=O)N" — primary amide
        let mol = parse("CC(=O)N").unwrap();

        let result = reionize(&mol);

        // Should NOT have a positively charged nitrogen
        let has_positive_nitrogen = result
            .atoms()
            .any(|(_, a)| a.element.atomic_number() == 7 && a.charge > 0);

        assert!(
            !has_positive_nitrogen,
            "reionize should NOT protonate amide nitrogen"
        );
    }

    #[test]
    fn reionize_protects_thioamide_nitrogen() {
        // BUG FIX #2: Thioamide nitrogen should also NOT be protonated
        // "CC(=S)N" — thioamide
        let mol = parse("CC(=S)N").unwrap();

        let result = reionize(&mol);

        // Should NOT have a positively charged nitrogen
        let has_positive_nitrogen = result
            .atoms()
            .any(|(_, a)| a.element.atomic_number() == 7 && a.charge > 0);

        assert!(
            !has_positive_nitrogen,
            "reionize should NOT protonate thioamide nitrogen (C=S conjugation)"
        );
    }

    // B2: normalize_groups expansion tests

    #[test]
    fn normalize_groups_nitro() {
        // Standard nitro group: [N+](=O)[O-]
        let mol = parse("C[N+](=O)[O-]").unwrap();
        let result = normalize_groups(&mol);

        // All atoms should be neutral after normalization
        let all_neutral = result.atoms().all(|(_, a)| a.charge == 0);
        assert!(all_neutral, "nitro group should be neutralized");

        // N-O bonds should be double
        let mut has_double_bond = false;
        for (_, bond) in result.bonds() {
            let a1 = result.atom(bond.atom1);
            let a2 = result.atom(bond.atom2);
            if ((a1.element.atomic_number() == 7 && a2.element.atomic_number() == 8)
                || (a1.element.atomic_number() == 8 && a2.element.atomic_number() == 7))
                && bond.order == chematic_core::BondOrder::Double
            {
                has_double_bond = true;
            }
        }
        assert!(has_double_bond, "nitro should have N=O double bond");
    }

    #[test]
    fn normalize_groups_azide() {
        // Azide: [N-][N+]#N
        let mol = parse("[N-][N+]#N").unwrap();
        let result = normalize_groups(&mol);

        // All atoms should be neutral
        let all_neutral = result.atoms().all(|(_, a)| a.charge == 0);
        assert!(all_neutral, "azide should be neutralized");

        // Check for N=N bonds (converted from single)
        let mut has_double_bond_count = 0;
        for (_, bond) in result.bonds() {
            let a1 = result.atom(bond.atom1);
            let a2 = result.atom(bond.atom2);
            if a1.element.atomic_number() == 7
                && a2.element.atomic_number() == 7
                && bond.order == chematic_core::BondOrder::Double
            {
                has_double_bond_count += 1;
            }
        }
        assert!(
            has_double_bond_count > 0,
            "azide should have N=N double bonds after normalization"
        );
    }

    #[test]
    fn normalize_groups_sulfoxide() {
        // Sulfoxide: S(=O)(C)(C)
        let mol = parse("C[S](=O)C").unwrap();
        let result = normalize_groups(&mol);

        // Sulfoxide structure should remain (S=O is already correct form)
        let mut has_s_double_o = false;
        for (_, bond) in result.bonds() {
            let a1 = result.atom(bond.atom1);
            let a2 = result.atom(bond.atom2);
            if ((a1.element.atomic_number() == 16 && a2.element.atomic_number() == 8)
                || (a1.element.atomic_number() == 8 && a2.element.atomic_number() == 16))
                && bond.order == chematic_core::BondOrder::Double
            {
                has_s_double_o = true;
            }
        }
        assert!(has_s_double_o, "sulfoxide should have S=O double bond");
    }

    // ── RDKit PR #9051: StereoGroup cleanup for non-chiral atoms ────────────

    #[test]
    fn remove_stereo_clears_stereo_groups() {
        // remove_stereo must produce a molecule with no stereo groups (RDKit PR #9051).
        // It rebuilds via MoleculeBuilder::new() so groups are already empty; this
        // test guards against regressions where from_molecule() is accidentally used.
        use chematic_core::{AtomIdx, StereoGroup, StereoGroupKind};
        let mut mol = parse("[C@@H](F)(Cl)Br").unwrap();
        mol.add_stereo_group(StereoGroup::new(
            StereoGroupKind::Absolute,
            vec![AtomIdx(0)],
        ));
        assert_eq!(
            mol.stereo_groups().len(),
            1,
            "precondition: group was added"
        );
        let stripped = remove_stereo(&mol);
        assert_eq!(
            stripped.stereo_groups().len(),
            0,
            "remove_stereo must clear stereo groups"
        );
    }

    #[test]
    fn clean_stereo_groups_drops_non_chiral_atoms() {
        // clean_stereo_groups filters out atoms with Chirality::None (RDKit PR #9051).
        use chematic_core::{AtomIdx, StereoGroup, StereoGroupKind};
        // atom 0 = CH3 (not chiral), atom 1 = @@ chiral center
        let mut mol = parse("C[C@@H](F)Cl").unwrap();
        mol.add_stereo_group(StereoGroup::new(
            StereoGroupKind::Absolute,
            vec![AtomIdx(0), AtomIdx(1)], // atom 0 is NOT chiral
        ));
        let cleaned = clean_stereo_groups(&mol);
        assert_eq!(
            cleaned.stereo_groups().len(),
            1,
            "group must survive with 1 atom"
        );
        assert_eq!(
            cleaned.stereo_groups()[0].atom_indices,
            vec![AtomIdx(1)],
            "only chiral atom must remain in group"
        );
    }

    #[test]
    fn clean_stereo_groups_drops_empty_groups() {
        // A group with no chiral atoms at all must be removed entirely.
        use chematic_core::{AtomIdx, StereoGroup, StereoGroupKind};
        let mut mol = parse("CC").unwrap(); // no chiral atoms
        mol.add_stereo_group(StereoGroup::new(
            StereoGroupKind::Absolute,
            vec![AtomIdx(0)],
        ));
        let cleaned = clean_stereo_groups(&mol);
        assert_eq!(
            cleaned.stereo_groups().len(),
            0,
            "empty group must be removed"
        );
    }

    #[test]
    fn clean_stereo_groups_preserves_valid_groups() {
        // A group whose atoms are all chiral must be kept intact.
        use chematic_core::{AtomIdx, StereoGroup, StereoGroupKind};
        let mut mol = parse("[C@@H](F)(Cl)Br").unwrap();
        mol.add_stereo_group(StereoGroup::new(
            StereoGroupKind::Absolute,
            vec![AtomIdx(0)],
        ));
        let cleaned = clean_stereo_groups(&mol);
        assert_eq!(
            cleaned.stereo_groups().len(),
            1,
            "valid group must be preserved"
        );
        assert_eq!(cleaned.stereo_groups()[0].atom_indices, vec![AtomIdx(0)]);
    }

    #[test]
    fn normalize_groups_mixed_nitro_and_azide() {
        // Molecule with both nitro and azide groups
        let mol = parse("C[N+](=O)[O-].N[N+](=O)[O-]").unwrap();
        let result = normalize_groups(&mol);

        // All atoms should be neutral
        let all_neutral = result.atoms().all(|(_, a)| a.charge == 0);
        assert!(all_neutral, "both nitro and azide should be neutralized");
    }

    // -- Phase 1 fragment-policy fixture tests --------------------------------
    // Mirrors validation/standardization_phase1_fixtures.jsonl and
    // validation/standardization_phase1_holdout.jsonl (ids referenced in each
    // test name/comment). See docs/rfcs/explainable_standardization_phase1_rfc.md.

    fn canon(s: &str) -> String {
        chematic_smiles::canonical_smiles(&parse(s).unwrap())
    }

    /// Assert `remove_salts(input)` keeps exactly the fragment `expected_kept`
    /// describes (compared by canonical SMILES, not string equality).
    fn assert_kept(id: &str, input: &str, expected_kept: &str) {
        let mol = parse(input).unwrap();
        let result = remove_salts(&mol);
        assert_eq!(
            chematic_smiles::canonical_smiles(&result),
            canon(expected_kept),
            "{id}: expected kept fragment {expected_kept:?} for input {input:?}"
        );
    }

    #[test]
    fn phase1_std_p1_01_sodium_acetate() {
        assert_kept("std-p1-01", "CC(=O)[O-].[Na+]", "CC(=O)[O-]");
    }

    #[test]
    fn phase1_std_p1_02_potassium_benzoate() {
        assert_kept("std-p1-02", "c1ccccc1C(=O)[O-].[K+]", "c1ccccc1C(=O)[O-]");
    }

    #[test]
    fn phase1_std_p1_03_amine_hcl_neutral_form() {
        assert_kept(
            "std-p1-03",
            "CC(C)NCC(O)COc1ccccc1CCOC.Cl",
            "CC(C)NCC(O)COc1ccccc1CCOC",
        );
    }

    #[test]
    fn phase1_std_p1_04_triethylamine_hbr() {
        assert_kept("std-p1-04", "CCN(CC)CC.Br", "CCN(CC)CC");
    }

    #[test]
    fn phase1_std_p1_05_sodium_octanoate_long_chain() {
        assert_kept("std-p1-05", "CCCCCCCC(=O)[O-].[Na+]", "CCCCCCCC(=O)[O-]");
    }

    #[test]
    fn phase1_std_p1_06_amine_hcl_hydrate_3_fragment() {
        assert_kept(
            "std-p1-06",
            "CC(C)NCC(O)c1ccc(O)c(CO)c1.Cl.O",
            "CC(C)NCC(O)c1ccc(O)c(CO)c1",
        );
    }

    #[test]
    fn phase1_std_p1_07_10_zwitterions_are_noop() {
        for (id, smi) in [
            ("std-p1-07", "NCC(=O)[O-]"),
            ("std-p1-08", "C[C@@H](N)C(=O)[O-]"),
            ("std-p1-09", "NCCS(=O)(=O)[O-]"),
            ("std-p1-10", "NCCCC(=O)[O-]"),
        ] {
            assert_kept(id, smi, smi);
        }
    }

    #[test]
    fn phase1_std_p1_11_13_hydrates() {
        assert_kept(
            "std-p1-11",
            "Cn1cnc2c1c(=O)n(C)c(=O)n2C.O",
            "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
        );
        assert_kept(
            "std-p1-12",
            "OCC(O)C(O)C(O)C(O)CO.O",
            "OCC(O)C(O)C(O)C(O)CO",
        );
        assert_kept("std-p1-13", "CC(=O)Nc1ccc(O)cc1.O", "CC(=O)Nc1ccc(O)cc1");
    }

    #[test]
    fn phase1_std_p1_14_anhydrous_negative_control() {
        // Same parent as std-p1-11 with no water fragment present: must be a
        // pure no-op, not strip anything from the one remaining fragment.
        assert_kept(
            "std-p1-14",
            "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
            "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
        );
    }

    #[test]
    fn phase1_std_p1_15_ferrocene_ionic_fragments() {
        // Fe2+ removed; either cyclopentadienide is an equally valid "kept"
        // choice since both are identical -- assert on the kept fragment's
        // OWN structure, not on which literal input position it came from.
        let mol = parse("[Fe+2].c1cc[cH-]c1.c1cc[cH-]c1").unwrap();
        let result = remove_salts(&mol);
        assert_eq!(
            chematic_smiles::canonical_smiles(&result),
            canon("c1cc[cH-]c1"),
            "std-p1-15: expected a cyclopentadienide fragment kept, Fe2+ removed"
        );
    }

    #[test]
    fn phase1_std_p1_16_ethylmagnesium_bromide() {
        assert_kept("std-p1-16", "CC[Mg+].[Br-]", "CC[Mg+]");
    }

    #[test]
    fn phase1_std_p1_17_cisplatin_single_fragment_noop() {
        // Bonded Cl/N ligands, not disconnected fragments -- must not be
        // misparsed as separate salt components.
        assert_kept("std-p1-17", "Cl[Pt](Cl)(N)N", "Cl[Pt](Cl)(N)N");
    }

    #[test]
    fn phase1_std_p1_18_amine_tartrate_salt() {
        assert_kept(
            "std-p1-18",
            "CN1CCC(CC1)Nc1ccccc1.OC(=O)C(O)C(O)C(=O)O",
            "CN1CCC(CC1)Nc1ccccc1",
        );
    }

    #[test]
    fn phase1_std_p1_19_two_api_cocrystal_ranked_by_size() {
        // Deliberate Phase 1 scope decision (see RFC section 6 / PR notes):
        // no ambiguity-margin detector is implemented -- a heavy-atom-count
        // margin alone cannot distinguish this genuinely-ambiguous cocrystal
        // from std-p1-holdout-10's unambiguous pentane/butane pair (both have
        // a 1-heavy-atom gap). This fixture's outcome is the same size-ranked
        // decision as any other multi-organic-fragment case: acetaminophen
        // (11 heavy atoms) outranks salicylic acid (10). The RFC's proposed
        // "flag close decisions" mechanism is a disclosed, not-yet-implemented
        // gap, not something this test asserts.
        assert_kept(
            "std-p1-19",
            "CC(=O)Nc1ccc(O)cc1.OC(=O)c1ccccc1O",
            "CC(=O)Nc1ccc(O)cc1",
        );
    }

    #[test]
    fn phase1_std_p1_20_22_isotope_containing() {
        assert_kept("std-p1-20", "[13CH3]C(=O)[O-].[Na+]", "[13CH3]C(=O)[O-]");
        assert_kept(
            "std-p1-21",
            "OC[C@H]1O[C@@H]([13OH])[C@H](O)[C@@H](O)[C@@H]1O.O",
            "OC[C@H]1O[C@@H]([13OH])[C@H](O)[C@@H](O)[C@@H]1O",
        );
        assert_kept("std-p1-22", "[13CH3][NH3+].[Cl-]", "[13CH3][NH3+]");
    }

    #[test]
    fn phase1_std_p1_23_tied_spelling_is_atom_order_invariant() {
        // REGRESSION for confirmed defect 1.1(a): pre-RFC largest_fragment()
        // picked a DIFFERENT fragment for these two spellings.
        let a = remove_salts(&parse("CCC.CCN").unwrap());
        let b = remove_salts(&parse("CCN.CCC").unwrap());
        assert_eq!(
            chematic_smiles::canonical_smiles(&a),
            chematic_smiles::canonical_smiles(&b),
            "std-p1-23: kept fragment must not depend on input spelling order"
        );
    }

    #[test]
    fn phase1_std_p1_24_heavy_atom_count_not_explicit_h_count() {
        // REGRESSION for confirmed defect 1.1(b).
        assert_kept("std-p1-24", "CCC.[H]C([H])([H])[H]", "CCC");
    }

    #[test]
    fn phase1_std_p1_25_benzene_hexane_tied_spelling_invariant() {
        let a = remove_salts(&parse("c1ccccc1.CCCCCC").unwrap());
        let b = remove_salts(&parse("CCCCCC.c1ccccc1").unwrap());
        assert_eq!(
            chematic_smiles::canonical_smiles(&a),
            chematic_smiles::canonical_smiles(&b),
            "std-p1-25: 6-heavy-atom tie must resolve the same way regardless of spelling"
        );
    }

    #[test]
    fn phase1_std_p1_26_27_charged_fragments() {
        assert_kept("std-p1-26", "CC(C)(C)[NH3+].[Br-]", "CC(C)(C)[NH3+]");
        assert_kept("std-p1-27", "c1cc[nH+]cc1.[Cl-]", "c1cc[nH+]cc1");
    }

    #[test]
    fn phase1_std_p1_28_choline_chloride_correct_for_the_right_reason() {
        // REGRESSION for confirmed defect 1.1(c): the legacy SaltCatalog's
        // "ammonium" SMARTS also classifies choline itself as salt (verified
        // empirically: SaltCatalog::default().is_salt() on choline alone
        // returns true), so the pre-RFC code only got this right via an
        // unrelated size-comparison fallback. The new structural policy
        // (no catalog consulted by default) must keep choline because it is
        // organic and heavy-atom-dominant, independent of any nitrogen-charge
        // pattern.
        assert_kept("std-p1-28", "C[N+](C)(C)CCO.[Cl-]", "C[N+](C)(C)CCO");
        assert!(
            SaltCatalog::default().is_salt(&parse("C[N+](C)(C)CCO").unwrap()),
            "sanity check: the legacy catalog's false positive on choline alone is still present \
             (expected -- SaltCatalog is intentionally untouched, kept as opt-in legacy behavior)"
        );
    }

    #[test]
    fn phase1_std_p1_29_30_no_organic_fragment_abstains_with_warning() {
        for (id, smi) in [("std-p1-29", "[Na+].[Cl-]"), ("std-p1-30", "[Cl-]")] {
            let mol = parse(smi).unwrap();
            let (_output, record) = select_fragment(&mol, &FragmentPolicy::default());
            assert!(
                record.abstained.is_some(),
                "{id}: expected an abstain reason when no fragment classifies as Kept"
            );
            assert!(
                !record.warnings.is_empty(),
                "{id}: expected a warning surfaced alongside the abstain reason"
            );
        }
    }

    #[test]
    fn phase1_std_p1_31_single_fragment_is_pure_noop() {
        let mol = parse("CCO").unwrap();
        let (output, record) = select_fragment(&mol, &FragmentPolicy::default());
        assert_eq!(chematic_smiles::canonical_smiles(&output), canon("CCO"));
        assert_eq!(record.fragments.len(), 1);
        assert!(record.abstained.is_none());
        assert!(matches!(
            record.fragments[0].decision,
            FragmentDecision::Kept { .. }
        ));
    }

    #[test]
    fn phase1_std_p1_32_empty_molecule_completes_not_errors() {
        let empty = MoleculeBuilder::new().build();
        let (output, record) = select_fragment(&empty, &FragmentPolicy::default());
        assert_eq!(output.atom_count(), 0);
        assert!(record.fragments.is_empty());
        assert!(record.abstained.is_none());
    }

    #[test]
    fn phase1_std_p1_33_duplicate_fragment_genuine_tie() {
        let mol = parse("CC(=O)[O-].[Na+].CC(=O)[O-].[Na+]").unwrap();
        let (output, record) = select_fragment(&mol, &FragmentPolicy::default());
        assert_eq!(
            chematic_smiles::canonical_smiles(&output),
            canon("CC(=O)[O-]")
        );
        let kept_count = record
            .fragments
            .iter()
            .filter(|f| matches!(f.decision, FragmentDecision::Kept { .. }))
            .count();
        assert_eq!(
            kept_count, 1,
            "exactly one fragment is reported Kept even though a duplicate exists"
        );
        let duplicate_marked = record.fragments.iter().any(|f| {
            matches!(
                &f.decision,
                FragmentDecision::Removed { rule_id, .. } if rule_id == "duplicate_of_kept_fragment"
            )
        });
        assert!(
            duplicate_marked,
            "the identical second acetate must be recorded as a duplicate, not a ranked-out loser"
        );
    }

    #[test]
    fn phase1_std_p1_34_carbon_containing_solvate_known_gap() {
        // Disclosed Phase 1 limitation: DMSO contains carbon, so it does NOT
        // match the no-carbon structural salt heuristic and is kept (and
        // must be recorded as Kept-by-ranking, not Removed-as-solvent).
        let mol = parse("c1ccc2c(c1)cccc2.CS(=O)C").unwrap();
        let (output, record) = select_fragment(&mol, &FragmentPolicy::default());
        assert_eq!(
            chematic_smiles::canonical_smiles(&output),
            canon("c1ccc2c(c1)cccc2")
        );
        let dmso_decision = record
            .fragments
            .iter()
            .find(|f| f.snapshot.heavy_atom_count == 4)
            .map(|f| f.decision.clone())
            .expect("DMSO fragment present");
        assert!(
            matches!(
                dmso_decision,
                FragmentDecision::Removed { ref rule_id, .. } if rule_id == "ranked_out_by_size"
            ),
            "DMSO must be recorded as ranked-out-by-size (Kept-class, just smaller), not as a matched salt/solvent rule"
        );
    }

    #[test]
    fn phase1_holdout_01_dopamine_wins_regardless_of_n_threshold() {
        // Dopamine (11 heavy atoms) vs H2SO4 (5 heavy atoms, S+4O -- corrected
        // from the RFC/holdout fixture's original "5 heavy atoms" note, which
        // is right) is not close enough for the open N threshold (RFC
        // section 6) to matter: dopamine wins by ranking regardless of
        // whether H2SO4 is classified StructuralSalt (excluded from
        // candidates) or Kept-but-outranked. This confirms the threshold
        // question only affects the audit log's rule_id here, not the
        // fragment-selection outcome.
        assert_kept(
            "std-p1-holdout-01",
            "CC(N)Cc1ccc(O)c(O)c1.OS(=O)(=O)O",
            "CC(N)Cc1ccc(O)c(O)c1",
        );
    }

    #[test]
    fn phase1_holdout_02_phosphoric_acid_outranks_small_amine_disclosed_gap() {
        // Correction to the original holdout fixture's arithmetic: H3PO4 has
        // 4 oxygens (P + 4*O = 5 heavy atoms), not 3 -- it is NOT a boundary
        // case at N=4, it is squarely over it, same bucket as H2SO4. Unlike
        // holdout-01, ethanolamine (NCCO, 4 heavy atoms) is SMALLER than
        // phosphoric acid (5 heavy atoms), so under a pure heavy-atom-count
        // policy with no named-acid recognition, phosphoric acid outranks the
        // amine and is kept. This is a disclosed limitation of a purely
        // structural, no-named-list policy -- not a bug -- and is exactly the
        // kind of case a future, carefully-scoped small always-strip
        // extension (common mineral acids) would need to address, per RFC
        // section 6.
        assert_kept("std-p1-holdout-02", "NCCO.OP(=O)(O)O", "OP(=O)(O)O");
    }

    #[test]
    fn phase1_holdout_03_calcium_diacetate_via_general_heuristic() {
        // Ca2+ is NOT on the tiny always-strip element list; it is still
        // removed, but via the general structural fallback rule, not
        // "monatomic_always_strip_ion" -- the audit log must say so.
        let mol = parse("CC(=O)[O-].CC(=O)[O-].[Ca+2]").unwrap();
        let (output, record) = select_fragment(&mol, &FragmentPolicy::default());
        assert_eq!(
            chematic_smiles::canonical_smiles(&output),
            canon("CC(=O)[O-]")
        );
        let ca_decision = record
            .fragments
            .iter()
            .find(|f| f.snapshot.heavy_atom_count == 1 && f.snapshot.formula.contains("Ca"))
            .map(|f| f.decision.clone())
            .expect("Ca2+ fragment present");
        assert!(
            matches!(
                ca_decision,
                FragmentDecision::Removed { ref rule_id, .. } if rule_id == "structural_no_carbon_small_fragment"
            ),
            "Ca2+ removal must be attributed to the general heuristic, not the always-strip-ion list"
        );
    }

    #[test]
    fn phase1_holdout_04_05_always_strip_coverage() {
        assert_kept("std-p1-holdout-04", "CC(=O)[O-].[Li+]", "CC(=O)[O-]");
        assert_kept("std-p1-holdout-05", "c1ccc(cc1)CCN.F", "c1ccc(cc1)CCN");
    }

    #[test]
    fn phase1_holdout_06_dihydrate_both_removed() {
        let mol = parse("Cn1cnc2c1c(=O)n(C)c(=O)n2C.O.O").unwrap();
        let (output, record) = select_fragment(&mol, &FragmentPolicy::default());
        assert_eq!(
            chematic_smiles::canonical_smiles(&output),
            canon("Cn1cnc2c1c(=O)n(C)c(=O)n2C")
        );
        let water_removed_count = record
            .fragments
            .iter()
            .filter(|f| matches!(&f.decision, FragmentDecision::Removed { rule_id, .. } if rule_id == "water_always_strip"))
            .count();
        assert_eq!(
            water_removed_count, 2,
            "both water fragments must be independently recorded"
        );
    }

    #[test]
    fn phase1_holdout_07_three_fragment_mixed_rationale() {
        let mol = parse("CN1CCC(CC1)Nc1ccccc1.OC(=O)c1ccccc1O.Cl").unwrap();
        let (output, record) = select_fragment(&mol, &FragmentPolicy::default());
        assert_eq!(
            chematic_smiles::canonical_smiles(&output),
            canon("CN1CCC(CC1)Nc1ccccc1")
        );
        assert_eq!(record.fragments.len(), 3);
        let rule_ids: Vec<&str> = record
            .fragments
            .iter()
            .filter_map(|f| match &f.decision {
                FragmentDecision::Removed { rule_id, .. } => Some(rule_id.as_str()),
                FragmentDecision::Kept { .. } => None,
            })
            .collect();
        assert!(rule_ids.contains(&"ranked_out_by_size"));
        assert!(rule_ids.contains(&"structural_no_carbon_small_fragment"));
    }

    #[test]
    fn phase1_holdout_08_isotope_sugar_hydrate_generalization() {
        assert_kept(
            "std-p1-holdout-08",
            "OC[C@H]1O[C@@H]([13OH])[C@H](O)[C@@H](O)[C@@H]1O.O",
            "OC[C@H]1O[C@@H]([13OH])[C@H](O)[C@@H](O)[C@@H]1O",
        );
    }

    #[test]
    fn phase1_holdout_09_10_non_tie_sanity_checks() {
        assert_kept("std-p1-holdout-09", "CC.CCCCCCCCCC", "CCCCCCCCCC");
        assert_kept("std-p1-holdout-10", "CCCCC.CCCC", "CCCCC");
    }

    // -- Phase 2 round-2B Parent-function fixture tests -----------------------
    // Mirrors validation/tautomer_parent_identity_phase2_fixtures.jsonl's
    // tp2-17..22. See docs/rfcs/tautomer_parent_identity_phase2_rfc.md
    // section 4.3.

    #[test]
    fn tp2_17_charge_parent_ammonium_acetate_single_fragment_result() {
        // charge_parent is NOT neutralize_charges: it selects the fragment
        // parent first (acetate: 4 heavy atoms, has carbon, over ammonium:
        // 1 heavy atom, no carbon), THEN neutralizes that one fragment.
        let mol = parse("CC(=O)[O-].[NH4+]").unwrap();
        let (result, record) = charge_parent(&mol);
        assert_eq!(chematic_smiles::canonical_smiles(&result), canon("CC(O)=O"));
        assert_eq!(record.rule_id, "charge_parent_v1");
        assert_eq!(
            record.fragments.len(),
            2,
            "inherited from fragment_parent's selection"
        );
    }

    #[test]
    fn tp2_18_charge_parent_zwitterion_amino_acid() {
        let mol = parse("[NH3+]CC(=O)[O-]").unwrap();
        let (result, _) = charge_parent(&mol);
        assert_eq!(
            chematic_smiles::canonical_smiles(&result),
            canon("NCC(O)=O")
        );
    }

    #[test]
    fn tp2_19_isotope_parent_deuterated_ethanol() {
        let mol = parse("[2H]C([2H])([2H])CO").unwrap();
        let (result, record) = isotope_parent(&mol);
        assert_eq!(
            chematic_smiles::canonical_smiles(&result),
            canon("[H]C(CO)([H])[H]")
        );
        assert!(record.fragments.is_empty());
    }

    #[test]
    fn tp2_20_isotope_parent_preserves_stereo() {
        // Expected value corrected under issue #399: `remove_isotopes`
        // rebuilt the molecule via a bare `MoleculeBuilder` without carrying
        // `stereo_neighbor_order` forward, so this test's original expected
        // string had the anomeric-adjacent ring stereocenter silently
        // flipped -- baked into the "golden" value instead of being caught.
        // Independently verified via RDKit: parse the input, zero every
        // atom's isotope, `AssignStereochemistry(cleanIt=True, force=True)`,
        // and compare `MolToInchi` -- only the value below (not the old one)
        // matches the isotope-free reference's InChI.
        let mol = parse("OC[C@H]1O[C@@H]([13CH2]O)[C@H](O)[C@@H](O)[C@@H]1O").unwrap();
        let (result, _) = isotope_parent(&mol);
        assert_eq!(
            chematic_smiles::canonical_smiles(&result),
            canon("[C@@H]1(CO)O[C@@H]([C@@H](O)[C@@H]([C@H]1O)O)CO")
        );
    }

    #[test]
    fn tp2_21_stereo_parent_alanine() {
        let mol = parse("C[C@@H](N)C(=O)O").unwrap();
        let (result, record) = stereo_parent(&mol);
        assert_eq!(
            chematic_smiles::canonical_smiles(&result),
            canon("C(C(O)=O)(C)N")
        );
        assert!(record.fragments.is_empty());
    }

    #[test]
    fn tp2_22_fragment_parent_hcl_salt() {
        let mol = parse("CN1CCC(CC1)Nc1ccccc1.Cl").unwrap();
        let (result, _) = fragment_parent(&mol);
        assert_eq!(
            chematic_smiles::canonical_smiles(&result),
            canon("N(c2ccccc2)C1CCN(CC1)C")
        );
    }

    #[test]
    fn holdout_04_parent_functions_noop_on_toluene() {
        let mol = parse("Cc1ccccc1").unwrap();
        let expected = chematic_smiles::canonical_smiles(&mol);
        assert_eq!(
            chematic_smiles::canonical_smiles(&fragment_parent(&mol).0),
            expected
        );
        assert_eq!(
            chematic_smiles::canonical_smiles(&charge_parent(&mol).0),
            expected
        );
        assert_eq!(
            chematic_smiles::canonical_smiles(&isotope_parent(&mol).0),
            expected
        );
        assert_eq!(
            chematic_smiles::canonical_smiles(&stereo_parent(&mol).0),
            expected
        );
    }

    // -----------------------------------------------------------------------
    // Issue #407: normalize_zwitterion must never invent a proton. Every
    // check here is a property (atom/element/H-count/charge conservation),
    // not a spot check against one hardcoded expected string -- per the
    // acceptance criteria, RDKit's MolStandardize output is deliberately
    // never used as a blind oracle for these.
    // -----------------------------------------------------------------------

    /// Per-element heavy-atom counts plus total H (implicit `hydrogen_count`
    /// summed across all atoms -- explicit H atoms, if any, are also heavy
    /// atoms already counted in `counts` and contribute 0 to this sum, so
    /// nothing is double-counted).
    fn atom_and_h_census(mol: &Molecule) -> (std::collections::BTreeMap<&'static str, u32>, u32) {
        let mut counts: std::collections::BTreeMap<&'static str, u32> =
            std::collections::BTreeMap::new();
        let mut total_h: u32 = 0;
        for (_, atom) in mol.atoms() {
            *counts.entry(atom.element.symbol()).or_insert(0) += 1;
            total_h += atom.hydrogen_count.unwrap_or(0) as u32;
        }
        (counts, total_h)
    }

    fn net_charge(mol: &Molecule) -> i32 {
        mol.atoms().map(|(_, a)| a.charge as i32).sum()
    }

    /// Rebuild `mol` with atom indices reversed (last atom becomes index 0,
    /// etc.), bonds remapped accordingly -- an isomorphic relabeling used to
    /// check that `normalize_zwitterion`'s output doesn't depend on which
    /// atom happens to be enumerated first.
    fn reverse_atom_order(mol: &Molecule) -> Molecule {
        let n = mol.atom_count();
        // old atom i -> new position (n - 1 - i): the last original atom
        // becomes index 0, and so on.
        let remap: HashMap<AtomIdx, AtomIdx> = (0..n)
            .map(|i| (AtomIdx(i as u32), AtomIdx((n - 1 - i) as u32)))
            .collect();
        let mut slots: Vec<Option<chematic_core::Atom>> = vec![None; n];
        for i in 0..n {
            let old_idx = AtomIdx(i as u32);
            slots[remap[&old_idx].0 as usize] = Some(mol.atom(old_idx).clone());
        }
        let mut builder = MoleculeBuilder::new();
        for slot in slots {
            builder.add_atom(slot.expect("every slot filled by construction"));
        }
        copy_bonds(mol, &mut builder, &remap);
        builder.build()
    }

    fn assert_conserves_atoms_and_protons(label: &str, smi: &str) {
        let mol = parse(smi).unwrap();
        let result = normalize_zwitterion(&mol);
        let (before_counts, before_h) = atom_and_h_census(&mol);
        let (after_counts, after_h) = atom_and_h_census(&result);
        assert_eq!(
            before_counts, after_counts,
            "{label}: per-element heavy-atom counts must be conserved"
        );
        assert_eq!(
            before_h, after_h,
            "{label}: total H count (implicit) must be conserved -- a proton must never be invented or destroyed"
        );
        assert_eq!(
            result.atom_count(),
            mol.atom_count(),
            "{label}: atom count must be conserved"
        );
    }

    fn assert_idempotent(label: &str, smi: &str) {
        let mol = parse(smi).unwrap();
        let once = normalize_zwitterion(&mol);
        let once_canon = chematic_smiles::canonical_smiles(&once);
        let reparsed = parse(&once_canon).unwrap();
        let twice = normalize_zwitterion(&reparsed);
        let twice_canon = chematic_smiles::canonical_smiles(&twice);
        assert_eq!(
            once_canon, twice_canon,
            "{label}: normalize_zwitterion must be idempotent (re-running on its own output is a no-op)"
        );
    }

    #[test]
    fn issue407_known_repro_molecules_conserve_atoms_and_are_idempotent() {
        // The exact 3 molecules from issue #407: diazo-N,N'-dioxide-like
        // fragments where neither charged nitrogen has an available proton.
        // Before the fix: the negative oxygen was unconditionally given a
        // proton it invented from nowhere, changing the formula and net
        // charge; the two nitrogens (which have no H to give) were untouched.
        for (label, smi) in [
            ("407-1", "CC12CCCCC1(Br)[N+]([O-])=[N+]2[O-]"),
            ("407-2", "CC12CCCCC1(Cl)[N+]([O-])=[N+]2[O-]"),
            ("407-3", "CC1(C)[N+]([O-])=[N+]([O-])C1(C)Br"),
        ] {
            assert_conserves_atoms_and_protons(label, smi);
            assert_idempotent(label, smi);
            let mol = parse(smi).unwrap();
            let result = normalize_zwitterion(&mol);
            assert_eq!(
                net_charge(&result),
                net_charge(&mol),
                "{label}: net charge must be unchanged when no atom has a transferable proton"
            );
            // Neither +N has an available H, so this pair must be left
            // completely untouched -- not just "charge-neutral overall".
            let canon_before = chematic_smiles::canonical_smiles(&mol);
            let canon_after = chematic_smiles::canonical_smiles(&result);
            assert_eq!(
                canon_before, canon_after,
                "{label}: with no transferable proton anywhere, normalize_zwitterion must be a full no-op"
            );
        }
    }

    #[test]
    fn issue407_genuine_amino_acid_zwitterion_still_normalizes() {
        // Positive control: a real protonation-state zwitterion (donor has
        // an available H) must still be neutralized as before the fix.
        let mol = parse("[NH3+]CC(=O)[O-]").unwrap();
        let result = normalize_zwitterion(&mol);
        assert_eq!(
            net_charge(&result),
            0,
            "genuine zwitterion must reach net charge 0"
        );
        assert_eq!(
            chematic_smiles::canonical_smiles(&result),
            chematic_smiles::canonical_smiles(&parse("NCC(=O)O").unwrap()),
            "genuine zwitterion must still normalize to the neutral amino acid form"
        );
        assert_conserves_atoms_and_protons("genuine-zwitterion", "[NH3+]CC(=O)[O-]");
        assert_idempotent("genuine-zwitterion", "[NH3+]CC(=O)[O-]");
    }

    #[test]
    fn issue407_negative_controls_no_regression() {
        // Nitro and pyridine-N-oxide: both have a +N/-O pair (has_zwitterion
        // fires) but the +N never has an available proton in either case --
        // must be a no-op, not a corruption, and must not regress into
        // "always neutralize nitro/N-oxide groups" (a different, wrong fix).
        for (label, smi) in [
            ("nitrobenzene", "c1ccccc1[N+](=O)[O-]"),
            ("pyridine-n-oxide", "[O-][n+]1ccccc1"),
        ] {
            assert_conserves_atoms_and_protons(label, smi);
            assert_idempotent(label, smi);
            let mol = parse(smi).unwrap();
            let result = normalize_zwitterion(&mol);
            assert_eq!(
                net_charge(&result),
                net_charge(&mol),
                "{label}: no transferable proton anywhere -- net charge must be unchanged"
            );
            assert_eq!(
                chematic_smiles::canonical_smiles(&mol),
                chematic_smiles::canonical_smiles(&result),
                "{label}: must be a full no-op, not a partial/silent mutation"
            );
        }
    }

    #[test]
    fn issue407_atom_permutation_does_not_change_outcome() {
        // A symmetric case with two equidistant candidate positive atoms
        // (both with an available proton) for one negative atom: whichever
        // one the BFS-tie-break picks, the OUTCOME (net charge, formula,
        // canonical structure after normalization) must not depend on which
        // atom happened to be enumerated first.
        let mol = parse("[NH3+]CC([O-])C[NH3+]").unwrap();
        let reversed = reverse_atom_order(&mol);
        let result = normalize_zwitterion(&mol);
        let result_reversed = normalize_zwitterion(&reversed);
        assert_eq!(
            net_charge(&result),
            net_charge(&result_reversed),
            "atom-order permutation must not change the net charge outcome"
        );
        assert_eq!(
            chematic_smiles::canonical_smiles(&result),
            chematic_smiles::canonical_smiles(&result_reversed),
            "atom-order permutation must not change the normalized structure \
             (whichever of the two equidistant +N atoms the BFS tie-break \
             picks, the two are chemically interchangeable here, so the \
             canonical result must converge either way)"
        );
    }
}
