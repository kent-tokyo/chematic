use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::VecDeque;

use chematic_core::{
    AtomIdx, BondIdx, BondOrder, Chirality, Molecule, MoleculeBuilder, STEREO_H_SENTINEL,
};
use chematic_perception::RingSet;
use chematic_smarts::{
    AtomPrimitive, AtomQuery, BondPrimitive, BondQuery, MatchConfig, QueryMolecule,
    find_matches_with_config, find_matches_with_rings_and_config,
};

use crate::reaction::{RxnError, parse_reaction_template};
use crate::requirements::ReactionRequirements;

mod diagnostics;
pub use diagnostics::{
    ProductAtomRejectionDiagnostic, ProductRejectionReason, ProductSanitizationDiagnostic,
    ReactionCompatibilityUnsupported, ReactionTransformDiagnostics, ReactionTransformLimits,
    ReactionVariantDiagnostics, TransformError,
};

/// Products and bounded per-match accounting from a reaction transformation.
pub struct ReactionTransformReport {
    pub products: Vec<Vec<Molecule>>,
    pub diagnostics: ReactionTransformDiagnostics,
}

/// The same bounded reaction accounting with per-product-atom provenance.
/// Product and match ordering matches [`ReactionTransformReport`].
pub struct TracedReactionTransformReport {
    pub products: Vec<Vec<TracedProduct>>,
    pub diagnostics: ReactionTransformDiagnostics,
}

/// Native-profile report with typed details for product sets rejected during
/// sanitization.
///
/// Rejected products are never returned in [`Self::report`]. The diagnostics
/// are ordered by rejected match, then by product position within that match.
/// This additive type keeps the existing aggregate report's struct-literal
/// shape stable for downstream Rust callers.
pub struct NativeDetailedReactionTransformReport {
    pub report: TracedReactionTransformReport,
    pub rejection_diagnostics: Vec<Vec<ProductSanitizationDiagnostic>>,
}

/// RDKit-profile report including product sets rejected during sanitization.
///
/// Kept separate from [`TracedReactionTransformReport`] so adding RDKit-only
/// diagnostics does not break existing Rust callers that construct the
/// general report with a struct literal.
pub struct RdkitTracedReactionTransformReport {
    pub report: TracedReactionTransformReport,
    /// RDKit profile only: the valence-rejected product sets as RDKit's
    /// `RunReactants` returns them (unsanitized), each product with whether
    /// RDKit's sanitization would accept it.
    pub rejected_products: Vec<Vec<(TracedProduct, bool)>>,
    /// Stable, bounded explanations aligned with [`Self::rejected_products`].
    /// Each outer entry is one rejected product set and each inner entry is
    /// the verdict for the product at the same position in that set.
    pub rejection_diagnostics: Vec<Vec<ProductSanitizationDiagnostic>>,
}

/// One match's product set: accepted, or rejected by product valence (with
/// each product's own verdict).
enum ProductSetOutcome {
    Valid(Vec<TracedProduct>),
    Rejected(Vec<TracedProduct>, Vec<ProductSanitizationDiagnostic>),
}

/// Apply a SMIRKS template to input reactant molecules.
///
/// Returns all combinations of product sets — one per unique match across all
/// reactant templates.  Each inner `Vec<Molecule>` contains one product per
/// product component in the SMIRKS right-hand side.
///
/// Returns `Ok(vec![])` when no match is found.
///
/// Unmapped atoms attached to a mapped core atom (substituents) are
/// automatically carried through to the matching product template.
/// Use [`run_reactants_strict`] to return only mapped atoms.
pub fn run_reactants(
    smirks: &str,
    reactants: &[&Molecule],
) -> Result<Vec<Vec<Molecule>>, TransformError> {
    run_reactants_with_limits(smirks, reactants, &ReactionTransformLimits::default())
}

/// Apply a SMIRKS template with an explicit accepted-match limit.
pub fn run_reactants_with_limits(
    smirks: &str,
    reactants: &[&Molecule],
    limits: &ReactionTransformLimits,
) -> Result<Vec<Vec<Molecule>>, TransformError> {
    run_reactants_impl(smirks, reactants, true, limits)
}

/// Apply a SMIRKS template and retain bounded accounting for filtered matches.
pub fn run_reactants_with_diagnostics(
    smirks: &str,
    reactants: &[&Molecule],
    limits: &ReactionTransformLimits,
) -> Result<ReactionTransformReport, TransformError> {
    PreparedReaction::shared(smirks)?.run_reactants_with_diagnostics(reactants, limits)
}

/// Apply a native-profile SMIRKS template and retain typed details for product
/// sets rejected during sanitization.
///
/// Accepted products retain atom provenance. Rejected products are not exposed
/// as successful reactions; only bounded, reliable rejection evidence is
/// returned alongside the aggregate accounting.
pub fn run_reactants_traced_with_rejection_diagnostics(
    smirks: &str,
    reactants: &[&Molecule],
    limits: &ReactionTransformLimits,
) -> Result<NativeDetailedReactionTransformReport, TransformError> {
    PreparedReaction::shared(smirks)?
        .run_reactants_traced_with_rejection_diagnostics(reactants, limits)
}

/// Like [`run_reactants`] but **does not carry through substituents**.
///
/// Only atoms that appear explicitly in the product template (via atom maps or
/// new template atoms) are included in each product.  Unmapped neighbors of
/// core atoms are **not** collected via BFS.
///
/// Useful when the SMIRKS describes a complete molecule transformation and
/// you do not want R-group carry-through behaviour.
pub fn run_reactants_strict(
    smirks: &str,
    reactants: &[&Molecule],
) -> Result<Vec<Vec<Molecule>>, TransformError> {
    run_reactants_strict_with_limits(smirks, reactants, &ReactionTransformLimits::default())
}

/// Apply a strict SMIRKS template with an explicit accepted-match limit.
pub fn run_reactants_strict_with_limits(
    smirks: &str,
    reactants: &[&Molecule],
    limits: &ReactionTransformLimits,
) -> Result<Vec<Vec<Molecule>>, TransformError> {
    run_reactants_impl(smirks, reactants, false, limits)
}

fn run_reactants_impl(
    smirks: &str,
    reactants: &[&Molecule],
    carry_substituents: bool,
    limits: &ReactionTransformLimits,
) -> Result<Vec<Vec<Molecule>>, TransformError> {
    let variants = crate::reaction::expand_atomic_number_primitives(smirks)?;
    let mut products = Vec::new();
    for variant in variants {
        products.extend(PreparedReaction::shared(&variant)?.run_reactants_impl(
            reactants,
            carry_substituents,
            limits,
        )?);
        if products.len() > limits.max_matches {
            return Err(TransformError::ResourceLimit {
                resource: "reaction products",
                actual: products.len(),
                limit: limits.max_matches,
            });
        }
    }
    Ok(products)
}

/// One accepted match of `smirks`'s reactant-side pattern(s) against a set of
/// input molecules — one map per reactant-template slot, keyed by the query
/// atom index (position within that reactant template) to the matched
/// [`AtomIdx`] in the corresponding input molecule. This is exactly the
/// per-combination mapping [`run_reactants_impl`] already builds internally,
/// given a name and a public seam between "enumerate matches" and
/// "apply one match" (issue #225).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReactionMatch {
    pub per_reactant: Vec<FxHashMap<usize, AtomIdx>>,
}

impl ReactionMatch {
    /// `atom_map` number → (reactant slot index, matched `AtomIdx`),
    /// resolved against `smirks`'s own atom-map annotations — so a caller
    /// can look up "where did mapped atom N end up in my input molecules"
    /// without re-deriving the reactant templates' atom maps itself.
    pub fn atom_map_positions(
        &self,
        smirks: &str,
    ) -> Result<FxHashMap<u16, (usize, AtomIdx)>, TransformError> {
        // Only the reactant side matters here; the product side may still
        // hold atomic-number primitives (`[#7:1]`) that only expansion reads.
        let Some(parts) = crate::reaction::split_reaction_parts(smirks) else {
            return Err(TransformError::SmirksParse(RxnError::MissingArrow));
        };
        let (_, queries) = parse_reactant_templates(parts[0], false)?;
        let n_templates = queries.len();
        if self.per_reactant.len() != n_templates {
            return Err(TransformError::ReactantCountMismatch {
                expected: n_templates,
                got: self.per_reactant.len(),
            });
        }
        Ok(global_map_of(
            &self.per_reactant,
            &template_atom_maps_of(&queries),
        ))
    }
}

/// Enumerate every match of `smirks`'s reactant-side pattern(s) against
/// `reactants`, without applying the transformation or building any
/// product. One entry per combination that passes the existing
/// chirality/E-Z stereo post-checks — i.e. exactly the matches
/// [`run_reactants`] would go on to build products for (issue #225).
pub fn find_reaction_matches(
    smirks: &str,
    reactants: &[&Molecule],
) -> Result<Vec<ReactionMatch>, TransformError> {
    find_reaction_matches_with_limits(smirks, reactants, &ReactionTransformLimits::default())
}

/// Enumerate reaction matches with an explicit accepted-match limit.
pub fn find_reaction_matches_with_limits(
    smirks: &str,
    reactants: &[&Molecule],
    limits: &ReactionTransformLimits,
) -> Result<Vec<ReactionMatch>, TransformError> {
    PreparedReaction::new(smirks)?.find_matches_with_limits(reactants, limits)
}

/// Apply the reaction for exactly one match (as returned by
/// [`find_reaction_matches`]), producing the product set for that match
/// alone. `Ok(None)` means this match's product set failed the existing
/// valence filter — the same case [`run_reactants`] silently drops today
/// (issue #225).
///
/// Does not re-run the chirality/E-Z stereo post-checks that
/// [`find_reaction_matches`] already applied — `m` is expected to be one
/// of the matches it returned (or otherwise already known to satisfy
/// them). Re-parses `smirks`, matching [`run_reactants`]'s own existing
/// per-call behavior. Use [`PreparedReaction::apply_match`] to reuse an
/// already compiled template.
pub fn apply_reaction_match(
    smirks: &str,
    reactants: &[&Molecule],
    m: &ReactionMatch,
    carry_substituents: bool,
) -> Result<Option<Vec<Molecule>>, TransformError> {
    PreparedReaction::new(smirks)?.apply_match(reactants, m, carry_substituents)
}

/// A reactant atom: slot `reactant` of the reactant list passed to the
/// reaction, atom `atom` of that molecule.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ReactantAtom {
    /// Index into the reactant slice the reaction was applied to.
    pub reactant: usize,
    /// Atom index within that reactant.
    pub atom: AtomIdx,
}

/// A reaction product together with the origin of each of its atoms.
///
/// `atom_sources[i]` is the reactant atom that product atom `AtomIdx(i)` was
/// copied from — an atom-mapped template atom or a carried substituent — or
/// `None` for an atom the product template creates. The molecule is exactly
/// the one the untraced API ([`apply_reaction_match`],
/// [`PreparedReaction::apply_match`]) returns.
#[derive(Clone)]
pub struct TracedProduct {
    /// The product molecule.
    pub molecule: Molecule,
    /// Per product atom, the reactant atom it came from (`None` = new atom).
    pub atom_sources: Vec<Option<ReactantAtom>>,
    /// Per product atom, the map label from the product template. Carried
    /// substituents have `None`; map labels are intentionally absent from
    /// the emitted molecule itself.
    pub template_maps: Vec<Option<u16>>,
}

/// [`apply_reaction_match`] with per-atom provenance: the same product set
/// (or `Ok(None)` for a valence-rejected match), each product paired with the
/// reactant atom every product atom was copied from (issue #650).
pub fn apply_reaction_match_traced(
    smirks: &str,
    reactants: &[&Molecule],
    m: &ReactionMatch,
    carry_substituents: bool,
) -> Result<Option<Vec<TracedProduct>>, TransformError> {
    PreparedReaction::new(smirks)?.apply_match_traced(reactants, m, carry_substituents)
}

/// Immutable target-side state for repeated prepared-reaction matching.
///
/// The context owns a snapshot of the target molecule and its ring index, so
/// callers can build it once and share it across threads and many query
/// templates without repeating ring perception.  Owning the snapshot also
/// prevents a caller from mutating or dropping the target while a batch of
/// matches is in flight.  Matching semantics and result ordering are the same
/// as [`PreparedReaction::find_matches_with_rings`].
#[derive(Clone)]
pub struct ReactionMatchContext {
    target: Molecule,
    rings: RingSet,
}

impl ReactionMatchContext {
    /// Build reusable target-side matching state.  Pass a precomputed ring
    /// set when one is already available; otherwise SSSR perception is done
    /// exactly once during construction.
    pub fn new(mol: &Molecule, rings: Option<&RingSet>) -> Self {
        Self {
            target: mol.clone(),
            rings: rings
                .cloned()
                .unwrap_or_else(|| chematic_perception::find_sssr(mol)),
        }
    }

    /// The immutable target snapshot held by this context.
    pub fn target(&self) -> &Molecule {
        &self.target
    }

    /// The precomputed ring index held by this context.
    pub fn rings(&self) -> &RingSet {
        &self.rings
    }
}

/// Parsed SMIRKS plus everything derived from it that matching needs —
/// shared by [`run_reactants_impl`], [`find_reaction_matches`], and
/// [`apply_reaction_match`] so the three can never compute it
/// inconsistently.
pub struct PreparedReaction {
    rxn: crate::reaction::Reaction,
    normalized_smirks: String,
    queries: Vec<QueryMolecule>,
    template_atom_maps: Vec<Vec<Option<u16>>>,
    has_stereo: bool,
    has_ez_stereo: bool,
    /// Compiled SMILES-compatible variants for supported atomic-number
    /// primitives. The first variant remains the compatibility view for
    /// existing match/apply helpers; application methods dispatch to all.
    variants: Option<Vec<PreparedReaction>>,
    requirements: ReactionRequirements,
    /// Per product template, per atom: what its SMILES spelling specifies.
    product_specs: Vec<Vec<ProductAtomSpec>>,
    /// Reactant-template atom facts per atom map, for `build_product`.
    reactant_map_atoms: FxHashMap<u16, ReactantMapAtom>,
    /// Atom-map pairs bonded in a reactant template (unmapped = 0).
    reactant_bond_maps: FxHashSet<(u16, u16)>,
    /// Reactant queries for the RDKit profile where they differ from
    /// `queries`: stereo templates matched by their SMARTS reading with
    /// `@`/`@@` dropped (RDKit does not match on them).
    rdkit_queries: Option<Vec<QueryMolecule>>,
    /// RDKit's reading of the template stereo, built on first use.
    rdkit_stereo: std::sync::OnceLock<Option<crate::rdkit_stereo::ReactionStereo>>,
}

/// Apply `smirks` with RDKit 2026.03.6 reaction semantics; see
/// [`PreparedReaction::run_reactants_traced_rdkit_2026_03_6`]. Unlike
/// [`PreparedReaction::new`], a SMARTS-only reactant template may carry
/// `@`/`@@` here, since RDKit reads them only for the product.
pub fn run_reactants_traced_rdkit_2026_03_6(
    smirks: &str,
    reactants: &[&Molecule],
    limits: &ReactionTransformLimits,
) -> Result<RdkitProfileOutcome, TransformError> {
    PreparedReaction::shared_with_reading(smirks, true)?
        .run_reactants_traced_rdkit_2026_03_6(reactants, limits)
}

/// Apply `smirks` with the pinned RDKit profile and retain rejected products.
///
/// This additive API carries RDKit-only sanitization diagnostics without
/// extending [`TracedReactionTransformReport`]'s stable struct-literal shape.
pub fn run_reactants_traced_rdkit_2026_03_6_detailed(
    smirks: &str,
    reactants: &[&Molecule],
    limits: &ReactionTransformLimits,
) -> Result<RdkitDetailedProfileOutcome, TransformError> {
    PreparedReaction::shared_with_reading(smirks, true)?
        .run_reactants_traced_rdkit_2026_03_6_detailed(reactants, limits)
}

/// Prepared templates by SMIRKS text, for the free functions that take a
/// SMIRKS string: a batch applies the same template to many molecules, and
/// preparing it (parsing and normalizing every component) took about a
/// third of a BioTransformer corpus run.
///
/// Bounded at [`PREPARED_CACHE_CAPACITY`] templates per reading (native and
/// RDKit) in [`PREPARED_CACHE_SHARDS`] shards (threads applying different
/// templates rarely share a lock), each evicting its least recently used
/// template when full. The cache never empties at once: it was cleared whole
/// when full, so every template was prepared again after each 2,048 new ones.
struct PreparedShard {
    entries: rustc_hash::FxHashMap<String, (std::sync::Arc<PreparedReaction>, u64)>,
    clock: u64,
}

const PREPARED_CACHE_CAPACITY: usize = 2048;
const PREPARED_CACHE_SHARDS: usize = 16;
const PREPARED_SHARD_CAPACITY: usize = PREPARED_CACHE_CAPACITY / PREPARED_CACHE_SHARDS;

/// Templates prepared because the cache did not hold them (for tests and
/// load measurements).
static PREPARED_CACHE_MISSES: std::sync::atomic::AtomicUsize =
    std::sync::atomic::AtomicUsize::new(0);

#[doc(hidden)]
pub fn prepared_cache_misses() -> usize {
    PREPARED_CACHE_MISSES.load(std::sync::atomic::Ordering::Relaxed)
}

impl PreparedShard {
    fn get(&mut self, smirks: &str) -> Option<std::sync::Arc<PreparedReaction>> {
        self.clock += 1;
        let clock = self.clock;
        self.entries.get_mut(smirks).map(|(hit, used)| {
            *used = clock;
            std::sync::Arc::clone(hit)
        })
    }

    fn insert(&mut self, smirks: String, prepared: std::sync::Arc<PreparedReaction>) {
        if self.entries.len() >= PREPARED_SHARD_CAPACITY && !self.entries.contains_key(&smirks) {
            // A miss already costs a template preparation; a scan of the
            // shard's entries is small beside it.
            if let Some(oldest) = self
                .entries
                .iter()
                .min_by_key(|(_, (_, used))| *used)
                .map(|(key, _)| key.clone())
            {
                self.entries.remove(&oldest);
            }
        }
        self.clock += 1;
        self.entries.insert(smirks, (prepared, self.clock));
    }
}

impl PreparedReaction {
    /// [`Self::new`], shared from a process-wide cache keyed by the SMIRKS.
    pub fn shared(smirks: &str) -> Result<std::sync::Arc<Self>, TransformError> {
        Self::shared_with_reading(smirks, false)
    }

    fn shared_with_reading(
        smirks: &str,
        rdkit_reading: bool,
    ) -> Result<std::sync::Arc<Self>, TransformError> {
        type Shards = [std::sync::Mutex<PreparedShard>; 2 * PREPARED_CACHE_SHARDS];
        static CACHE: std::sync::OnceLock<Shards> = std::sync::OnceLock::new();
        let shards = CACHE.get_or_init(|| {
            std::array::from_fn(|_| {
                std::sync::Mutex::new(PreparedShard {
                    entries: Default::default(),
                    clock: 0,
                })
            })
        });
        let hash = {
            use std::hash::BuildHasher;
            rustc_hash::FxBuildHasher.hash_one(smirks)
        };
        // One half of the shards per reading.
        let shard = &shards[(hash as usize % PREPARED_CACHE_SHARDS)
            + usize::from(rdkit_reading) * PREPARED_CACHE_SHARDS];
        if let Some(hit) = shard.lock().ok().and_then(|mut s| s.get(smirks)) {
            return Ok(hit);
        }
        PREPARED_CACHE_MISSES.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let prepared = std::sync::Arc::new(Self::new_with_reading(smirks, rdkit_reading)?);
        if let Ok(mut s) = shard.lock() {
            s.insert(smirks.to_string(), std::sync::Arc::clone(&prepared));
        }
        Ok(prepared)
    }
}

/// Outcome of applying a reaction under the pinned RDKit 2026.03.6 profile.
pub enum RdkitProfileOutcome {
    /// Products, accounting and provenance as RDKit would give them.
    Report(TracedReactionTransformReport),
    /// The profile cannot claim a result for these inputs.
    Unsupported(ReactionCompatibilityUnsupported),
}

/// Detailed outcome of the pinned RDKit profile, including rejected products.
pub enum RdkitDetailedProfileOutcome {
    Report(RdkitTracedReactionTransformReport),
    Unsupported(ReactionCompatibilityUnsupported),
}

/// Which reaction semantics a run follows.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Profile {
    Native,
    Rdkit,
}

/// What `build_product` needs to know about a mapped reactant-template atom:
/// its element (or wildcard) for the element-change rule and its template
/// degree for the hydrogen rule.
#[derive(Clone, Copy, Debug)]
struct ReactantMapAtom {
    element: chematic_core::Element,
    wildcard: bool,
    degree: usize,
}

mod prepared;
mod templates;
use templates::*;

mod matching;
use matching::*;

mod sanitization;
use sanitization::*;
mod stereo;
use stereo::*;
mod product;
use product::*;
/// Standard Cartesian product: given `sets[0], sets[1], …`, return all
/// ordered selections of one element from each set.
fn cartesian_product<T: Clone>(sets: &[Vec<T>]) -> Vec<Vec<T>> {
    let mut result: Vec<Vec<T>> = vec![vec![]];
    for set in sets {
        result = result
            .into_iter()
            .flat_map(|combo| {
                set.iter().map(move |item| {
                    let mut new_combo = combo.clone();
                    new_combo.push(item.clone());
                    new_combo
                })
            })
            .collect();
    }
    result
}

#[cfg(test)]
#[path = "transform/tests.rs"]
mod tests;

#[cfg(test)]
#[path = "transform/input_boundary_contract_tests.rs"]
mod input_boundary_contract_tests;
