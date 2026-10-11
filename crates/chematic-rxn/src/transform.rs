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

/// Error type for SMIRKS transformation.
#[derive(Debug)]
pub enum TransformError {
    /// A configured match-enumeration limit was exceeded.
    ResourceLimit {
        resource: &'static str,
        actual: usize,
        limit: usize,
    },
    SmirksParse(RxnError),
    ReactantCountMismatch {
        expected: usize,
        got: usize,
    },
}

impl core::fmt::Display for TransformError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::ResourceLimit {
                resource,
                actual,
                limit,
            } => write!(
                f,
                "reaction {resource} exceeds limit {limit} (got {actual})"
            ),
            Self::SmirksParse(e) => write!(f, "SMIRKS parse error: {e}"),
            Self::ReactantCountMismatch { expected, got } => {
                write!(f, "reactant count mismatch: expected {expected}, got {got}")
            }
        }
    }
}

/// Resource limits for reaction-template matching and product generation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReactionTransformLimits {
    /// Maximum number of accepted match combinations.
    pub max_matches: usize,
}

/// Bounded accounting for one reaction-template application.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReactionTransformDiagnostics {
    /// Number of matches accepted by the matcher and stereo post-checks.
    pub accepted_matches: usize,
    /// Number of accepted matches that produced a valence-valid product set.
    pub applied_products: usize,
    /// Number of accepted matches rejected by product valence validation.
    pub valence_rejected_matches: usize,
    /// Always false for the current fail-on-limit policy; reserved for a future
    /// truncating mode so consumers can distinguish partial reports.
    pub truncated_matches: bool,
}

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

/// Why RDKit-compatible product sanitization rejected a generated product.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProductRejectionReason {
    /// An aromatic atom or bond is inconsistent with ring membership.
    Aromaticity,
    /// The aromatic bond system cannot be assigned a valid Kekule form.
    Kekulization,
    /// One or more atoms exceed RDKit's allowed explicit valence.
    Valence,
    /// The product was rejected but no reliable narrower cause was found.
    Unknown,
}

impl ProductRejectionReason {
    /// Stable machine-readable code shared by Rust, Python and WASM.
    pub const fn reason_code(self) -> &'static str {
        match self {
            Self::Aromaticity => "aromaticity",
            Self::Kekulization => "kekulization",
            Self::Valence => "valence",
            Self::Unknown => "unknown",
        }
    }
}

/// Reliable atom-level context for one rejected product.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProductAtomRejectionDiagnostic {
    /// Atom index in the generated product before any serialization reorder.
    pub atom: AtomIdx,
    /// Product-template map label, when this atom came from a mapped slot.
    pub template_map: Option<u16>,
    /// Element of the generated product atom.
    pub element: chematic_core::Element,
    /// Explicit hydrogen count recorded on the generated atom.
    pub explicit_hydrogen_count: Option<u8>,
    /// Kekule explicit valence after RDKit cleanup rules, when applicable.
    pub observed_valence: Option<i16>,
    /// Largest valence RDKit permits after cleanup, when bounded.
    pub max_allowed_valence: Option<i16>,
}

/// Sanitization verdict for one generated product in a rejected product set.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProductSanitizationDiagnostic {
    /// Zero-based product position within the generated product set.
    pub product_index: usize,
    /// Whether this product itself is sanitizable. A set is rejected when at
    /// least one product is not, so accepted siblings remain visible here.
    pub accepted: bool,
    /// Typed rejection reason (`None` for an accepted sibling product).
    pub reason: Option<ProductRejectionReason>,
    /// At most 32 implicated atoms, in ascending product atom order.
    pub atoms: Vec<ProductAtomRejectionDiagnostic>,
}

/// One match's product set: accepted, or rejected by product valence (with
/// each product's own verdict).
enum ProductSetOutcome {
    Valid(Vec<TracedProduct>),
    Rejected(Vec<TracedProduct>, Vec<ProductSanitizationDiagnostic>),
}

/// Why the pinned RDKit 2026.03.6 compatibility profile cannot safely claim
/// a reaction result. Native CheMatic reaction semantics are unchanged.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactionCompatibilityUnsupported {
    /// A stereo template uses syntax the RDKit reading here does not follow
    /// (component grouping, `@?`, non-tetrahedral chirality classes).
    ChiralReactantTemplateSemantics,
    /// A product's tetrahedral stereo depends on a reactant bond order the
    /// input molecule does not record (RDKit orders an atom's ring-closure
    /// bonds by ring number).
    AmbiguousStereoBondOrder,
    /// RDKit ignores reactant-template `/` `\\` when matching; this match
    /// was rejected by the native E/Z check.
    EzReactantTemplateSemantics,
    /// A reactant template asks for a ring count or ring size (`R<n>`,
    /// `r<n>`) and RDKit's symmetrized ring set for this reactant could not
    /// be settled (the ring model is ambiguous or over its budget).
    RingModelUnresolved,
}

impl ReactionCompatibilityUnsupported {
    /// Stable machine-readable refusal code shared by bindings and audit gates.
    pub const fn reason_code(self) -> &'static str {
        match self {
            Self::ChiralReactantTemplateSemantics => "chiral_reactant_template_semantics",
            Self::AmbiguousStereoBondOrder => "ambiguous_stereo_bond_order",
            Self::EzReactantTemplateSemantics => "ez_reactant_template_semantics",
            Self::RingModelUnresolved => "ring_model_unresolved",
        }
    }
}

/// Diagnostics for one normalized atomic-number variant.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReactionVariantDiagnostics {
    /// Stable zero-based order in the deterministic expansion.
    pub variant_index: usize,
    /// SMILES-compatible reaction template used for this variant.
    pub normalized_smirks: String,
    pub diagnostics: ReactionTransformDiagnostics,
}

impl Default for ReactionTransformLimits {
    fn default() -> Self {
        Self {
            max_matches: 100_000,
        }
    }
}

impl std::error::Error for TransformError {}

impl From<RxnError> for TransformError {
    fn from(e: RxnError) -> Self {
        Self::SmirksParse(e)
    }
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

impl PreparedReaction {
    /// Return a typed reason when the pinned RDKit 2026.03.6 profile cannot
    /// read this template's stereo at all (`@?`, `@TH1`, component grouping).
    /// Ordinary `@`/`@@` templates are supported; inputs whose result the
    /// profile cannot reproduce are reported per call by
    /// [`Self::run_reactants_traced_rdkit_2026_03_6`].
    pub fn rdkit_2026_03_6_unsupported_reason(&self) -> Option<ReactionCompatibilityUnsupported> {
        let unreadable =
            |p: &PreparedReaction| p.rdkit_stereo().is_none() && p.has_any_stereo_text();
        let any = match &self.variants {
            Some(variants) => variants.iter().any(unreadable),
            None => unreadable(self),
        };
        any.then_some(ReactionCompatibilityUnsupported::ChiralReactantTemplateSemantics)
    }

    /// RDKit's reading of this template's stereo (`None`: not followed).
    fn rdkit_stereo(&self) -> Option<&crate::rdkit_stereo::ReactionStereo> {
        self.rdkit_stereo
            .get_or_init(|| {
                let parts = crate::reaction::split_reaction_parts(&self.normalized_smirks)?;
                crate::rdkit_stereo::ReactionStereo::new(parts[0], parts[2])
            })
            .as_ref()
    }

    /// Whether a reactant template spells a hydrogen atom (`[H]`, `#1`).
    fn reactant_side_spells_hydrogen(&self) -> bool {
        let reactants = self
            .normalized_smirks
            .split('>')
            .next()
            .unwrap_or_default()
            .as_bytes();
        reactants.windows(2).any(|w| w == b"[H")
            || reactants
                .windows(3)
                .any(|w| w[0] == b'#' && w[1] == b'1' && !w[2].is_ascii_digit())
            || reactants.ends_with(b"#1")
    }

    /// Whether the template text spells any tetrahedral stereo.
    fn has_any_stereo_text(&self) -> bool {
        self.normalized_smirks.contains('@')
    }

    /// Apply this template with RDKit 2026.03.6 reaction semantics: reactant
    /// `@`/`@@` do not restrict matching, and product tetrahedral tags follow
    /// RDKit's inversion flags and bond orders (see `rdkit_stereo`). Products,
    /// accounting and provenance are otherwise those of
    /// [`Self::run_reactants_traced_with_diagnostics`].
    pub fn run_reactants_traced_rdkit_2026_03_6(
        &self,
        reactants: &[&Molecule],
        limits: &ReactionTransformLimits,
    ) -> Result<RdkitProfileOutcome, TransformError> {
        Ok(
            match self.run_reactants_traced_rdkit_2026_03_6_detailed(reactants, limits)? {
                RdkitDetailedProfileOutcome::Report(report) => {
                    RdkitProfileOutcome::Report(report.report)
                }
                RdkitDetailedProfileOutcome::Unsupported(reason) => {
                    RdkitProfileOutcome::Unsupported(reason)
                }
            },
        )
    }

    /// Detailed RDKit-profile application including valence-rejected products.
    pub fn run_reactants_traced_rdkit_2026_03_6_detailed(
        &self,
        reactants: &[&Molecule],
        limits: &ReactionTransformLimits,
    ) -> Result<RdkitDetailedProfileOutcome, TransformError> {
        crate::perf_counters::record_run_reactants_call();
        if let Some(reason) = self.rdkit_2026_03_6_unsupported_reason() {
            return Ok(RdkitDetailedProfileOutcome::Unsupported(reason));
        }
        let variants: Vec<&PreparedReaction> = match &self.variants {
            Some(variants) => variants.iter().collect(),
            None => vec![self],
        };
        let mut products = Vec::new();
        let mut rejected_products = Vec::new();
        let mut rejection_diagnostics = Vec::new();
        let mut diagnostics = ReactionTransformDiagnostics {
            accepted_matches: 0,
            applied_products: 0,
            valence_rejected_matches: 0,
            truncated_matches: false,
        };
        for variant in variants {
            let report =
                match variant.run_traced_profile(reactants, true, limits, Profile::Rdkit)? {
                    Ok(report) => report,
                    Err(reason) => return Ok(RdkitDetailedProfileOutcome::Unsupported(reason)),
                };
            diagnostics.accepted_matches += report.report.diagnostics.accepted_matches;
            diagnostics.applied_products += report.report.diagnostics.applied_products;
            diagnostics.valence_rejected_matches +=
                report.report.diagnostics.valence_rejected_matches;
            diagnostics.truncated_matches |= report.report.diagnostics.truncated_matches;
            products.extend(report.report.products);
            rejected_products.extend(report.rejected_products);
            rejection_diagnostics.extend(report.rejection_diagnostics);
        }
        Ok(RdkitDetailedProfileOutcome::Report(
            RdkitTracedReactionTransformReport {
                report: TracedReactionTransformReport {
                    products,
                    diagnostics,
                },
                rejected_products,
                rejection_diagnostics,
            },
        ))
    }

    /// Whether the reactant pattern has tetrahedral `@`/`@@` constraints.
    ///
    /// This is deliberately separate from E/Z bond stereo. Callers that need
    /// RDKit-compatible reaction semantics can refuse this profile explicitly
    /// rather than returning a confident, differently interpreted product.
    pub fn has_tetrahedral_reactant_stereo(&self) -> bool {
        self.has_stereo
            || self
                .variants
                .as_ref()
                .is_some_and(|variants| variants.iter().any(|variant| variant.has_stereo))
    }

    /// Parse and compile one SMIRKS template for repeated matching and
    /// application. The returned value owns all query state, is safe to share
    /// between threads, and never reparses the template during its methods.
    pub fn new(smirks: &str) -> Result<Self, TransformError> {
        Self::new_with_reading(smirks, false)
    }

    /// `rdkit_reading`: SMARTS-only reactant templates may carry `@`/`@@`,
    /// which then do not restrict matching (RDKit's reading). Only the RDKit
    /// profile entry points build with it.
    fn new_with_reading(smirks: &str, rdkit_reading: bool) -> Result<Self, TransformError> {
        // Surrounding whitespace is not part of the template (RDKit ignores
        // it; the stereo tokenizer refused the template).
        let smirks = smirks.trim();
        let variants = if rdkit_reading {
            crate::reaction::expand_atomic_number_primitives_rdkit(smirks)?
        } else {
            crate::reaction::expand_atomic_number_primitives(smirks)?
        };
        if variants.len() > 1 || variants.first().is_none_or(|variant| variant != smirks) {
            let mut compiled = Vec::with_capacity(variants.len());
            for variant in &variants {
                compiled.push(Self::new_normalized(variant, rdkit_reading)?);
            }
            let mut primary = Self::new_normalized(&variants[0], rdkit_reading)?;
            primary.variants = Some(compiled);
            primary.requirements.element_lower_bounds.clear();
            primary.requirements.aromatic_element_lower_bounds.clear();
            primary.requirements.aliphatic_element_lower_bounds.clear();
            primary.requirements.bond_lower_bounds.clear();
            return Ok(primary);
        }
        Self::new_normalized(smirks, rdkit_reading)
    }

    fn new_normalized(smirks: &str, rdkit_reading: bool) -> Result<Self, TransformError> {
        crate::perf_counters::record_reaction_parse_call();
        let (rxn, queries, product_specs, rdkit_queries) =
            parse_smirks_templates(smirks, rdkit_reading)?;
        let mut reactant_map_atoms: FxHashMap<u16, ReactantMapAtom> = FxHashMap::default();
        for template in &rxn.reactants {
            for (idx, atom) in template.atoms() {
                if let Some(map) = atom.atom_map {
                    reactant_map_atoms.entry(map).or_insert(ReactantMapAtom {
                        element: atom.element,
                        wildcard: atom.wildcard,
                        degree: template.degree(idx),
                    });
                }
            }
        }

        let mut reactant_bond_maps: FxHashSet<(u16, u16)> = FxHashSet::default();
        for template in &rxn.reactants {
            for (_, bond) in template.bonds() {
                let a = template.atom(bond.atom1).atom_map.unwrap_or(0);
                let b = template.atom(bond.atom2).atom_map.unwrap_or(0);
                reactant_bond_maps.insert((a.min(b), a.max(b)));
            }
        }

        // Record the atom-map number for each query atom index.
        let requirements = ReactionRequirements::from_queries(&queries);
        let template_atom_maps = template_atom_maps_of(&queries);

        // Detect whether any reactant template carries @/@@ stereo, so we can apply
        // the parity-aware post-check after VF2 completes.  Chirality is NOT encoded
        // into the VF2 query because the raw flag comparison in eval_chirality is
        // SMILES-write-order-dependent; the correct check requires the full mapping
        // (see smirks_chirality_ok below).
        let has_stereo = rxn
            .reactants
            .iter()
            .any(|r| r.atoms().any(|(_, a)| a.chirality != Chirality::None));
        // Similarly, E/Z double-bond stereo (/ and \) is NOT encoded into the VF2
        // query; it is checked post-VF2 via smirks_ez_stereo_ok.
        let has_ez_stereo = rxn.reactants.iter().any(|r| {
            r.bonds()
                .any(|(_, b)| matches!(b.order, BondOrder::Up | BondOrder::Down))
        });

        Ok(Self {
            rxn,
            normalized_smirks: smirks.to_string(),
            queries,
            template_atom_maps,
            has_stereo,
            has_ez_stereo,
            variants: None,
            requirements,
            product_specs,
            reactant_map_atoms,
            reactant_bond_maps,
            rdkit_queries,
            rdkit_stereo: std::sync::OnceLock::new(),
        })
    }

    /// Conservative lower bounds derived from the compiled reactant queries.
    pub fn requirements(&self) -> &ReactionRequirements {
        &self.requirements
    }

    /// Apply the cheap requirements-only prefilter without running VF2.
    pub fn could_match(&self, reactants: &[&Molecule]) -> bool {
        self.requirements.could_match(reactants)
    }

    /// Apply this compiled template, carrying unmapped substituents through.
    pub fn run_reactants(
        &self,
        reactants: &[&Molecule],
    ) -> Result<Vec<Vec<Molecule>>, TransformError> {
        self.run_reactants_with_limits(reactants, &ReactionTransformLimits::default())
    }

    /// Apply this compiled template with an explicit accepted-match limit.
    pub fn run_reactants_with_limits(
        &self,
        reactants: &[&Molecule],
        limits: &ReactionTransformLimits,
    ) -> Result<Vec<Vec<Molecule>>, TransformError> {
        crate::perf_counters::record_run_reactants_call();
        if let Some(variants) = &self.variants {
            let mut products = Vec::new();
            for variant in variants {
                products.extend(variant.run_reactants_impl(reactants, true, limits)?);
            }
            return Ok(products);
        }
        self.run_reactants_impl(reactants, true, limits)
    }

    /// Apply this compiled template without carrying unmapped substituents.
    pub fn run_reactants_strict(
        &self,
        reactants: &[&Molecule],
    ) -> Result<Vec<Vec<Molecule>>, TransformError> {
        self.run_reactants_strict_with_limits(reactants, &ReactionTransformLimits::default())
    }

    /// Strict application with an explicit accepted-match limit.
    pub fn run_reactants_strict_with_limits(
        &self,
        reactants: &[&Molecule],
        limits: &ReactionTransformLimits,
    ) -> Result<Vec<Vec<Molecule>>, TransformError> {
        crate::perf_counters::record_run_reactants_call();
        if let Some(variants) = &self.variants {
            let mut products = Vec::new();
            for variant in variants {
                products.extend(variant.run_reactants_impl(reactants, false, limits)?);
            }
            return Ok(products);
        }
        self.run_reactants_impl(reactants, false, limits)
    }

    fn run_reactants_impl(
        &self,
        reactants: &[&Molecule],
        carry_substituents: bool,
        limits: &ReactionTransformLimits,
    ) -> Result<Vec<Vec<Molecule>>, TransformError> {
        self.run_reactants_impl_with_rings(reactants, carry_substituents, limits, None)
    }

    fn run_reactants_impl_with_rings(
        &self,
        reactants: &[&Molecule],
        carry_substituents: bool,
        limits: &ReactionTransformLimits,
        rings: Option<&[&RingSet]>,
    ) -> Result<Vec<Vec<Molecule>>, TransformError> {
        Ok(self
            .run_reactants_with_diagnostics_impl(reactants, carry_substituents, limits, rings)?
            .products)
    }

    fn run_reactants_with_diagnostics_impl(
        &self,
        reactants: &[&Molecule],
        carry_substituents: bool,
        limits: &ReactionTransformLimits,
        rings: Option<&[&RingSet]>,
    ) -> Result<ReactionTransformReport, TransformError> {
        let matches = find_matches_impl(self, reactants, limits, rings)?;
        let accepted_matches = matches.len();
        let mut products = Vec::with_capacity(accepted_matches);
        let mut valence_rejected_matches = 0;
        for m in &matches {
            match apply_match_impl(self, reactants, m, carry_substituents) {
                Some(product_set) => products.push(product_set),
                None => valence_rejected_matches += 1,
            }
        }
        Ok(ReactionTransformReport {
            diagnostics: ReactionTransformDiagnostics {
                accepted_matches,
                applied_products: products.len(),
                valence_rejected_matches,
                truncated_matches: false,
            },
            products,
        })
    }

    fn run_reactants_traced_with_diagnostics_impl(
        &self,
        reactants: &[&Molecule],
        carry_substituents: bool,
        limits: &ReactionTransformLimits,
    ) -> Result<TracedReactionTransformReport, TransformError> {
        Ok(self
            .run_traced_profile(reactants, carry_substituents, limits, Profile::Native)?
            .unwrap_or_else(|_| unreachable!("the native profile never declines"))
            .report)
    }

    fn run_traced_profile(
        &self,
        reactants: &[&Molecule],
        carry_substituents: bool,
        limits: &ReactionTransformLimits,
        profile: Profile,
    ) -> Result<
        Result<RdkitTracedReactionTransformReport, ReactionCompatibilityUnsupported>,
        TransformError,
    > {
        // RDKit's SMILES parser removes hydrogen atoms (`removeHs`): the
        // profile runs on that molecule and maps atom sources back.
        // A template that spells hydrogen atoms is for explicit-H reactants
        // (`AddHs`, as BioTransformer runs RDKit): those keep them.
        let suppressed: Vec<Option<(Molecule, Vec<usize>)>> = match profile {
            Profile::Rdkit if self.reactant_side_spells_hydrogen() => Vec::new(),
            Profile::Native => Vec::new(),
            Profile::Rdkit => reactants
                .iter()
                .map(|m| {
                    // An explicit H atom is part of the caller's molecule,
                    // as with RDKit `AddHs`; keep it. Suppression is only for
                    // bracket-H counts that RDKit's SMILES parser folds into
                    // their heavy atom.
                    if m.atoms()
                        .any(|(_, atom)| atom.element == chematic_core::Element::H)
                    {
                        None
                    } else {
                        chematic_smiles::rdkit_hydrogen_suppressed_with_map(m)
                    }
                })
                .collect(),
        };
        let suppressed_refs: Vec<&Molecule>;
        let reactants: &[&Molecule] = if suppressed.iter().any(Option::is_some) {
            suppressed_refs = reactants
                .iter()
                .zip(&suppressed)
                .map(|(&m, s)| s.as_ref().map_or(m, |(s, _)| s))
                .collect();
            &suppressed_refs
        } else {
            reactants
        };
        let restore_sources = |product: &mut TracedProduct| {
            for source in product.atom_sources.iter_mut().flatten() {
                if let Some(Some((_, kept))) = suppressed.get(source.reactant) {
                    source.atom = AtomIdx(kept[source.atom.0 as usize] as u32);
                }
            }
        };
        // RDKit's SMILES parser drops tags on atoms that cannot be centres.
        let cleaned: Vec<std::sync::Arc<Option<Molecule>>> = match profile {
            Profile::Native => Vec::new(),
            Profile::Rdkit => reactants
                .iter()
                .map(|m| crate::rdkit_stereo::rdkit_parse_cleanup(m))
                .collect(),
        };
        let cleaned_refs: Vec<&Molecule>;
        let reactants: &[&Molecule] = if cleaned.iter().any(|c| c.is_some()) {
            cleaned_refs = reactants
                .iter()
                .zip(&cleaned)
                .map(|(&m, c)| c.as_ref().as_ref().unwrap_or(m))
                .collect();
            &cleaned_refs
        } else {
            reactants
        };
        let (matches, unsupported) = find_matches_profile(self, reactants, limits, None, profile)?;
        if let Some(reason) = unsupported {
            return Ok(Err(reason));
        }
        let accepted_matches = matches.len();
        let mut products = Vec::with_capacity(accepted_matches);
        let mut valence_rejected_matches = 0;
        let mut rejected_products = Vec::new();
        let mut rejection_diagnostics = Vec::new();
        for m in &matches {
            match apply_match_profile(self, reactants, m, carry_substituents, profile) {
                Ok(ProductSetOutcome::Valid(mut product_set)) => {
                    product_set.iter_mut().for_each(restore_sources);
                    products.push(product_set);
                }
                Ok(ProductSetOutcome::Rejected(mut set, set_diagnostics)) => {
                    valence_rejected_matches += 1;
                    if profile == Profile::Rdkit {
                        set.iter_mut().for_each(restore_sources);
                        rejected_products.push(
                            set.into_iter()
                                .zip(set_diagnostics.iter().map(|item| item.accepted))
                                .collect(),
                        );
                        rejection_diagnostics.push(set_diagnostics);
                    }
                }
                Err(reason) => return Ok(Err(reason)),
            }
        }
        Ok(Ok(RdkitTracedReactionTransformReport {
            report: TracedReactionTransformReport {
                diagnostics: ReactionTransformDiagnostics {
                    accepted_matches,
                    applied_products: products.len(),
                    valence_rejected_matches,
                    truncated_matches: false,
                },
                products,
            },
            rejected_products,
            rejection_diagnostics,
        }))
    }

    /// Apply this compiled template and retain bounded accounting for filtered matches.
    pub fn run_reactants_with_diagnostics(
        &self,
        reactants: &[&Molecule],
        limits: &ReactionTransformLimits,
    ) -> Result<ReactionTransformReport, TransformError> {
        crate::perf_counters::record_run_reactants_call();
        if let Some(variants) = &self.variants {
            return aggregate_variant_reports(variants, reactants, limits, true, None);
        }
        self.run_reactants_with_diagnostics_impl(reactants, true, limits, None)
    }

    /// Apply this template with the same bounded accounting and product order
    /// as [`Self::run_reactants_with_diagnostics`], retaining atom origins.
    pub fn run_reactants_traced_with_diagnostics(
        &self,
        reactants: &[&Molecule],
        limits: &ReactionTransformLimits,
    ) -> Result<TracedReactionTransformReport, TransformError> {
        crate::perf_counters::record_run_reactants_call();
        if let Some(variants) = &self.variants {
            let mut products = Vec::new();
            let mut diagnostics = ReactionTransformDiagnostics {
                accepted_matches: 0,
                applied_products: 0,
                valence_rejected_matches: 0,
                truncated_matches: false,
            };
            for variant in variants {
                let report =
                    variant.run_reactants_traced_with_diagnostics_impl(reactants, true, limits)?;
                diagnostics.accepted_matches += report.diagnostics.accepted_matches;
                diagnostics.applied_products += report.diagnostics.applied_products;
                diagnostics.valence_rejected_matches += report.diagnostics.valence_rejected_matches;
                diagnostics.truncated_matches |= report.diagnostics.truncated_matches;
                products.extend(report.products);
            }
            return Ok(TracedReactionTransformReport {
                products,
                diagnostics,
            });
        }
        self.run_reactants_traced_with_diagnostics_impl(reactants, true, limits)
    }

    /// Apply this compiled template and return diagnostics separately for
    /// every normalized atomic-number variant. This is additive to the
    /// aggregate report and keeps variant identity explicit for callers that
    /// need to explain aromatic/aliphatic alternatives.
    pub fn run_reactants_with_variant_diagnostics(
        &self,
        reactants: &[&Molecule],
        limits: &ReactionTransformLimits,
    ) -> Result<Vec<ReactionVariantDiagnostics>, TransformError> {
        crate::perf_counters::record_run_reactants_call();
        let mut reports = Vec::new();
        if let Some(variants) = &self.variants {
            for (variant_index, variant) in variants.iter().enumerate() {
                reports.push(variant_diagnostics(
                    variant,
                    variant_index,
                    reactants,
                    limits,
                    None,
                )?);
            }
        } else {
            reports.push(variant_diagnostics(self, 0, reactants, limits, None)?);
        }
        Ok(reports)
    }

    /// Enumerate accepted matches without applying the transformation.
    pub fn find_matches(
        &self,
        reactants: &[&Molecule],
    ) -> Result<Vec<ReactionMatch>, TransformError> {
        self.find_matches_with_limits(reactants, &ReactionTransformLimits::default())
    }

    /// Enumerate accepted matches with an explicit accepted-match limit.
    pub fn find_matches_with_limits(
        &self,
        reactants: &[&Molecule],
        limits: &ReactionTransformLimits,
    ) -> Result<Vec<ReactionMatch>, TransformError> {
        find_matches_impl(self, reactants, limits, None)
    }

    /// Enumerate matches against one reusable target-side context.
    ///
    /// This is the convenient repeated-query path for single-reactant
    /// templates. Multi-reactant templates should continue to use
    /// [`Self::find_matches_with_rings_and_limits`] with one context per
    /// reactant slot.
    pub fn find_matches_with_context(
        &self,
        context: &ReactionMatchContext,
    ) -> Result<Vec<ReactionMatch>, TransformError> {
        self.find_matches_with_context_and_limits(context, &ReactionTransformLimits::default())
    }

    /// Context-backed matching with an explicit accepted-match limit.
    pub fn find_matches_with_context_and_limits(
        &self,
        context: &ReactionMatchContext,
        limits: &ReactionTransformLimits,
    ) -> Result<Vec<ReactionMatch>, TransformError> {
        let rings = [&context.rings];
        find_matches_impl(self, &[&context.target], limits, Some(&rings))
    }

    /// Apply this compiled template while reusing ring perception already
    /// computed for each reactant molecule.
    pub fn run_reactants_with_rings(
        &self,
        reactants: &[&Molecule],
        rings: &[&RingSet],
    ) -> Result<Vec<Vec<Molecule>>, TransformError> {
        self.run_reactants_with_rings_and_limits(
            reactants,
            rings,
            &ReactionTransformLimits::default(),
        )
    }

    /// Apply this compiled template with caller-provided ring perception and
    /// an explicit accepted-match limit.
    pub fn run_reactants_with_rings_and_limits(
        &self,
        reactants: &[&Molecule],
        rings: &[&RingSet],
        limits: &ReactionTransformLimits,
    ) -> Result<Vec<Vec<Molecule>>, TransformError> {
        crate::perf_counters::record_run_reactants_call();
        if let Some(variants) = &self.variants {
            let mut products = Vec::new();
            for variant in variants {
                products.extend(variant.run_reactants_impl_with_rings(
                    reactants,
                    true,
                    limits,
                    Some(rings),
                )?);
            }
            return Ok(products);
        }
        self.run_reactants_impl_with_rings(reactants, true, limits, Some(rings))
    }

    /// Apply this compiled template with caller-provided rings and diagnostics.
    pub fn run_reactants_with_rings_and_limits_with_diagnostics(
        &self,
        reactants: &[&Molecule],
        rings: &[&RingSet],
        limits: &ReactionTransformLimits,
    ) -> Result<ReactionTransformReport, TransformError> {
        crate::perf_counters::record_run_reactants_call();
        if let Some(variants) = &self.variants {
            return aggregate_variant_reports(variants, reactants, limits, true, Some(rings));
        }
        self.run_reactants_with_diagnostics_impl(reactants, true, limits, Some(rings))
    }

    /// Apply this compiled template with caller-provided rings and return
    /// diagnostics separately for every normalized atomic-number variant.
    /// This preserves the same variant identity contract as
    /// [`Self::run_reactants_with_variant_diagnostics`] while reusing the
    /// caller's ring perception.
    pub fn run_reactants_with_rings_and_limits_with_variant_diagnostics(
        &self,
        reactants: &[&Molecule],
        rings: &[&RingSet],
        limits: &ReactionTransformLimits,
    ) -> Result<Vec<ReactionVariantDiagnostics>, TransformError> {
        crate::perf_counters::record_run_reactants_call();
        let mut reports = Vec::new();
        if let Some(variants) = &self.variants {
            for (variant_index, variant) in variants.iter().enumerate() {
                reports.push(variant_diagnostics(
                    variant,
                    variant_index,
                    reactants,
                    limits,
                    Some(rings),
                )?);
            }
        } else {
            reports.push(variant_diagnostics(
                self,
                0,
                reactants,
                limits,
                Some(rings),
            )?);
        }
        Ok(reports)
    }

    /// Apply this compiled template in strict mode with caller-provided ring
    /// perception.
    pub fn run_reactants_strict_with_rings(
        &self,
        reactants: &[&Molecule],
        rings: &[&RingSet],
    ) -> Result<Vec<Vec<Molecule>>, TransformError> {
        self.run_reactants_strict_with_rings_and_limits(
            reactants,
            rings,
            &ReactionTransformLimits::default(),
        )
    }

    /// Strict application with caller-provided ring perception and an
    /// explicit accepted-match limit.
    pub fn run_reactants_strict_with_rings_and_limits(
        &self,
        reactants: &[&Molecule],
        rings: &[&RingSet],
        limits: &ReactionTransformLimits,
    ) -> Result<Vec<Vec<Molecule>>, TransformError> {
        crate::perf_counters::record_run_reactants_call();
        if let Some(variants) = &self.variants {
            let mut products = Vec::new();
            for variant in variants {
                products.extend(variant.run_reactants_impl_with_rings(
                    reactants,
                    false,
                    limits,
                    Some(rings),
                )?);
            }
            return Ok(products);
        }
        self.run_reactants_impl_with_rings(reactants, false, limits, Some(rings))
    }

    /// Enumerate accepted matches while reusing precomputed ring perception.
    pub fn find_matches_with_rings(
        &self,
        reactants: &[&Molecule],
        rings: &[&RingSet],
    ) -> Result<Vec<ReactionMatch>, TransformError> {
        self.find_matches_with_rings_and_limits(
            reactants,
            rings,
            &ReactionTransformLimits::default(),
        )
    }

    /// Enumerate accepted matches using caller-provided ring perception and an
    /// explicit accepted-match limit.
    pub fn find_matches_with_rings_and_limits(
        &self,
        reactants: &[&Molecule],
        rings: &[&RingSet],
        limits: &ReactionTransformLimits,
    ) -> Result<Vec<ReactionMatch>, TransformError> {
        find_matches_impl(self, reactants, limits, Some(rings))
    }

    /// Apply this template for exactly one previously accepted match.
    pub fn apply_match(
        &self,
        reactants: &[&Molecule],
        m: &ReactionMatch,
        carry_substituents: bool,
    ) -> Result<Option<Vec<Molecule>>, TransformError> {
        let n_templates = self.rxn.reactants.len();
        if reactants.len() != n_templates {
            return Err(TransformError::ReactantCountMismatch {
                expected: n_templates,
                got: reactants.len(),
            });
        }
        if m.per_reactant.len() != n_templates {
            return Err(TransformError::ReactantCountMismatch {
                expected: n_templates,
                got: m.per_reactant.len(),
            });
        }
        Ok(apply_match_impl(self, reactants, m, carry_substituents))
    }

    /// [`Self::apply_match`] with per-atom provenance (see [`TracedProduct`]).
    /// The molecules are identical to [`Self::apply_match`]'s.
    pub fn apply_match_traced(
        &self,
        reactants: &[&Molecule],
        m: &ReactionMatch,
        carry_substituents: bool,
    ) -> Result<Option<Vec<TracedProduct>>, TransformError> {
        let n_templates = self.rxn.reactants.len();
        if reactants.len() != n_templates || m.per_reactant.len() != n_templates {
            return Err(TransformError::ReactantCountMismatch {
                expected: n_templates,
                got: if reactants.len() != n_templates {
                    reactants.len()
                } else {
                    m.per_reactant.len()
                },
            });
        }
        Ok(apply_match_traced_impl(
            self,
            reactants,
            m,
            carry_substituents,
        ))
    }
}

fn aggregate_variant_reports(
    variants: &[PreparedReaction],
    reactants: &[&Molecule],
    limits: &ReactionTransformLimits,
    carry_substituents: bool,
    rings: Option<&[&RingSet]>,
) -> Result<ReactionTransformReport, TransformError> {
    let mut products = Vec::new();
    let mut diagnostics = ReactionTransformDiagnostics {
        accepted_matches: 0,
        applied_products: 0,
        valence_rejected_matches: 0,
        truncated_matches: false,
    };
    for variant in variants {
        let report = variant.run_reactants_with_diagnostics_impl(
            reactants,
            carry_substituents,
            limits,
            rings,
        )?;
        diagnostics.accepted_matches += report.diagnostics.accepted_matches;
        diagnostics.applied_products += report.diagnostics.applied_products;
        diagnostics.valence_rejected_matches += report.diagnostics.valence_rejected_matches;
        diagnostics.truncated_matches |= report.diagnostics.truncated_matches;
        products.extend(report.products);
    }
    Ok(ReactionTransformReport {
        products,
        diagnostics,
    })
}

fn variant_diagnostics(
    variant: &PreparedReaction,
    variant_index: usize,
    reactants: &[&Molecule],
    limits: &ReactionTransformLimits,
    rings: Option<&[&RingSet]>,
) -> Result<ReactionVariantDiagnostics, TransformError> {
    let report = variant.run_reactants_with_diagnostics_impl(reactants, true, limits, rings)?;
    Ok(ReactionVariantDiagnostics {
        variant_index,
        normalized_smirks: variant.normalized_smirks.clone(),
        diagnostics: report.diagnostics,
    })
}

fn template_atom_maps_of(queries: &[QueryMolecule]) -> Vec<Vec<Option<u16>>> {
    queries
        .iter()
        .map(|q| q.atoms.iter().map(|a| a.atom_map).collect())
        .collect()
}

/// Parse a SMIRKS into its reaction templates and the compiled reactant
/// queries (issue #734).
///
/// Reactant templates are SMARTS. A component that is also valid SMILES
/// (`[C:1]`, `[nH:2]`, `[O-:3]`, `[C@H:4]`…) goes through the SMILES parser
/// and [`mol_to_query`] exactly as before, which keeps the tetrahedral and
/// E/Z post-checks that read the template molecule. Any other component
/// (`[CX4:1]`, `[C;H3:1]`, `[C,N:1]`, `[!O:1]`, `[$([OH]):1]`…) is parsed
/// with the full SMARTS parser; its query is matched directly and a shadow
/// molecule (same atoms and atom maps, implied elements where the query pins
/// one) stands in for the template wherever the engine only needs atom
/// counts and maps. Stereo primitives in such a component are refused with a
/// typed error rather than silently ignored: the engine checks `@`/`@@` and
/// `/`/`\` against a template *molecule*, which a general query does not
/// provide.
///
/// Product templates stay SMILES (they are specifications, not queries); the
/// SMARTS spelling `[C;H3:1]` of an H count is accepted and rewritten to
/// `[CH3:1]`, any other `;`/`,`/`!`/`$()` product primitive is a typed error.
type SmirksTemplates = (
    crate::reaction::Reaction,
    Vec<QueryMolecule>,
    Vec<Vec<ProductAtomSpec>>,
    Option<Vec<QueryMolecule>>,
);

fn parse_smirks_templates(
    smirks: &str,
    rdkit_reading: bool,
) -> Result<SmirksTemplates, TransformError> {
    let limits = crate::reaction::ReactionParseLimits::default();
    if smirks.len() > limits.max_input_bytes {
        return Err(TransformError::SmirksParse(RxnError::ResourceLimit {
            resource: "input bytes",
            actual: smirks.len(),
            limit: limits.max_input_bytes,
        }));
    }
    // Reaction SMILES/SMIRKS: reactants>agents>products; a dative `->` is
    // a bond, not an arrow.
    let Some(parts) = crate::reaction::split_reaction_parts(smirks) else {
        return Err(TransformError::SmirksParse(RxnError::MissingArrow));
    };
    // Product query atoms (`[OX2H1:2]`) are reduced to what a product can
    // apply; callers reaching here without the atomic-number expansion
    // (e.g. `ReactionMatch::atom_map_positions`) get the same reading.
    let products = crate::reaction::normalize_product_query_atoms_reading(
        &format!(">>{}", parts[2]),
        rdkit_reading,
    )
    .map_err(TransformError::SmirksParse)?;
    let written = normalize_product_templates(&products[2..])?;
    let products = rdkit_product_bonds(&written);
    // Product component grouping, `([C:1].[O:2])`: as in RDKit, a group is
    // one product object (one molecule with disconnected components).
    let grouped = grouped_product_components(&products)?;
    let written_components: Vec<String> = match grouped_product_components(&written)? {
        Some(components) => components,
        None => written
            .split('.')
            .filter(|p| !p.is_empty())
            .map(str::to_string)
            .collect(),
    };
    let ungrouped = grouped.as_ref().map_or(products.clone(), |c| c.join("."));
    // Agents and products through the ordinary reaction parser (with its
    // limits); the reactant slot is filled below.
    // Agents take no part in template application (RDKit ignores them too)
    // and may be SMARTS (`>[O;X2]>`); when they are not SMILES they are
    // dropped rather than failing the whole template.
    let mut rxn = parse_reaction_template(&format!(">{}>{}", parts[1], ungrouped))
        .or_else(|_| parse_reaction_template(&format!(">>{ungrouped}")))?;
    let components: Vec<String> = match grouped {
        Some(components) => {
            rxn.products = components
                .iter()
                .map(|c| {
                    chematic_smiles::parse_template(c).map_err(|e| {
                        TransformError::SmirksParse(RxnError::SmilesParse {
                            part: c.clone(),
                            source: e.to_string(),
                        })
                    })
                })
                .collect::<Result<_, _>>()?;
            components
        }
        None => products
            .split('.')
            .filter(|p| !p.is_empty())
            .map(str::to_string)
            .collect(),
    };
    if written_components.len() == rxn.products.len() {
        for (component, original) in rxn.products.iter_mut().zip(&written_components) {
            if original.contains(['@', ',', ';', '&', '!', '~']) {
                zero_order_product_bonds(component, original);
            }
        }
    }
    let product_specs = components.iter().map(|c| product_atom_specs(c)).collect();
    let (reactants, queries) = parse_reactant_templates(parts[0], rdkit_reading)?;
    let rdkit_queries = rdkit_reactant_queries(parts[0], &reactants, &queries);
    rxn.reactants = reactants;
    Ok((rxn, queries, product_specs, rdkit_queries))
}

/// What a product template atom's SMILES spelling specifies. A bare bracket
/// atom (`[O:1]`) and an explicit `[OH0:1]` / `[O+0:1]` parse to the same
/// atom, but RDKit keeps the matched reactant atom's charge and hydrogens for
/// the former and applies the template's for the latter.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct ProductAtomSpec {
    charge: bool,
    hcount: bool,
}

/// [`ProductAtomSpec`] for each atom of one SMILES product component, in
/// parse order.
fn product_atom_specs(component: &str) -> Vec<ProductAtomSpec> {
    let b = component.as_bytes();
    let mut specs = Vec::new();
    let mut i = 0;
    while i < b.len() {
        match b[i] {
            b'[' => {
                let close = component[i..].find(']').map_or(b.len(), |c| i + c);
                specs.push(bracket_atom_spec(&component[i + 1..close.min(b.len())]));
                i = close + 1;
                continue;
            }
            b'%' => i += 2,
            b'B' | b'C' | b'N' | b'O' | b'P' | b'S' | b'F' | b'I' | b'b' | b'c' | b'n' | b'o'
            | b'p' | b's' | b'*' => specs.push(ProductAtomSpec::default()),
            _ => {}
        }
        i += 1;
    }
    specs
}

/// [`ProductAtomSpec`] of a SMILES bracket atom's text (between `[` and `]`).
fn bracket_atom_spec(inner: &str) -> ProductAtomSpec {
    let b = inner.as_bytes();
    let mut i = 0;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1; // isotope
    }
    // Element symbol (one or two letters) or `*`.
    if i < b.len() {
        let two = inner.get(i..i + 2).filter(|t| {
            t.as_bytes()[1].is_ascii_lowercase()
                && (chematic_core::Element::from_symbol(t).is_some()
                    || matches!(*t, "se" | "as" | "te"))
        });
        i += two.map_or(1, |_| 2);
    }
    // Chirality: `@`, `@@`, `@TH1`, `@SP2`…
    if b.get(i) == Some(&b'@') {
        i += 1;
        if b.get(i) == Some(&b'@') {
            i += 1;
        } else if inner
            .get(i..i + 2)
            .is_some_and(|t| matches!(t, "TH" | "AL" | "SP" | "TB" | "OH"))
        {
            i += 2;
            while b.get(i).is_some_and(u8::is_ascii_digit) {
                i += 1;
            }
        }
    }
    let hcount = b.get(i) == Some(&b'H');
    let rest = &inner[i.min(inner.len())..];
    let before_map = rest.split(':').next().unwrap_or("");
    ProductAtomSpec {
        charge: before_map.contains(['+', '-']),
        hcount,
    }
}

/// Reactant templates of one SMIRKS reactant side and their compiled
/// queries (see [`parse_smirks_templates`]).
fn parse_reactant_templates(
    side: &str,
    rdkit_reading: bool,
) -> Result<(Vec<Molecule>, Vec<QueryMolecule>), TransformError> {
    parse_reactant_templates_with_limits(
        side,
        rdkit_reading,
        &crate::reaction::ReactionParseLimits::default(),
    )
}

fn parse_reactant_templates_with_limits(
    side: &str,
    rdkit_reading: bool,
    limits: &crate::reaction::ReactionParseLimits,
) -> Result<(Vec<Molecule>, Vec<QueryMolecule>), TransformError> {
    let components: Vec<&str> = split_components(side);
    if components.len() > limits.max_components_per_side {
        return Err(TransformError::SmirksParse(RxnError::ResourceLimit {
            resource: "reactants",
            actual: components.len(),
            limit: limits.max_components_per_side,
        }));
    }
    let mut reactants = Vec::with_capacity(components.len());
    let mut queries = Vec::with_capacity(components.len());
    for part in components {
        let (mol, query) = match chematic_smiles::parse_template(part) {
            // Stereo templates keep the SMILES reading: the `@`/`@@` and
            // `/`/`\` post-checks read the template molecule.
            Ok(mol) if has_stereo(&mol) => {
                let query = mol_to_query(&mol);
                (mol, query)
            }
            // Otherwise the SMARTS query is what matches, as in RDKit: SMILES
            // and SMARTS readings differ for `[*:1]` (any atom, not carbon),
            // `[N+0:1]`/`[CH0:1]` (explicit zero charge / H count), `->`, and
            // the implicit bond between aromatic atoms (single or aromatic).
            // The SMILES molecule (same atoms in the same order) stays the
            // template molecule.
            Ok(mol) => match chematic_smarts::parse_smarts(part) {
                Ok(mut query)
                    if query.atoms.len() == mol.atom_count()
                        && query.bonds.len() == mol.bond_count() =>
                {
                    simplify_reactant_query(&mut query);
                    (mol, query)
                }
                _ => {
                    let query = mol_to_query(&mol);
                    (mol, query)
                }
            },
            Err(smiles_err) => {
                let mut query = chematic_smarts::parse_smarts(part).map_err(|smarts_err| {
                    TransformError::SmirksParse(RxnError::SmartsParse {
                        part: part.to_string(),
                        source: format!("{smarts_err} (as SMILES: {smiles_err})"),
                    })
                })?;
                if rdkit_reading {
                    // RDKit does not match on `@`/`@@`; their effect on the
                    // product is read from the template text.
                    for atom in &mut query.atoms {
                        if let Some(stripped) = strip_chirality(&atom.query) {
                            atom.query = stripped;
                        }
                    }
                }
                // RDKit does not match on `/`/`\` either (they match a
                // single or aromatic bond, as in plain SMARTS): the RDKit
                // profile accepts them; the native profile, which enforces
                // E/Z on SMILES templates, refuses what it cannot check.
                let shadow = shadow_molecule_of(&query, rdkit_reading).map_err(|reason| {
                    TransformError::SmirksParse(RxnError::UnsupportedReactantTemplate {
                        part: part.to_string(),
                        reason,
                    })
                })?;
                (shadow, query)
            }
        };
        if mol.atom_count() > limits.max_atoms_per_molecule {
            return Err(TransformError::SmirksParse(RxnError::ResourceLimit {
                resource: "atoms per molecule",
                actual: mol.atom_count(),
                limit: limits.max_atoms_per_molecule,
            }));
        }
        if mol.bond_count() > limits.max_bonds_per_molecule {
            return Err(TransformError::SmirksParse(RxnError::ResourceLimit {
                resource: "bonds per molecule",
                actual: mol.bond_count(),
                limit: limits.max_bonds_per_molecule,
            }));
        }
        reactants.push(mol);
        queries.push(query);
    }
    Ok((reactants, queries))
}

/// Reactant queries for the RDKit profile, where a stereo template's
/// SMARTS reading (with `@`/`@@` dropped: RDKit does not match on them)
/// differs from the SMILES reading the native profile matches. `None` when
/// no template has tetrahedral stereo.
fn rdkit_reactant_queries(
    side: &str,
    reactants: &[Molecule],
    queries: &[QueryMolecule],
) -> Option<Vec<QueryMolecule>> {
    if !reactants
        .iter()
        .any(|t| t.atoms().any(|(_, a)| a.chirality != Chirality::None))
    {
        return None;
    }
    let components = split_components(side);
    let mut any = false;
    let mut out = Vec::with_capacity(queries.len());
    for (k, query) in queries.iter().enumerate() {
        let template = reactants.get(k)?;
        let tetrahedral = template
            .atoms()
            .any(|(_, a)| a.chirality != Chirality::None);
        let replacement = tetrahedral
            .then(|| components.get(k))
            .flatten()
            .and_then(|part| chematic_smarts::parse_smarts(part).ok())
            .and_then(|mut q| {
                for atom in &mut q.atoms {
                    atom.query = strip_chirality(&atom.query)?;
                }
                for bond in &mut q.bonds {
                    if matches!(
                        bond.query,
                        BondQuery::Primitive(BondPrimitive::Up | BondPrimitive::Down)
                    ) {
                        bond.query = BondQuery::Primitive(BondPrimitive::Single);
                    }
                }
                (q.atoms.len() == template.atom_count() && q.bonds.len() == template.bond_count())
                    .then(|| {
                        simplify_reactant_query(&mut q);
                        q
                    })
            });
        any |= replacement.is_some();
        out.push(replacement.unwrap_or_else(|| query.clone()));
    }
    any.then_some(out)
}

/// `q` without its chirality primitives; `None` when one sits under `!` or
/// `,` (not a plain conjunction) or nothing else is left.
fn strip_chirality(q: &AtomQuery) -> Option<AtomQuery> {
    fn has(q: &AtomQuery) -> bool {
        match q {
            AtomQuery::Primitive(AtomPrimitive::Chirality(_)) => true,
            AtomQuery::Primitive(_) => false,
            AtomQuery::And(a, b) | AtomQuery::Or(a, b) => has(a) || has(b),
            AtomQuery::Not(a) => has(a),
        }
    }
    fn strip(q: &AtomQuery) -> Result<Option<AtomQuery>, ()> {
        match q {
            AtomQuery::Primitive(AtomPrimitive::Chirality(_)) => Ok(None),
            AtomQuery::And(a, b) => Ok(match (strip(a)?, strip(b)?) {
                (Some(a), Some(b)) => Some(AtomQuery::And(Box::new(a), Box::new(b))),
                (Some(x), None) | (None, Some(x)) => Some(x),
                (None, None) => None,
            }),
            other if has(other) => Err(()),
            other => Ok(Some(other.clone())),
        }
    }
    strip(q).ok().flatten()
}

/// Rewrite a SMARTS reactant query into the equivalent form the engine's
/// prefilter and matcher handle fastest (the form [`mol_to_query`] builds):
/// an element symbol becomes its atomic number, and an implicit bond with an
/// aliphatic endpoint — which can only match a single bond — becomes `-`.
fn simplify_reactant_query(query: &mut QueryMolecule) {
    fn symbols_to_numbers(q: &mut AtomQuery) {
        match q {
            AtomQuery::Primitive(AtomPrimitive::Symbol(sym)) => {
                if let Some(e) =
                    chematic_core::Element::from_symbol(sym).filter(|e| e.symbol() == sym.as_str())
                {
                    *q = AtomQuery::Primitive(AtomPrimitive::AtomicNum(e.atomic_number()));
                }
            }
            AtomQuery::Primitive(_) => {}
            AtomQuery::And(a, b) | AtomQuery::Or(a, b) => {
                symbols_to_numbers(a);
                symbols_to_numbers(b);
            }
            AtomQuery::Not(a) => symbols_to_numbers(a),
        }
    }
    fn aliphatic(q: &AtomQuery) -> bool {
        match q {
            AtomQuery::Primitive(AtomPrimitive::Aromatic(false)) => true,
            AtomQuery::And(a, b) => aliphatic(a) || aliphatic(b),
            _ => false,
        }
    }
    for atom in &mut query.atoms {
        symbols_to_numbers(&mut atom.query);
    }
    let aliphatic_atoms: Vec<bool> = query.atoms.iter().map(|a| aliphatic(&a.query)).collect();
    for bond in &mut query.bonds {
        if bond.query == BondQuery::Any
            && (aliphatic_atoms[bond.atom1] || aliphatic_atoms[bond.atom2])
        {
            bond.query = BondQuery::Primitive(BondPrimitive::Single);
        }
    }
}

/// Whether a SMILES-parsed reactant template carries tetrahedral or
/// double-bond stereo.
fn has_stereo(mol: &Molecule) -> bool {
    mol.atoms().any(|(_, a)| a.chirality != Chirality::None)
        || mol
            .bonds()
            .any(|(_, b)| matches!(b.order, BondOrder::Up | BondOrder::Down))
}

/// Dot-separated components of one reaction side, ignoring dots inside
/// bracket atoms and recursive SMARTS parentheses.
fn split_components(side: &str) -> Vec<&str> {
    let mut out = Vec::new();
    let (mut bracket, mut paren) = (0usize, 0usize);
    let mut start = 0usize;
    // Component grouping, `(A.B)`: as in RDKit's reaction SMARTS, the
    // group is one template whose disconnected parts match one molecule.
    let mut group: Option<usize> = None;
    let mut group_depth = 0usize;
    for (i, b) in side.bytes().enumerate() {
        if group.is_some() {
            match b {
                b'[' => bracket += 1,
                b']' => bracket = bracket.saturating_sub(1),
                b'(' if bracket == 0 => group_depth += 1,
                b')' if bracket == 0 => {
                    group_depth -= 1;
                    if group_depth == 0 {
                        let g = group.take().expect("group");
                        if i > g + 1 {
                            out.push(&side[g + 1..i]);
                        }
                        start = i + 1;
                    }
                }
                _ => {}
            }
            continue;
        }
        if b == b'(' && bracket == 0 && i == start {
            group = Some(i);
            group_depth = 1;
            continue;
        }
        match b {
            b'[' => bracket += 1,
            b']' => bracket = bracket.saturating_sub(1),
            b'(' if bracket > 0 => paren += 1,
            b')' if bracket > 0 => paren = paren.saturating_sub(1),
            b'.' if bracket == 0 && paren == 0 => {
                if i > start {
                    out.push(&side[start..i]);
                }
                start = i + 1;
            }
            _ => {}
        }
    }
    if let Some(g) = group {
        // Unclosed group: left to the template parser to refuse.
        out.push(&side[g..]);
    } else if side.len() > start {
        out.push(&side[start..]);
    }
    out
}

/// The template molecule standing in for a SMARTS-only reactant query: one
/// atom per query atom (the implied element when the query pins one, else a
/// wildcard) carrying the atom map, and one single bond per query bond.
/// `Err` names the stereo feature that makes the query unsupported.
fn shadow_molecule_of(
    query: &QueryMolecule,
    allow_direction: bool,
) -> Result<Molecule, &'static str> {
    fn has_chirality(q: &AtomQuery) -> bool {
        match q {
            AtomQuery::Primitive(AtomPrimitive::Chirality(_)) => true,
            AtomQuery::Primitive(AtomPrimitive::Recursive(sub)) => {
                sub.atoms.iter().any(|a| has_chirality(&a.query))
            }
            AtomQuery::Primitive(_) => false,
            AtomQuery::And(a, b) | AtomQuery::Or(a, b) => has_chirality(a) || has_chirality(b),
            AtomQuery::Not(a) => has_chirality(a),
        }
    }
    fn has_direction(q: &BondQuery) -> bool {
        match q {
            BondQuery::Primitive(BondPrimitive::Up | BondPrimitive::Down) => true,
            BondQuery::Primitive(_) | BondQuery::Any => false,
            BondQuery::And(a, b) | BondQuery::Or(a, b) => has_direction(a) || has_direction(b),
            BondQuery::Not(a) => has_direction(a),
        }
    }
    fn implied_element(q: &AtomQuery) -> Option<chematic_core::Element> {
        match q {
            AtomQuery::Primitive(AtomPrimitive::AtomicNum(n)) => {
                chematic_core::Element::from_atomic_number(*n)
            }
            AtomQuery::Primitive(AtomPrimitive::Symbol(sym)) => {
                chematic_core::Element::from_symbol(sym).filter(|e| e.symbol() == sym.as_str())
            }
            // The reaction engine compares this element with the product
            // template's to decide whether a mapped atom changes element
            // (issue #734). As in RDKit, an alternative (`[C,N]`), a negation
            // (`[!O]`) or a recursive query pins no element, so any concrete
            // product element is applied.
            AtomQuery::And(a, b) => implied_element(a).or_else(|| implied_element(b)),
            _ => None,
        }
    }
    if query.atoms.iter().any(|a| has_chirality(&a.query)) {
        return Err(
            "tetrahedral stereo (@/@@) in a SMARTS-only reactant template is not \
                    interpreted; spell the atom as a SMILES bracket atom or drop the stereo",
        );
    }
    if !allow_direction && query.bonds.iter().any(|b| has_direction(&b.query)) {
        return Err(
            "double-bond stereo (/ or \\) in a SMARTS-only reactant template is not \
                    interpreted; spell the template with SMILES bonds or drop the stereo",
        );
    }
    let mut builder = MoleculeBuilder::new();
    for qa in &query.atoms {
        let mut atom = match implied_element(&qa.query) {
            Some(element) => chematic_core::Atom::new(element),
            None => chematic_core::Atom::wildcard(),
        };
        atom.atom_map = qa.atom_map;
        builder.add_atom(atom);
    }
    for qb in &query.bonds {
        let _ = builder.add_bond(
            AtomIdx(qb.atom1 as u32),
            AtomIdx(qb.atom2 as u32),
            BondOrder::Single,
        );
    }
    Ok(builder.build())
}

/// A product bond spelled as a SMARTS expression takes the bond type RDKit's
/// SMARTS parser gives it: the first primitive's order (`=,:` is double,
/// `-;!@` single, `!=` double: a leading `!` is skipped as RDKit does), or
/// no order when the expression starts with a ring primitive (`!@-`, `@=`,
/// `@`). An order-less bond is written `~` here and becomes a zero-order
/// bond in [`zero_order_product_bonds`]. Only text outside bracket atoms is
/// rewritten; a dative `->`/`<-` is left alone.
fn rdkit_product_bonds(products: &str) -> String {
    const EXPR: &[u8] = b"-=#:$/\\~@!,;&";
    let b = products.as_bytes();
    let mut out = String::with_capacity(products.len());
    let mut depth = 0usize;
    let mut i = 0;
    while i < b.len() {
        let c = b[i];
        if c == b'[' {
            depth += 1;
        } else if c == b']' {
            depth = depth.saturating_sub(1);
        } else if depth == 0 && EXPR.contains(&c) {
            let mut j = i;
            while j < b.len() && EXPR.contains(&b[j]) {
                j += 1;
            }
            let run = &products[i..j];
            let dative = b.get(j) == Some(&b'>') || (i > 0 && b[i - 1] == b'<');
            if (run.len() > 1 || run == "@") && !dative {
                out.push(rdkit_bond_symbol(run));
                i = j;
                continue;
            }
            out.push_str(run);
            i = j;
            continue;
        }
        out.push(c as char);
        i += 1;
    }
    out
}

/// The single bond symbol for a product bond expression (see
/// [`rdkit_product_bonds`]): its first primitive's order, or `~`.
fn rdkit_bond_symbol(expr: &str) -> char {
    expr.trim_start_matches('!')
        .chars()
        .next()
        .filter(|c| "-=#:$/\\".contains(*c))
        .unwrap_or('~')
}

/// Whether RDKit gives the product bond `q` no order: its first primitive
/// is a ring primitive or `~` inside a larger expression (a lone `~` is
/// RDKit's "any" bond, which keeps a reactant bond).
fn zero_order_product_bond(q: &BondQuery) -> bool {
    fn first(q: &BondQuery) -> &BondQuery {
        match q {
            BondQuery::And(a, _) | BondQuery::Or(a, _) => first(a),
            BondQuery::Not(a) => first(a),
            other => other,
        }
    }
    let lone_any = matches!(q, BondQuery::Primitive(BondPrimitive::Any));
    match first(q) {
        BondQuery::Primitive(BondPrimitive::Ring) => true,
        BondQuery::Primitive(BondPrimitive::Any) => !lone_any,
        _ => false,
    }
}

/// Set the bonds of a parsed product component that RDKit builds without an
/// order (see [`rdkit_product_bonds`]) to [`BondOrder::Zero`]. `original` is
/// the component as written, read as SMARTS for its bond expressions; atoms
/// are numbered in text order by both parsers, so bonds pair by endpoints.
fn zero_order_product_bonds(component: &mut Molecule, original: &str) {
    let Ok(query) = chematic_smarts::parse_smarts(original) else {
        return;
    };
    if query.atoms.len() != component.atom_count() {
        return;
    }
    for qb in &query.bonds {
        if !zero_order_product_bond(&qb.query) {
            continue;
        }
        let (a, b) = (AtomIdx(qb.atom1 as u32), AtomIdx(qb.atom2 as u32));
        if let Some((bidx, _)) = component.bond_between(a, b) {
            component.set_bond_order(bidx, BondOrder::Zero);
        }
    }
}

/// The product components when the product side uses component grouping
/// (`([C:1].[O:2]).[N:3]`), with each group's parentheses removed; `None`
/// when it does not.
fn grouped_product_components(products: &str) -> Result<Option<Vec<String>>, TransformError> {
    let b = products.as_bytes();
    let (mut components, mut any_group) = (Vec::new(), false);
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'(' {
            // A component starting with `(` is a group: up to the matching `)`.
            let (mut depth, mut bracket, mut close) = (0usize, 0usize, None);
            for (k, &c) in b.iter().enumerate().skip(i) {
                match c {
                    b'[' => bracket += 1,
                    b']' => bracket = bracket.saturating_sub(1),
                    b'(' if bracket == 0 => depth += 1,
                    b')' if bracket == 0 => {
                        depth -= 1;
                        if depth == 0 {
                            close = Some(k);
                            break;
                        }
                    }
                    _ => {}
                }
            }
            let close = close.ok_or_else(|| {
                TransformError::SmirksParse(RxnError::SmilesParse {
                    part: products.to_string(),
                    source: "unclosed product component group".to_string(),
                })
            })?;
            if close + 1 < b.len() && b[close + 1] != b'.' {
                return Err(TransformError::SmirksParse(RxnError::SmilesParse {
                    part: products.to_string(),
                    source: "a product component group must be followed by `.` or end".into(),
                }));
            }
            components.push(products[i + 1..close].to_string());
            any_group = true;
            i = close + 2;
        } else {
            let mut bracket = 0usize;
            let end = b[i..]
                .iter()
                .position(|&c| {
                    match c {
                        b'[' => bracket += 1,
                        b']' => bracket = bracket.saturating_sub(1),
                        _ => {}
                    }
                    c == b'.' && bracket == 0
                })
                .map_or(b.len(), |p| i + p);
            if end > i {
                components.push(products[i..end].to_string());
            }
            i = end + 1;
        }
    }
    Ok(any_group.then_some(components))
}

/// Product templates are SMILES specifications. Accept the SMARTS spelling of
/// an explicit H count, `[<atom>;H<n>[:map]]` → `[<atom>H<n>[:map]]`; refuse
/// any other `;`, `,`, `!`, `&` or `$(` inside a product bracket atom with a
/// typed error instead of the SMILES parser's position message.
fn normalize_product_templates(products: &str) -> Result<String, TransformError> {
    let mut out = String::with_capacity(products.len());
    let mut rest = products;
    while let Some(open) = rest.find('[') {
        out.push_str(&rest[..open]);
        let Some(close_rel) = rest[open + 1..].find(']') else {
            out.push_str(&rest[open..]);
            return Ok(out);
        };
        let inner = &rest[open + 1..open + 1 + close_rel];
        rest = &rest[open + 1 + close_rel + 1..];
        let rewritten = if inner.starts_with('#') {
            // Atomic-number primitives are handled by
            // `expand_atomic_number_primitives` before this point.
            inner.to_string()
        } else if let Some((atom, tail)) = inner.split_once(';') {
            let (hcount, map) = match tail.split_once(':') {
                Some((h, map)) => (h, Some(map)),
                None => (tail, None),
            };
            let valid_h = hcount.strip_prefix('H').is_some_and(|n| {
                n.is_empty() || (n.len() <= 2 && n.bytes().all(|b| b.is_ascii_digit()))
            });
            let valid_map =
                map.is_none_or(|m| !m.is_empty() && m.bytes().all(|b| b.is_ascii_digit()));
            let valid_atom = !atom.is_empty() && !atom.contains([';', ',', '!', '&', '$', 'H']);
            if !(valid_h && valid_map && valid_atom) {
                return Err(TransformError::SmirksParse(
                    RxnError::UnsupportedProductPrimitive {
                        primitive: format!("[{inner}]"),
                    },
                ));
            }
            match map {
                Some(map) => format!("{atom}{hcount}:{map}"),
                None => format!("{atom}{hcount}"),
            }
        } else if inner.contains([',', '!', '&']) || inner.contains("$(") {
            return Err(TransformError::SmirksParse(
                RxnError::UnsupportedProductPrimitive {
                    primitive: format!("[{inner}]"),
                },
            ));
        } else {
            inner.to_string()
        };
        out.push('[');
        out.push_str(&rewritten);
        out.push(']');
    }
    out.push_str(rest);
    Ok(out)
}

/// `atom_map` number → (reactant slot index, matched `AtomIdx`), built from
/// one match's per-reactant maps and the reactant templates' own atom-map
/// annotations. Shared by [`ReactionMatch::atom_map_positions`] and
/// [`apply_match_impl`].
fn global_map_of(
    per_reactant: &[FxHashMap<usize, AtomIdx>],
    template_atom_maps: &[Vec<Option<u16>>],
) -> FxHashMap<u16, (usize, AtomIdx)> {
    let mut global_map: FxHashMap<u16, (usize, AtomIdx)> = FxHashMap::default();
    for (ri, match_map) in per_reactant.iter().enumerate() {
        for (&qi, &t_idx) in match_map {
            if let Some(am) = template_atom_maps[ri][qi] {
                global_map.insert(am, (ri, t_idx));
            }
        }
    }
    global_map
}

/// Steps 1–2 of the original `run_reactants_impl`: VF2-match every reactant
/// template against its input molecule, take the cartesian product across
/// template slots, and keep only combinations that survive the
/// chirality/E-Z stereo post-checks. Does not build any product.
fn find_matches_impl(
    prepared: &PreparedReaction,
    reactants: &[&Molecule],
    limits: &ReactionTransformLimits,
    rings: Option<&[&RingSet]>,
) -> Result<Vec<ReactionMatch>, TransformError> {
    Ok(find_matches_profile(prepared, reactants, limits, rings, Profile::Native)?.0)
}

/// [`find_matches_impl`] for either profile. Under RDKit's, reactant `@`/`@@`
/// do not restrict matching; the second value reports whether the native E/Z
/// check rejected a combination (RDKit would not have).
fn find_matches_profile(
    prepared: &PreparedReaction,
    reactants: &[&Molecule],
    limits: &ReactionTransformLimits,
    rings: Option<&[&RingSet]>,
    profile: Profile,
) -> Result<(Vec<ReactionMatch>, Option<ReactionCompatibilityUnsupported>), TransformError> {
    let views = aromatic_views(reactants);
    let perceived = perceived_refs(reactants, &views);
    let reactants: &[&Molecule] = &perceived;
    let n_templates = prepared.rxn.reactants.len();
    if reactants.len() != n_templates {
        return Err(TransformError::ReactantCountMismatch {
            expected: n_templates,
            got: reactants.len(),
        });
    }
    if let Some(rings) = rings
        && rings.len() != n_templates
    {
        return Err(TransformError::ReactantCountMismatch {
            expected: n_templates,
            got: rings.len(),
        });
    }

    // VF2 match: for each (template_query, input_mol) pair.
    // A substructure API may deduplicate embeddings by target atom *set*.
    // Reactions cannot: swapping two template atoms mapped to the same set
    // may produce distinct products (e.g. asymmetric ether cleavage).
    let match_config = MatchConfig {
        uniquify: false,
        // A reaction reactant template is an exact molecular specification,
        // unlike the general substructure API's isotope-agnostic default.
        use_isotopes: true,
        ..MatchConfig::default()
    };
    let queries = match (profile, &prepared.rdkit_queries) {
        (Profile::Rdkit, Some(rdkit)) => rdkit,
        _ => &prepared.queries,
    };
    let mut ring_model_unresolved = false;
    let all_match_sets: Vec<Vec<FxHashMap<usize, AtomIdx>>> = queries
        .iter()
        .zip(reactants.iter())
        .enumerate()
        .map(|(index, (q, mol))| {
            let matches = match rings {
                // RDKit counts `R<n>` and ring sizes on its symmetrized ring
                // set (a quinuclidine N is in three rings there, two in the
                // SSSR); the RDKit profile uses the SMARTS opt-in matcher
                // that reads them that way.
                None if profile == Profile::Rdkit
                    && query_uses_ring_model(q)
                    && !rdkit_ring_model_is_sssr(mol) =>
                {
                    let config = chematic_smarts::RdkitParityConfig {
                        base: match_config.clone(),
                        ..chematic_smarts::RdkitParityConfig::default()
                    };
                    match chematic_smarts::find_matches_rdkit_parity(q, mol, &config) {
                        Ok((matches, _)) => matches,
                        Err(_) => {
                            ring_model_unresolved = true;
                            Vec::new()
                        }
                    }
                }
                Some(rings) => {
                    find_matches_with_rings_and_config(q, mol, rings[index], &match_config)
                }
                None => find_matches_with_config(q, mol, &match_config),
            };
            crate::perf_counters::record_reactant_query_match_call(matches.len());
            if matches.len() > limits.max_matches {
                return Err(TransformError::ResourceLimit {
                    resource: "reaction matches",
                    actual: matches.len(),
                    limit: limits.max_matches,
                });
            }
            Ok(matches)
        })
        .collect::<Result<_, _>>()?;

    if ring_model_unresolved {
        return Ok((
            vec![],
            Some(ReactionCompatibilityUnsupported::RingModelUnresolved),
        ));
    }
    // No matches when any template has no match.
    if all_match_sets.iter().any(|ms| ms.is_empty()) {
        return Ok((vec![], None));
    }

    let total_combinations = all_match_sets
        .iter()
        .map(|matches| matches.len())
        .try_fold(1usize, |total, count| total.checked_mul(count))
        .unwrap_or(usize::MAX);
    if total_combinations > limits.max_matches {
        return Err(TransformError::ResourceLimit {
            resource: "reaction match combinations",
            actual: total_combinations,
            limit: limits.max_matches,
        });
    }

    let mut matches: Vec<ReactionMatch> = Vec::new();
    let mut ez_rejected = false;

    for combo in cartesian_product(&all_match_sets) {
        crate::perf_counters::record_match_combination();

        // Parity-aware chirality post-check.  Runs only when the SMIRKS has @/@@.
        // This must happen after the complete VF2 mapping is known, because
        // correct chirality comparison requires the full neighbor permutation.
        if prepared.has_stereo && profile == Profile::Native {
            let ok = (0..prepared.rxn.reactants.len()).all(|ri| {
                smirks_chirality_ok(&prepared.rxn.reactants[ri], reactants[ri], &combo[ri])
            });
            if !ok {
                continue;
            }
        }
        // E/Z double-bond stereo post-check.  Runs only when the SMIRKS has /\.
        if prepared.has_ez_stereo {
            let ok = (0..prepared.rxn.reactants.len()).all(|ri| {
                smirks_ez_stereo_ok(&prepared.rxn.reactants[ri], reactants[ri], &combo[ri])
            });
            if !ok {
                ez_rejected = true;
                continue;
            }
        }

        matches.push(ReactionMatch {
            per_reactant: combo,
        });
    }

    Ok((
        matches,
        ez_rejected.then_some(ReactionCompatibilityUnsupported::EzReactantTemplateSemantics),
    ))
}

/// Whether RDKit's ring models give the same rings as the SSSR for `mol`
/// (no extra symmetrized rings, no organometallic view), so the ordinary
/// matcher answers `R<n>` and ring-size queries as RDKit does. Memoized on
/// the (cached) aromatic view; only cages such as quinuclidine or adamantane
/// need the slower RDKit-parity matcher.
fn rdkit_ring_model_is_sssr(mol: &Molecule) -> bool {
    *mol.derived(chematic_core::DerivedSlot::RdkitRingModelIsSssr, || {
        if mol
            .atoms()
            .any(|(_, a)| a.element == chematic_core::Element::FE)
        {
            return false;
        }
        let rings = chematic_perception::find_sssr(mol);
        let budget = chematic_smarts::RdkitRingModelBudget::default();
        let Ok(model) = chematic_smarts::build_rdkit_parity_ring_model(mol, &rings, &budget) else {
            return false;
        };
        if model.extra_ring_count() > 0 {
            return false;
        }
        let symmetrized = chematic_perception::find_symmetrized_sssr_with_diagnostics_bounded(
            mol,
            Some(budget.max_candidates),
        );
        symmetrized.status() != chematic_perception::SymmetrizedSssrStatus::CapExhausted
            && symmetrized.rings().rings().len() == rings.rings().len()
    })
}

/// Whether a reactant template asks for a positive ring count or a ring size,
/// the primitives RDKit reads on its symmetrized ring set.
fn query_uses_ring_model(query: &chematic_smarts::QueryMolecule) -> bool {
    use chematic_smarts::{AtomPrimitive, AtomQuery};
    fn atom(q: &AtomQuery) -> bool {
        match q {
            AtomQuery::Primitive(AtomPrimitive::RingCount(n)) => *n > 0,
            AtomQuery::Primitive(AtomPrimitive::RingSize(_)) => true,
            AtomQuery::Primitive(AtomPrimitive::Recursive(sub)) => query_uses_ring_model(sub),
            AtomQuery::Primitive(_) => false,
            AtomQuery::And(a, b) | AtomQuery::Or(a, b) => atom(a) || atom(b),
            AtomQuery::Not(a) => atom(a),
        }
    }
    query.atoms.iter().any(|a| atom(&a.query))
}

/// The RDKit-parity aromatic view of each reactant, where it differs from
/// the reactant (`None`: the reactant is its own view). RDKit's SMILES parser
/// aromatizes reactants before a reaction sees them, so a Kekulé benzene or
/// a 4-pyranone matches `[c]`, not `[C]` (#754); `smarts_find` already
/// matches on this view. Atom and bond indices and stereo records are
/// unchanged.
type AromaticView = Option<std::sync::Arc<Result<Molecule, chematic_perception::AromaticityError>>>;

fn aromatic_views(reactants: &[&Molecule]) -> Vec<AromaticView> {
    reactants
        .iter()
        .map(|m| {
            (!chematic_perception::rdkit_parity_view_is_identity(m))
                .then(|| chematic_perception::apply_aromaticity_rdkit_parity_shared(m))
        })
        .collect()
}

fn perceived_refs<'a>(reactants: &[&'a Molecule], views: &'a [AromaticView]) -> Vec<&'a Molecule> {
    reactants
        .iter()
        .zip(views)
        .map(|(&m, view)| match view {
            Some(view) => view.as_ref().as_ref().unwrap_or(m),
            None => m,
        })
        .collect()
}

/// Step 3 of the original `run_reactants_impl`: build the product set for
/// one already-accepted match and apply the valence filter. `None` means
/// the product set contained an over-valenced atom.
fn apply_match_impl(
    prepared: &PreparedReaction,
    reactants: &[&Molecule],
    m: &ReactionMatch,
    carry_substituents: bool,
) -> Option<Vec<Molecule>> {
    apply_match_traced_impl(prepared, reactants, m, carry_substituents)
        .map(|products| products.into_iter().map(|p| p.molecule).collect())
}

fn apply_match_traced_impl(
    prepared: &PreparedReaction,
    reactants: &[&Molecule],
    m: &ReactionMatch,
    carry_substituents: bool,
) -> Option<Vec<TracedProduct>> {
    match apply_match_profile(prepared, reactants, m, carry_substituents, Profile::Native) {
        Ok(ProductSetOutcome::Valid(set)) => Some(set),
        Ok(ProductSetOutcome::Rejected(..)) => None,
        Err(_) => unreachable!("the native profile never declines"),
    }
}

fn apply_match_profile(
    prepared: &PreparedReaction,
    reactants: &[&Molecule],
    m: &ReactionMatch,
    carry_substituents: bool,
    profile: Profile,
) -> Result<ProductSetOutcome, ReactionCompatibilityUnsupported> {
    let originals = reactants;
    let views = aromatic_views(reactants);
    let perceived = perceived_refs(reactants, &views);
    let reactants: &[&Molecule] = &perceived;
    let stereo = match profile {
        Profile::Native => None,
        Profile::Rdkit => Some(
            prepared
                .rdkit_stereo()
                .ok_or(ReactionCompatibilityUnsupported::ChiralReactantTemplateSemantics)?,
        ),
    };
    // global_map: atom_map_number → (reactant_mol_idx, matched_AtomIdx)
    let global_map = global_map_of(&m.per_reactant, &prepared.template_atom_maps);

    // all_template_atoms: every (mol_idx, AtomIdx) matched by any reactant template atom.
    // Used as BFS walls to prevent substituent collection from crossing into the
    // template region, and to identify bonds that the product template replaces.
    let mut all_template_atoms: FxHashSet<(usize, AtomIdx)> = FxHashSet::default();
    for (ri, match_map) in m.per_reactant.iter().enumerate() {
        for &t_idx in match_map.values() {
            all_template_atoms.insert((ri, t_idx));
        }
    }

    let products: Vec<TracedProduct> = prepared
        .rxn
        .products
        .iter()
        .enumerate()
        .map(|(pi, pt)| {
            let (mut product, template_idx_to_new) = build_product(
                pt,
                &global_map,
                &prepared.reactant_map_atoms,
                &prepared.reactant_bond_maps,
                reactants,
                &all_template_atoms,
                carry_substituents,
                prepared.product_specs.get(pi).map_or(&[], Vec::as_slice),
                profile == Profile::Rdkit,
            );
            if let Some(stereo) = stereo {
                let src_to_new: FxHashMap<(usize, AtomIdx), AtomIdx> = product
                    .atom_sources
                    .iter()
                    .enumerate()
                    .filter_map(|(i, s)| s.map(|s| ((s.reactant, s.atom), AtomIdx(i as u32))))
                    .collect();
                crate::rdkit_stereo::apply_product_stereo(
                    stereo,
                    reactants,
                    &m.per_reactant,
                    pi,
                    &template_idx_to_new,
                    &src_to_new,
                    &mut product.molecule,
                )
                .map_err(|_| ReactionCompatibilityUnsupported::AmbiguousStereoBondOrder)?;
                rdkit_carried_ez(&mut product, reactants, pt, &template_idx_to_new);
            }
            // The same provenance map serves traced and ordinary reaction
            // APIs. Born atoms remain untagged; caller labels are not remapped
            // or made unique across reactants.
            for (i, source) in product.atom_sources.iter().enumerate() {
                if let Some(source) = source
                    && let Some(tag) = originals[source.reactant].atom_tag(source.atom)
                {
                    product.molecule.set_tag(AtomIdx(i as u32), Some(tag.get()));
                }
            }
            crate::perf_counters::record_build_product_call(
                product.molecule.atom_count(),
                product.molecule.bond_count(),
            );
            Ok(product)
        })
        .collect::<Result<_, ReactionCompatibilityUnsupported>>()?;

    // A template that writes Kekulé bonds into an aromatic ring (`=,:` takes
    // its first alternative) leaves atoms flagged aromatic beside plain
    // double bonds, as in RDKit's raw product; RDKit's sanitize kekulizes
    // what is still aromatic and perceives aromaticity again. Do the same
    // before judging the product.
    let mut products = products;
    for product in &mut products {
        if mixed_aromatic_state(&product.molecule)
            && mixed_state_introduced(product, reactants)
            && let Some(view) = rdkit_sanitized_view(&product.molecule)
        {
            product.molecule = view;
        }
    }

    // Skip product sets RDKit's sanitize step would reject (issue #734).
    let sanitization: Vec<ProductSanitizationDiagnostic> = products
        .iter()
        .enumerate()
        .map(|(product_index, product)| product_sanitization_diagnostic(product, product_index))
        .collect();
    if sanitization.iter().all(|item| item.accepted) {
        crate::perf_counters::record_product_set();
        Ok(ProductSetOutcome::Valid(products))
    } else {
        Ok(ProductSetOutcome::Rejected(products, sanitization))
    }
}

/// RDKit's legacy CIP ranks (`assignStereochemistry`'s atom invariants,
/// which pick a double bond's stereo atoms): atomic number and isotope,
/// refined by the sorted ranks of the neighbours, each counted by bond order
/// (twice the order, so a double bond's partner counts twice), and implicit
/// H as the lowest rank. Explicit H atoms are ordinary neighbours, so an
/// explicit-H molecule ranks as its implicit-H form.
fn rdkit_legacy_cip_ranks(mol: &Molecule) -> Vec<u32> {
    let n = mol.atom_count();
    let key0 = |i: usize| {
        let a = mol.atom(AtomIdx(i as u32));
        (
            u32::from(a.element.atomic_number()),
            u32::from(a.isotope.unwrap_or(0)),
        )
    };
    let mut keys0: Vec<(u32, u32)> = (0..n).map(key0).collect();
    keys0.sort_unstable();
    keys0.dedup();
    let mut ranks: Vec<u32> = (0..n)
        .map(|i| keys0.binary_search(&key0(i)).unwrap_or(0) as u32 + 1)
        .collect();
    let mut distinct = keys0.len();
    for _ in 0..n {
        let keys: Vec<(u32, Vec<u32>)> = (0..n)
            .map(|i| {
                let idx = AtomIdx(i as u32);
                let mut nbrs: Vec<u32> = Vec::new();
                for (nb, b) in mol.neighbors(idx) {
                    let count = match mol.bond(b).order {
                        BondOrder::Double => 4,
                        BondOrder::Triple => 6,
                        BondOrder::Aromatic => 3,
                        BondOrder::Zero => 0,
                        _ => 2,
                    };
                    nbrs.extend(std::iter::repeat_n(ranks[nb.0 as usize], count));
                }
                // Implicit H: rank 0, below every atom.
                let h = chematic_core::implicit_hcount(mol, idx) as usize;
                nbrs.extend(std::iter::repeat_n(0, 2 * h));
                nbrs.sort_unstable_by(|a, b| b.cmp(a));
                (ranks[i], nbrs)
            })
            .collect();
        let mut sorted: Vec<&(u32, Vec<u32>)> = keys.iter().collect();
        sorted.sort_unstable();
        sorted.dedup();
        let next: Vec<u32> = keys
            .iter()
            .map(|k| sorted.binary_search(&k).unwrap_or(0) as u32 + 1)
            .collect();
        let now = sorted.len();
        ranks = next;
        if now == distinct {
            break;
        }
        distinct = now;
    }
    ranks
}

/// RDKit's reading of a reactant double bond's E/Z at one end: its stereo
/// atom (the substituent of higher legacy CIP rank, which RDKit records
/// whatever the markers), that atom's side seen from the end, and the end's
/// other substituent. `None` when the end has no stereo atom.
fn rdkit_reactant_ez_end(
    mol: &Molecule,
    ranks: &[u32],
    double: BondIdx,
    end: AtomIdx,
) -> Option<(AtomIdx, bool, Option<AtomIdx>)> {
    use chematic_core::ez_markers::marker_side;
    let subs: Vec<(AtomIdx, BondIdx)> = mol
        .neighbors(end)
        .filter(|&(_, b)| b != double && mol.bond(b).order != BondOrder::Double)
        .collect();
    let (stereo, other) = match subs.as_slice() {
        [one] => (*one, None),
        [a, b] => match ranks[a.0.0 as usize].cmp(&ranks[b.0.0 as usize]) {
            std::cmp::Ordering::Less => (*b, Some(*a)),
            std::cmp::Ordering::Greater => (*a, Some(*b)),
            std::cmp::Ordering::Equal => return None,
        },
        _ => return None,
    };
    let side = marker_side(mol, stereo.1, end)
        .or_else(|| other.and_then(|o| marker_side(mol, o.1, end)).map(|up| !up))?;
    Some((stereo.0, side, other.map(|o| o.0)))
}

/// E/Z of product double bonds copied from a reactant double bond, as
/// RDKit's reaction runner carries it: the reactant bond's stereo atoms are
/// followed into the product. A stereo atom left out of the product is
/// replaced by the end's other substituent, on the other side (the E/Z
/// survives losing the bond that carried its marker); a stereo atom still in
/// the product but no longer bonded to its end drops the E/Z (a ring
/// expansion inserting an atom on that side). A template double bond with
/// template neighbours at both ends takes the template's stereo: the
/// template's own when it spells one, else none.
fn rdkit_carried_ez(
    product: &mut TracedProduct,
    reactants: &[&Molecule],
    template: &Molecule,
    template_idx_to_new: &[Option<AtomIdx>],
) {
    use chematic_core::ez_markers::{
        EzFact, EzFactEnd, raw_direction, read_ez_facts, write_ez_facts,
    };
    if !reactants
        .iter()
        .any(|r| r.bonds().any(|(b, _)| raw_direction(r, b).is_some()))
    {
        return;
    }
    let mol = &product.molecule;
    let template_atom_of: FxHashMap<AtomIdx, AtomIdx> = template_idx_to_new
        .iter()
        .enumerate()
        .filter_map(|(t, p)| p.map(|p| (p, AtomIdx(t as u32))))
        .collect();
    // The template double bond between two product atoms, with whether each
    // template end has another template neighbour and a marked one.
    let template_double = |a: AtomIdx, b: AtomIdx| -> Option<(bool, bool)> {
        let (&ta, &tb) = (template_atom_of.get(&a)?, template_atom_of.get(&b)?);
        let (tbond, entry) = template.bond_between(ta, tb)?;
        (entry.order == BondOrder::Double).then(|| {
            let ends = [ta, tb];
            let neighbours = ends
                .iter()
                .all(|&end| template.neighbors(end).any(|(_, nb)| nb != tbond));
            let marked = ends.iter().all(|&end| {
                template
                    .neighbors(end)
                    .any(|(_, nb)| nb != tbond && raw_direction(template, nb).is_some())
            });
            (neighbours, marked)
        })
    };
    let product_of = |ri: usize, atom: AtomIdx| -> Option<AtomIdx> {
        product
            .atom_sources
            .iter()
            .position(|s| s.is_some_and(|s| s.reactant == ri && s.atom == atom))
            .map(|i| AtomIdx(i as u32))
    };
    let mut dropped: Vec<BondIdx> = Vec::new();
    let mut derived: Vec<EzFact> = Vec::new();
    let mut rank_cache: FxHashMap<usize, Vec<u32>> = FxHashMap::default();
    for (pb, bond) in mol.bonds() {
        if bond.order != BondOrder::Double {
            continue;
        }
        let (pa, pz) = (bond.atom1, bond.atom2);
        let (Some(sa), Some(sz)) = (
            product.atom_sources[pa.0 as usize],
            product.atom_sources[pz.0 as usize],
        ) else {
            continue;
        };
        if sa.reactant != sz.reactant {
            continue;
        }
        let r = reactants[sa.reactant];
        let Some((rb, rbond)) = r.bond_between(sa.atom, sz.atom) else {
            continue;
        };
        if rbond.order != BondOrder::Double || r.atom(sa.atom).aromatic && r.atom(sz.atom).aromatic
        {
            continue;
        }
        let ranks = rank_cache
            .entry(sa.reactant)
            .or_insert_with(|| rdkit_legacy_cip_ranks(r));
        let (Some(ea), Some(ez)) = (
            rdkit_reactant_ez_end(r, ranks, rb, sa.atom),
            rdkit_reactant_ez_end(r, ranks, rb, sz.atom),
        ) else {
            continue;
        };
        match template_double(pa, pz) {
            // The template spells this double bond's stereo: keep it.
            Some((_, true)) => continue,
            // `[C:4][C:1]=[C:2][C:3]>>[C:4][C:1]=[C:2][C:3]` on an alkene,
            // BTMR0072 on a bis-allylic CH2: no E/Z.
            Some((true, false)) => {
                dropped.push(pb);
                continue;
            }
            _ => {}
        }
        let mut ends: Vec<EzFactEnd> = Vec::with_capacity(2);
        let mut drop = false;
        for (p_end, (stereo, side, other)) in [(pa, ea), (pz, ez)] {
            let reference = match product_of(sa.reactant, stereo) {
                Some(ps) => match mol.bond_between(p_end, ps) {
                    Some((b, _)) => Some((b, side)),
                    None => {
                        drop = true;
                        None
                    }
                },
                // Only the end's other reactant substituent replaces it; an
                // atom the template adds does not (a vinyl iodide reduced to
                // `C=C[H]` has no E/Z, BTMR0671).
                None => other
                    .and_then(|o| product_of(sa.reactant, o))
                    .and_then(|po| mol.bond_between(p_end, po))
                    .map(|(b, _)| (b, !side)),
            };
            match reference {
                Some((reference, up)) => ends.push(EzFactEnd {
                    atom: p_end,
                    reference,
                    up,
                }),
                None => {
                    drop = true;
                    break;
                }
            }
        }
        if drop {
            dropped.push(pb);
        } else if let [a, b] = ends.as_slice() {
            derived.push(EzFact {
                double: pb,
                ends: [*a, *b],
            });
        }
    }
    if dropped.is_empty() && derived.is_empty() {
        return;
    }
    // Markers of dropped double bonds go first, so they cannot be read as
    // another double bond's (or conflict with its).
    let mut trial = product.molecule.clone();
    let keep: FxHashSet<BondIdx> = derived
        .iter()
        .flat_map(|f| f.ends.iter().map(|e| e.reference))
        .collect();
    for &pb in &dropped {
        let bond = trial.bond(pb);
        for end in [bond.atom1, bond.atom2] {
            let subs: Vec<BondIdx> = trial
                .neighbors(end)
                .filter(|&(_, b)| b != pb && !keep.contains(&b))
                .map(|(_, b)| b)
                .collect();
            for b in subs {
                if matches!(trial.bond(b).order, BondOrder::Up | BondOrder::Down) {
                    trial.set_bond_order(b, BondOrder::Single);
                }
            }
        }
    }
    let mut facts = read_ez_facts(&trial).map(|(f, _)| f).unwrap_or_default();
    facts.retain(|f| !dropped.contains(&f.double) && !derived.iter().any(|d| d.double == f.double));
    facts.extend(derived);
    if write_ez_facts(&mut trial, &facts) {
        product.molecule = trial;
    }
}

/// Atoms flagged aromatic with a plain double bond between two of them or
/// with fewer than two aromatic bonds, or an aromatic bond on an atom not
/// flagged aromatic: a template rewrote part of an aromatic system.
fn mixed_aromatic_state(mol: &Molecule) -> bool {
    mol.bonds().any(|(_, b)| {
        b.order == BondOrder::Double && mol.atom(b.atom1).aromatic && mol.atom(b.atom2).aromatic
    }) || mol.atoms().any(|(idx, a)| {
        let aromatic_bonds = mol
            .neighbors(idx)
            .filter(|&(_, b)| mol.bond(b).order == BondOrder::Aromatic)
            .count();
        (a.aromatic && aromatic_bonds < 2) || (!a.aromatic && aromatic_bonds > 0)
    })
}

/// Per atom, whether it is in a mixed aromatic state (see
/// [`mixed_aromatic_state`]).
fn mixed_aromatic_atoms(mol: &Molecule) -> Vec<bool> {
    let mut out = vec![false; mol.atom_count()];
    for (_, b) in mol.bonds() {
        if b.order == BondOrder::Double && mol.atom(b.atom1).aromatic && mol.atom(b.atom2).aromatic
        {
            out[b.atom1.0 as usize] = true;
            out[b.atom2.0 as usize] = true;
        }
    }
    for (idx, a) in mol.atoms() {
        let aromatic_bonds = mol
            .neighbors(idx)
            .filter(|&(_, b)| mol.bond(b).order == BondOrder::Aromatic)
            .count();
        if (a.aromatic && aromatic_bonds < 2) || (!a.aromatic && aromatic_bonds > 0) {
            out[idx.0 as usize] = true;
        }
    }
    out
}

/// Whether the template made a product's mixed aromatic state: some mixed
/// atom is new or was not mixed in its reactant. A state carried over from
/// the reactant unchanged (an aromatic ring written with `#`) stays as
/// RDKit's product keeps it.
fn mixed_state_introduced(product: &TracedProduct, reactants: &[&Molecule]) -> bool {
    let mixed = mixed_aromatic_atoms(&product.molecule);
    let mut reactant_mixed: Vec<Option<Vec<bool>>> = vec![None; reactants.len()];
    mixed.iter().enumerate().any(|(i, &m)| {
        if !m {
            return false;
        }
        let Some(src) = product.atom_sources.get(i).copied().flatten() else {
            return true;
        };
        let rm = reactant_mixed[src.reactant]
            .get_or_insert_with(|| mixed_aromatic_atoms(reactants[src.reactant]));
        !rm[src.atom.0 as usize]
    })
}

/// `mol` after RDKit's sanitize steps for a mixed aromatic state: the
/// remaining aromatic bonds kekulized, aromatic flags cleared, aromaticity
/// perceived again with RDKit's model. `None` when the aromatic bonds do
/// not kekulize (RDKit's sanitize fails too).
fn rdkit_sanitized_view(mol: &Molecule) -> Option<Molecule> {
    // RDKit's kekulization refuses an aromatic atom outside a ring.
    let in_ring = atoms_in_rings(mol);
    if mol
        .atoms()
        .any(|(idx, a)| a.aromatic && !in_ring[idx.0 as usize])
    {
        return None;
    }
    // RDKit kekulizes the aromatic bonds; see [`rdkit_kekule_candidate`].
    let mut kekule = if mol.bonds().any(|(_, b)| b.order == BondOrder::Aromatic) {
        let k = chematic_core::kekulize_with(mol, |idx| rdkit_kekule_candidate(mol, idx)).ok()?;
        chematic_core::apply_kekule(mol, &k)
    } else {
        mol.clone()
    };
    let flagged: Vec<AtomIdx> = kekule
        .atoms()
        .filter(|(_, a)| a.aromatic)
        .map(|(i, _)| i)
        .collect();
    for idx in flagged {
        kekule.set_atom_aromatic(idx, false);
    }
    let mut view = match &*chematic_perception::apply_aromaticity_rdkit_parity_shared(&kekule) {
        Ok(view) => view.clone(),
        Err(_) => kekule,
    };
    // RDKit recomputes the H count of a carbon read without brackets: one
    // left non-aromatic takes the hydrogens its Kekulé bonds leave room for
    // (an explicit-H ring CH next to a carbon the template made sp3 becomes
    // CH2, BTMR0655).
    for (idx, atom) in mol.atoms() {
        let now = view.atom(idx);
        if atom.aromatic
            && atom.element == chematic_core::Element::C
            && atom.charge == 0
            && !now.aromatic
            && now.hydrogen_count.is_some()
        {
            view.set_hydrogen_count(idx, None);
        }
    }
    Some(view)
}

/// Whether RDKit's kekulization would give atom `idx` of `mol` a double
/// bond. RDKit makes an atom a candidate when its bond orders (aromatic
/// bonds as 1) and explicit H come to one less than its default valence,
/// `sbo == dv - 1`. That rule decides for an atom that already has a
/// non-aromatic double or triple bond (a benzene carbon whose ring bond a
/// template rewrote as `=` still is a candidate, RDKit's `C1=C=CC=CC=1`; a
/// ring carbon with a substituent or a second ring bond rewritten is not)
/// for a chalcogen cation (`[o+]`/`[s+]` with two ring bonds, as in
/// pyrylium, is one; a thiophene S-oxide `[O-][s+]1cccc1` is not) and for a
/// neutral N or P (pyridine's N is one; a ring-fusion N with three aromatic
/// bonds is not, also when its third ring neighbour is a `c(=O)` carbon, as
/// in imidazo[1,2-a]pyrimidin-5-one, where chematic's bridgehead fallback
/// does not apply). Other atoms follow chematic's ordinary rule.
fn rdkit_kekule_candidate(mol: &Molecule, idx: AtomIdx) -> bool {
    let atom = mol.atom(idx);
    let pi = mol
        .neighbors(idx)
        .any(|(_, b)| matches!(mol.bond(b).order, BondOrder::Double | BondOrder::Triple));
    let chalcogen_cation =
        matches!(atom.element.atomic_number(), 8 | 16 | 34 | 52) && atom.charge == 1;
    let neutral_pnictogen = matches!(atom.element.atomic_number(), 7 | 15)
        && atom.charge == 0
        && !mol
            .neighbors(idx)
            .any(|(nb, _)| mol.atom(nb).element.atomic_number() == 1);
    // A zero-order bond (RDKit's unspecified bond from a product `~` or
    // ring-led `@`) is not counted for a neutral N or P: an [nH] of an
    // explicit-H indole left with one aromatic bond and a zero-order bond
    // to its fused ring is then a candidate, as in RDKit, and the product
    // fails to kekulize (BTMR1031). Carbon keeps the ordinary rule there
    // (RDKit kekulizes `c1ccc~cc1`).
    let pnictogen_beside_zero = matches!(atom.element.atomic_number(), 7 | 15)
        && atom.charge == 0
        && mol
            .neighbors(idx)
            .any(|(_, b)| mol.bond(b).order == BondOrder::Zero);
    if !pi && !chalcogen_cation && !neutral_pnictogen && !pnictogen_beside_zero {
        // A three-connected aromatic P, As, S, Se or Te is judged by RDKit's
        // own rule; chematic's never makes it a candidate.
        if matches!(atom.element.atomic_number(), 15 | 16 | 33 | 34 | 52)
            && atom.charge == 0
            && mol.degree(idx) >= 3
            && let Some(candidate) = rdkit_mark_dbond_cand(mol, idx)
        {
            return candidate;
        }
        return chematic_core::atom_must_be_matched(mol, idx);
    }
    let sbo: i16 = mol
        .neighbors(idx)
        .map(|(_, b)| match mol.bond(b).order {
            BondOrder::Aromatic => 1,
            order => i16::from(order.order_int()),
        })
        .sum::<i16>()
        + i16::from(atom.hydrogen_count.unwrap_or(0));
    crate::rdkit_valence::default_valence(atom.element.atomic_number(), atom.charge)
        .is_some_and(|dv| sbo == dv - 1)
}

/// RDKit's `markDbondCands` decision for an aromatic atom with no double
/// or triple bond (`Kekulize.cpp`, RDKit 2026.03): aromatic bonds count 1
/// toward `sbo`, other bonds their order, plus all hydrogens; the default
/// valence, shifted by the charge (negated for a carbocation and for the
/// early elements), is raised through the element's valence list while the
/// total valence (aromatic bonds as 1.5) exceeds it; an atom whose degree
/// plus implicit hydrogens reaches that valence is no candidate; otherwise
/// it is one when the valence is `sbo + 1`, or `sbo + 2` for a bracket
/// atom. `None` for an element RDKit leaves unrestricted or a non-aromatic
/// atom. [`rdkit_kekule_candidate`] uses it for three-connected aromatic
/// P, As, S, Se and Te: a ring carbon a template rewrites to `[S:1]` keeps
/// the reactant's aromatic flag in RDKit, and 2-chloropyridine gives the
/// hypervalent `ClS1=CC=CC=N1` (#754 fuzz lane). Applied to every aromatic
/// atom it accepted products RDKit's sanitize rejects, so the other atoms
/// keep chematic's rule.
fn rdkit_mark_dbond_cand(mol: &Molecule, idx: AtomIdx) -> Option<bool> {
    let atom = mol.atom(idx);
    if !atom.aromatic || atom.wildcard {
        return None;
    }
    let z = atom.element.atomic_number();
    let list = crate::rdkit_valence::valence_list(z)?;
    let early = matches!(z, 1..=5 | 11..=13 | 19..=31 | 37..=49 | 55..=81 | 87..=112);
    let mut chrg = i16::from(atom.charge);
    if early || (z == 6 && chrg > 0) {
        chrg = -chrg;
    }
    let hydrogens = i16::from(chematic_core::implicit_hcount(mol, idx));
    let mut sbo = hydrogens;
    let mut total2 = 2 * hydrogens;
    let mut ignored = 0i16;
    for (_, b) in mol.neighbors(idx) {
        let bond = mol.bond(b);
        match bond.order {
            BondOrder::Aromatic => {
                sbo += 1;
                total2 += 3;
            }
            BondOrder::Zero | BondOrder::QueryAny => ignored += 1,
            BondOrder::Dative if bond.atom1 == idx => ignored += 1,
            order => {
                let o = i16::from(order.order_int());
                sbo += o;
                total2 += 2 * o;
            }
        }
    }
    let tbo = (total2 + 1) / 2;
    let mut dv = i16::from(list[0]) + chrg;
    for &v in &list[1..] {
        if tbo <= dv {
            break;
        }
        dv = i16::from(v) + chrg;
    }
    let implicit = if atom.hydrogen_count.is_none() {
        hydrogens
    } else {
        0
    };
    let total_degree = mol.degree(idx) as i16 + implicit - ignored;
    if total_degree >= dv {
        return Some(false);
    }
    Some(dv == sbo + 1 || (atom.hydrogen_count.is_some() && dv == sbo + 2))
}

/// Whether RDKit's sanitize step would accept `mol` (issue #734).
///
/// Follows RDKit's `SanitizeMol`: aromatic atoms must lie in rings and
/// aromatic systems must kekulize; RDKit's clean-up spellings are accepted
/// (a neutral five-valent N with `=O` or `#N` is read as `[N+][O-]` /
/// `[N+]=[N-]`, a neutral Cl/Br/I of valence 3, 5 or 7 bonded only to O as
/// `[X+k]([O-])…`); and each atom's explicit valence — Kekulé bond orders
/// (a dative bond counts for its acceptor only) plus its explicit H count —
/// must not exceed RDKit's largest allowed valence for its element and
/// charge (see [`crate::rdkit_valence`]).
fn product_sanitization_diagnostic(
    product: &TracedProduct,
    product_index: usize,
) -> ProductSanitizationDiagnostic {
    use crate::rdkit_valence::max_valence;
    const MAX_DIAGNOSTIC_ATOMS: usize = 32;

    let mol = &product.molecule;
    let atom_context =
        |idx: AtomIdx, observed_valence: Option<i16>, max_allowed_valence: Option<i16>| {
            let atom = mol.atom(idx);
            ProductAtomRejectionDiagnostic {
                atom: idx,
                template_map: product.template_maps.get(idx.0 as usize).copied().flatten(),
                element: atom.element,
                explicit_hydrogen_count: atom.hydrogen_count,
                observed_valence,
                max_allowed_valence,
            }
        };
    let rejected = |reason: ProductRejectionReason,
                    mut atoms: Vec<ProductAtomRejectionDiagnostic>| {
        atoms.sort_unstable_by_key(|atom| atom.atom.0);
        atoms.dedup_by_key(|atom| atom.atom.0);
        atoms.truncate(MAX_DIAGNOSTIC_ATOMS);
        ProductSanitizationDiagnostic {
            product_index,
            accepted: false,
            reason: Some(reason),
            atoms,
        }
    };

    if mol.atoms().any(|(_, a)| a.aromatic) {
        let in_ring = atoms_in_rings(mol);
        let outside_ring: Vec<_> = mol
            .atoms()
            .filter(|(idx, atom)| atom.aromatic && !in_ring[idx.0 as usize])
            .take(MAX_DIAGNOSTIC_ATOMS)
            .map(|(idx, _)| atom_context(idx, None, None))
            .collect();
        if !outside_ring.is_empty() {
            return rejected(ProductRejectionReason::Aromaticity, outside_ring);
        }
        // A plain double bond between two aromatic ring atoms (an aromatic
        // ring bond a template rewrote as `=`, `c1cc=ncc1`) is not a form
        // RDKit's kekulization accepts; the written SMILES would not parse.
        if let Some((_, bond)) = mol.bonds().find(|(_, bond)| {
            bond.order == BondOrder::Double
                && mol.atom(bond.atom1).aromatic
                && mol.atom(bond.atom2).aromatic
                && in_ring[bond.atom1.0 as usize]
                && in_ring[bond.atom2.0 as usize]
        }) {
            return rejected(
                ProductRejectionReason::Aromaticity,
                [bond.atom1, bond.atom2]
                    .into_iter()
                    .map(|idx| atom_context(idx, None, None))
                    .collect(),
            );
        }
    }
    // Kekulé orders of the aromatic bonds, read in place (the product is
    // not copied).
    let mut orders: Vec<BondOrder> = mol.bonds().map(|(_, b)| b.order).collect();
    if orders.contains(&BondOrder::Aromatic) {
        match chematic_core::kekulize_with(mol, |idx| rdkit_kekule_candidate(mol, idx)) {
            Ok(k) => {
                for (b, order) in k {
                    orders[b.0 as usize] = order;
                }
            }
            Err(_) => {
                let atoms = mol
                    .atoms()
                    .filter(|(_, atom)| atom.aromatic)
                    .take(MAX_DIAGNOSTIC_ATOMS)
                    .map(|(idx, _)| atom_context(idx, None, None))
                    .collect();
                return rejected(ProductRejectionReason::Kekulization, atoms);
            }
        }
    }
    let order_of = |b: BondIdx| orders[b.0 as usize];
    let explicit_valence = |idx: AtomIdx| -> i16 {
        let bonds: i16 = mol
            .neighbors(idx)
            .map(|(_, b)| {
                let bond = mol.bond(b);
                match order_of(b) {
                    BondOrder::Dative if bond.atom1 == idx => 0,
                    order => i16::from(order.order_int()),
                }
            })
            .sum();
        bonds + i16::from(mol.atom(idx).hydrogen_count.unwrap_or(0))
    };
    let mut invalid_atoms = Vec::new();
    for (idx, atom) in mol.atoms() {
        if atom.wildcard {
            continue;
        }
        let z = atom.element.atomic_number();
        let mut charge = atom.charge;
        let mut used = explicit_valence(idx);
        // RDKit's cleanUp, which runs before the valence check.
        if charge == 0 {
            let bond_to = |order: BondOrder, nz: u8| {
                mol.neighbors(idx).any(|(nb, b)| {
                    order_of(b) == order && {
                        let n = mol.atom(nb);
                        n.element.atomic_number() == nz && n.charge == 0
                    }
                })
            };
            if z == 7
                && used == 5
                && (bond_to(BondOrder::Double, 8) || bond_to(BondOrder::Triple, 7))
            {
                charge = 1;
                used = 4;
            } else if matches!(z, 17 | 35 | 53)
                && matches!(used, 3 | 5 | 7)
                && mol
                    .neighbors(idx)
                    .all(|(nb, _)| mol.atom(nb).element.atomic_number() == 8)
            {
                let doubles = mol
                    .neighbors(idx)
                    .filter(|&(_, b)| order_of(b) == BondOrder::Double)
                    .count() as i16;
                charge = doubles as i8;
                used -= doubles;
            }
        }
        let max = max_valence(z, charge);
        if max.is_some_and(|max| used > max) && invalid_atoms.len() < MAX_DIAGNOSTIC_ATOMS {
            invalid_atoms.push(atom_context(idx, Some(used), max));
        }
    }
    if !invalid_atoms.is_empty() {
        return rejected(ProductRejectionReason::Valence, invalid_atoms);
    }
    ProductSanitizationDiagnostic {
        product_index,
        accepted: true,
        reason: None,
        atoms: Vec::new(),
    }
}

/// Copy each reactant bond's stashed direction ([`Molecule::bond_direction`]
/// with its anchor) onto the product bond between the same two atoms, when
/// that bond is still single or aromatic and one of its atoms keeps a double
/// bond to a third atom (#734: `[Cl:1]>>[Cl:1]` dropped the E/Z of a
/// squaraine imine).
fn carry_direction_stash(
    product: &mut Molecule,
    input_mols: &[&Molecule],
    src_to_new: &FxHashMap<(usize, AtomIdx), AtomIdx>,
) {
    if !input_mols.iter().any(|m| m.has_bond_directions()) {
        return;
    }
    let mut source_of: FxHashMap<AtomIdx, (usize, AtomIdx)> =
        FxHashMap::with_capacity_and_hasher(src_to_new.len(), Default::default());
    for (&key, &new_idx) in src_to_new {
        source_of.insert(new_idx, key);
    }
    let stashed: Vec<(BondIdx, BondOrder, AtomIdx)> = product
        .bonds()
        .filter_map(|(pb, bond)| {
            if !matches!(bond.order, BondOrder::Single | BondOrder::Aromatic)
                || product.bond_direction(pb).is_some()
            {
                return None;
            }
            let &(ma, sa) = source_of.get(&bond.atom1)?;
            let &(mb, sb) = source_of.get(&bond.atom2)?;
            if ma != mb {
                return None;
            }
            let src = input_mols[ma];
            let (sbid, sbond) = src.bond_between(sa, sb)?;
            let dir = src.bond_direction(sbid)?;
            let anchor = src.bond_direction_anchor(sbid).unwrap_or(sbond.atom1);
            let new_anchor = *src_to_new.get(&(ma, anchor))?;
            let flanks_double = [bond.atom1, bond.atom2].into_iter().any(|end| {
                product
                    .neighbors(end)
                    .any(|(_, b)| b != pb && product.bond(b).order == BondOrder::Double)
            });
            flanks_double.then_some((pb, dir, new_anchor))
        })
        .collect();
    for (pb, dir, anchor) in stashed {
        product.set_bond_direction(pb, dir);
        product.set_bond_direction_anchor(pb, anchor);
    }
}

/// Per atom, whether it lies on a ring (has an incident non-bridge bond).
fn atoms_in_rings(mol: &Molecule) -> Vec<bool> {
    let n = mol.atom_count();
    let mut order = vec![usize::MAX; n];
    let mut low = vec![0usize; n];
    let mut in_ring = vec![false; n];
    let mut counter = 0usize;
    // (atom, bond used to reach it, index of its next neighbour)
    let mut stack: Vec<(usize, Option<BondIdx>, usize)> = Vec::with_capacity(n);
    for root in 0..n {
        if order[root] != usize::MAX {
            continue;
        }
        stack.clear();
        stack.push((root, None, 0));
        order[root] = counter;
        low[root] = counter;
        counter += 1;
        while let Some(top) = stack.last_mut() {
            let (v, via, k) = *top;
            if let Some((w, b)) = mol.neighbors(AtomIdx(v as u32)).nth(k) {
                top.2 += 1;
                if Some(b) == via {
                    continue;
                }
                let w = w.0 as usize;
                if order[w] == usize::MAX {
                    order[w] = counter;
                    low[w] = counter;
                    counter += 1;
                    stack.push((w, Some(b), 0));
                } else {
                    low[v] = low[v].min(order[w]);
                    // A back edge closes a cycle through both ends.
                    in_ring[v] = true;
                    in_ring[w] = true;
                }
            } else {
                stack.pop();
                if let Some(&(p, _, _)) = stack.last() {
                    low[p] = low[p].min(low[v]);
                    if low[v] <= order[p] {
                        // The tree edge p–v is not a bridge.
                        in_ring[v] = true;
                        in_ring[p] = true;
                    }
                }
            }
        }
    }
    in_ring
}

// ---------------------------------------------------------------------------
// SMIRKS chirality post-check (parity-aware)
// ---------------------------------------------------------------------------

/// Returns the parity of the permutation P that maps `from_seq` to `to_seq`:
/// `Some(true)` = even (same chirality sense), `Some(false)` = odd (inverted).
/// Returns `None` when the sequences differ in length or contain elements that
/// cannot be aligned (e.g. an unmapped neighbour).
fn permutation_parity(from_seq: &[u32], to_seq: &[u32]) -> Option<bool> {
    let n = from_seq.len();
    if n != to_seq.len() {
        return None;
    }
    // Build perm[j] = index i in from_seq where from_seq[i] == to_seq[j].
    let mut perm = Vec::with_capacity(n);
    for &t in to_seq {
        let pos = from_seq.iter().position(|&f| f == t)?;
        perm.push(pos);
    }
    // Count inversions to determine parity.
    let mut inv = 0usize;
    for i in 0..n {
        for j in (i + 1)..n {
            if perm[i] > perm[j] {
                inv += 1;
            }
        }
    }
    Some(inv.is_multiple_of(2)) // true = even = chirality flags must agree for same config
}

/// Parity-aware stereo check for one (template, reactant, mapping) triple.
///
/// For each chiral atom in `tmpl`, maps the template's recorded SMILES stereo
/// neighbor order through the VF2 `match_map` into the reactant atom-index space,
/// then computes the parity of the permutation relative to the reactant's recorded
/// SMILES stereo neighbor order.
///
/// - Even parity → chirality flags must agree for same absolute configuration.
/// - Odd parity  → chirality flags must differ for same absolute configuration.
///
/// Returns `true` if all chiral centres are consistent, `false` if any mismatch.
fn smirks_chirality_ok(
    tmpl: &Molecule,
    reactant: &Molecule,
    match_map: &FxHashMap<usize, AtomIdx>,
) -> bool {
    for i in 0..tmpl.atom_count() {
        let tmpl_atom = tmpl.atom(AtomIdx(i as u32));
        if tmpl_atom.chirality == Chirality::None {
            continue;
        }

        // Template atom's SMILES stereo neighbour order (template atom indices).
        let Some(tmpl_order) = tmpl.stereo_neighbor_order(AtomIdx(i as u32)) else {
            continue; // No recorded order — skip this centre.
        };

        // Corresponding matched reactant atom.
        let Some(&react_idx) = match_map.get(&i) else {
            continue; // Template atom not in mapping (shouldn't happen for complete match).
        };

        let react_atom = reactant.atom(react_idx);
        if react_atom.chirality == Chirality::None {
            return false; // Template requires stereo; reactant atom has none.
        }

        // Map each template stereo-neighbour index to the corresponding reactant atom index.
        let mut mapped: Vec<u32> = Vec::with_capacity(tmpl_order.len());
        let mut all_mapped = true;
        for &t in tmpl_order {
            if t == STEREO_H_SENTINEL {
                mapped.push(STEREO_H_SENTINEL);
            } else {
                match match_map.get(&(t as usize)) {
                    Some(ri) => mapped.push(ri.0),
                    None => {
                        all_mapped = false;
                        break;
                    }
                }
            }
        }
        if !all_mapped {
            continue; // Partial substructure match — cannot verify chirality.
        }

        // Reactant atom's SMILES stereo neighbour order (reactant atom indices).
        let Some(react_order) = reactant.stereo_neighbor_order(react_idx) else {
            // No recorded order in reactant — fall back to raw flag comparison.
            if react_atom.chirality != tmpl_atom.chirality {
                return false;
            }
            continue;
        };

        let Some(even_parity) = permutation_parity(&mapped, react_order) else {
            continue; // Alignment failed — skip.
        };

        // Check: even parity → flags must agree; odd parity → flags must differ.
        let same_flag = tmpl_atom.chirality == react_atom.chirality;
        if same_flag != even_parity {
            return false;
        }
    }
    true
}

// ---------------------------------------------------------------------------
// SMIRKS E/Z double-bond stereo post-check
// ---------------------------------------------------------------------------

/// Returns the "outward" stereo-bond direction from `atom` (one endpoint of a
/// double bond) toward its Up/Down-annotated substituent, or `None` if no
/// such bond exists.
///
/// E/Z stereo is encoded on the *substituent* bonds adjacent to a C=C, not on
/// the double bond itself.  Whether the stored bond goes *into* or *out of*
/// `atom` determines how to read its direction:
/// - Outgoing (`atom` is `bond.atom1`): direction is as stored.
/// - Incoming (`atom` is `bond.atom2`): direction is flipped (Up ↔ Down).
///
/// The `other` parameter is the opposite endpoint of the double bond; that
/// bond is skipped so we look only at substituents.
fn ez_stereo_outward(mol: &Molecule, atom: AtomIdx, other: AtomIdx) -> Option<BondOrder> {
    for (nb, bidx) in mol.neighbors(atom) {
        if nb == other {
            continue; // skip the double bond itself
        }
        let bond = mol.bond(bidx);
        match bond.order {
            BondOrder::Up | BondOrder::Down => {
                let outward = if bond.atom1 == atom {
                    // bond goes FROM atom outward → direction as stored
                    bond.order
                } else {
                    // bond comes INTO atom → flip to get outward direction
                    match bond.order {
                        BondOrder::Up => BondOrder::Down,
                        _ => BondOrder::Up,
                    }
                };
                return Some(outward);
            }
            _ => {}
        }
    }
    None
}

/// E/Z stereo post-check for one (template, reactant, mapping) triple.
///
/// For each double bond in `tmpl` that has Up/Down substituent bonds on **both**
/// sides, verify that the corresponding double bond in `reactant` (found via
/// the VF2 `match_map`) encodes the same E/Z parity.
///
/// Parity is determined by comparing the "outward direction" from each end of
/// the double bond (see [`ez_stereo_outward`]).  Two outward directions that are
/// equal → same-side (Z/cis); directions that differ → opposite-side (E/trans).
///
/// If only one side of a template double bond has a stereo bond (or the
/// reactant doesn't annotate stereo at all), the constraint is skipped — this
/// matches the behaviour of SMIRKS templates extracted from rdchiral where a
/// single-sided annotation is common.
///
/// Returns `true` if all constrained double bonds are consistent.
fn smirks_ez_stereo_ok(
    tmpl: &Molecule,
    reactant: &Molecule,
    match_map: &FxHashMap<usize, AtomIdx>,
) -> bool {
    for (_, bond) in tmpl.bonds() {
        if bond.order != BondOrder::Double {
            continue;
        }
        let ta = bond.atom1;
        let tb = bond.atom2;

        // Outward directions from each end of the template double bond.
        let sa = ez_stereo_outward(tmpl, ta, tb);
        let sb = ez_stereo_outward(tmpl, tb, ta);

        // Both sides must be specified to establish an E/Z constraint.
        let (sa, sb) = match (sa, sb) {
            (Some(a), Some(b)) => (a, b),
            _ => continue, // no constraint on this double bond
        };

        // Map template atoms to reactant atoms.
        let Some(&ra) = match_map.get(&(ta.0 as usize)) else {
            continue;
        };
        let Some(&rb) = match_map.get(&(tb.0 as usize)) else {
            continue;
        };

        // Outward directions from the corresponding reactant double bond ends.
        let ma = ez_stereo_outward(reactant, ra, rb);
        let mb = ez_stereo_outward(reactant, rb, ra);

        // If the reactant doesn't encode stereo on either end, skip (don't reject).
        let (ma, mb) = match (ma, mb) {
            (Some(a), Some(b)) => (a, b),
            _ => continue,
        };

        // Compare E/Z parity: same outward direction = Z, different = E.
        // If template and reactant disagree, reject this mapping.
        if (sa == sb) != (ma == mb) {
            return false;
        }
    }
    true
}

/// Convert a SMIRKS reactant-template `Molecule` to a `QueryMolecule` for VF2.
///
/// Constraints included:
/// - `AtomicNum` and `Aromatic` (always)
/// - `Charge` when non-zero
/// - `HCount` when a bracket atom specifies H > 0 (e.g. `[NH2:1]`)
///   Zero-H bracket atoms (`[N:1]`) are treated as "any H count" because
///   the parser returns 0 for both unspecified and explicit-zero H.
///
/// Chirality (`@`/`@@`) is NOT encoded into the query here; it is checked after
/// VF2 matching via `smirks_chirality_ok`, which uses a permutation-parity
/// comparison so that the same absolute configuration is recognised regardless
/// of how the reactant molecule was written as SMILES.
fn mol_to_query(mol: &Molecule) -> QueryMolecule {
    let mut qmol = QueryMolecule::new();

    for (_, atom) in mol.atoms() {
        let mut q = AtomQuery::And(
            Box::new(AtomQuery::Primitive(AtomPrimitive::AtomicNum(
                atom.element.atomic_number(),
            ))),
            Box::new(AtomQuery::Primitive(AtomPrimitive::Aromatic(atom.aromatic))),
        );

        if atom.charge != 0 {
            q = AtomQuery::And(
                Box::new(q),
                Box::new(AtomQuery::Primitive(AtomPrimitive::Charge(atom.charge))),
            );
        }

        if let Some(mass) = atom.isotope {
            q = AtomQuery::And(
                Box::new(q),
                Box::new(AtomQuery::Primitive(AtomPrimitive::Isotope(mass))),
            );
        }

        if let Some(h) = atom.hydrogen_count
            && h > 0
        {
            q = AtomQuery::And(
                Box::new(q),
                Box::new(AtomQuery::Primitive(AtomPrimitive::HCount(h))),
            );
        }

        qmol.add_atom_with_map(q, atom.atom_map);
    }

    for (_bidx, bond) in mol.bonds() {
        let bq = match bond.order {
            BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::Dative => {
                BondQuery::Primitive(BondPrimitive::Single)
            }
            BondOrder::Double => BondQuery::Primitive(BondPrimitive::Double),
            BondOrder::Triple => BondQuery::Primitive(BondPrimitive::Triple),
            BondOrder::Aromatic => BondQuery::Primitive(BondPrimitive::Aromatic),
            BondOrder::QuerySingleOrDouble => BondQuery::Or(
                Box::new(BondQuery::Primitive(BondPrimitive::Single)),
                Box::new(BondQuery::Primitive(BondPrimitive::Double)),
            ),
            BondOrder::QuerySingleOrAromatic => BondQuery::Or(
                Box::new(BondQuery::Primitive(BondPrimitive::Single)),
                Box::new(BondQuery::Primitive(BondPrimitive::Aromatic)),
            ),
            BondOrder::QueryDoubleOrAromatic => BondQuery::Or(
                Box::new(BondQuery::Primitive(BondPrimitive::Double)),
                Box::new(BondQuery::Primitive(BondPrimitive::Aromatic)),
            ),
            BondOrder::Quadruple => BondQuery::Primitive(BondPrimitive::Quadruple),
            BondOrder::Zero | BondOrder::QueryAny => BondQuery::Primitive(BondPrimitive::Any),
        };
        qmol.add_bond(bond.atom1.0 as usize, bond.atom2.0 as usize, bq);
    }

    qmol
}

/// Clear Up/Down stereo markers from bonds that have no adjacent double bond.
///
/// After a SMIRKS reaction, C=C → C=O conversions can leave stale Up/Down
/// markers (E/Z direction indicators) on single bonds that are no longer
/// adjacent to any double bond.  Such orphaned markers produce invalid SMILES
/// (`/C=O` is nonsensical) and must be demoted to plain Single bonds (RDKit #9339).
fn clear_orphaned_stereo_bonds(mol: Molecule) -> Molecule {
    let orphaned: Vec<BondIdx> = mol
        .bonds()
        .filter_map(|(bidx, bond)| {
            if bond.order != BondOrder::Up && bond.order != BondOrder::Down {
                return None;
            }
            // Up/Down is valid only when at least one endpoint has an adjacent
            // double bond (the one that the Up/Down bond specifies direction for).
            let has_double = [bond.atom1, bond.atom2].iter().any(|&a| {
                mol.neighbors(a)
                    .any(|(_, nb_bidx)| mol.bond(nb_bidx).order == BondOrder::Double)
            });
            if has_double { None } else { Some(bidx) }
        })
        .collect();

    if orphaned.is_empty() {
        return mol;
    }

    let mut builder = chematic_core::MoleculeBuilder::new();
    for (_, atom) in mol.atoms() {
        builder.add_atom(atom.clone());
    }
    for (bidx, bond) in mol.bonds() {
        let order = if orphaned.contains(&bidx) {
            BondOrder::Single
        } else {
            bond.order
        };
        let _ = builder.add_bond(bond.atom1, bond.atom2, order);
    }
    // copy_stereo_from copies stereo_neighbor_order but NOT stereo_groups.
    // Preserve both by applying each separately.
    builder.copy_stereo_from(&mol);
    let mut result = builder.build();
    // Restore enhanced stereo groups (ABS/OR/AND) that copy_stereo_from omits.
    result.set_stereo_groups(mol.stereo_groups().to_vec());
    result
}

// ---------------------------------------------------------------------------
// Product-side chirality correction (parity-aware)
// ---------------------------------------------------------------------------
//
// `build_product`'s Step 1 leaves every mapped core atom's chirality as
// whatever `src_atom.clone()` produced -- the reactant's own flag, if any.
// Two things can make that flag wrong or meaningless by the time the
// product molecule's bonds are fully built (Steps 3-4):
//
//   - An EXPLICIT product-template flag (`[C@:1]`/`[C@@:1]`) describes an
//     absolute configuration relative to the TEMPLATE's own neighbour
//     write-order, not the reactant's. Naively copying the symbol without
//     accounting for a reordered template ignores exactly the kind of
//     reorder that inverts configuration (mirroring the parity math
//     `smirks_chirality_ok`/`permutation_parity` already does for
//     REACTANT-side match validation, just never applied to the product).
//   - An INHERITED flag (no explicit template chirality) can survive on an
//     atom whose real bonding topology changed within the reaction's own
//     mapped/core region (e.g. a ring bond broken by the reaction) --
//     `build_product`'s old invalidation heuristic only ever compared
//     *unmapped* substituent element sets, never mapped-neighbour identity.
//
// Both are handled here, in a single post-build pass (run after all of
// `build_product`'s bonds exist, since validating either case requires the
// atom's REAL final neighbour set): an explicit template flag is
// transcribed together with its own neighbour order (validated against the
// atom's real adjacency, not re-derived); an inherited flag is kept only if
// both an order exists to trust it by and the core neighbourhood provably
// didn't change, and cleared otherwise. Two mechanisms, not one unified
// permutation-parity re-expression, because an unmapped/template-literal
// substituent has no reactant identity to map an order back to at all.

/// Map a template-space stereo order into product-index space via
/// `template_idx_to_new`. [`STEREO_H_SENTINEL`] passes through unchanged.
/// `None` if any real (non-sentinel) template atom didn't make it into the
/// product (shouldn't happen for a well-formed template, but fails closed).
fn map_stereo_order(order: &[u32], template_idx_to_new: &[Option<AtomIdx>]) -> Option<Vec<u32>> {
    order
        .iter()
        .map(|&t| {
            if t == STEREO_H_SENTINEL {
                Some(STEREO_H_SENTINEL)
            } else {
                template_idx_to_new
                    .get(t as usize)
                    .copied()
                    .flatten()
                    .map(|a| a.0)
            }
        })
        .collect()
}

/// True when `order` (already in product-index space) matches `new_idx`'s
/// REAL final neighbour set in `product`: at most one H-sentinel, no
/// duplicate real entries, right length, and the real entries are exactly
/// `new_idx`'s actual adjacency. Deliberately a set/degree check, not a
/// permutation-parity one -- `order` is stored verbatim once validated;
/// `corrected_chirality` (the SMILES writer) and `chematic-cip` already do
/// the write-order-independent parity math whenever they consume it.
fn order_matches_final_topology(product: &Molecule, new_idx: AtomIdx, order: &[u32]) -> bool {
    let sentinels = order.iter().filter(|&&x| x == STEREO_H_SENTINEL).count();
    if sentinels > 1 {
        return false;
    }
    let real: Vec<u32> = order
        .iter()
        .copied()
        .filter(|&x| x != STEREO_H_SENTINEL)
        .collect();
    let real_set: FxHashSet<u32> = real.iter().copied().collect();
    if real.len() != real_set.len() {
        return false; // Duplicate entry -- malformed order.
    }
    if order.len() != product.degree(new_idx) + sentinels {
        return false;
    }
    let actual: FxHashSet<u32> = product.neighbors(new_idx).map(|(nb, _)| nb.0).collect();
    real_set == actual
}

/// Bug-B validity gate for an atom that inherited its chirality flag from
/// the reactant (no explicit product-template chirality): the original
/// unmapped-substituent element-multiset comparison, PLUS a symmetric
/// mapped-neighbour atom-map-number-set comparison -- closing the gap where
/// a core/mapped neighbour's topology changes (e.g. a ring bond broken by
/// the reaction) invisibly to the element-multiset-only check.
fn bug_b_topology_unchanged(
    product_template: &Molecule,
    ti: AtomIdx,
    reactant: &Molecule,
    src_idx: AtomIdx,
    all_template_atoms: &FxHashSet<(usize, AtomIdx)>,
    mol_idx: usize,
    atom_to_map: &FxHashMap<(usize, AtomIdx), u16>,
) -> bool {
    let mut prod_elems: FxHashMap<u8, usize> = FxHashMap::default();
    let mut prod_mapped: FxHashSet<u16> = FxHashSet::default();
    for (nb, _) in product_template.neighbors(ti) {
        match product_template.atom(nb).atom_map {
            None => {
                *prod_elems
                    .entry(product_template.atom(nb).element.atomic_number())
                    .or_insert(0) += 1;
            }
            Some(am) => {
                prod_mapped.insert(am);
            }
        }
    }

    let mut rxn_elems: FxHashMap<u8, usize> = FxHashMap::default();
    let mut rxn_mapped: FxHashSet<u16> = FxHashSet::default();
    for (nb, _) in reactant.neighbors(src_idx) {
        if !all_template_atoms.contains(&(mol_idx, nb)) {
            continue;
        }
        match atom_to_map.get(&(mol_idx, nb)) {
            Some(&am) => {
                rxn_mapped.insert(am);
            }
            None => {
                *rxn_elems
                    .entry(reactant.atom(nb).element.atomic_number())
                    .or_insert(0) += 1;
            }
        }
    }

    // NOTE: `prod_elems.is_empty()` is deliberately NOT special-cased here.
    // An empty `prod_elems` legitimately means "product-template atom has no
    // unmapped neighbours" -- which is also true when the reactant's
    // template-matched unmapped substituents (F/Cl/Br etc., matched by the
    // reactant SMARTS itself, not carried-through remote substituents --
    // those never enter `rxn_elems` at all, since `rxn_elems` only counts
    // neighbours inside `all_template_atoms`) were silently DELETED by the
    // product template. `prod_elems == rxn_elems` alone correctly requires
    // both to be empty together, or both to hold the same multiset --
    // exactly what "unchanged" means; special-casing empty-on-one-side
    // would let a genuine deletion through undetected.
    prod_elems == rxn_elems && prod_mapped == rxn_mapped
}

/// Map a REACTANT-space stereo order (atom indices into
/// `input_mols[mol_idx]`) into product-index space via `src_to_new`.
/// [`STEREO_H_SENTINEL`] passes through unchanged. `None` if any real
/// (non-sentinel) reactant neighbour never made it into the product (e.g. it
/// was part of a leaving group the reaction consumed) -- the correspondence
/// is not derivable, so the caller must fail closed rather than guess.
fn remap_reactant_stereo_order(
    order: &[u32],
    mol_idx: usize,
    src_to_new: &FxHashMap<(usize, AtomIdx), AtomIdx>,
) -> Option<Vec<u32>> {
    order
        .iter()
        .map(|&t| {
            if t == STEREO_H_SENTINEL {
                Some(STEREO_H_SENTINEL)
            } else {
                src_to_new.get(&(mol_idx, AtomIdx(t))).map(|a| a.0)
            }
        })
        .collect()
}

/// The reactant stereo order with its implicit-H slot given to the one
/// neighbour the product atom gained, when that is the only change: every
/// other recorded neighbour maps into the product and stays bonded, and the
/// product atom has no H left. `None` otherwise.
fn gained_neighbour_order(
    order: &[u32],
    mol_idx: usize,
    src_to_new: &FxHashMap<(usize, AtomIdx), AtomIdx>,
    product: &Molecule,
    new_idx: AtomIdx,
) -> Option<Vec<u32>> {
    if order.iter().filter(|&&t| t == STEREO_H_SENTINEL).count() != 1 {
        return None;
    }
    let remapped = remap_reactant_stereo_order(order, mol_idx, src_to_new)?;
    let kept: FxHashSet<u32> = remapped
        .iter()
        .copied()
        .filter(|&t| t != STEREO_H_SENTINEL)
        .collect();
    let neighbours: Vec<AtomIdx> = product.neighbors(new_idx).map(|(nb, _)| nb).collect();
    if neighbours.len() != kept.len() + 1
        || !kept.iter().all(|&t| neighbours.contains(&AtomIdx(t)))
        || chematic_core::implicit_hcount(product, new_idx) != 0
    {
        return None;
    }
    let gained = neighbours.into_iter().find(|nb| !kept.contains(&nb.0))?;
    if product.atom(gained).element == chematic_core::Element::H {
        return None;
    }
    Some(
        remapped
            .into_iter()
            .map(|t| if t == STEREO_H_SENTINEL { gained.0 } else { t })
            .collect(),
    )
}

/// Post-build chirality correction. Must run after `build_product`'s Steps
/// 3-4 (all bonds added), since validation needs each atom's real final
/// degree/neighbour set. See the module-level doc above for why this is two
/// mechanisms (transcribe-and-validate for an explicit template flag,
/// gate-and-preserve-or-clear for an inherited one), not one.
fn correct_product_stereo(
    mut product: Molecule,
    product_template: &Molecule,
    template_idx_to_new: &[Option<AtomIdx>],
    global_map: &FxHashMap<u16, (usize, AtomIdx)>,
    all_template_atoms: &FxHashSet<(usize, AtomIdx)>,
    input_mols: &[&Molecule],
    src_to_new: &FxHashMap<(usize, AtomIdx), AtomIdx>,
) -> Molecule {
    let atom_to_map: FxHashMap<(usize, AtomIdx), u16> =
        global_map.iter().map(|(&am, &k)| (k, am)).collect();

    for i in 0..product_template.atom_count() {
        let ti = AtomIdx(i as u32);
        let Some(new_idx) = template_idx_to_new[i] else {
            continue;
        };
        let tmpl_atom = product_template.atom(ti);

        if tmpl_atom.chirality != Chirality::None {
            // Explicit product-template chirality: the template is the sole
            // source of truth here regardless of whether this atom is a
            // matched core atom, an unmatched map number, or fully new --
            // no global_map lookup needed on this branch.
            let resolved = product_template
                .stereo_neighbor_order(ti)
                .and_then(|order| map_stereo_order(order, template_idx_to_new))
                .filter(|order| order_matches_final_topology(&product, new_idx, order));
            match resolved {
                Some(order) => {
                    product.set_chirality(new_idx, tmpl_atom.chirality);
                    product.set_stereo_neighbor_order(new_idx, order);
                }
                None => product.set_chirality(new_idx, Chirality::None),
            }
        } else if let Some(am) = tmpl_atom.atom_map
            && let Some(&(mol_idx, src_idx)) = global_map.get(&am)
        {
            // Matched core atom, no explicit template chirality: keep the
            // flag `src_atom.clone()` gave it in Step 1 only if it's
            // actually usable and provably still valid. Three-way rule,
            // all conditions required:
            //   1. bug_b_topology_unchanged: the TEMPLATE-level neighbour
            //      composition (unmapped element multiset + mapped atom-map
            //      set) is unchanged reactant-template -> product-template.
            //   2. remap_reactant_stereo_order succeeds: every atom the
            //      reactant's own recorded order references is uniquely
            //      identifiable in the product (via src_to_new) -- fails
            //      closed (None) for a leaving-group neighbour the reaction
            //      consumed, where no correspondence is derivable at all.
            //   3. order_matches_final_topology: the remapped order's real
            //      entries are EXACTLY this atom's real final adjacency in
            //      the built product molecule -- catches a neighbour that
            //      survived into the product (so (2) succeeds) but whose
            //      specific BOND to this atom did not (e.g. reattached
            //      elsewhere by the reaction), which a template-level
            //      heuristic alone cannot see.
            // Never keep a raw flag without a validated order alongside it.
            let src_atom = input_mols[mol_idx].atom(src_idx);
            if src_atom.chirality != Chirality::None {
                let topology_unchanged = bug_b_topology_unchanged(
                    product_template,
                    ti,
                    input_mols[mol_idx],
                    src_idx,
                    all_template_atoms,
                    mol_idx,
                    &atom_to_map,
                );
                let remapped_order = topology_unchanged
                    .then(|| input_mols[mol_idx].stereo_neighbor_order(src_idx))
                    .flatten()
                    .and_then(|order| remap_reactant_stereo_order(order, mol_idx, src_to_new))
                    .filter(|order| order_matches_final_topology(&product, new_idx, order));
                // Case B: the centre lost its one implicit H and gained one
                // new neighbour, every other neighbour kept. The new bond
                // takes the H's place, as RDKit's reaction runner copies the
                // tag (C-H functionalization keeps the configuration).
                let remapped_order = remapped_order.or_else(|| {
                    input_mols[mol_idx]
                        .stereo_neighbor_order(src_idx)
                        .and_then(|order| {
                            gained_neighbour_order(order, mol_idx, src_to_new, &product, new_idx)
                        })
                });
                match remapped_order {
                    Some(order) => {
                        product.set_chirality(new_idx, src_atom.chirality);
                        product.set_stereo_neighbor_order(new_idx, order);
                    }
                    None => product.set_chirality(new_idx, Chirality::None),
                }
            }
        }
        // else: no template chirality, and not a matched core atom --
        // new_atom.chirality is already Chirality::None from Step 1's clone
        // of a template-literal atom.
    }

    product
}

/// RDKit's bond order around a carried, non-ring atom of the product: the
/// bond to `discoverer` (the neighbour the walk reached it from) first, then
/// its other bonds in reactant order. `None` for an atom in a ring (other
/// neighbours may add their bonds earlier) or with an implicit-H slot.
fn rdkit_carried_bond_order(
    mol: &Molecule,
    atom: AtomIdx,
    discoverer: Option<AtomIdx>,
) -> Option<Vec<u32>> {
    let d = discoverer?;
    if atoms_in_rings(mol)[atom.0 as usize] {
        return None;
    }
    let reactant: Vec<u32> = match mol.stereo_neighbor_order(atom) {
        Some(o) => o.to_vec(),
        None => mol.neighbors(atom).map(|(n, _)| n.0).collect(),
    };
    if reactant.contains(&chematic_core::STEREO_H_SENTINEL) || !reactant.contains(&d.0) {
        return None;
    }
    let mut order = vec![d.0];
    order.extend(reactant.into_iter().filter(|&n| n != d.0));
    Some(order)
}

fn add_product_template_bonds(
    builder: &mut MoleculeBuilder,
    product_template: &Molecule,
    template_idx_to_new: &[Option<AtomIdx>],
    source_of: &[Option<(usize, AtomIdx)>],
    input_mols: &[&Molecule],
) -> Vec<AtomIdx> {
    // Atoms of a product `:` bond are flagged aromatic, as RDKit does.
    let mut colon_atoms = Vec::new();
    for (_bidx, bond) in product_template.bonds() {
        let a_new = template_idx_to_new[bond.atom1.0 as usize].unwrap();
        let b_new = template_idx_to_new[bond.atom2.0 as usize].unwrap();
        if bond.order == BondOrder::Aromatic {
            colon_atoms.extend([a_new, b_new]);
        }
        let mut order = bond.order;
        let (mut from, mut to) = (a_new, b_new);
        // A template single bond over a reactant `/`/`\` bond keeps the
        // direction: it carries the E/Z of an adjacent double bond the
        // template leaves alone.
        if order == BondOrder::Single
            && let (Some((ma, sa)), Some((mb, sb))) =
                (source_of[a_new.0 as usize], source_of[b_new.0 as usize])
            && ma == mb
            && let Some((_, ob)) = input_mols[ma].bond_between(sa, sb)
            && matches!(ob.order, BondOrder::Up | BondOrder::Down)
        {
            order = ob.order;
            if ob.atom1 != sa {
                (from, to) = (b_new, a_new);
            }
        }
        // A product `~` keeps the reactant bond between mapped atoms and is
        // a zero-order bond otherwise.
        if order == BondOrder::QueryAny {
            order = BondOrder::Zero;
            if let (Some((ma, sa)), Some((mb, sb))) =
                (source_of[a_new.0 as usize], source_of[b_new.0 as usize])
                && ma == mb
                && let Some((_, ob)) = input_mols[ma].bond_between(sa, sb)
            {
                order = ob.order;
                if ob.atom1 != sa {
                    (from, to) = (b_new, a_new);
                }
            }
        }
        let _ = builder.add_bond(from, to, order);
    }
    colon_atoms
}

#[allow(clippy::too_many_arguments)]
fn carry_source_bonds(
    builder: &mut MoleculeBuilder,
    source_of: &[Option<(usize, AtomIdx)>],
    new_of: &[Vec<Option<AtomIdx>>],
    input_mols: &[&Molecule],
    all_template_atoms: &FxHashSet<(usize, AtomIdx)>,
    global_map: &FxHashMap<u16, (usize, AtomIdx)>,
    reactant_bond_maps: &FxHashSet<(u16, u16)>,
) {
    let map_of = |key: (usize, AtomIdx)| {
        global_map
            .iter()
            .find(|&(_, &mapped)| mapped == key)
            .map(|(&atom_map, _)| atom_map)
    };
    let carried = source_of
        .iter()
        .enumerate()
        .filter_map(|(new, src)| src.map(|src| (src, AtomIdx(new as u32))));
    for ((mol_idx, src_idx), a_new) in carried {
        for (nb_idx, bond_idx) in input_mols[mol_idx].neighbors(src_idx) {
            let nb_key = (mol_idx, nb_idx);
            let Some(b_new) = new_of[mol_idx][nb_idx.0 as usize] else {
                continue;
            };
            if all_template_atoms.contains(&(mol_idx, src_idx))
                && all_template_atoms.contains(&nb_key)
            {
                // Both ends matched: templates decide this bond, except for
                // a reactant bond between mapped atoms the template does not
                // itself connect (RDKit, github #1387).
                let kept = match (map_of((mol_idx, src_idx)), map_of(nb_key)) {
                    (Some(a), Some(b)) => !reactant_bond_maps.contains(&(a.min(b), a.max(b))),
                    _ => false,
                };
                if !kept {
                    continue;
                }
            }
            let source_bond = input_mols[mol_idx].bond(bond_idx);
            let (a, b) = if source_bond.atom1 == src_idx {
                (a_new, b_new)
            } else {
                (b_new, a_new)
            };
            let _ = builder.add_bond(a, b, source_bond.order);
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn restore_carried_atom_stereo(
    product: &mut Molecule,
    product_template: &Molecule,
    template_idx_to_new: &[Option<AtomIdx>],
    input_mols: &[&Molecule],
    source_of: &[Option<(usize, AtomIdx)>],
    discoverer: &FxHashMap<(usize, AtomIdx), AtomIdx>,
    src_to_new: &FxHashMap<(usize, AtomIdx), AtomIdx>,
    keep_carried_maps: bool,
) -> Vec<((usize, AtomIdx), AtomIdx)> {
    let sourced_atoms: Vec<_> = source_of
        .iter()
        .enumerate()
        .filter_map(|(new, src)| src.map(|src| (src, AtomIdx(new as u32))))
        .collect();
    for &((mol_idx, src_idx), new_idx) in &sourced_atoms {
        let src = input_mols[mol_idx];
        if keep_carried_maps
            && let Some((class, _)) = src.atom(src_idx).chirality.nontetrahedral()
            && let Some(template_idx) = template_idx_to_new
                .iter()
                .position(|idx| *idx == Some(new_idx))
            && product_template
                .atom(AtomIdx(template_idx as u32))
                .chirality
                == Chirality::None
        {
            product.set_chirality(new_idx, Chirality::from_nontetrahedral(class, 0).unwrap());
            continue;
        }
        if keep_carried_maps
            && src.atom(src_idx).chirality.nontetrahedral().is_some()
            && !template_idx_to_new.contains(&Some(new_idx))
        {
            let order = rdkit_carried_bond_order(
                src,
                src_idx,
                discoverer.get(&(mol_idx, src_idx)).copied(),
            )
            .and_then(|order| remap_reactant_stereo_order(&order, mol_idx, src_to_new))
            .filter(|order| order_matches_final_topology(product, new_idx, order));
            match order {
                Some(order) => {
                    product.set_chirality(new_idx, src.atom(src_idx).chirality);
                    product.set_stereo_neighbor_order(new_idx, order);
                }
                None => product.set_chirality(new_idx, Chirality::None),
            }
            continue;
        }
        if !src.atom(src_idx).chirality.is_tetrahedral()
            || template_idx_to_new.contains(&Some(new_idx))
        {
            continue;
        }
        match src
            .stereo_neighbor_order(src_idx)
            .and_then(|order| remap_reactant_stereo_order(order, mol_idx, src_to_new))
            .filter(|order| order_matches_final_topology(product, new_idx, order))
        {
            Some(order) => product.set_stereo_neighbor_order(new_idx, order),
            None => product.set_chirality(new_idx, Chirality::None),
        }
    }
    sourced_atoms
}

fn transfer_folded_h_stereo(
    product: &mut Molecule,
    folded_h: &FxHashSet<(usize, AtomIdx)>,
    input_mols: &[&Molecule],
    new_of: &[Vec<Option<AtomIdx>>],
) {
    let from_core = |order: BondOrder, atom1: AtomIdx, at: AtomIdx| match (order, atom1 == at) {
        (o, true) => o,
        (BondOrder::Up, false) => BondOrder::Down,
        (BondOrder::Down, false) => BondOrder::Up,
        (o, false) => o,
    };
    for &(mol_idx, h) in folded_h {
        let src = input_mols[mol_idx];
        let Some((core, h_bond_idx)) = src.neighbors(h).next() else {
            continue;
        };
        let h_bond = src.bond(h_bond_idx);
        let h_dir = from_core(h_bond.order, h_bond.atom1, core);
        if !matches!(h_dir, BondOrder::Up | BondOrder::Down) {
            continue;
        }
        let Some(core_new) = new_of[mol_idx][core.0 as usize] else {
            continue;
        };
        let on_double = product
            .neighbors(core_new)
            .any(|(_, bond)| product.bond(bond).order == BondOrder::Double);
        let marked = product
            .neighbors(core_new)
            .any(|(_, bond)| matches!(product.bond(bond).order, BondOrder::Up | BondOrder::Down));
        if !on_double || marked {
            continue;
        }
        let other = if h_dir == BondOrder::Up {
            BondOrder::Down
        } else {
            BondOrder::Up
        };
        let carrier = product
            .neighbors(core_new)
            .find(|&(_, bond)| product.bond(bond).order == BondOrder::Single)
            .map(|(_, bond)| bond);
        if let Some(bond) = carrier {
            let atom1 = product.bond(bond).atom1;
            product.set_bond_order(bond, from_core(other, atom1, core_new));
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn carry_substituent_atoms(
    builder: &mut MoleculeBuilder,
    input_mols: &[&Molecule],
    all_template_atoms: &FxHashSet<(usize, AtomIdx)>,
    folded_h: &FxHashSet<(usize, AtomIdx)>,
    core_keys: Vec<(usize, AtomIdx)>,
    carry_substituents: bool,
    keep_carried_maps: bool,
    src_to_new: &mut FxHashMap<(usize, AtomIdx), AtomIdx>,
) -> FxHashMap<(usize, AtomIdx), AtomIdx> {
    let mut visited: Vec<Vec<bool>> = input_mols
        .iter()
        .map(|molecule| vec![false; molecule.atom_count()])
        .collect();
    for &(mol_idx, atom) in all_template_atoms.iter().chain(folded_h.iter()) {
        visited[mol_idx][atom.0 as usize] = true;
    }

    // Each carried atom's neighbour from which the walk reached it. RDKit's
    // product bond order around that atom starts with this bond.
    let mut discoverer = FxHashMap::default();
    if !carry_substituents {
        return discoverer;
    }

    // The seeds are in reactant atom order so product atom order never
    // depends on hash-table iteration order.
    let mut queue: VecDeque<(usize, AtomIdx)> = core_keys.into();
    while let Some((mol_idx, current)) = queue.pop_front() {
        for (neighbor, _bond) in input_mols[mol_idx].neighbors(current) {
            let key = (mol_idx, neighbor);
            if std::mem::replace(&mut visited[mol_idx][neighbor.0 as usize], true) {
                continue;
            }
            discoverer.insert(key, current);
            let mut atom = input_mols[mol_idx].atom(neighbor).clone();
            if !keep_carried_maps {
                atom.atom_map = None;
            }
            let product_atom = builder.add_atom(atom);
            src_to_new.insert(key, product_atom);
            queue.push_back(key);
        }
    }
    discoverer
}

struct DenseAtomMaps {
    source_of: Vec<Option<(usize, AtomIdx)>>,
    new_of: Vec<Vec<Option<AtomIdx>>>,
}

fn dense_atom_maps(
    builder_atom_count: usize,
    input_mols: &[&Molecule],
    src_to_new: &FxHashMap<(usize, AtomIdx), AtomIdx>,
) -> DenseAtomMaps {
    let mut source_of = vec![None; builder_atom_count];
    let mut new_of: Vec<Vec<Option<AtomIdx>>> = input_mols
        .iter()
        .map(|molecule| vec![None; molecule.atom_count()])
        .collect();
    for (&source, &product) in src_to_new {
        source_of[product.0 as usize] = Some(source);
        new_of[source.0][source.1.0 as usize] = Some(product);
    }
    DenseAtomMaps { source_of, new_of }
}

fn restore_product_aromaticity(
    molecule: &mut Molecule,
    colon_atoms: &[AtomIdx],
    dearomatized: &[AtomIdx],
) {
    for &atom in colon_atoms {
        if !molecule.atom(atom).aromatic {
            molecule.set_atom_aromatic(atom, true);
        }
    }
    for &atom in dearomatized {
        let aromatic_bonds = molecule
            .neighbors(atom)
            .filter(|&(_, bond)| molecule.bond(bond).order == BondOrder::Aromatic)
            .count();
        if aromatic_bonds >= 2 {
            molecule.set_atom_aromatic(atom, true);
        }
    }
}

fn reactant_has_explicit_h(
    input_mols: &[&Molecule],
    cache: &[std::cell::OnceCell<bool>],
    mol_idx: usize,
) -> bool {
    *cache[mol_idx].get_or_init(|| {
        let molecule = input_mols[mol_idx];
        molecule.atoms().any(|(atom_idx, atom)| {
            atom.element == chematic_core::Element::H && molecule.degree(atom_idx) == 1
        })
    })
}

#[allow(clippy::too_many_arguments)]
fn refill_product_hydrogens(
    molecule: &mut Molecule,
    sourced_atoms: &[((usize, AtomIdx), AtomIdx)],
    input_mols: &[&Molecule],
    template_idx_to_new: &[Option<AtomIdx>],
    explicit_h_cache: &[std::cell::OnceCell<bool>],
    mut folded_cores: Vec<AtomIdx>,
    source_of: &[Option<(usize, AtomIdx)>],
    folded_h: &FxHashSet<(usize, AtomIdx)>,
) {
    // A carried atom of an explicit-H reactant that lost a bond gets its
    // hydrogens back, matching RDKit's refill after an edit.
    for &((mol_idx, source_idx), product_idx) in sourced_atoms {
        let source = input_mols[mol_idx];
        let atom = source.atom(source_idx);
        if atom.hydrogen_count != Some(0)
            || atom.charge != 0
            || atom.aromatic
            || !atom.element.is_organic_subset()
            || atom.element == chematic_core::Element::H
            || molecule.degree(product_idx) >= source.degree(source_idx)
            || template_idx_to_new.contains(&Some(product_idx))
        {
            continue;
        }
        if reactant_has_explicit_h(input_mols, explicit_h_cache, mol_idx) {
            molecule.set_hydrogen_count(product_idx, None);
            folded_cores.push(product_idx);
        }
    }

    for core in folded_cores {
        let atom = molecule.atom(core);
        let mixed_aromatic = !atom.aromatic
            && molecule
                .neighbors(core)
                .any(|(_, bond)| molecule.bond(bond).order == BondOrder::Aromatic);
        let element_changed = source_of[core.0 as usize].is_some_and(|(mol_idx, source_idx)| {
            input_mols[mol_idx].atom(source_idx).element != atom.element
        });
        if mixed_aromatic && element_changed {
            molecule.set_hydrogen_count(core, None);
            continue;
        }
        let count = if mixed_aromatic {
            source_of[core.0 as usize].map_or(0, |(mol_idx, source_idx)| {
                input_mols[mol_idx]
                    .neighbors(source_idx)
                    .filter(|&(neighbor, _)| folded_h.contains(&(mol_idx, neighbor)))
                    .count() as u8
            })
        } else if atom.charge != 0 && !atom.aromatic {
            let bond_order_sum: i16 = molecule
                .neighbors(core)
                .map(|(_, bond)| i16::from(molecule.bond(bond).order.order_int()))
                .sum();
            crate::rdkit_valence::implicit_hydrogens(
                atom.element.atomic_number(),
                atom.charge,
                bond_order_sum,
            )
        } else {
            chematic_core::implicit_hcount(molecule, core)
        };
        molecule.set_hydrogen_count(core, Some(0));
        for _ in 0..count {
            let mut hydrogen = chematic_core::Atom::new(chematic_core::Element::H);
            hydrogen.hydrogen_count = Some(0);
            let hydrogen_idx = molecule.add_atom(hydrogen);
            let _ = molecule.add_bond(core, hydrogen_idx, BondOrder::Single);
        }
    }
}

fn apply_rdkit_template_hydrogens(
    molecule: &mut Molecule,
    template_idx_to_new: &[Option<AtomIdx>],
) {
    for &idx in template_idx_to_new.iter().flatten() {
        let atom = molecule.atom(idx);
        let aromatic_bonds = molecule
            .neighbors(idx)
            .filter(|&(_, bond)| molecule.bond(bond).order == BondOrder::Aromatic)
            .count();
        if atom.aromatic
            && atom.charge != 0
            && atom.hydrogen_count.is_none()
            && !atom.wildcard
            && aromatic_bonds == 2
        {
            let other_bonds: i16 = molecule
                .neighbors(idx)
                .filter(|&(_, bond)| molecule.bond(bond).order != BondOrder::Aromatic)
                .map(|(_, bond)| i16::from(molecule.bond(bond).order.order_int()))
                .sum();
            let rdkit_h = crate::rdkit_valence::implicit_hydrogens(
                atom.element.atomic_number(),
                atom.charge,
                3 + other_bonds,
            );
            if rdkit_h != chematic_core::implicit_hcount(molecule, idx) {
                molecule.set_hydrogen_count(idx, Some(rdkit_h));
            }
            continue;
        }
        if (atom.charge == 0 && atom.element.is_organic_subset())
            || atom.hydrogen_count.is_some()
            || atom.wildcard
            || atom.aromatic
            || molecule
                .neighbors(idx)
                .any(|(_, bond)| molecule.bond(bond).order == BondOrder::Aromatic)
        {
            continue;
        }
        let bond_order_sum: i16 = molecule
            .neighbors(idx)
            .map(|(_, bond)| {
                let bond = molecule.bond(bond);
                match bond.order {
                    BondOrder::Dative if bond.atom1 == idx => 0,
                    order => i16::from(order.order_int()),
                }
            })
            .sum();
        let rdkit_h = crate::rdkit_valence::implicit_hydrogens(
            atom.element.atomic_number(),
            atom.charge,
            bond_order_sum,
        );
        if rdkit_h != chematic_core::implicit_hcount(molecule, idx) {
            molecule.set_hydrogen_count(idx, Some(rdkit_h));
        }
    }
}

fn traced_product(
    molecule: Molecule,
    product_template: &Molecule,
    template_idx_to_new: Vec<Option<AtomIdx>>,
    sourced_atoms: &[((usize, AtomIdx), AtomIdx)],
) -> (TracedProduct, Vec<Option<AtomIdx>>) {
    let mut atom_sources = vec![None; molecule.atom_count()];
    for &((reactant, atom), product_atom) in sourced_atoms {
        atom_sources[product_atom.0 as usize] = Some(ReactantAtom { reactant, atom });
    }
    let mut template_maps = vec![None; molecule.atom_count()];
    for (template_atom, product_atom) in template_idx_to_new.iter().enumerate() {
        if let Some(product_atom) = product_atom {
            template_maps[product_atom.0 as usize] = product_template
                .atom(AtomIdx(template_atom as u32))
                .atom_map;
        }
    }
    (
        TracedProduct {
            molecule,
            atom_sources,
            template_maps,
        },
        template_idx_to_new,
    )
}

/// Build one product molecule applying full SMIRKS semantics.
///
/// Atom-mapped and new template atoms are assembled first. Substituents and
/// source bonds are then carried into the product before stereo, explicit-H,
/// aromaticity, and trace metadata are finalized.
#[allow(clippy::too_many_arguments)]
fn build_product(
    product_template: &Molecule,
    global_map: &FxHashMap<u16, (usize, AtomIdx)>,
    reactant_template_atoms: &FxHashMap<u16, ReactantMapAtom>,
    reactant_bond_maps: &FxHashSet<(u16, u16)>,
    input_mols: &[&Molecule],
    all_template_atoms: &FxHashSet<(usize, AtomIdx)>,
    carry_substituents: bool,
    specs: &[ProductAtomSpec],
    keep_carried_maps: bool,
) -> (TracedProduct, Vec<Option<AtomIdx>>) {
    // A product holds at most the reactants' atoms and bonds plus the
    // template's: reserve them once instead of growing per atom and bond.
    let mut builder = MoleculeBuilder::with_capacity(
        input_mols.iter().map(|m| m.atom_count()).sum::<usize>() + product_template.atom_count(),
        input_mols.iter().map(|m| m.bond_count()).sum::<usize>() + product_template.bond_count(),
    );
    // What each template atom's spelling specifies; without spelling
    // information a charge always applies and only `H<n>` with n > 0 pins
    // the H count (a bare bracket atom's `Some(0)` is "unspecified", #18).
    let specs = (specs.len() == product_template.atom_count()).then_some(specs);
    let spec_of = |i: usize| -> ProductAtomSpec {
        let atom = product_template.atom(AtomIdx(i as u32));
        specs.map_or(
            ProductAtomSpec {
                charge: true,
                hcount: atom.hydrogen_count.is_some_and(|h| h > 0),
            },
            |s| s[i],
        )
    };
    let template_h = |i: usize| -> Option<u8> {
        let atom = product_template.atom(AtomIdx(i as u32));
        spec_of(i).hcount.then(|| atom.hydrogen_count.unwrap_or(0))
    };

    // template_idx_to_new[i]: new AtomIdx for product template atom i.
    let mut template_idx_to_new: Vec<Option<AtomIdx>> = vec![None; product_template.atom_count()];
    // src_to_new: (mol_idx, src_AtomIdx) → new AtomIdx in the product.
    let mut src_to_new: FxHashMap<(usize, AtomIdx), AtomIdx> = FxHashMap::with_capacity_and_hasher(
        input_mols.iter().map(|m| m.atom_count()).sum(),
        Default::default(),
    );

    // --- Step 1: add product template atoms ---
    // core_keys: only source atoms that are mapped by THIS product template.
    // Using global_map.values() (all matched atoms across all templates) would
    // seed the BFS in Step 2 from atoms belonging to *other* product templates,
    // causing their substituents to leak into this product (issue #13).
    let mut core_keys: Vec<(usize, AtomIdx)> = product_template
        .atoms()
        .filter_map(|(_, atom)| atom.atom_map.and_then(|am| global_map.get(&am).copied()))
        .collect();
    core_keys.sort_unstable();
    core_keys.dedup();

    // Explicit hydrogen atoms (issue #734): a mapped, non-stereo atom's plain
    // H-atom neighbours that the template did not match are not carried as
    // substituents. The atom's hydrogens are re-derived from its product
    // valence and added back as explicit H atoms after assembly, so an
    // `add_hydrogens` reactant gives the same products as its implicit-H
    // form instead of keeping every H and failing the valence check
    // (`[C:1]-[C:2]>>[C:1]=[C:2]`, `[O:1]>>[O:1]C`, `[CH3:1]>>[CH3:1]`).
    let mut folded_h: FxHashSet<(usize, AtomIdx)> = FxHashSet::default();
    // Mapped aromatic atoms the template spells aliphatic (see below).
    let mut dearomatized: Vec<AtomIdx> = Vec::new();
    let mut folded_cores: Vec<AtomIdx> = Vec::new();
    let plain_h = |mol: &Molecule, idx: AtomIdx| {
        let a = mol.atom(idx);
        a.element.atomic_number() == 1
            && !a.wildcard
            && a.isotope.is_none()
            && a.charge == 0
            && a.atom_map.is_none()
            && mol.degree(idx) == 1
    };

    // Reactants given with explicit H atoms (`add_hydrogens`).
    let explicit_h_cache: Vec<std::cell::OnceCell<bool>> = input_mols
        .iter()
        .map(|_| std::cell::OnceCell::new())
        .collect();
    for (i, slot) in template_idx_to_new.iter_mut().enumerate() {
        let tmpl_atom = product_template.atom(AtomIdx(i as u32));
        let new_idx = if let Some(am) = tmpl_atom.atom_map {
            if let Some(&(mol_idx, src_idx)) = global_map.get(&am) {
                // Core atom: copy source, then override electronic state from template.
                let src_atom = input_mols[mol_idx].atom(src_idx);
                let mut new_atom = src_atom.clone();
                // Element change (issue #734): as in RDKit, a mapped atom takes
                // the product template's element when that element differs
                // from the reactant template's (`[C:1]>>[N:1]`). A wildcard
                // product atom (`[*:1]`) or an unchanged element keeps the
                // matched atom's element; a reactant query naming no single
                // element (`[C,N:1]`, `[*:1]`) takes any concrete product one.
                let reactant_tmpl = reactant_template_atoms.get(&am);
                if !tmpl_atom.wildcard {
                    let changes =
                        reactant_tmpl.is_none_or(|r| r.wildcard || r.element != tmpl_atom.element);
                    if changes && new_atom.element != tmpl_atom.element {
                        new_atom.element = tmpl_atom.element;
                        new_atom.isotope = tmpl_atom.isotope;
                    }
                }
                // A wildcard product atom (`[*:1]`) keeps the matched atom's
                // aromaticity along with its element.
                if !tmpl_atom.wildcard {
                    new_atom.aromatic = tmpl_atom.aromatic;
                }
                // Unspecified charge keeps the reactant atom's (RDKit:
                // `[O-:1]>>[O:1]` leaves the alkoxide charged), unless the
                // element changes: RDKit then takes the template atom's
                // charge too (`[N:1]>>[S:1]` on a quaternary ammonium gives
                // a neutral sulfur, not `[SH+]`; found by the SMIRKS fuzz
                // lane, #754).
                if spec_of(i).charge || new_atom.element != src_atom.element {
                    new_atom.charge = tmpl_atom.charge;
                }
                // An isotope in the template applies; none keeps the
                // reactant's (`[C:1]>>[13C:1]` labels the atom).
                if tmpl_atom.isotope.is_some() {
                    new_atom.isotope = tmpl_atom.isotope;
                }
                // H count: a template `H<n>` with n > 0 pins it (`[NH2:1]`); a
                // bare bracket atom's `Some(0)` means "unspecified" (issue #18).
                // Unspecified, RDKit keeps the reactant atom's explicit H count
                // when the atom's template degree is unchanged (`[n:1]>>[n:1]`
                // keeps pyrrole's `[nH]`) and re-derives it from valence when
                // bonds were added or removed (`[n:1]>>[n:1]C`).
                let degree_unchanged = reactant_tmpl
                    .is_some_and(|r| r.degree == product_template.degree(AtomIdx(i as u32)));
                // An element change re-derives H as well (`[C:1]>>[N:1]` on
                // `[13CH3]C` gives `CN`, not a five-valent N).
                let element_unchanged = new_atom.element == src_atom.element;
                new_atom.hydrogen_count = match template_h(i) {
                    Some(h) => Some(h),
                    None if degree_unchanged && element_unchanged => src_atom.hydrogen_count,
                    None => None,
                };
                // Chirality is intentionally left as whatever src_atom.clone()
                // produced above (i.e. the reactant's own flag, if any) --
                // correct_product_stereo, run after all bonds are added
                // below, is the sole authority on the final chirality/order
                // for every atom, whether it inherits from the reactant or
                // carries an explicit product-template @/@@.
                new_atom.atom_map = None;
                // Only organic-subset atoms: their H count can be re-derived
                // from valence.
                let explicit_h: Vec<AtomIdx> = if carry_substituents
                    && src_atom.chirality == Chirality::None
                    && src_atom.element.is_organic_subset()
                    && new_atom.element.is_organic_subset()
                {
                    input_mols[mol_idx]
                        .neighbors(src_idx)
                        .map(|(nb, _)| nb)
                        .filter(|&nb| {
                            plain_h(input_mols[mol_idx], nb)
                                && !all_template_atoms.contains(&(mol_idx, nb))
                        })
                        .collect()
                } else {
                    Vec::new()
                };
                // An H-free mapped atom of an explicit-H reactant has its
                // count pinned to 0 by `add_hydrogens`; RDKit refills it when
                // the edit lowers its valence (`[C:1]=[O:2]>>[C:1]-[O:2]` on
                // explicit-H acetic acid gives `CC(O)O`, #754).
                let refill = explicit_h.is_empty()
                    && carry_substituents
                    && src_atom.chirality == Chirality::None
                    && src_atom.hydrogen_count == Some(0)
                    && src_atom.charge == 0
                    // An aromatic carbon the template spells aliphatic
                    // (`[#6;a:1]` in RDKit's reading) is refilled too.
                    && (!src_atom.aromatic || src_atom.element == chematic_core::Element::C)
                    && !new_atom.aromatic
                    && src_atom.element.is_organic_subset()
                    && new_atom.element.is_organic_subset()
                    && template_h(i).is_none()
                    && reactant_has_explicit_h(input_mols, &explicit_h_cache, mol_idx);
                if !explicit_h.is_empty() || refill {
                    // The folded H atoms are re-derived from valence below.
                    new_atom.hydrogen_count = template_h(i);
                    // An aromatic atom keeps its H count when the template
                    // leaves its degree, element and charge alone, as the
                    // implicit form keeps `[nH]`: valence inference cannot
                    // tell pyrrole's N from pyridine's.
                    if new_atom.hydrogen_count.is_none()
                        && new_atom.aromatic
                        && src_atom.aromatic
                        && degree_unchanged
                        && element_unchanged
                        && new_atom.charge == src_atom.charge
                    {
                        new_atom.hydrogen_count = Some(explicit_h.len() as u8);
                    }
                }
                let dearomatize = src_atom.aromatic && !new_atom.aromatic;
                let idx = builder.add_atom(new_atom);
                if dearomatize {
                    dearomatized.push(idx);
                }
                if refill {
                    folded_cores.push(idx);
                }
                if !explicit_h.is_empty() {
                    folded_h.extend(explicit_h.into_iter().map(|h| (mol_idx, h)));
                    folded_cores.push(idx);
                }
                src_to_new.insert((mol_idx, src_idx), idx);
                idx
            } else {
                // Map number not in reactants — new atom from template.
                let mut new_atom = tmpl_atom.clone();
                new_atom.atom_map = None;
                // RDKit does not pin H0 on a product-only atom: it uses
                // ordinary valence inference, as for a bare `[C]`/`[O]`.
                // A mapped reactant atom above still honours explicit H0.
                new_atom.hydrogen_count = template_h(i).filter(|&h| h != 0);
                builder.add_atom(new_atom)
            }
        } else {
            // No atom_map — entirely new atom from template. As for core atoms,
            // a bare bracket atom's `Some(0)` means "unspecified" in a product
            // template (`[*:1][C](=[O])[C]` adds an acetyl group, not three
            // radicals — issue #679); only an explicit `H<n>` with n > 0 pins
            // the count, as in SMIRKS. RDKit also leaves an explicit H0 on
            // a new atom unspecified; it is not a radical marker.
            let mut new_atom = tmpl_atom.clone();
            new_atom.atom_map = None;
            new_atom.hydrogen_count = template_h(i).filter(|&h| h != 0);
            builder.add_atom(new_atom)
        };
        *slot = Some(new_idx);
    }

    // --- Step 2: collect substituents and index source/product atoms ---
    let discoverer = carry_substituent_atoms(
        &mut builder,
        input_mols,
        all_template_atoms,
        &folded_h,
        core_keys,
        carry_substituents,
        keep_carried_maps,
        &mut src_to_new,
    );
    let DenseAtomMaps { source_of, new_of } =
        dense_atom_maps(builder.atom_count(), input_mols, &src_to_new);

    // --- Steps 3 and 4: template bonds, then carried source bonds ---
    let colon_atoms = add_product_template_bonds(
        &mut builder,
        product_template,
        &template_idx_to_new,
        &source_of,
        input_mols,
    );
    carry_source_bonds(
        &mut builder,
        &source_of,
        &new_of,
        input_mols,
        all_template_atoms,
        global_map,
        reactant_bond_maps,
    );

    // Parity-aware atom chirality correction (see the doc above
    // correct_product_stereo) -- must run after all bonds above are added,
    // since it validates against each atom's real final neighbour set.
    let product = correct_product_stereo(
        builder.build(),
        product_template,
        &template_idx_to_new,
        global_map,
        all_template_atoms,
        input_mols,
        &src_to_new,
    );

    // Carried stereocentres keep their geometry: re-express the reactant's
    // stereo record in product atoms, since the product's bond order around
    // them can differ from the reactant's (a ring re-entered from another
    // side).
    let mut product = product;
    let sourced_atoms = restore_carried_atom_stereo(
        &mut product,
        product_template,
        &template_idx_to_new,
        input_mols,
        &source_of,
        &discoverer,
        &src_to_new,
        keep_carried_maps,
    );

    // A folded H atom (explicit-H reactant, see above) can be the only
    // `/`/`\` carrier of its double-bond end (`[H]/C(=C(/[H])C)C`). Move
    // its direction to another substituent of that end, on the other side,
    // so the E/Z survives as it does in RDKit.
    transfer_folded_h_stereo(&mut product, &folded_h, input_mols, &new_of);

    // A `/`/`\` read on an aromatic ring bond is stashed beside the bond
    // (`Molecule::bond_direction`): it carries the E/Z of an exocyclic
    // double bond (`C/N=c1\c(O)c(O)c1`). Carry it onto the product bond
    // built from that ring bond while the double bond is still there.
    carry_direction_stash(&mut product, input_mols, &src_to_new);

    // Template markers on a double bond the template creates can sit next
    // to a carried marker of a conjugated double bond, asking the shared
    // end for the same side twice; readers then drop or misread one E/Z
    // (BioTransformer BTMR0709 on an enyne gave both E and Z). Rewrite the
    // markers consistently, each double bond keeping its own geometry.
    chematic_core::ez_markers::reconcile_ez_markers(&mut product);

    // Clear any Up/Down stereo markers left on bonds that are no longer adjacent
    // to a double bond (e.g. after C=C → C=O conversion via SMIRKS).
    let mut molecule = clear_orphaned_stereo_bonds(product);

    refill_product_hydrogens(
        &mut molecule,
        &sourced_atoms,
        input_mols,
        &template_idx_to_new,
        &explicit_h_cache,
        folded_cores,
        &source_of,
        &folded_h,
    );

    restore_product_aromaticity(&mut molecule, &colon_atoms, &dearomatized);
    apply_rdkit_template_hydrogens(&mut molecule, &template_idx_to_new);
    traced_product(
        molecule,
        product_template,
        template_idx_to_new,
        &sourced_atoms,
    )
}

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
