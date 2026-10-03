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

use crate::reaction::{RxnError, parse_reaction};
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

/// Why the pinned RDKit 2026.03.6 compatibility profile cannot safely claim
/// a reaction result. Native CheMatic reaction semantics are unchanged.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactionCompatibilityUnsupported {
    /// Reactant-side tetrahedral constraints have different matching semantics.
    ChiralReactantTemplateSemantics,
}

impl ReactionCompatibilityUnsupported {
    /// Stable machine-readable refusal code shared by bindings and audit gates.
    pub const fn reason_code(self) -> &'static str {
        match self {
            Self::ChiralReactantTemplateSemantics => "chiral_reactant_template_semantics",
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
    PreparedReaction::new(smirks)?.run_reactants_with_diagnostics(reactants, limits)
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
        products.extend(PreparedReaction::new(&variant)?.run_reactants_impl(
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
        let (_, queries) = parse_reactant_templates(parts[0])?;
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
    /// safely compare this template. Product-side stereo and E/Z are not
    /// blanket-refused; only reactant-side tetrahedral matching is affected.
    pub fn rdkit_2026_03_6_unsupported_reason(&self) -> Option<ReactionCompatibilityUnsupported> {
        self.has_tetrahedral_reactant_stereo()
            .then_some(ReactionCompatibilityUnsupported::ChiralReactantTemplateSemantics)
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
        let variants = crate::reaction::expand_atomic_number_primitives(smirks)?;
        if variants.len() > 1 || variants.first().is_none_or(|variant| variant != smirks) {
            let mut compiled = Vec::with_capacity(variants.len());
            for variant in &variants {
                compiled.push(Self::new_normalized(variant)?);
            }
            let mut primary = Self::new_normalized(&variants[0])?;
            primary.variants = Some(compiled);
            primary.requirements.element_lower_bounds.clear();
            primary.requirements.aromatic_element_lower_bounds.clear();
            primary.requirements.aliphatic_element_lower_bounds.clear();
            primary.requirements.bond_lower_bounds.clear();
            return Ok(primary);
        }
        Self::new_normalized(smirks)
    }

    fn new_normalized(smirks: &str) -> Result<Self, TransformError> {
        crate::perf_counters::record_reaction_parse_call();
        let (rxn, queries, product_specs) = parse_smirks_templates(smirks)?;
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
        let matches = find_matches_impl(self, reactants, limits, None)?;
        let accepted_matches = matches.len();
        let mut products = Vec::with_capacity(accepted_matches);
        let mut valence_rejected_matches = 0;
        for m in &matches {
            match apply_match_traced_impl(self, reactants, m, carry_substituents) {
                Some(product_set) => products.push(product_set),
                None => valence_rejected_matches += 1,
            }
        }
        Ok(TracedReactionTransformReport {
            diagnostics: ReactionTransformDiagnostics {
                accepted_matches,
                applied_products: products.len(),
                valence_rejected_matches,
                truncated_matches: false,
            },
            products,
        })
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
);

fn parse_smirks_templates(smirks: &str) -> Result<SmirksTemplates, TransformError> {
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
    let products = crate::reaction::normalize_product_query_atoms(&format!(">>{}", parts[2]))
        .map_err(TransformError::SmirksParse)?;
    let products = normalize_product_templates(&products[2..])?;
    // Agents and products through the ordinary reaction parser (with its
    // limits); the reactant slot is filled below.
    // Agents take no part in template application (RDKit ignores them too)
    // and may be SMARTS (`>[O;X2]>`); when they are not SMILES they are
    // dropped rather than failing the whole template.
    let mut rxn = parse_reaction(&format!(">{}>{}", parts[1], products))
        .or_else(|_| parse_reaction(&format!(">>{products}")))?;
    let product_specs = products
        .split('.')
        .filter(|p| !p.is_empty())
        .map(product_atom_specs)
        .collect();
    let (reactants, queries) = parse_reactant_templates(parts[0])?;
    rxn.reactants = reactants;
    Ok((rxn, queries, product_specs))
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
) -> Result<(Vec<Molecule>, Vec<QueryMolecule>), TransformError> {
    let limits = crate::reaction::ReactionParseLimits::default();
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
        let (mol, query) = match chematic_smiles::parse(part) {
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
                let query = chematic_smarts::parse_smarts(part).map_err(|smarts_err| {
                    TransformError::SmirksParse(RxnError::SmartsParse {
                        part: part.to_string(),
                        source: format!("{smarts_err} (as SMILES: {smiles_err})"),
                    })
                })?;
                let shadow = shadow_molecule_of(&query).map_err(|reason| {
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
    for (i, b) in side.bytes().enumerate() {
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
    if side.len() > start {
        out.push(&side[start..]);
    }
    out
}

/// The template molecule standing in for a SMARTS-only reactant query: one
/// atom per query atom (the implied element when the query pins one, else a
/// wildcard) carrying the atom map, and one single bond per query bond.
/// `Err` names the stereo feature that makes the query unsupported.
fn shadow_molecule_of(query: &QueryMolecule) -> Result<Molecule, &'static str> {
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
    if query.bonds.iter().any(|b| has_direction(&b.query)) {
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
                n.is_empty() || (n.len() == 1 && n.bytes().all(|b| b.is_ascii_digit()))
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
    let all_match_sets: Vec<Vec<FxHashMap<usize, AtomIdx>>> = prepared
        .queries
        .iter()
        .zip(reactants.iter())
        .enumerate()
        .map(|(index, (q, mol))| {
            let matches = match rings {
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

    // No matches when any template has no match.
    if all_match_sets.iter().any(|ms| ms.is_empty()) {
        return Ok(vec![]);
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

    for combo in cartesian_product(&all_match_sets) {
        crate::perf_counters::record_match_combination();

        // Parity-aware chirality post-check.  Runs only when the SMIRKS has @/@@.
        // This must happen after the complete VF2 mapping is known, because
        // correct chirality comparison requires the full neighbor permutation.
        if prepared.has_stereo {
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
                continue;
            }
        }

        matches.push(ReactionMatch {
            per_reactant: combo,
        });
    }

    Ok(matches)
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
            let mut product = build_product(
                pt,
                &global_map,
                &prepared.reactant_map_atoms,
                reactants,
                &all_template_atoms,
                carry_substituents,
                prepared.product_specs.get(pi).map_or(&[], Vec::as_slice),
            );
            // The same provenance map serves traced and ordinary reaction
            // APIs. Born atoms remain untagged; caller labels are not remapped
            // or made unique across reactants.
            for (i, source) in product.atom_sources.iter().enumerate() {
                if let Some(source) = source
                    && let Some(tag) = reactants[source.reactant].atom_tag(source.atom)
                {
                    product.molecule.set_tag(AtomIdx(i as u32), Some(tag.get()));
                }
            }
            crate::perf_counters::record_build_product_call(
                product.molecule.atom_count(),
                product.molecule.bond_count(),
            );
            product
        })
        .collect();

    // Skip product sets RDKit's sanitize step would reject (issue #734).
    if products.iter().all(|p| sanitizable_product(&p.molecule)) {
        crate::perf_counters::record_product_set();
        Some(products)
    } else {
        None
    }
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
fn sanitizable_product(mol: &Molecule) -> bool {
    use crate::rdkit_valence::max_valence;
    if mol.atoms().any(|(_, a)| a.aromatic) {
        let in_ring = atoms_in_rings(mol);
        if mol
            .atoms()
            .any(|(idx, a)| a.aromatic && !in_ring[idx.0 as usize])
        {
            return false;
        }
    }
    let kekule;
    let mol = if mol.bonds().any(|(_, b)| b.order == BondOrder::Aromatic) {
        match chematic_core::kekulize(mol) {
            Ok(k) => {
                kekule = chematic_core::apply_kekule(mol, &k);
                &kekule
            }
            Err(_) => return false,
        }
    } else {
        mol
    };
    let explicit_valence = |idx: AtomIdx| -> i16 {
        let bonds: i16 = mol
            .neighbors(idx)
            .map(|(_, b)| {
                let bond = mol.bond(b);
                match bond.order {
                    BondOrder::Dative if bond.atom1 == idx => 0,
                    order => i16::from(order.order_int()),
                }
            })
            .sum();
        bonds + i16::from(mol.atom(idx).hydrogen_count.unwrap_or(0))
    };
    mol.atoms().all(|(idx, atom)| {
        if atom.wildcard {
            return true;
        }
        let z = atom.element.atomic_number();
        let mut charge = atom.charge;
        let mut used = explicit_valence(idx);
        // RDKit's cleanUp, which runs before the valence check.
        if charge == 0 {
            let bond_to = |order: BondOrder, nz: u8| {
                mol.neighbors(idx).any(|(nb, b)| {
                    mol.bond(b).order == order && {
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
                    .filter(|&(_, b)| mol.bond(b).order == BondOrder::Double)
                    .count() as i16;
                charge = doubles as i8;
                used -= doubles;
            }
        }
        max_valence(z, charge).is_none_or(|max| used <= max)
    })
}

/// Per atom, whether it lies on a ring (has an incident non-bridge bond).
fn atoms_in_rings(mol: &Molecule) -> Vec<bool> {
    let n = mol.atom_count();
    let mut order = vec![usize::MAX; n];
    let mut low = vec![0usize; n];
    let mut in_ring = vec![false; n];
    let mut counter = 0usize;
    for root in 0..n {
        if order[root] != usize::MAX {
            continue;
        }
        // (atom, bond used to reach it, index of its next neighbour)
        let mut stack: Vec<(usize, Option<BondIdx>, usize)> = vec![(root, None, 0)];
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
            BondOrder::Quadruple | BondOrder::Zero | BondOrder::QueryAny => {
                BondQuery::Primitive(BondPrimitive::Any)
            }
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

/// Build one product molecule applying full SMIRKS semantics.
///
/// 1. Atom-mapped product atoms: copy source atom + override aromatic/charge/H from template.
/// 2. New product atoms (no map): clone from template.
/// 3. BFS from core (mapped) atoms through input molecules, collecting substituents
///    (non-template atoms reachable without crossing template-atom walls).
/// 4. Add product-template bonds (new/changed bonds).
/// 5. Carry through bonds from source molecules where at least one endpoint is a substituent.
fn build_product(
    product_template: &Molecule,
    global_map: &FxHashMap<u16, (usize, AtomIdx)>,
    reactant_template_atoms: &FxHashMap<u16, ReactantMapAtom>,
    input_mols: &[&Molecule],
    all_template_atoms: &FxHashSet<(usize, AtomIdx)>,
    carry_substituents: bool,
    specs: &[ProductAtomSpec],
) -> TracedProduct {
    let mut builder = MoleculeBuilder::new();
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
    let mut src_to_new: FxHashMap<(usize, AtomIdx), AtomIdx> = FxHashMap::default();

    // --- Step 1: add product template atoms ---
    // core_keys: only source atoms that are mapped by THIS product template.
    // Using global_map.values() (all matched atoms across all templates) would
    // seed the BFS in Step 2 from atoms belonging to *other* product templates,
    // causing their substituents to leak into this product (issue #13).
    let product_maps: FxHashSet<u16> = (0..product_template.atom_count())
        .filter_map(|i| product_template.atom(AtomIdx(i as u32)).atom_map)
        .collect();
    let core_keys: FxHashSet<(usize, AtomIdx)> = global_map
        .iter()
        .filter(|(am, _)| product_maps.contains(am))
        .map(|(_, &src)| src)
        .collect();

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
                // `[O-:1]>>[O:1]` leaves the alkoxide charged).
                if spec_of(i).charge {
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
                new_atom.hydrogen_count = match template_h(i) {
                    Some(h) => Some(h),
                    None if degree_unchanged => src_atom.hydrogen_count,
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
                if !explicit_h.is_empty() {
                    // The folded H atoms are re-derived from valence below.
                    new_atom.hydrogen_count = template_h(i);
                }
                let dearomatize = src_atom.aromatic && !new_atom.aromatic;
                let idx = builder.add_atom(new_atom);
                if dearomatize {
                    dearomatized.push(idx);
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

    // --- Step 2: BFS from core atoms to collect substituents ---
    // Skipped when carry_substituents = false (run_reactants_strict mode).
    // Seed visited with all template atoms so BFS cannot cross into the template region.
    let mut visited: FxHashSet<(usize, AtomIdx)> = all_template_atoms.clone();
    if !folded_h.is_empty() {
        visited.extend(folded_h.iter().copied());
    }
    if carry_substituents {
        let mut queue: VecDeque<(usize, AtomIdx)> = core_keys.iter().cloned().collect();

        while let Some((mol_idx, cur_idx)) = queue.pop_front() {
            for (nb_idx, _bond_idx) in input_mols[mol_idx].neighbors(cur_idx) {
                let key = (mol_idx, nb_idx);
                if visited.contains(&key) {
                    continue;
                }
                visited.insert(key);
                let src_atom = input_mols[mol_idx].atom(nb_idx);
                let mut new_atom = src_atom.clone();
                new_atom.atom_map = None;
                let new_idx = builder.add_atom(new_atom);
                src_to_new.insert(key, new_idx);
                queue.push_back(key);
            }
        }
    }

    // --- Step 3: add product template bonds ---
    let mut added_bond_pairs: FxHashSet<(AtomIdx, AtomIdx)> = FxHashSet::default();

    for (_bidx, bond) in product_template.bonds() {
        let a_new = template_idx_to_new[bond.atom1.0 as usize].unwrap();
        let b_new = template_idx_to_new[bond.atom2.0 as usize].unwrap();
        let _ = builder.add_bond(a_new, b_new, bond.order);
        added_bond_pairs.insert((a_new.min(b_new), a_new.max(b_new)));
    }

    // --- Step 4: carry-through bonds from source molecules ---
    // Bonds where both endpoints are template atoms are replaced or broken by the template;
    // bonds where at least one endpoint is a substituent are carried through.
    for (&(mol_idx, src_idx), &a_new) in &src_to_new {
        for (nb_idx, bond_idx) in input_mols[mol_idx].neighbors(src_idx) {
            let nb_key = (mol_idx, nb_idx);
            let Some(&b_new) = src_to_new.get(&nb_key) else {
                continue;
            };
            if all_template_atoms.contains(&(mol_idx, src_idx))
                && all_template_atoms.contains(&nb_key)
            {
                continue;
            }
            let pair = (a_new.min(b_new), a_new.max(b_new));
            if added_bond_pairs.contains(&pair) {
                continue;
            }
            added_bond_pairs.insert(pair);
            let ob = input_mols[mol_idx].bond(bond_idx);
            // Preserve the original atom1→atom2 orientation. Up/Down (E/Z)
            // bond semantics are direction-dependent, so adding the bond with
            // endpoints swapped relative to the source would flip the geometry.
            let (a, b) = if ob.atom1 == src_idx {
                (a_new, b_new)
            } else {
                (b_new, a_new)
            };
            let _ = builder.add_bond(a, b, ob.order);
        }
    }

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

    // Clear any Up/Down stereo markers left on bonds that are no longer adjacent
    // to a double bond (e.g. after C=C → C=O conversion via SMIRKS).
    let mut molecule = clear_orphaned_stereo_bonds(product);

    // Re-add the folded explicit hydrogens (see Step 1): as many H atoms as
    // the product valence implies, after which the atom's count is pinned to
    // its explicit H atoms, as `add_hydrogens` leaves it.
    for &core in &folded_cores {
        let n = chematic_core::implicit_hcount(&molecule, core);
        molecule.set_hydrogen_count(core, Some(0));
        for _ in 0..n {
            let mut h = chematic_core::Atom::new(chematic_core::Element::H);
            h.hydrogen_count = Some(0);
            let h_idx = molecule.add_atom(h);
            let _ = molecule.add_bond(core, h_idx, BondOrder::Single);
        }
    }

    // A mapped aromatic atom spelled aliphatic in the product (`[#6:1]` from
    // a SMARTS reactant expands to `[C:1]`) that still sits in its aromatic
    // ring stays aromatic, as RDKit's sanitize re-perceives it.
    for &idx in &dearomatized {
        if molecule
            .neighbors(idx)
            .filter(|&(_, b)| molecule.bond(b).order == BondOrder::Aromatic)
            .count()
            >= 2
        {
            molecule.set_atom_aromatic(idx, true);
        }
    }

    // Template atoms whose H count is left to valence get RDKit's implicit
    // count where the native inference differs: chematic infers hydrogens
    // for the organic subset only and with its own valence lists, RDKit for
    // every element (`[C:1]>>[Se:1]` gives `C[SeH]`, `[Cl+](C)(C)C` gains an
    // H). Atoms in aromatic systems keep the native, kekulization-aware
    // count.
    // Neutral organic-subset atoms are skipped: the native inference already
    // agrees with RDKit wherever the product can be valid.
    for &idx in template_idx_to_new.iter().flatten() {
        let atom = molecule.atom(idx);
        if (atom.charge == 0 && atom.element.is_organic_subset())
            || atom.hydrogen_count.is_some()
            || atom.wildcard
            || atom.aromatic
            || molecule
                .neighbors(idx)
                .any(|(_, b)| molecule.bond(b).order == BondOrder::Aromatic)
        {
            continue;
        }
        let bonds: i16 = molecule
            .neighbors(idx)
            .map(|(_, b)| {
                let bond = molecule.bond(b);
                match bond.order {
                    BondOrder::Dative if bond.atom1 == idx => 0,
                    order => i16::from(order.order_int()),
                }
            })
            .sum();
        let rdkit_h = crate::rdkit_valence::implicit_hydrogens(
            atom.element.atomic_number(),
            atom.charge,
            bonds,
        );
        if rdkit_h != chematic_core::implicit_hcount(&molecule, idx) {
            molecule.set_hydrogen_count(idx, Some(rdkit_h));
        }
    }

    // Both passes above keep atom indices, so `src_to_new` indexes the final
    // product. Atoms absent from it were created from the product template.
    let mut atom_sources = vec![None; molecule.atom_count()];
    for (&(reactant, atom), &new_idx) in &src_to_new {
        atom_sources[new_idx.0 as usize] = Some(ReactantAtom { reactant, atom });
    }
    let mut template_maps = vec![None; molecule.atom_count()];
    for (template_idx, product_idx) in template_idx_to_new.iter().enumerate() {
        if let Some(product_idx) = product_idx {
            template_maps[product_idx.0 as usize] =
                product_template.atom(AtomIdx(template_idx as u32)).atom_map;
        }
    }
    TracedProduct {
        molecule,
        atom_sources,
        template_maps,
    }
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
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn canonical(mol: &Molecule) -> String {
        chematic_smiles::canonical_smiles(mol)
    }

    fn canonical_set(set: Vec<Molecule>) -> Vec<String> {
        set.into_iter().map(|mol| canonical(&mol)).collect()
    }

    #[test]
    fn identity_single_atom() {
        let mol = parse("C").unwrap();
        let results = run_reactants("[C:1]>>[C:1]", &[&mol]).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].len(), 1);
        assert_eq!(results[0][0].atom_count(), 1);
    }

    #[test]
    fn mapped_asymmetric_ether_keeps_both_embeddings_on_one_atom_set() {
        // The C-O-C query covers the same three target atoms in both
        // orientations, but the mapped cleavage yields two different product
        // tuples. Target-set uniquification must not discard either one.
        let mol = parse("COCC").unwrap();
        let results = run_reactants("[C:1][O:2][C:3]>>[C:1][O:2].[C:3]", &[&mol]).unwrap();
        let outcomes: std::collections::BTreeSet<Vec<String>> = results
            .into_iter()
            .map(|set| {
                let mut products = canonical_set(set);
                products.sort();
                products
            })
            .collect();
        let expected: std::collections::BTreeSet<Vec<String>> =
            [vec!["C", "CCO"], vec!["CC", "CO"]]
                .into_iter()
                .map(|set| {
                    let mut products: Vec<String> = set
                        .into_iter()
                        .map(|smiles| canonical(&parse(smiles).unwrap()))
                        .collect();
                    products.sort();
                    products
                })
                .collect();
        assert_eq!(outcomes, expected);
    }

    #[test]
    fn applies_atomic_number_smirks_and_keeps_atom_maps() {
        let reactant = parse("NC=O").unwrap();
        let smirks = "[#7:1][C:2](=[O:3])>>[#7:1][C:2](=[O:3])";
        let normalized = crate::reaction::expand_atomic_number_primitives(smirks)
            .unwrap()
            .into_iter()
            .next()
            .unwrap();
        let matches = find_reaction_matches(&normalized, &[&reactant]).unwrap();
        assert_eq!(matches.len(), 1);
        assert_eq!(
            matches[0].atom_map_positions(&normalized).unwrap().len(),
            3,
            "map identity must survive query matching"
        );
        let results = run_reactants(smirks, &[&reactant]).unwrap();
        assert!(!results.is_empty(), "atomic-number SMIRKS must apply");
        let product = &results[0][0];
        assert_eq!(product.atom_count(), 3);
    }

    #[test]
    fn styrene_hydrogenation_keeps_both_template_map_orientations() {
        let reactant = parse("C=Cc1ccccc1").unwrap();
        let smirks = "[C:1]=[C:2]>>[C:1][C:2]";
        let prepared = PreparedReaction::new(smirks).unwrap();
        let matches = prepared.find_matches(&[&reactant]).unwrap();
        let maps: Vec<_> = matches
            .iter()
            .map(|m| m.atom_map_positions(smirks).unwrap())
            .collect();
        assert_eq!(
            maps.len(),
            2,
            "each atom-map orientation is observable: {maps:?}"
        );
        let orientations: std::collections::BTreeSet<_> = maps
            .iter()
            .map(|positions| (positions[&1].1.0, positions[&2].1.0))
            .collect();
        assert_eq!(orientations, [(0, 1), (1, 0)].into());
    }

    #[test]
    fn reactant_side_smarts_queries_apply(/* issue #734 */) {
        // The SMIRKS reactant side now parses as full SMARTS: query-only
        // features (D, X, ;, ',', !, $()) match and the identity transform
        // returns the input, where before every one of these failed to parse.
        let canon = |m: &Molecule| chematic_smiles::canonical_smiles(m);
        let ethane = parse("CC").unwrap();
        let ethanol = parse("CCO").unwrap();
        for (smirks, reactant) in [
            ("[CD1:1]>>[C:1]", &ethane),
            ("[CX4:1]>>[C:1]", &ethane),
            ("[C;H3:1]>>[C:1]", &ethane),
            ("[C,N:1]>>[C:1]", &ethane),
            ("[!O:1]>>[C:1]", &ethane),
            ("[$([OH]):1]>>[O:1]", &ethanol),
        ] {
            let products = run_reactants(smirks, &[reactant]).unwrap();
            assert!(!products.is_empty(), "{smirks}: expected a product");
            assert!(
                products.iter().all(|set| canon(&set[0]) == canon(reactant)),
                "{smirks}: identity transform must return the input"
            );
        }
        // Product-side query features are specifications in RDKit's sense:
        // query-only parts are dropped (`[OX2H1]` applies H1), a mapped atom
        // without a single element keeps the reactant's element, and an
        // unmapped one cannot be built.
        let anisole = parse("COc1ccccc1").unwrap();
        let out = run_reactants("[c:1][OX2:2][CH3:3]>>[c:1][OX2H1:2]", &[&anisole]).unwrap();
        assert_eq!(out.len(), 1);
        assert_eq!(canon(&out[0][0]), canon(&parse("Oc1ccccc1").unwrap()));
        let out = run_reactants("[C:1]>>[C,N:1]", &[&ethane]).unwrap();
        assert!(out.iter().all(|set| canon(&set[0]) == canon(&ethane)));
        assert!(matches!(
            run_reactants("[C:1]>>[C:1][C,N]", &[&ethane]),
            Err(TransformError::SmirksParse(
                crate::reaction::RxnError::UnsupportedProductPrimitive { .. }
            ))
        ));
        // `;H<n>` is the SMARTS spelling of an H count.
        let out = run_reactants("[C:1]>>[C;H3:1]", &[&ethane]).unwrap();
        assert!(!out.is_empty());
    }

    #[test]
    fn mapped_atom_takes_a_changed_product_element(/* issue #734 item 2 */) {
        // RDKit 2026.03.6 reference products for each case.
        let canon = |m: &Molecule| chematic_smiles::canonical_smiles(m);
        for (smirks, reactant, expected) in [
            ("[C:1]>>[N:1]", "CC", vec!["CN"]),
            ("[#6:1]>>[#7:1]", "CC", vec!["CN"]),
            ("[CH3:1]>>[N:1]", "CC", vec!["CN"]),
            ("[C:1][O:2]>>[C:1][S:2]", "CO", vec!["CS"]),
            ("[c:1]>>[n:1]", "c1ccccc1", vec!["c1ccncc1"]),
            ("[*:1]>>[N:1]", "C", vec!["N"]),
            // No change: wildcard product, or the product element equals the
            // reactant query's (first) element.
            ("[C:1]>>[*:1]", "CO", vec!["CO"]),
            ("[C,N:1]>>[C:1]", "CN", vec!["CC", "CN"]),
            ("[C,N:1]>>[N:1]", "CC", vec!["CN"]),
        ] {
            let mol = parse(reactant).unwrap();
            let mut got: Vec<String> = run_reactants(smirks, &[&mol])
                .unwrap()
                .iter()
                .map(|set| canon(&set[0]))
                .collect();
            got.sort();
            got.dedup();
            let mut want: Vec<String> =
                expected.iter().map(|s| canon(&parse(s).unwrap())).collect();
            want.sort();
            assert_eq!(got, want, "{smirks} on {reactant}");
        }
    }

    #[test]
    fn charge_change_keeps_rdkit_hydrogen_count(/* issue #734 item 4 */) {
        // RDKit 2026.03.6 reference products.
        let canon = |m: &Molecule| chematic_smiles::canonical_smiles(m);
        for (smirks, reactant, expected) in [
            ("[N:1]>>[N+2:1]", "CN", "C[NH2+2]"),
            ("[N:1]>>[N++:1]", "CN", "C[NH2+2]"),
            ("[N:1]>>[N+:1]", "CN", "C[NH3+]"),
            ("[N:1]>>[N-:1]", "CN", "C[NH-]"),
            ("[O:1]>>[O-:1]", "CO", "C[O-]"),
            ("[C:1]>>[C-:1]", "CC", "[CH2-]C"),
            ("[C:1]>>[C+:1]", "CC", "[CH2+]C"),
        ] {
            let mol = parse(reactant).unwrap();
            let out = run_reactants(smirks, &[&mol]).unwrap();
            assert!(!out.is_empty(), "{smirks}");
            assert!(
                out.iter()
                    .all(|set| canon(&set[0]) == canon(&parse(expected).unwrap())),
                "{smirks} on {reactant}: {:?}",
                out.iter().map(|s| canon(&s[0])).collect::<Vec<_>>()
            );
        }
    }

    #[test]
    fn explicit_hydrogen_reactants_give_the_implicit_products(/* issue #734 item 5 */) {
        // Every case the issue lists: the explicit-H reactant must give the
        // same products (compared without explicit H) as its implicit form.
        let canon = |m: &Molecule| chematic_smiles::canonical_smiles(m);
        let strip = |m: &Molecule| canon(&chematic_chem::remove_hydrogens(m));
        for (smirks, reactant) in [
            ("[C:1]-[C:2]>>[C:1]=[C:2]", "CC"),
            ("[O:1]>>[O:1]C(=O)C", "CCO"),
            ("[C:1]>>[C:1]Cl", "CC"),
            ("[O:1]>>[O:1]C", "CCO"),
            ("[C:1][O:2]>>[C:1]=[O:2]", "CO"),
            ("[C:1][O:2]>>[C:1]=[O:2]", "CCO"),
            ("[C:1][N:2]>>[C:1]=[N+:2]", "CN(C)C"),
            ("[CH3:1]>>[CH3:1]", "CC"),
            ("[C:1]>>[C:1]", "CC"),
            ("[N:1]>>[N+:1]", "CN"),
            ("[C:1][O:2]>>[C:1].[O:2]", "CCO"),
        ] {
            let implicit = parse(reactant).unwrap();
            let explicit = chematic_chem::add_hydrogens(&implicit);
            let mut want: Vec<String> = run_reactants(smirks, &[&implicit])
                .unwrap()
                .iter()
                .map(|set| set.iter().map(canon).collect::<Vec<_>>().join("."))
                .collect();
            let mut got: Vec<String> = run_reactants(smirks, &[&explicit])
                .unwrap()
                .iter()
                .map(|set| set.iter().map(&strip).collect::<Vec<_>>().join("."))
                .collect();
            want.sort();
            want.dedup();
            got.sort();
            got.dedup();
            assert!(!want.is_empty(), "{smirks}: implicit form must react");
            assert_eq!(got, want, "{smirks} on explicit-H {reactant}");
        }
        // The edited atoms keep the input's explicit-H mode: ethane with
        // explicit H gives ethylene as C2 plus four H atoms.
        let ethane = chematic_chem::add_hydrogens(&parse("CC").unwrap());
        let out = run_reactants("[C:1]-[C:2]>>[C:1]=[C:2]", &[&ethane]).unwrap();
        assert!(!out.is_empty());
        for set in out {
            assert_eq!(set[0].atom_count(), 6);
            assert!(
                set[0]
                    .atoms()
                    .all(|(i, _)| chematic_core::implicit_hcount(&set[0], i) == 0)
            );
        }
        // A template that matches explicit hydrogens still sees them.
        let methane = chematic_chem::add_hydrogens(&parse("C").unwrap());
        assert!(
            !run_reactants("[C:1][H]>>[C:1]O", &[&methane])
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn bug_check_products_match_rdkit(/* issue #734 bug check */) {
        // Each row: SMIRKS, reactant, RDKit 2026.03.6's sanitized products
        // (empty: no product survives sanitize).
        let canon = |s: &str| chematic_smiles::canonical_smiles(&parse(s).unwrap());
        for (smirks, reactant, expected) in [
            // A dative `->` in a template is a bond, not a reaction arrow.
            ("[N:1]->[Pt:2]>>[N:1].[Pt:2]", "CN->[Pt]", vec!["CN.[Pt]"]),
            // RDKit's clean-up spellings are valid products.
            (
                "[C:1][Br:2]>>[C:1][Cl:2]",
                "BrCCN(=O)=O",
                vec!["ClCCN(=O)=O"],
            ),
            (
                "[c:1][Cl:2]>>[c:1]N(=O)=O",
                "Clc1ccccc1",
                vec!["O=N(=O)c1ccccc1"],
            ),
            ("[O:1]>>[O:1]Cl(=O)(=O)=O", "CO", vec!["COCl(=O)(=O)=O"]),
            ("[C:1][Br:2]>>[C:1]N=N#N", "CBr", vec!["CN=N#N"]),
            // RDKit's valence model: Mg is unrestricted, Cl+ is S-like.
            ("[C:1]>>[Mg:1](C)(C)C", "CC", vec!["C[Mg](C)(C)C"]),
            ("[Cl:1]>>[Cl+:1](C)C", "CCl", vec!["C[ClH+](C)C"]),
            ("[S:1]>>[S+:1](C)(C)(C)(C)C", "CS", vec![]),
            // Unchanged template degree keeps the reactant's explicit H.
            ("[n:1]>>[n:1]", "c1cc[nH]c1", vec!["c1cc[nH]c1"]),
            ("[n:1]>>[n:1]C", "c1cc[nH]c1", vec!["Cn1cccc1"]),
            // RDKit implicit hydrogens beyond the organic subset.
            ("[C:1]>>[Se:1]", "CC", vec!["C[SeH]"]),
            ("[C:1]>>[Si:1]", "CC", vec!["C[SiH3]"]),
            // SMARTS reactant semantics.
            ("[*:1]>>[*:1]O", "CN", vec!["CNO", "NCO"]),
            ("[#7;H0:1]>>[#7:1]C", "CNC", vec![]),
            ("[NH0:1]>>[N:1]C", "CNC", vec![]),
            ("[N+0:1]>>[N:1]C", "C[NH3+]", vec![]),
            // Aromatic atoms outside rings are rejected.
            ("[C:1]>>[c:1]", "CC", vec![]),
            // Unspecified product charge and isotope keep the reactant's;
            // an explicit H0 applies.
            ("[O-:1]>>[O:1]", "C[O-]", vec!["C[O-]"]),
            ("[N:1]>>[N:1]C", "C[NH3+]", vec!["C[NH2+]C"]),
            ("[C:1]>>[13C:1]", "CC", vec!["C[13CH3]"]),
            ("[C:1]>>[CH0:1]", "CC", vec!["[C]C"]),
            // RDKit 2026.03.6 treats H0 on a newly created, unmapped atom
            // as unspecified, unlike H0 on a mapped product atom above.
            ("[C:1]>>[C:1][CH0]", "C", vec!["CC"]),
            ("[C:1]>>[C:1][CH0:2]", "C", vec!["CC"]),
            ("[C:1]>>[C:1][OH0]", "C", vec!["CO"]),
            ("[C:1]>>[C:1][NH0]", "C", vec!["CN"]),
        ] {
            let mol = parse(reactant).unwrap();
            let mut got: Vec<String> = run_reactants(smirks, &[&mol])
                .unwrap_or_else(|e| panic!("{smirks}: {e}"))
                .iter()
                .map(|set| {
                    let mut parts: Vec<String> = set.iter().map(canonical).collect();
                    parts.sort();
                    parts.join(".")
                })
                .collect();
            got.sort();
            got.dedup();
            let mut want: Vec<String> = expected
                .iter()
                .map(|e| {
                    let mut parts: Vec<String> = e.split('.').map(canon).collect();
                    parts.sort();
                    parts.join(".")
                })
                .collect();
            want.sort();
            assert_eq!(got, want, "{smirks} on {reactant}");
        }
        // Atom-map positions read only the reactant side, so a product
        // atomic-number primitive is fine.
        let mol = parse("CN").unwrap();
        let smirks = "[N:1]>>[#7:1]C";
        let matches = find_reaction_matches(smirks, &[&mol]).unwrap();
        assert_eq!(matches[0].atom_map_positions(smirks).unwrap().len(), 1);
    }

    #[test]
    fn products_rdkit_sanitize_rejects_are_dropped(/* issue #734 item 3 */) {
        // RDKit 2026.03.6 rejects each of these raw products at sanitize; the
        // native model kept the last three (N valence 3 or 5, aromatic bond
        // counted as one).
        for (smirks, reactant) in [
            ("[C:1][O:2]>>[C:1]=[O:2]", "COC"),
            ("[C:1][N:2]>>[C:1]=[N:2]", "CN(C)C"),
            ("[N:1]>>[N:1]=O", "CNC"),
        ] {
            let mol = parse(reactant).unwrap();
            assert!(
                run_reactants(smirks, &[&mol]).unwrap().is_empty(),
                "{smirks} on {reactant}"
            );
        }
        // Ipso attack on toluene's methyl-bearing carbon is dropped; the
        // ortho/meta/para products remain (RDKit: the same three).
        let canon = |m: &Molecule| chematic_smiles::canonical_smiles(m);
        let toluene = parse("Cc1ccccc1").unwrap();
        let mut got: Vec<String> = run_reactants("[c:1]>>[c:1]O", &[&toluene])
            .unwrap()
            .iter()
            .map(|set| canon(&set[0]))
            .collect();
        got.sort();
        got.dedup();
        let mut want: Vec<String> = ["Cc1ccccc1O", "Cc1cccc(O)c1", "Cc1ccc(O)cc1"]
            .iter()
            .map(|s| canon(&parse(s).unwrap()))
            .collect();
        want.sort();
        assert_eq!(got, want);
        // The charged spellings are valid and kept.
        let tma = parse("CN(C)C").unwrap();
        assert!(
            !run_reactants("[C:1][N:2]>>[C:1]=[N+:2]", &[&tma])
                .unwrap()
                .is_empty()
        );
        assert!(
            !run_reactants("[N:1]>>[N+:1][O-]", &[&tma])
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn reactant_smarts_query_selects_the_matching_atom(/* issue #734 */) {
        // A reactant query that only some atoms satisfy must transform only
        // those atoms, proving the query is actually evaluated (not ignored).
        // Primary alcohol carbon oxidation: only the CH2 next to OH matches.
        let canon = |m: &Molecule| chematic_smiles::canonical_smiles(m);
        let ethanol = parse("CCO").unwrap();
        let products = run_reactants("[CX4;H2:1][OX2H:2]>>[C:1]=[O:2]", &[&ethanol]).unwrap();
        assert_eq!(products.len(), 1);
        assert_eq!(canon(&products[0][0]), canon(&parse("CC=O").unwrap()));
        // The isopropanol carbon is CH, so the H2 query must not match it.
        let isopropanol = parse("CC(C)O").unwrap();
        assert!(
            run_reactants("[CX4;H2:1][OX2H:2]>>[C:1]=[O:2]", &[&isopropanol])
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn reactant_smarts_only_stereo_is_a_typed_error(/* issue #734 */) {
        // A SMARTS-only reactant template (one that is not valid SMILES) that
        // carries tetrahedral or E/Z stereo is refused with a typed error
        // rather than silently dropping the constraint.
        let subject = parse("F[C@H](Cl)Br").unwrap();
        assert!(matches!(
            run_reactants("[C@;X4:1](F)(Cl)Br>>[C:1](F)(Cl)I", &[&subject]),
            Err(TransformError::SmirksParse(
                crate::reaction::RxnError::UnsupportedReactantTemplate { .. }
            ))
        ));
        // A purely SMILES reactant template with the same stereo still works
        // (it goes through the SMILES parser and the existing post-check).
        assert!(run_reactants("[C@H:1](F)(Cl)Br>>[C:1](F)(Cl)I", &[&subject]).is_ok());
    }

    #[test]
    fn bracket_hydrogen_atom_parses_as_hydrogen(/* issue #734 */) {
        // `[H]` inside a SMARTS reactant template is a hydrogen atom, not an
        // H-count primitive: it matches an explicit hydrogen neighbour.
        let mol = parse("C").unwrap();
        let mol = chematic_chem::add_hydrogens(&mol);
        let products = run_reactants("[C:1][H]>>[C:1]", &[&mol]).unwrap();
        assert!(!products.is_empty(), "[H] must match an explicit hydrogen");
    }

    #[test]
    fn product_atomic_number_primitives_are_literal_atoms_with_implicit_h() {
        // Issue #679: an unmapped product-side `[#6](=[#8])[#6]` adds an acetyl
        // group exactly like `C(=O)C`; no aromatic-flag variants, no `[C]`
        // radicals, one product.
        let ethanol = parse("CCO").unwrap();
        let canon = |m: &Molecule| chematic_smiles::canonical_smiles(m);
        let organic = run_reactants("[O:1]>>[*:1]C(=O)C", &[&ethanol]).unwrap();
        let atomic = run_reactants("[O:1]>>[*:1][#6](=[#8])[#6]", &[&ethanol]).unwrap();
        assert_eq!(organic.len(), 1);
        assert_eq!(atomic.len(), 1);
        assert_eq!(canon(&atomic[0][0]), canon(&organic[0][0]));
        assert_eq!(canon(&atomic[0][0]), "CC(=O)OCC");
        // Bracket product atoms behave the same: `[C]` is an unspecified
        // hydrogen count in a product template, as for mapped atoms (#18).
        let bracket = run_reactants("[O:1]>>[*:1][C](=[O])[C]", &[&ethanol]).unwrap();
        assert_eq!(canon(&bracket[0][0]), "CC(=O)OCC");
        // `;H1` pins the count; non-organic-subset elements stay bracketed.
        let vinyl = run_reactants("[O:1]>>[*:1][#6;H1]=[#6]", &[&ethanol]).unwrap();
        assert_eq!(canon(&vinyl[0][0]), canon(&parse("C=COCC").unwrap()));
        let silyl = run_reactants("[O:1]>>[*:1][#14]", &[&ethanol]).unwrap();
        assert_eq!(silyl.len(), 1);
        assert_eq!(
            silyl[0][0]
                .atoms()
                .filter(|(_, a)| a.element.atomic_number() == 14)
                .count(),
            1
        );
    }

    #[test]
    fn mapped_product_atomic_number_primitive_keeps_reactant_aromaticity() {
        // `[#7:1]` on both sides shares one spelling per variant, so an
        // aliphatic amine gives one aliphatic product and pyridine one aromatic
        // product (previously each gave an extra product with the flag flipped).
        let canon = |m: &Molecule| chematic_smiles::canonical_smiles(m);
        let amine = parse("CCN").unwrap();
        let out = run_reactants("[#7:1]>>[#7:1]C", &[&amine]).unwrap();
        assert_eq!(out.len(), 1);
        assert_eq!(canon(&out[0][0]), canon(&parse("CCNC").unwrap()));
        // Pyridine needs the charge (a neutral N-methylpyridine has a
        // four-valent N, rejected like RDKit's sanitize does, #734).
        let pyridine = parse("c1ccncc1").unwrap();
        assert!(
            run_reactants("[#7:1]>>[#7:1]C", &[&pyridine])
                .unwrap()
                .is_empty()
        );
        let out = run_reactants("[#7:1]>>[#7+:1]C", &[&pyridine]).unwrap();
        assert_eq!(out.len(), 1);
        assert!(
            out[0][0]
                .atoms()
                .any(|(_, a)| a.aromatic && a.element.atomic_number() == 7)
        );
        // A product map spelled with a plain symbol on the reactant side
        // follows that spelling.
        let variants = crate::reaction::expand_atomic_number_primitives("[n:1]>>[#7:1]C").unwrap();
        assert_eq!(variants, vec!["[n:1]>>[n:1]C".to_string()]);
        let variants = crate::reaction::expand_atomic_number_primitives("[N:1]>>[#7:1]C").unwrap();
        assert_eq!(variants, vec!["[N:1]>>[N:1]C".to_string()]);
    }

    #[test]
    fn prepared_reaction_applies_atomic_number_smirks_variants() {
        let reactant = parse("NC=O").unwrap();
        let prepared = PreparedReaction::new("[#7:1][C:2](=[O:3])>>[#7:1][C:2](=[O:3])").unwrap();
        let products = prepared.run_reactants(&[&reactant]).unwrap();
        // Two variants (`[N:1]`/`[n:1]` on both sides together); only the
        // aliphatic one matches formamide, and the product keeps the reactant's
        // aliphatic nitrogen (issue #679: no aromatic-flag duplicate).
        assert_eq!(products.len(), 1);
        assert!(products.iter().all(|set| set.len() == 1));
        assert!(!products[0][0].atoms().any(|(_, a)| a.aromatic));
    }

    #[test]
    fn prepared_reaction_reports_each_atomic_number_variant() {
        let reactant = parse("NC=O").unwrap();
        let prepared = PreparedReaction::new("[#7:1][C:2](=[O:3])>>[#7:1][C:2](=[O:3])").unwrap();
        let reports = prepared
            .run_reactants_with_variant_diagnostics(
                &[&reactant],
                &ReactionTransformLimits::default(),
            )
            .unwrap();
        assert_eq!(reports.len(), 2);
        assert_eq!(reports[0].variant_index, 0);
        assert_eq!(reports[1].variant_index, 1);
        assert_eq!(
            reports[0].normalized_smirks,
            "[N:1][C:2](=[O:3])>>[N:1][C:2](=[O:3])"
        );
        assert_eq!(
            reports[1].normalized_smirks,
            "[n:1][C:2](=[O:3])>>[n:1][C:2](=[O:3])"
        );
        assert!(reports[0].diagnostics.accepted_matches > 0);
        assert!(reports[0].diagnostics.applied_products > 0);
        assert!(reports[1].diagnostics.accepted_matches == 0);
    }

    #[test]
    fn atomic_number_smirks_rejects_unsupported_compound_product_primitive() {
        let reactant = parse("N").unwrap();
        // An unmapped product atom naming no single element cannot be built
        // and is rejected (#734); `;X3` alone would just be dropped.
        let err = run_reactants("[#7:1]>>[#7:1][#6,#7]", &[&reactant]);
        assert!(matches!(
            err,
            Err(TransformError::SmirksParse(
                crate::reaction::RxnError::UnsupportedProductPrimitive { .. }
            ))
        ));
        assert!(run_reactants("[#7:1]>>[#7;X3:1]", &[&reactant]).is_ok());
    }

    #[test]
    fn prepared_reaction_is_send_sync_and_matches_legacy_output() {
        fn assert_send_sync<T: Send + Sync>() {}
        assert_send_sync::<PreparedReaction>();
        assert_send_sync::<ReactionMatchContext>();

        let mol = parse("CCOC(=O)C").unwrap();
        let smirks = "[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])O.[O:3][C:4]";
        let prepared = PreparedReaction::new(smirks).unwrap();
        let legacy = run_reactants(smirks, &[&mol]).unwrap();
        let reused = prepared.run_reactants(&[&mol]).unwrap();
        let rings = chematic_perception::find_sssr(&mol);
        let reused_with_rings = prepared
            .run_reactants_with_rings(&[&mol], &[&rings])
            .unwrap();
        let canonical_sets = |sets: Vec<Vec<Molecule>>| {
            sets.into_iter()
                .map(|set| set.into_iter().map(|m| canonical(&m)).collect::<Vec<_>>())
                .collect::<Vec<_>>()
        };
        let legacy = canonical_sets(legacy);
        assert_eq!(canonical_sets(reused), legacy);
        assert_eq!(canonical_sets(reused_with_rings), legacy);

        let strict = prepared
            .run_reactants_strict_with_rings(&[&mol], &[&rings])
            .unwrap();
        let strict_limited = prepared
            .run_reactants_strict_with_rings_and_limits(
                &[&mol],
                &[&rings],
                &ReactionTransformLimits::default(),
            )
            .unwrap();
        assert_eq!(canonical_sets(strict), canonical_sets(strict_limited));
        assert_eq!(
            prepared
                .find_matches_with_rings_and_limits(
                    &[&mol],
                    &[&rings],
                    &ReactionTransformLimits::default(),
                )
                .unwrap(),
            prepared
                .find_matches_with_rings(&[&mol], &[&rings])
                .unwrap()
        );
        let context = ReactionMatchContext::new(&mol, Some(&rings));
        assert_eq!(
            prepared.find_matches_with_context(&context).unwrap(),
            prepared
                .find_matches_with_rings(&[context.target()], &[context.rings()])
                .unwrap()
        );
    }

    #[test]
    fn prepared_reaction_validates_once_and_preserves_limit_errors() {
        assert!(PreparedReaction::new("not-smirks").is_err());
        let mol = parse("CC").unwrap();
        let prepared = PreparedReaction::new("[C:1]>>[C:1]").unwrap();
        let err = match prepared
            .run_reactants_with_limits(&[&mol], &ReactionTransformLimits { max_matches: 0 })
        {
            Ok(_) => panic!("zero match limit must reject the first match set"),
            Err(err) => err,
        };
        assert!(matches!(
            err,
            TransformError::ResourceLimit {
                resource: "reaction matches",
                ..
            }
        ));
    }

    #[test]
    fn prepared_find_and_apply_match_legacy_split_api() {
        let mol = parse("NCCN").unwrap();
        let smirks = "[N:1]>>[N:1]";
        let prepared = PreparedReaction::new(smirks).unwrap();

        let legacy_matches = find_reaction_matches(smirks, &[&mol]).unwrap();
        let prepared_matches = prepared.find_matches(&[&mol]).unwrap();
        assert_eq!(prepared_matches, legacy_matches);
        for m in &prepared_matches {
            let legacy = apply_reaction_match(smirks, &[&mol], m, true).unwrap();
            let reused = prepared.apply_match(&[&mol], m, true).unwrap();
            assert_eq!(
                reused.map(canonical_set),
                legacy.map(canonical_set),
                "prepared match application must preserve legacy output"
            );
        }
    }

    #[test]
    fn reaction_transform_limits_reject_matches_before_cartesian_product() {
        let mol = parse("CC").unwrap();
        let err = find_reaction_matches_with_limits(
            "[C:1]>>[C:1]",
            &[&mol],
            &ReactionTransformLimits { max_matches: 0 },
        );
        assert!(matches!(
            err,
            Err(TransformError::ResourceLimit {
                resource: "reaction matches",
                ..
            })
        ));
    }

    #[test]
    fn no_match_returns_empty() {
        let mol = parse("C").unwrap();
        let results = run_reactants("[N:1]>>[N:1]", &[&mol]).unwrap();
        assert!(
            results.is_empty(),
            "nitrogen template must not match methane"
        );
    }

    #[test]
    fn multiple_matches_in_single_mol() {
        let mol = parse("NCCN").unwrap();
        let results = run_reactants("[N:1]>>[N:1]", &[&mol]).unwrap();
        assert_eq!(results.len(), 2, "two N atoms in NCCN → two product sets");
    }

    #[test]
    fn bond_formation_two_mols() {
        let n_mol = parse("N").unwrap();
        let c_mol = parse("C").unwrap();
        let results = run_reactants("[N:1].[C:2]>>[N:1][C:2]", &[&n_mol, &c_mol]).unwrap();
        assert!(!results.is_empty());
        let prod = &results[0][0];
        assert_eq!(prod.atom_count(), 2, "product must have 2 atoms");
        assert_eq!(prod.bonds().count(), 1, "product must have 1 bond");
    }

    #[test]
    fn bond_cleavage_two_products() {
        let mol = parse("CC").unwrap();
        let results = run_reactants("[C:1][C:2]>>[C:1].[C:2]", &[&mol]).unwrap();
        assert!(!results.is_empty());
        let products = &results[0];
        assert_eq!(products.len(), 2, "two product templates → two products");
        assert_eq!(products[0].atom_count(), 1);
        assert_eq!(products[1].atom_count(), 1);
    }

    #[test]
    fn reactant_count_mismatch_error() {
        let mol = parse("C").unwrap();
        let err = run_reactants("[N:1].[C:2]>>[N:1][C:2]", &[&mol]);
        assert!(
            matches!(
                err,
                Err(TransformError::ReactantCountMismatch {
                    expected: 2,
                    got: 1
                })
            ),
            "two-template SMIRKS with one reactant must error"
        );
    }

    #[test]
    fn invalid_smirks_error() {
        let mol = parse("C").unwrap();
        let err = run_reactants("[X]>>[X]", &[&mol]);
        assert!(
            matches!(err, Err(TransformError::SmirksParse(_))),
            "unknown element must yield SmirksParse error"
        );
    }

    #[test]
    fn overvalent_product_filtered_oxygen() {
        // O normally has max valence 2.
        // SMIRKS adds two carbons to an oxygen that already has one bond → 3 bonds on O → invalid.
        // CCO: the O is bonded to 1 C (bond_sum=1). Template [O:1]>>[O:1](C)C adds 2 more.
        let ethanol = parse("CCO").unwrap();
        let results = run_reactants("[O:1]>>[O:1](C)C", &[&ethanol]).unwrap();
        // The O that already had 1 bond would get 3 → over-valenced → filtered out.
        // The only match is the terminal O (1 bond → +2 = 3 bonds, invalid).
        assert!(
            results.is_empty(),
            "product with O having 3 bonds must be filtered out, got {} sets",
            results.len()
        );
    }

    #[test]
    fn valid_charged_product_kept() {
        // N with charge +1 can have up to 4 bonds (normal valences [3,5], +1 allows 4).
        // trimethylamine N(C)(C)C has N with bond_sum=3, charge=0.
        // Template [N:1]>>[N+:1] just changes charge, keeps 3 bonds → valid.
        let tma = parse("N(C)(C)C").unwrap();
        let results = run_reactants("[N:1]>>[N+:1]", &[&tma]).unwrap();
        assert!(
            !results.is_empty(),
            "N+ with 3 bonds must be valid and kept"
        );
    }

    #[test]
    fn new_atom_in_product() {
        let mol = parse("C").unwrap();
        let results = run_reactants("[C:1]>>[C:1]=O", &[&mol]).unwrap();
        assert!(!results.is_empty());
        let prod = &results[0][0];
        assert_eq!(prod.atom_count(), 2, "C + new O = 2 atoms");
    }

    #[test]
    fn amide_bond_formation() {
        // NH3 + H-C(=O)-Cl → H-C(=O)-NH2 (formamide)
        let nh3 = parse("N").unwrap();
        let hcocl = parse("C(=O)Cl").unwrap();
        let results = run_reactants("[N:1].[C:2](=O)Cl>>[C:2](=O)[N:1]", &[&nh3, &hcocl]).unwrap();
        assert!(!results.is_empty());
        let prod = &results[0][0];
        assert_eq!(prod.atom_count(), 3, "C + O(new) + N = 3 atoms");
    }

    #[test]
    fn double_bond_product() {
        let mol = parse("CC").unwrap();
        let results = run_reactants("[C:1][C:2]>>[C:1]=[C:2]", &[&mol]).unwrap();
        assert!(!results.is_empty());
        let prod = &results[0][0];
        assert_eq!(prod.atom_count(), 2);
        let bond_orders: Vec<BondOrder> = prod.bonds().map(|(_, b)| b.order).collect();
        assert!(
            bond_orders.contains(&BondOrder::Double),
            "product must contain a double bond"
        );
    }

    #[test]
    fn substituent_carry_through() {
        // Methylamine + acetyl chloride → N-methylacetamide (5 heavy atoms)
        // CH3-NH2 + CH3-C(=O)-Cl → CH3-C(=O)-NH-CH3
        let methylamine = parse("NC").unwrap();
        let acetyl_cl = parse("CC(=O)Cl").unwrap();
        let results = run_reactants(
            "[N:1].[C:2](=O)Cl>>[C:2](=O)[N:1]",
            &[&methylamine, &acetyl_cl],
        )
        .unwrap();
        assert!(!results.is_empty(), "must produce at least one product set");
        let prod = &results[0][0];
        assert_eq!(
            prod.atom_count(),
            5,
            "N-methylacetamide has 5 heavy atoms, got {}",
            prod.atom_count()
        );
    }

    #[test]
    fn bfs_no_leakage_into_other_product_template_atoms() {
        // Issue #13: in diethylamine (CCNCC), the SMIRKS [N:1][C:2]>>[N:1].[C:2]
        // should cleave the N-C bond and produce:
        //   product1 [N:1] = N + right ethyl chain  (3 atoms: N, C, C)
        //   product2 [C:2] = left ethyl fragment     (2 atoms: C, C)
        //
        // Before the #13 fix, the BFS for product2 was seeded from BOTH N and C:2,
        // causing the right ethyl chain (atoms beyond N) to leak into product2
        // → product2 would have 4 atoms instead of 2.
        let diethylamine = parse("CCNCC").unwrap(); // C-C-N-C-C, 5 heavy atoms
        let results = run_reactants("[N:1][C:2]>>[N:1].[C:2]", &[&diethylamine]).unwrap();
        assert!(
            !results.is_empty(),
            "should find at least one N-C bond match"
        );

        // Find a result where product2 ([C:2]) has exactly 2 atoms (ethyl fragment)
        // — this is only possible when BFS does NOT leak the other ethyl chain.
        let clean_cleavage = results.iter().find(|ps| {
            ps.len() == 2
                && ((ps[0].atom_count() == 3 && ps[1].atom_count() == 2)
                    || (ps[0].atom_count() == 2 && ps[1].atom_count() == 3))
        });
        assert!(
            clean_cleavage.is_some(),
            "expected at least one product set with sizes {{3, 2}} (N+ethyl, ethyl); \
             all sets: {:?}",
            results
                .iter()
                .map(|ps| ps.iter().map(|p| p.atom_count()).collect::<Vec<_>>())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn single_product_no_leakage_from_other_template_core() {
        // Ethane cleavage: each product should be a single carbon atom.
        let ethane = parse("CC").unwrap();
        let results = run_reactants("[C:1][C:2]>>[C:1].[C:2]", &[&ethane]).unwrap();
        assert!(!results.is_empty());
        for ps in &results {
            assert_eq!(ps.len(), 2, "two product templates → two products");
            assert_eq!(ps[0].atom_count(), 1, "each product is a single carbon");
            assert_eq!(ps[1].atom_count(), 1, "each product is a single carbon");
        }
    }

    // ── Stereo SMIRKS tests ───────────────────────────────────────────────────

    #[test]
    fn stereo_cleared_when_all_neighbors_are_unmapped_template_literals() {
        // Product template [C:1] has no chirality, and all three substituents
        // (F, Cl, Br) are unmapped literal atoms in BOTH templates -- SMIRKS
        // gives no per-atom reactant<->product correspondence for an unmapped
        // atom (only `:n` atom-map numbers establish identity), so there is no
        // derivable stereo_neighbor_order to carry the inherited flag with.
        // Fail-closed: clear rather than keep a flag with an undefined order
        // (chematic-cip already treats order-less chirality as unassigned, so
        // keeping the raw flag here was never actually meaningful downstream).
        let mol = parse("[C@@H](F)(Cl)Br").unwrap();
        let results = run_reactants("[C@@H:1](F)(Cl)Br>>[C:1](F)(Cl)Br", &[&mol]).unwrap();
        assert!(!results.is_empty(), "should match and produce a product");
        let prod = &results[0][0];
        // The core C atom is first in the builder (index 0).
        let core_chirality = prod.atom(AtomIdx(0)).chirality;
        assert_eq!(
            core_chirality,
            Chirality::None,
            "no derivable stereo_neighbor_order for an all-unmapped-substituent \
             inherited flag must clear chirality, not keep an order-less flag"
        );
    }

    #[test]
    fn stereo_preserved_when_all_neighbors_are_mapped_and_remappable() {
        // Every substituent around the stereocenter is atom-mapped, so each
        // has a real, unique reactant->product correspondence via src_to_new.
        // This is the case the remap mechanism *can* and must handle: keep
        // both the flag and a validated, remapped stereo_neighbor_order.
        let mol = parse("[C@@H](F)(Cl)Br").unwrap();
        let results = run_reactants(
            "[C@@H:1]([F:2])([Cl:3])[Br:4]>>[C:1]([F:2])([Cl:3])[Br:4]",
            &[&mol],
        )
        .unwrap();
        assert!(!results.is_empty(), "should match and produce a product");
        let prod = &results[0][0];
        let core_chirality = prod.atom(AtomIdx(0)).chirality;
        assert_eq!(
            core_chirality,
            Chirality::Clockwise,
            "source @@ chirality must be preserved when every neighbor is \
             mapped and remappable, with an identity-order template"
        );
        assert!(
            prod.stereo_neighbor_order(AtomIdx(0)).is_some(),
            "a kept inherited flag must always carry a validated stereo_neighbor_order"
        );
    }

    #[test]
    fn stereo_unmapped_leaving_groups_fully_removed_clears_chirality() {
        // Blocker #1 from review: the reactant's literal unmapped
        // substituents (F, Cl, Br) are removed entirely by the product
        // template (bare [C:1], zero neighbours) -- bug_b_topology_unchanged
        // must not special-case an empty product-side unmapped set as
        // "unchanged" when the reactant side was non-empty. In the current
        // architecture this is enforced twice over: the topology gate itself
        // (fixed here) AND independently by remap_reactant_stereo_order,
        // since a literal/unmapped reactant-template neighbour never has a
        // src_to_new entry regardless of what the product does with it. The
        // observable behavior this test pins is the one the reviewer asked
        // for either way: chirality must clear, not survive.
        let mol = parse("[C@@H](F)(Cl)Br").unwrap();
        let results = run_reactants("[C@@H:1](F)(Cl)Br>>[C:1]", &[&mol]).unwrap();
        assert!(!results.is_empty(), "should match and produce a product");
        let prod = &results[0][0];
        assert_eq!(
            prod.atom(AtomIdx(0)).chirality,
            Chirality::None,
            "deleting all of the stereocenter's unmapped substituents must clear chirality"
        );
    }

    #[test]
    fn stereo_remote_reaction_preserves_stereocenter() {
        // The bond change (Br leaving, N arriving) happens two bonds away
        // from the stereocenter, through a mapped, unchanged carrier chain
        // (F, Cl, C4 all keep the same atom-map numbers and the same
        // relative template order on both sides) -- the stereocenter's own
        // immediate neighbor set and order are completely untouched by the
        // remote reaction. Configuration must survive.
        let mol = parse("[C@@H](F)(Cl)CCBr").unwrap();
        let amine = parse("N").unwrap();
        let results = run_reactants(
            "[C@@H:1]([F:2])([Cl:3])[C:4][C:5][Br:6].[N:7]\
             >>[C:1]([F:2])([Cl:3])[C:4][C:5][N:7]",
            &[&mol, &amine],
        )
        .unwrap();
        assert!(!results.is_empty(), "should match and produce a product");
        let prod = &results[0][0];
        assert_eq!(
            prod.atom(AtomIdx(0)).chirality,
            Chirality::Clockwise,
            "a remote bond change 2 bonds away must not disturb the \
             stereocenter's own unchanged, fully-mapped neighbor set"
        );
        assert!(
            prod.stereo_neighbor_order(AtomIdx(0)).is_some(),
            "preserved chirality must carry a validated stereo_neighbor_order"
        );
    }

    #[test]
    fn stereo_survives_16_plus_atom_map_relabelings_smiles_and_cip_invariant() {
        // The chemistry must not depend on which integers a SMIRKS author
        // picks for `:n` map numbers -- only the correspondence they encode.
        // Re-run the same reaction with 20 distinct map-number assignments
        // (same template text shape and order each time, only the four
        // integers change) and assert canonical SMILES and accurate CIP
        // agree with the first run every time.
        let mol = parse("[C@@H](F)(Cl)Br").unwrap();
        // Central atom is fixed at map `:1`; offset the other three ranges
        // so none of them ever collides with it or each other.
        let labelings: Vec<[u32; 3]> = (2..=21).map(|i| [i, i + 100, i + 200]).collect();
        assert!(labelings.len() >= 16, "need at least 16 relabelings");

        let mut baseline_smiles: Option<String> = None;
        let mut baseline_cip: Option<Option<chematic_core::CipCode>> = None;
        for [b, c, d] in labelings {
            let smirks =
                format!("[C@@H:1]([F:{b}])([Cl:{c}])[Br:{d}]>>[C:1]([F:{b}])([Cl:{c}])[Br:{d}]");
            let results = run_reactants(&smirks, &[&mol]).unwrap();
            assert!(!results.is_empty(), "labeling {b},{c},{d} must still match");
            let prod = &results[0][0];
            let smi = canonical(prod);
            let cip = chematic_chem::assign_cip_with_mode(prod, chematic_chem::CipMode::Accurate)
                .unwrap()
                .get(AtomIdx(0));

            match (&baseline_smiles, &baseline_cip) {
                (None, None) => {
                    baseline_smiles = Some(smi);
                    baseline_cip = Some(cip);
                }
                (Some(base_smi), Some(base_cip)) => {
                    assert_eq!(
                        &smi, base_smi,
                        "canonical SMILES must be invariant to atom-map relabeling \
                         (labeling {b},{c},{d})"
                    );
                    assert_eq!(
                        &cip, base_cip,
                        "accurate CIP code must be invariant to atom-map relabeling \
                         (labeling {b},{c},{d})"
                    );
                }
                _ => unreachable!(),
            }
        }
    }

    #[test]
    fn stereo_chirality_implies_valid_stereo_neighbor_order_invariant() {
        // Cross-cutting invariant the reviewer required: whenever an atom
        // carries a non-None chirality flag, it must always have a valid
        // stereo_neighbor_order alongside it -- never a raw flag with no
        // order (which chematic-cip and the SMILES writer cannot interpret
        // meaningfully). Checked across every atom of a representative
        // spread of already-covered product-generating scenarios.
        let assert_invariant_holds = |mol: &Molecule, label: &str| {
            for i in 0..mol.atom_count() {
                let idx = AtomIdx(i as u32);
                if mol.atom(idx).chirality != Chirality::None {
                    assert!(
                        mol.stereo_neighbor_order(idx).is_some(),
                        "{label}: atom {i} has chirality {:?} but no stereo_neighbor_order",
                        mol.atom(idx).chirality
                    );
                }
            }
        };

        let mol = parse("[C@@H](F)(Cl)Br").unwrap();
        let chain_mol = parse("[C@@H](F)(Cl)CCBr").unwrap();
        let amine = parse("N").unwrap();

        let scenarios: Vec<(&str, Vec<Molecule>)> = vec![
            (
                "identity template order",
                vec![
                    run_reactants("[C@@H:1](F)(Cl)Br>>[C@@H:1](F)(Cl)Br", &[&mol]).unwrap()[0][0]
                        .clone(),
                ],
            ),
            (
                "reordered template inverts",
                vec![
                    run_reactants("[C@@H:1](F)(Cl)Br>>[C@@H:1](Cl)(F)Br", &[&mol]).unwrap()[0][0]
                        .clone(),
                ],
            ),
            (
                "all neighbors mapped and remappable",
                vec![
                    run_reactants(
                        "[C@@H:1]([F:2])([Cl:3])[Br:4]>>[C:1]([F:2])([Cl:3])[Br:4]",
                        &[&mol],
                    )
                    .unwrap()[0][0]
                        .clone(),
                ],
            ),
            (
                "all unmapped substituents deleted",
                vec![run_reactants("[C@@H:1](F)(Cl)Br>>[C:1]", &[&mol]).unwrap()[0][0].clone()],
            ),
            (
                "substituent replacement",
                vec![
                    run_reactants("[C@@H:1](F)(Cl)Br>>[C:1](F)(Cl)I", &[&mol]).unwrap()[0][0]
                        .clone(),
                ],
            ),
            (
                "remote reaction",
                vec![
                    run_reactants(
                        "[C@@H:1]([F:2])([Cl:3])[C:4][C:5][Br:6].[N:7]\
                         >>[C:1]([F:2])([Cl:3])[C:4][C:5][N:7]",
                        &[&chain_mol, &amine],
                    )
                    .unwrap()[0][0]
                        .clone(),
                ],
            ),
        ];

        for (label, products) in &scenarios {
            for prod in products {
                assert_invariant_holds(prod, label);
            }
        }
    }

    #[test]
    fn stereo_inverted_by_template() {
        // Product template [C@H:1] has @ (CounterClockwise) → overrides source @@ (Clockwise).
        let mol = parse("[C@@H](F)(Cl)Br").unwrap();
        let results = run_reactants("[C@@H:1](F)(Cl)Br>>[C@H:1](F)(Cl)Br", &[&mol]).unwrap();
        assert!(!results.is_empty(), "should match and produce a product");
        let prod = &results[0][0];
        let core_chirality = prod.atom(AtomIdx(0)).chirality;
        assert_eq!(
            core_chirality,
            Chirality::CounterClockwise,
            "product template @ must override source @@ → CounterClockwise"
        );
    }

    #[test]
    fn stereo_identity_template_order_preserves_configuration() {
        // Product template repeats the exact same neighbour order as the
        // reactant template (F, Cl, Br) with an explicit @@ -- the simplest
        // "nothing changed" case for the parity-aware correction.
        let mol = parse("[C@@H](F)(Cl)Br").unwrap();
        let results = run_reactants("[C@@H:1](F)(Cl)Br>>[C@@H:1](F)(Cl)Br", &[&mol]).unwrap();
        let prod = &results[0][0];
        let smi = canonical(prod);
        let input_smi = canonical(&mol);
        assert_eq!(
            smi, input_smi,
            "identity template order should reproduce the input unchanged"
        );
    }

    #[test]
    fn stereo_reordered_template_inverts_absolute_configuration() {
        // Issue found while surveying RDKit's open issues (analogous to
        // RDKit #9257): reordering two substituents in the product template
        // while keeping the SAME @@ symbol must invert the absolute
        // configuration -- the symbol describes an order-relative sense,
        // not an absolute one. Cross-checked against a live RDKit oracle
        // (`rdkit.Chem.rdChemReactions`): RDKit's `reordered` output equals
        // its own `inverted` output, both differing from `identity`.
        let mol = parse("[C@@H](F)(Cl)Br").unwrap();
        let identity = canonical(
            &run_reactants("[C@@H:1](F)(Cl)Br>>[C@@H:1](F)(Cl)Br", &[&mol]).unwrap()[0][0],
        );
        let reordered = canonical(
            &run_reactants("[C@@H:1](F)(Cl)Br>>[C@@H:1](Cl)(F)Br", &[&mol]).unwrap()[0][0],
        );
        let explicit_invert = canonical(
            &run_reactants("[C@@H:1](F)(Cl)Br>>[C@H:1](F)(Cl)Br", &[&mol]).unwrap()[0][0],
        );
        assert_ne!(
            reordered, identity,
            "reordering two substituents under the same @@ symbol must change \
             the absolute configuration"
        );
        assert_eq!(
            reordered, explicit_invert,
            "reordering two substituents under @@ must match the explicit-@ \
             (same order, opposite symbol) result -- both express the same \
             inverted configuration"
        );
    }

    #[test]
    fn stereo_substituent_replacement_clears_chirality() {
        // Neighbour SET changes (Br -> I): the old flag can no longer mean
        // anything -- must clear to Chirality::None, not carry a stale @@
        // onto a differently-substituted center. Pre-existing behavior,
        // now driven by correct_product_stereo's Bug-B gate instead of the
        // deleted ad hoc heuristic; must not regress.
        let mol = parse("[C@@H](F)(Cl)Br").unwrap();
        let results = run_reactants("[C@@H:1](F)(Cl)Br>>[C:1](F)(Cl)I", &[&mol]).unwrap();
        let prod = &results[0][0];
        let core_chirality = prod.atom(AtomIdx(0)).chirality;
        assert_eq!(
            core_chirality,
            Chirality::None,
            "substituent replacement (Br->I) must clear chirality, not preserve a stale flag"
        );
    }

    #[test]
    fn stereo_lost_mapped_neighbor_clears_spurious_chirality() {
        // Issue found while surveying RDKit's open issues: a mapped/core
        // neighbour (N) is dropped by the product template, while the
        // template's own C:1 atom carries no explicit chirality. The old
        // invalidation heuristic only ever compared *unmapped* substituent
        // element sets (both empty here, since every neighbour of C:1 is
        // mapped) and missed this case entirely, producing a spurious
        // chirality tag on a carbon with 3 identical implicit hydrogens.
        let mol = parse("[C@H](N)(O)Cl").unwrap();
        let results = run_reactants("[C@H:1]([N:2])([O:3])Cl>>[C:1][O:3]", &[&mol]).unwrap();
        let prod = &results[0][0];
        let smi = canonical(prod);
        assert!(
            !smi.contains('@'),
            "carbon losing its N neighbour is no longer a stereocenter (3 identical \
             implicit H's) -- product must carry no chirality symbol, got {smi}"
        );
    }

    #[test]
    fn stereo_polyol_ring_contraction_no_spurious_ch2oh_chirality() {
        // The originally reported repro (analogous to RDKit #9257), minimized
        // versions of which are the two tests directly above. A pyranose ->
        // furanose ring contraction where one ring carbon (mapped :2) loses
        // its own explicit chirality in the product template and becomes a
        // plain CH2OH (degree 2, bonded to a degree-1 O) -- pre-fix this
        // atom kept a spurious inherited chirality tag ([C@H2], a carbon
        // with 3 identical implicit hydrogens). Does NOT assert a specific
        // sign at the two atoms whose absolute configuration this fix does
        // not by itself resolve either way (mapped :11, unchanged reactant->
        // product template and therefore architecturally untouched by this
        // fix; and :15, the debated center in the original RDKit report) --
        // only that the reaction produces exactly 4 real stereocenters, not
        // 5, and that the correct atom (the one that structurally lost its
        // stereocenter) is the one that lost its tag.
        let smirks = "[O:1][C@H:2]1[O:3][C@H:4]([C:5][C:6])[C@@H:11]([O:12])[C@H:13]([O:14])\
                       [C@H:15]1[O:16]>>[O:1][C:2][C@:15]([O:16])1[O:3][C@H:4]([C:5][C:6])\
                       [C@@H:11]([O:12])[C@H:13]1[O:14]";
        let mol = parse("CC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O").unwrap();
        let results = run_reactants(smirks, &[&mol]).unwrap();
        let prod = &results[0][0];

        let n = prod.atom_count();
        let chiral_atoms: Vec<AtomIdx> = (0..n)
            .map(|i| AtomIdx(i as u32))
            .filter(|&a| prod.atom(a).chirality != Chirality::None)
            .collect();
        assert_eq!(
            chiral_atoms.len(),
            4,
            "expected exactly 4 real stereocenters in the product, got {}: {:?}",
            chiral_atoms.len(),
            chiral_atoms
        );

        // Structurally identify the CH2OH carbon: degree 2, bonded to a
        // degree-1 oxygen -- not by a hardcoded index, since atom order
        // depends on build_product's own construction order.
        let ch2oh = (0..n).map(|i| AtomIdx(i as u32)).find(|&a| {
            let atom = prod.atom(a);
            atom.element == chematic_core::Element::C
                && prod.degree(a) == 2
                && prod.neighbors(a).any(|(nb, _)| {
                    prod.atom(nb).element == chematic_core::Element::O && prod.degree(nb) == 1
                })
        });
        let ch2oh = ch2oh.expect("product must contain a CH2OH carbon");
        assert_eq!(
            prod.atom(ch2oh).chirality,
            Chirality::None,
            "the CH2OH carbon (lost its own chirality in the product template, degree 2, \
             no longer has 4 distinct substituents) must not carry a stereo tag"
        );
    }

    // ── run_reactants_strict tests ────────────────────────────────────────────

    #[test]
    fn strict_mode_excludes_substituents() {
        // Methylamine (NC): in normal mode [N:1]>>[N:1] carries C through as substituent.
        // In strict mode only N is returned (no C).
        let mol = parse("NC").unwrap();
        let normal = run_reactants("[N:1]>>[N:1]", &[&mol]).unwrap();
        let strict = run_reactants_strict("[N:1]>>[N:1]", &[&mol]).unwrap();
        assert!(!normal.is_empty());
        assert!(!strict.is_empty());
        let normal_atoms = normal[0][0].atom_count();
        let strict_atoms = strict[0][0].atom_count();
        assert!(
            normal_atoms > strict_atoms,
            "normal mode carries substituent C (got {normal_atoms}), \
             strict mode only mapped N (got {strict_atoms})"
        );
        assert_eq!(strict_atoms, 1, "strict mode: only the mapped N atom");
    }

    #[test]
    fn strict_mode_bond_cleavage() {
        // Ethane cleavage: strict mode gives 1-atom products, same as normal here
        // (no unmapped substituents on either C).
        let ethane = parse("CC").unwrap();
        let results = run_reactants_strict("[C:1][C:2]>>[C:1].[C:2]", &[&ethane]).unwrap();
        assert!(!results.is_empty());
        for ps in &results {
            assert_eq!(ps[0].atom_count(), 1);
            assert_eq!(ps[1].atom_count(), 1);
        }
    }

    // ── Issue #18: product bracket notation cleanup ───────────────────────────

    #[test]
    fn product_removes_bracket_from_bare_bracket_atoms() {
        // Issue #18: [O:1] in product template (hydrogen_count=Some(0)) must produce
        // clean `O` SMILES, not `[O]`.
        use chematic_smiles::canonical_smiles;
        let mol = parse("OCC").unwrap();
        let results = run_reactants("[OH:1]>>[O:1]", &[&mol]).unwrap();
        assert!(!results.is_empty(), "should match hydroxyl");
        let prod_smi = canonical_smiles(&results[0][0]);
        assert!(
            !prod_smi.contains("[O]"),
            "bare [O:1] product must write as O, not [O], got: {prod_smi}"
        );
    }

    #[test]
    fn product_preserves_explicit_h_from_template() {
        // [NH2:1] in product template specifies 2H explicitly — the built
        // atom's `hydrogen_count` must be `Some(2)`, not silently dropped or
        // recomputed.
        //
        // Checks the data-level field directly rather than the canonical
        // SMILES string's bracket notation (issue #205): after
        // `initial_invariant`/`emit_atom`'s explicit/implicit-H-count
        // unification fix, `canonical_smiles` correctly stops forcing
        // brackets for an atom whose explicit H count merely repeats what
        // organic-subset valence inference would already give -- this atom
        // (N bonded to one carbon) infers to 2 implicit H anyway, so its
        // canonical form is now the fully standard, bracket-free "CN"
        // (methylamine), not "C[NH2]". The template's explicit
        // specification is still honored -- it is what the "2" came from --
        // just no longer forced into visible bracket notation once it's
        // redundant with inference.
        let mol = parse("NC").unwrap();
        let results = run_reactants("[N:1]>>[NH2:1]", &[&mol]).unwrap();
        assert!(!results.is_empty(), "should match amine N");
        let product = &results[0][0];
        let n_atom = product
            .atoms()
            .find(|(_, a)| a.element == chematic_core::Element::N)
            .map(|(_, a)| a)
            .expect("product must contain the templated N atom");
        assert_eq!(
            n_atom.hydrogen_count,
            Some(2),
            "explicit [NH2:1] in product must set hydrogen_count = Some(2)"
        );
    }

    #[test]
    fn reaction_derived_matches_direct_parse_chlorobenzene() {
        // The exact fixture from kent-tokyo/renkin PR #65: a reaction-
        // derived molecule (built via `run_reactants`, whose product-
        // template atom comes from bracket notation `[Cl]` in the SMIRKS)
        // must canonicalize identically to a directly-parsed organic-subset
        // "Clc1ccccc1" of the same compound (issue #205).
        use chematic_smiles::canonical_smiles;
        let fwd = "[c:1][Br]>>[c:1][Cl]";
        let known = parse("Brc1ccccc1").unwrap();
        let results = run_reactants(fwd, &[&known]).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].len(), 1);
        let reaction_derived = &results[0][0];
        let direct = parse("Clc1ccccc1").unwrap();

        // Structural sanity first: same atom count and element multiset,
        // so canonical-string equality below is a real invariance proof,
        // not an assumption that the two are the same molecule.
        assert_eq!(reaction_derived.atom_count(), direct.atom_count());
        let mut els_a: Vec<_> = reaction_derived.atoms().map(|(_, a)| a.element).collect();
        let mut els_b: Vec<_> = direct.atoms().map(|(_, a)| a.element).collect();
        els_a.sort();
        els_b.sort();
        assert_eq!(els_a, els_b, "same element multiset required");

        assert_eq!(
            canonical_smiles(reaction_derived),
            canonical_smiles(&direct),
            "reaction-derived and directly-parsed chlorobenzene must canonicalize identically"
        );
    }

    // ── Issue #20: SMIRKS stereo filtering ───────────────────────────────────

    #[test]
    fn stereo_filter_rejects_wrong_enantiomer() {
        // Issue #20: SMIRKS with @@ reactant template must match @@ but NOT @ reactant.
        let l_ala = parse("N[C@@H](C)C(=O)O").unwrap(); // L-alanine (@@)
        let d_ala = parse("N[C@H](C)C(=O)O").unwrap(); // D-alanine (@)

        let smirks = "[N:1][C@@H:2](C)C(=O)O>>[N:1][C@@H:2](C)C(=O)O";
        let results_l = run_reactants(smirks, &[&l_ala]).unwrap();
        let results_d = run_reactants(smirks, &[&d_ala]).unwrap();

        assert!(
            !results_l.is_empty(),
            "L-alanine (@@) must match @@ template"
        );
        assert!(
            results_d.is_empty(),
            "D-alanine (@) must NOT match @@ template (stereo filter, issue #20)"
        );
    }

    #[test]
    fn stereo_neutral_smirks_matches_both_enantiomers() {
        // SMIRKS without @/@@ must still match both enantiomers (backward compat).
        let l_ala = parse("N[C@@H](C)C(=O)O").unwrap();
        let d_ala = parse("N[C@H](C)C(=O)O").unwrap();
        let smirks = "[N:1][CH:2](C)C(=O)O>>[N:1][CH:2](C)C(=O)O";
        let r_l = run_reactants(smirks, &[&l_ala]).unwrap();
        let r_d = run_reactants(smirks, &[&d_ala]).unwrap();
        assert!(!r_l.is_empty(), "L-alanine must match non-stereo template");
        assert!(!r_d.is_empty(), "D-alanine must match non-stereo template");
    }

    #[test]
    fn pinned_rdkit_profile_refuses_only_reactant_tetrahedral_templates() {
        let chiral =
            PreparedReaction::new("[N:1][C@@H:2](C)C(=O)O>>[N:1][C@@H:2](C)C(=O)O").unwrap();
        assert_eq!(
            chiral.rdkit_2026_03_6_unsupported_reason(),
            Some(ReactionCompatibilityUnsupported::ChiralReactantTemplateSemantics)
        );
        assert_eq!(
            chiral
                .rdkit_2026_03_6_unsupported_reason()
                .unwrap()
                .reason_code(),
            "chiral_reactant_template_semantics"
        );
        let product_only = PreparedReaction::new("[C:1][C:2]>>[C@H:1][C:2]").unwrap();
        assert_eq!(product_only.rdkit_2026_03_6_unsupported_reason(), None);
        let ez = PreparedReaction::new("[C:1]/[C:2]=[C:3]/[C:4]>>[C:1]/[C:2]=[C:3]/[C:4]").unwrap();
        assert_eq!(ez.rdkit_2026_03_6_unsupported_reason(), None);
    }

    #[test]
    fn stereo_filter_same_config_different_write_order() {
        // Parity-aware matching must accept the same absolute configuration
        // regardless of SMILES atom write order (the confirmed bug in raw flag
        // comparison). Both molecules ARE L-alanine.
        //   Form A: N[C@@H](C)C(=O)O  — N written first, stored as Clockwise
        //   Form B: C[C@H](N)C(=O)O   — C_methyl first, stored as CounterClockwise
        // A raw flag comparison would reject Form B against an @@ template.
        let l_form_a = parse("N[C@@H](C)C(=O)O").unwrap();
        let l_form_b = parse("C[C@H](N)C(=O)O").unwrap(); // same absolute config, diff write order
        let d_form = parse("N[C@H](C)C(=O)O").unwrap(); // D-alanine (opposite config)

        let smirks = "[N:1][C@@H:2](C)C(=O)O>>[N:1][C@@H:2](C)C(=O)O";

        let r_a = run_reactants(smirks, &[&l_form_a]).unwrap();
        let r_b = run_reactants(smirks, &[&l_form_b]).unwrap();
        let r_d = run_reactants(smirks, &[&d_form]).unwrap();

        assert!(!r_a.is_empty(), "L-alanine form A (N-first @@) must match");
        assert!(
            !r_b.is_empty(),
            "L-alanine form B (C-first @, same absolute config) must also match \
             — parity-aware comparison required"
        );
        assert!(r_d.is_empty(), "D-alanine must still be rejected");
    }

    // ── RDKit #9339: orphaned stereo bonds cleared from products ─────────────

    #[test]
    fn smirks_reaction_clears_orphaned_stereo_bonds() {
        // (E)-2-butene C/C=C/C has Up/Down bonds adjacent to the C=C double bond.
        // After SMIRKS [C:1]=[C:2]>>[C:1][C:2] the double bond is reduced to a
        // single bond.  The Up/Down single bonds on C0-C1 and C2-C3 are no longer
        // adjacent to ANY double bond → orphaned → must be cleared (RDKit PR #9339).
        let mol = parse("C/C=C/C").unwrap(); // (E)-2-butene
        let results = run_reactants("[C:1]=[C:2]>>[C:1][C:2]", &[&mol]).unwrap();
        assert!(!results.is_empty(), "should produce at least one product");
        for prod_set in &results {
            for prod in prod_set {
                for (_, bond) in prod.bonds() {
                    assert_ne!(
                        bond.order,
                        BondOrder::Up,
                        "stray Up bond in product after C=C→C-C (RDKit #9339)"
                    );
                    assert_ne!(
                        bond.order,
                        BondOrder::Down,
                        "stray Down bond in product after C=C→C-C (RDKit #9339)"
                    );
                }
            }
        }
    }

    #[test]
    fn smirks_preserves_stereo_bonds_adjacent_to_remaining_double() {
        // If the double bond is kept unchanged, the *exact* E/Z geometry must be
        // preserved — not merely "some Up/Down bond survives" (which passed even
        // while the geometry was flipping E↔Z, issue #50). Verify by comparing the
        // product's canonical SMILES to the canonical SMILES of the known input.
        use chematic_smiles::canonical_smiles;
        for input in ["C/C=C/C", "C/C=C\\C"] {
            let mol = parse(input).unwrap();
            let results = run_reactants("[C:1]=[C:2]>>[C:1]=[C:2]", &[&mol]).unwrap();
            assert!(!results.is_empty());
            let expected = canonical_smiles(&mol);
            let got = canonical_smiles(&results[0][0]);
            assert_eq!(
                got, expected,
                "identity SMIRKS must preserve exact E/Z geometry for {input}"
            );
        }
    }

    // -----------------------------------------------------------------------
    // E/Z double-bond stereo filtering (issue #21)
    // -----------------------------------------------------------------------

    #[test]
    fn ez_stereo_e_template_matches_e_alkene() {
        // Template specifies E: [C:1]/[C:2]=[C:3]/[C:4]
        // E-2-butene (C/C=C/C) should produce 1 result.
        let e_alkene = parse("C/C=C/C").unwrap();
        let smirks = "[C:1]/[C:2]=[C:3]/[C:4]>>[C:1][C:2][C:3][C:4]";
        let results = run_reactants(smirks, &[&e_alkene]).unwrap();
        assert!(!results.is_empty(), "E-template must match E-alkene");
    }

    #[test]
    fn ez_stereo_e_template_rejects_z_alkene() {
        // Template specifies E: [C:1]/[C:2]=[C:3]/[C:4]
        // Z-2-butene (C/C=C\C) should produce 0 results.
        let z_alkene = parse("C/C=C\\C").unwrap();
        let smirks = "[C:1]/[C:2]=[C:3]/[C:4]>>[C:1][C:2][C:3][C:4]";
        let results = run_reactants(smirks, &[&z_alkene]).unwrap();
        assert!(results.is_empty(), "E-template must reject Z-alkene");
    }

    #[test]
    fn ez_stereo_neutral_template_matches_both_geometries() {
        // Template without stereo: [C:1][C:2]=[C:3][C:4]>>[C:1]
        // Both E and Z alkenes should match.
        let e_alkene = parse("C/C=C/C").unwrap();
        let z_alkene = parse("C/C=C\\C").unwrap();
        let smirks = "[C:1][C:2]=[C:3][C:4]>>[C:1]";
        assert!(
            !run_reactants(smirks, &[&e_alkene]).unwrap().is_empty(),
            "neutral template must match E-alkene"
        );
        assert!(
            !run_reactants(smirks, &[&z_alkene]).unwrap().is_empty(),
            "neutral template must match Z-alkene"
        );
    }

    #[test]
    fn ez_stereo_one_sided_template_matches_both_geometries() {
        // Single-sided stereo bond in template: [C:1]/[C:2]=[C:3][C:4]
        // Without both sides specified, E/Z is ambiguous → no filtering.
        let e_alkene = parse("C/C=C/C").unwrap();
        let z_alkene = parse("C/C=C\\C").unwrap();
        let smirks = "[C:1]/[C:2]=[C:3][C:4]>>[C:1]";
        assert!(
            !run_reactants(smirks, &[&e_alkene]).unwrap().is_empty(),
            "one-sided template must match E-alkene"
        );
        assert!(
            !run_reactants(smirks, &[&z_alkene]).unwrap().is_empty(),
            "one-sided template must match Z-alkene"
        );
    }

    #[test]
    fn ez_stereo_retro_wittig_z_matches_z_hexene() {
        // Retro-Wittig (Z-alkene → two carbonyls).
        // SMIRKS: [C:1]/[C:2]=[C:3]\[C:4]>>[C:1][C:2]=O.[O:3]=[C:4]
        //   reads: C:2 and C:3 on OPPOSITE sides (E/trans for those two)
        //   but the substituents C:1 and C:4 are on the SAME side (Z-selectivity)
        //
        // Z-3-hexene (CC/C=C\CC) should match; E-3-hexene (CC/C=C/CC) should not.
        let z_hexene = parse("CC/C=C\\CC").unwrap();
        let e_hexene = parse("CC/C=C/CC").unwrap();
        let smirks = "[C:1]/[C:2]=[C:3]\\[C:4]>>[C:1][C:2]=O.[O:3]=[C:4]";
        assert!(
            !run_reactants(smirks, &[&z_hexene]).unwrap().is_empty(),
            "Z-template must match Z-3-hexene"
        );
        assert!(
            run_reactants(smirks, &[&e_hexene]).unwrap().is_empty(),
            "Z-template must reject E-3-hexene"
        );
    }

    #[test]
    fn ez_stereo_z_template_matches_z_alkene() {
        // Template specifies Z: [C:1]/[C:2]=[C:3]\[C:4]
        // Z-2-butene (C/C=C\C) should match.
        let z_alkene = parse("C/C=C\\C").unwrap();
        let e_alkene = parse("C/C=C/C").unwrap();
        let smirks = "[C:1]/[C:2]=[C:3]\\[C:4]>>[C:1][C:2][C:3][C:4]";
        assert!(
            !run_reactants(smirks, &[&z_alkene]).unwrap().is_empty(),
            "Z-template must match Z-alkene"
        );
        assert!(
            run_reactants(smirks, &[&e_alkene]).unwrap().is_empty(),
            "Z-template must reject E-alkene"
        );
    }

    // -----------------------------------------------------------------------
    // E/Z stereo transfer & creation in products (issue #50)
    //
    // Geometry is verified by comparing the product's canonical SMILES to the
    // canonical SMILES of a reference molecule of known E/Z — exact, not "some
    // Up/Down survives". (CipCode-based verification lives in the Python tests;
    // chematic-rxn cannot depend on chematic-chem without a dependency cycle.)
    // -----------------------------------------------------------------------

    /// Canonical SMILES of the single product of `smirks` applied to `inputs`.
    fn product_canon(smirks: &str, inputs: &[&str]) -> String {
        use chematic_smiles::canonical_smiles;
        let mols: Vec<Molecule> = inputs.iter().map(|s| parse(s).unwrap()).collect();
        let refs: Vec<&Molecule> = mols.iter().collect();
        let results = run_reactants(smirks, &refs).unwrap();
        assert!(!results.is_empty(), "no product for {smirks} on {inputs:?}");
        canonical_smiles(&results[0][0])
    }

    fn canon(smiles: &str) -> String {
        chematic_smiles::canonical_smiles(&parse(smiles).unwrap())
    }

    /// Writer-invariant E/Z of the first C=C double bond: `Some(true)` = E,
    /// `Some(false)` = Z, `None` = no specified geometry. Reuses the crate's
    /// own `ez_stereo_outward` (same convention as the #21 filter): equal
    /// outward directions → Z, opposite → E.
    fn double_bond_is_e(smiles: &str) -> Option<bool> {
        let mol = parse(smiles).unwrap();
        let (a1, a2) = mol
            .bonds()
            .find(|(_, b)| b.order == BondOrder::Double)
            .map(|(_, b)| (b.atom1, b.atom2))?;
        let sa = ez_stereo_outward(&mol, a1, a2)?;
        let sb = ez_stereo_outward(&mol, a2, a1)?;
        Some(sa != sb)
    }

    #[test]
    fn issue50_transfer_identity_preserves_e() {
        // Identity SMIRKS on E-2-butene must yield an E product (was Z before Fix A).
        assert_eq!(
            product_canon("[C:1]=[C:2]>>[C:1]=[C:2]", &["C/C=C/C"]),
            canon("C/C=C/C"),
        );
    }

    #[test]
    fn issue50_transfer_identity_preserves_z() {
        assert_eq!(
            product_canon("[C:1]=[C:2]>>[C:1]=[C:2]", &["C/C=C\\C"]),
            canon("C/C=C\\C"),
        );
    }

    #[test]
    fn issue50_create_e_from_template() {
        // Product template introduces an E double bond from a saturated chain.
        assert_eq!(
            product_canon("[C:1][C:2][C:3][C:4]>>[C:1]/[C:2]=[C:3]/[C:4]", &["CCCC"]),
            canon("C/C=C/C"),
        );
    }

    #[test]
    fn issue50_create_z_from_template() {
        assert_eq!(
            product_canon("[C:1][C:2][C:3][C:4]>>[C:1]/[C:2]=[C:3]\\[C:4]", &["CCCC"]),
            canon("C/C=C\\C"),
        );
    }

    #[test]
    fn issue50_transfer_remote_reaction_keeps_e() {
        // Reaction at a remote site (aldehyde→alcohol) must not disturb the
        // E geometry of a carried-through alkene. The canonical writer may pick
        // `/C=C/` or `\C=C\` (both E) depending on traversal, so assert geometry
        // directly rather than the exact string.
        let got = product_canon("[CH:1]=O>>[C:1]O", &["CC/C=C/CC=O"]);
        assert_eq!(
            double_bond_is_e(&got),
            Some(true),
            "E geometry must survive a remote edit"
        );
        // And a Z input stays Z.
        let got_z = product_canon("[CH:1]=O>>[C:1]O", &["CC/C=C\\CC=O"]);
        assert_eq!(
            double_bond_is_e(&got_z),
            Some(false),
            "Z geometry must survive a remote edit"
        );
    }

    #[test]
    fn issue50_geometry_is_deterministic() {
        // The pre-fix bug was nondeterministic (FxHashMap iteration order).
        // The same transform must give the same geometry on every run.
        let first = product_canon("[C:1]=[C:2]>>[C:1]=[C:2]", &["CC/C=C/CC"]);
        for _ in 0..6 {
            assert_eq!(
                product_canon("[C:1]=[C:2]>>[C:1]=[C:2]", &["CC/C=C/CC"]),
                first,
                "product geometry must be deterministic across runs"
            );
        }
    }

    // -----------------------------------------------------------------------
    // Reaction-transform performance regression witnesses (see
    // docs/rfcs/reaction_transform_perf.md). Root cause: `chematic-smiles`'s
    // `canonical_smiles()` wrote the winning individualize-refine branch's
    // string, threw it away, and had `winning_individualized_ranks`'s caller
    // write it a *second* time -- one fully redundant DFS-and-format pass on
    // every single call, tied or not. Fixed by returning the already-written
    // string instead of recomputing it. These three cases mirror the ones
    // used to characterize and fix the regression:
    // (a) a highly symmetric molecule (many individualize-refine branches --
    //     this is the case that actually reproduces a large, measured slowdown
    //     between chematic 0.4.25 and 0.4.30, NOT run_reactants match/product
    //     volume, which stayed flat across versions);
    // (b) an E/Z stereo control, since the fix touches the same
    //     canonical-writer code path `resolve_ez_markers` depends on;
    // (c) a negative control (asymmetric, no ties) that should show only the
    //     universal (small, single-redundant-write) improvement, not the
    //     symmetric-molecule-specific one.
    // -----------------------------------------------------------------------

    #[test]
    fn perf_witness_a_symmetric_molecule_product_is_correct() {
        // Positive witness: adamantane (Td cage symmetry, 24
        // individualize-refine branches at time of writing) run through a
        // simple ring-opening SMIRKS. The fix must not change *which* string
        // wins -- only how many times it gets written -- so the product's
        // canonical SMILES must still round-trip to the same structure.
        let mol = parse("C1C2CC3CC1CC(C2)C3").unwrap(); // adamantane
        let results = run_reactants("[C:1][C:2]>>[C:1][C:2]", &[&mol]).unwrap();
        assert!(!results.is_empty(), "expected at least one C-C bond match");
        let canon = chematic_smiles::canonical_smiles(&results[0][0]);
        // Adamantane's own canonical form is a fixed point of this
        // identity-shaped SMIRKS: it must reparse to the exact same molecule
        // (same atom/bond count -- the transform doesn't add/remove atoms).
        let reparsed = parse(&canon).unwrap();
        assert_eq!(reparsed.atom_count(), mol.atom_count());
        assert_eq!(reparsed.bond_count(), mol.bond_count());
    }

    #[test]
    fn perf_witness_b_ez_stereo_control_survives_symmetric_fix() {
        // Stereo control: identity/remote transforms on the E/Z pair used in
        // the perf investigation's own fixture
        // (crates/chematic-rxn/fixtures/witness_molecules.smi) must still
        // preserve exact geometry -- this is the non-negotiable issue #50
        // gate, re-run here against the specific molecules this perf fix
        // touched.
        assert_eq!(
            product_canon("[C:1]=[C:2]>>[C:1]=[C:2]", &["CC/C=C/CC(=O)O"]),
            canon("CC/C=C/CC(=O)O"),
            "(E)-hex-3-enoic acid must keep its E geometry"
        );
        assert_eq!(
            product_canon("[C:1]=[C:2]>>[C:1]=[C:2]", &["CC/C=C\\CC(=O)O"]),
            canon("CC/C=C\\CC(=O)O"),
            "(Z)-hex-3-enoic acid must keep its Z geometry"
        );
    }

    #[test]
    fn perf_witness_c_negative_control_asymmetric_molecule() {
        // Negative control: aspirin has no automorphism ties (every ring
        // carbon has a distinct substitution environment), so
        // `winning_individualized_ranks` takes exactly one branch. This
        // exercises the same code path as (a) but should show only the
        // universal single-redundant-write saving, not a
        // symmetric-molecule-specific one -- included so a future reader can
        // tell the two effects apart empirically, not just by reasoning.
        let mol = parse("CC(=O)OC1=CC=CC=C1C(=O)O").unwrap(); // aspirin
        let results = run_reactants("[OH:1]-[C:2]=[O:3]>>C-[O:1]-[C:2]=[O:3]", &[&mol]).unwrap();
        assert!(!results.is_empty(), "expected the carboxylic acid to match");
        let canon = chematic_smiles::canonical_smiles(&results[0][0]);
        let reparsed = parse(&canon).unwrap();
        assert_eq!(reparsed.atom_count(), mol.atom_count() + 1); // +1 methyl carbon
    }

    // -------------------------------------------------------------------
    // Match-level reaction application (issue #225)
    // -------------------------------------------------------------------

    /// `run_reactants(smirks, reactants)` must equal
    /// `find_reaction_matches(...).filter_map(|m| apply_reaction_match(...))`
    /// -- the exact equivalence issue #225's proposed API is built on.
    /// Checked against several existing SMIRKS/reactant pairs already used
    /// elsewhere in this file, not just one.
    #[test]
    fn find_and_apply_match_equals_run_reactants() {
        let mut cases: Vec<(&str, Vec<chematic_core::Molecule>)> = vec![
            ("[N:1]>>[N:1]", vec![parse("NCCN").unwrap()]),
            (
                "[N:1].[C:2]>>[N:1][C:2]",
                vec![parse("N").unwrap(), parse("C").unwrap()],
            ),
            ("[C:1][C:2]>>[C:1].[C:2]", vec![parse("CC").unwrap()]),
            (
                "[OH:1]-[C:2]=[O:3]>>C-[O:1]-[C:2]=[O:3]",
                vec![parse("CC(=O)OC1=CC=CC=C1C(=O)O").unwrap()],
            ),
        ];

        for (_, mols) in &mut cases {
            for mol in mols {
                for i in 0..mol.atom_count() {
                    // Intentionally reused across reactants: tags are labels,
                    // while atom_sources disambiguates reactant identity.
                    mol.set_tag(AtomIdx(i as u32), Some(i as u16 + 1));
                }
            }
        }
        for (smirks, mols) in &cases {
            let reactants: Vec<&Molecule> = mols.iter().collect();
            let direct = run_reactants(smirks, &reactants).unwrap();

            let matches = find_reaction_matches(smirks, &reactants).unwrap();
            let via_matches: Vec<Vec<Molecule>> = matches
                .iter()
                .filter_map(|m| apply_reaction_match(smirks, &reactants, m, true).unwrap())
                .collect();

            assert_eq!(
                direct.len(),
                via_matches.len(),
                "{smirks}: product-set count must match"
            );
            for (d, v) in direct.iter().zip(via_matches.iter()) {
                assert_eq!(
                    d.len(),
                    v.len(),
                    "{smirks}: product count per set must match"
                );
                for (dp, vp) in d.iter().zip(v.iter()) {
                    for (idx, _) in dp.atoms() {
                        assert_eq!(dp.atom_tag(idx), vp.atom_tag(idx), "{smirks}: product tags");
                    }
                    assert_eq!(
                        chematic_smiles::canonical_smiles(dp),
                        chematic_smiles::canonical_smiles(vp),
                        "{smirks}: product molecule must match run_reactants exactly"
                    );
                }
            }
        }
    }

    /// Issue #650: the traced application returns the same molecules as the
    /// untraced one, and every traced product atom is a copy of the reactant
    /// atom it names; mapped template atoms trace to their matched atoms and
    /// template-created atoms trace to `None`.
    #[test]
    fn traced_apply_match_reports_atom_sources() {
        let mut cases: Vec<(&str, Vec<chematic_core::Molecule>)> = vec![
            ("[N:1]>>[N:1]", vec![parse("NCCN").unwrap()]),
            (
                "[N:1].[C:2]>>[N:1][C:2]",
                vec![parse("NCC").unwrap(), parse("CO").unwrap()],
            ),
            ("[C:1][C:2]>>[C:1].[C:2]", vec![parse("CCO").unwrap()]),
            // A product-only map labels a newly created atom; it has no
            // reactant origin, including when emitted as a second product.
            ("[C:1]>>[C:1][O:2]", vec![parse("C").unwrap()]),
            ("[C:1]>>[C:1].[O:2]", vec![parse("C").unwrap()]),
            (
                "[OH:1]-[C:2]=[O:3]>>C-[O:1]-[C:2]=[O:3]",
                vec![parse("CC(=O)OC1=CC=CC=C1C(=O)O").unwrap()],
            ),
            (
                "[C:1](=[O:2])[OH:3].[OH:4][C:5]>>[C:1](=[O:2])[O:4][C:5]",
                vec![parse("CC(=O)O").unwrap(), parse("OCC").unwrap()],
            ),
        ];
        for (_, mols) in &mut cases {
            for mol in mols {
                for i in 0..mol.atom_count() {
                    mol.set_tag(AtomIdx(i as u32), Some(i as u16 + 1));
                }
            }
        }
        for (smirks, mols) in &cases {
            let reactants: Vec<&Molecule> = mols.iter().collect();
            let prepared = PreparedReaction::new(smirks).unwrap();
            let matches = find_reaction_matches(smirks, &reactants).unwrap();
            assert!(!matches.is_empty(), "{smirks}: expected a match");
            for m in &matches {
                for carry in [true, false] {
                    let plain = prepared.apply_match(&reactants, m, carry).unwrap();
                    let traced = prepared.apply_match_traced(&reactants, m, carry).unwrap();
                    let free = apply_reaction_match_traced(smirks, &reactants, m, carry).unwrap();
                    let (Some(plain), Some(traced), Some(free)) = (plain, traced, free) else {
                        continue;
                    };
                    assert_eq!(plain.len(), traced.len());
                    assert_eq!(traced.len(), free.len());
                    let global = m.atom_map_positions(smirks).unwrap();
                    for (((p, t), f), template) in plain
                        .iter()
                        .zip(&traced)
                        .zip(&free)
                        .zip(&prepared.rxn.products)
                    {
                        assert_eq!(
                            chematic_smiles::write(p),
                            chematic_smiles::write(&t.molecule),
                            "{smirks}: traced molecule differs"
                        );
                        assert_eq!(t.atom_sources, f.atom_sources);
                        assert_eq!(t.atom_sources.len(), t.molecule.atom_count());
                        for (i, src) in t.atom_sources.iter().enumerate() {
                            let expected_tag =
                                src.and_then(|src| reactants[src.reactant].atom_tag(src.atom));
                            for product in [p, &t.molecule, &f.molecule] {
                                assert_eq!(
                                    product.atom_tag(AtomIdx(i as u32)),
                                    expected_tag,
                                    "{smirks}: tag at product atom {i}, carry={carry}"
                                );
                            }
                            if let Some(src) = src {
                                assert_eq!(
                                    t.molecule.atom(AtomIdx(i as u32)).element,
                                    reactants[src.reactant].atom(src.atom).element,
                                    "{smirks}: atom {i} element"
                                );
                            }
                        }
                        // Sources are injective.
                        let mut seen = FxHashSet::default();
                        for src in t.atom_sources.iter().flatten() {
                            assert!(seen.insert(*src), "{smirks}: duplicate source {src:?}");
                        }
                        // Product template atom i is product atom i: a mapped
                        // one traces to where the match put that map number,
                        // an unmapped one is new.
                        for i in 0..template.atom_count() {
                            let expected = template
                                .atom(AtomIdx(i as u32))
                                .atom_map
                                .and_then(|am| global.get(&am))
                                .map(|&(reactant, atom)| ReactantAtom { reactant, atom });
                            assert_eq!(t.atom_sources[i], expected, "{smirks}: template atom {i}");
                        }
                    }
                }
            }
        }
        // Template-created atom: the methyl C in the ester SMIRKS is new.
        let mol = parse("OC(=O)c1ccccc1").unwrap();
        let smirks = "[OH:1]-[C:2]=[O:3]>>C-[O:1]-[C:2]=[O:3]";
        let m = &find_reaction_matches(smirks, &[&mol]).unwrap()[0];
        let product = &apply_reaction_match_traced(smirks, &[&mol], m, true)
            .unwrap()
            .unwrap()[0];
        assert_eq!(
            product.atom_sources.iter().filter(|s| s.is_none()).count(),
            1
        );
        assert_eq!(
            product.atom_sources.iter().flatten().count(),
            mol.atom_count(),
            "every reactant atom is carried"
        );
    }

    #[test]
    fn traced_report_preserves_product_order_and_diagnostics() {
        let cases = [
            ("[C:1]>>[C:1].[O:2]", vec![parse("C").unwrap()]),
            ("[#6:1]>>[#6:1]", vec![parse("CC").unwrap()]),
            (
                "[N:1].[C:2]>>[N:1][C:2]",
                vec![parse("N").unwrap(), parse("C").unwrap()],
            ),
            (
                "[13CH3:1][O:2]>>[13CH3:1].[O:2]",
                vec![parse("[13CH3]O").unwrap()],
            ),
        ];
        for (smirks, molecules) in cases {
            let refs = molecules.iter().collect::<Vec<_>>();
            let prepared = PreparedReaction::new(smirks).unwrap();
            let limits = ReactionTransformLimits::default();
            let plain = prepared
                .run_reactants_with_diagnostics(&refs, &limits)
                .unwrap();
            let traced = prepared
                .run_reactants_traced_with_diagnostics(&refs, &limits)
                .unwrap();
            assert_eq!(plain.diagnostics, traced.diagnostics, "{smirks}");
            assert_eq!(plain.products.len(), traced.products.len(), "{smirks}");
            for (plain_set, traced_set) in plain.products.iter().zip(&traced.products) {
                assert_eq!(plain_set.len(), traced_set.len(), "{smirks}");
                for (plain_product, traced_product) in plain_set.iter().zip(traced_set) {
                    assert_eq!(
                        chematic_smiles::write(plain_product),
                        chematic_smiles::write(&traced_product.molecule),
                        "{smirks}"
                    );
                    assert_eq!(
                        traced_product.atom_sources.len(),
                        traced_product.molecule.atom_count(),
                        "{smirks}"
                    );
                    assert_eq!(
                        traced_product.template_maps.len(),
                        traced_product.molecule.atom_count(),
                        "{smirks}"
                    );
                }
            }
            if smirks == "[C:1]>>[C:1].[O:2]" {
                assert_eq!(traced.products[0][0].template_maps, vec![Some(1)]);
                assert_eq!(traced.products[0][1].template_maps, vec![Some(2)]);
                assert_eq!(traced.products[0][1].atom_sources, vec![None]);
            }
        }
    }

    #[test]
    fn symmetric_ether_cleavage_retains_both_source_assignments() {
        let mut ether = parse("COC").unwrap();
        for i in 0..ether.atom_count() {
            ether.set_tag(AtomIdx(i as u32), Some(i as u16 + 1));
        }
        let smirks = "[C:1][O:2][C:3]>>[C:1][O:2].[C:3]";
        let products = run_reactants(smirks, &[&ether]).unwrap();
        let mut assignments: Vec<_> = products
            .iter()
            .map(|set| {
                assert_eq!(set.len(), 2);
                (
                    set[0].atom_tag(AtomIdx(0)).unwrap().get(),
                    set[1].atom_tag(AtomIdx(0)).unwrap().get(),
                )
            })
            .collect();
        assignments.sort_unstable();
        assert_eq!(assignments, [(1, 3), (3, 1)]);
    }

    /// `run_reactants_strict` (carry_substituents=false) must also compose
    /// the same way as the `carry_substituents=true` case above.
    #[test]
    fn find_and_apply_match_equals_run_reactants_strict() {
        let mol = parse("CC(=O)OC1=CC=CC=C1C(=O)O").unwrap();
        let reactants = [&mol];
        let smirks = "[OH:1]-[C:2]=[O:3]>>C-[O:1]-[C:2]=[O:3]";

        let direct = run_reactants_strict(smirks, &reactants).unwrap();
        let matches = find_reaction_matches(smirks, &reactants).unwrap();
        let via_matches: Vec<Vec<Molecule>> = matches
            .iter()
            .filter_map(|m| apply_reaction_match(smirks, &reactants, m, false).unwrap())
            .collect();

        assert_eq!(direct.len(), via_matches.len());
        for (d, v) in direct.iter().zip(via_matches.iter()) {
            for (dp, vp) in d.iter().zip(v.iter()) {
                assert_eq!(
                    chematic_smiles::canonical_smiles(dp),
                    chematic_smiles::canonical_smiles(vp)
                );
            }
        }
    }

    /// The core motivating use case from issue #225: enumerate matches
    /// independently of applying them, reject some based on a property of
    /// the match itself, and apply only the accepted ones -- without
    /// discarding the legitimate matches along with the rejected one.
    #[test]
    fn selective_match_application() {
        let mol = parse("NCCN").unwrap();
        let reactants = [&mol];
        let smirks = "[N:1]>>[N:1]";

        let matches = find_reaction_matches(smirks, &reactants).unwrap();
        assert_eq!(matches.len(), 2, "two N atoms in NCCN → two matches");

        // Reject the match touching the higher-numbered atom (arbitrary
        // match-specific policy standing in for RENKIN's real ring-bond
        // rejection rule), keep the other.
        let positions: Vec<_> = matches
            .iter()
            .map(|m| m.atom_map_positions(smirks).unwrap()[&1])
            .collect();
        let keep_idx = if positions[0].1.0 < positions[1].1.0 {
            0
        } else {
            1
        };

        let applied = apply_reaction_match(smirks, &reactants, &matches[keep_idx], true).unwrap();
        assert!(applied.is_some(), "the accepted match must still apply");

        // Applying only one match must yield exactly one of the two
        // products `run_reactants` would have returned for the whole call,
        // not both and not neither.
        let full = run_reactants(smirks, &reactants).unwrap();
        assert_eq!(full.len(), 2);
        let applied_canon = chematic_smiles::canonical_smiles(&applied.unwrap()[0]);
        assert!(
            full.iter()
                .any(|ps| chematic_smiles::canonical_smiles(&ps[0]) == applied_canon),
            "the selectively-applied product must be one of run_reactants's own outputs"
        );
    }

    /// `ReactionMatch::atom_map_positions` must resolve atom_map:1 to the
    /// actual matched N atom, for each of the two NCCN matches separately.
    #[test]
    fn atom_map_positions_resolves_matched_atom() {
        let mol = parse("NCCN").unwrap();
        let reactants = [&mol];
        let smirks = "[N:1]>>[N:1]";

        let matches = find_reaction_matches(smirks, &reactants).unwrap();
        let mut matched_atoms: Vec<AtomIdx> = matches
            .iter()
            .map(|m| {
                let positions = m.atom_map_positions(smirks).unwrap();
                let (reactant_slot, atom_idx) = positions[&1];
                assert_eq!(reactant_slot, 0, "single-reactant SMIRKS: slot must be 0");
                assert_eq!(mol.atom(atom_idx).element.symbol(), "N");
                atom_idx
            })
            .collect();
        matched_atoms.sort_by_key(|a| a.0);
        assert_eq!(
            matched_atoms,
            vec![AtomIdx(0), AtomIdx(3)],
            "NCCN's two N atoms are at index 0 and 3"
        );
    }

    /// `apply_reaction_match` must return `Ok(None)` -- not an error and not
    /// a product -- for a match whose product set fails the existing
    /// valence filter, matching the case [`run_reactants`] silently drops
    /// (`overvalent_product_filtered_oxygen` above).
    #[test]
    fn apply_reaction_match_none_on_valence_violation() {
        let ethanol = parse("CCO").unwrap();
        let reactants = [&ethanol];
        let smirks = "[O:1]>>[O:1](C)C";

        let matches = find_reaction_matches(smirks, &reactants).unwrap();
        assert_eq!(matches.len(), 1, "exactly one O in ethanol");

        let applied = apply_reaction_match(smirks, &reactants, &matches[0], true).unwrap();
        assert!(
            applied.is_none(),
            "over-valenced product must come back as Ok(None), not Some(..) or Err(..)"
        );
    }

    #[test]
    fn diagnostics_classify_valence_filtered_match_without_changing_products() {
        let ethanol = parse("CCO").unwrap();
        let smirks = "[O:1]>>[O:1](C)C";
        let report = run_reactants_with_diagnostics(
            smirks,
            &[&ethanol],
            &ReactionTransformLimits::default(),
        )
        .unwrap();
        assert_eq!(report.products.len(), 0);
        assert_eq!(report.diagnostics.accepted_matches, 1);
        assert_eq!(report.diagnostics.applied_products, 0);
        assert_eq!(report.diagnostics.valence_rejected_matches, 1);
        assert!(!report.diagnostics.truncated_matches);
        assert!(run_reactants(smirks, &[&ethanol]).unwrap().is_empty());
    }

    #[test]
    fn mapped_isotope_template_matches_only_the_specified_nuclide() {
        let labeled = parse("[13CH3]O").unwrap();
        let unlabeled = parse("CO").unwrap();
        let smirks = "[13CH3:1][O:2]>>[13CH3:1][O:2]";

        let products = run_reactants(smirks, &[&labeled]).unwrap();
        assert_eq!(products.len(), 1);
        assert_eq!(products[0].len(), 1);
        assert!(
            products[0]
                .iter()
                .any(|m| m.atoms().any(|(_, a)| a.isotope == Some(13)))
        );
        assert!(run_reactants(smirks, &[&unlabeled]).unwrap().is_empty());
    }

    #[test]
    fn diagnostics_match_ring_aware_products_and_counts() {
        let mol = parse("C1CCCCC1").unwrap();
        let rings = chematic_perception::find_sssr(&mol);
        let prepared = PreparedReaction::new("[C:1]>>[C:1]").unwrap();
        let report = prepared
            .run_reactants_with_rings_and_limits_with_diagnostics(
                &[&mol],
                &[&rings],
                &ReactionTransformLimits::default(),
            )
            .unwrap();
        assert_eq!(
            report.diagnostics.accepted_matches,
            report.diagnostics.applied_products
        );
        assert_eq!(report.diagnostics.valence_rejected_matches, 0);
        assert_eq!(report.products.len(), report.diagnostics.applied_products);
    }

    #[test]
    fn variant_diagnostics_match_ring_aware_application() {
        let mol = parse("C1CCCCC1").unwrap();
        let rings = chematic_perception::find_sssr(&mol);
        let prepared = PreparedReaction::new("[#6:1]>>[#6:1]").unwrap();
        let reports = prepared
            .run_reactants_with_rings_and_limits_with_variant_diagnostics(
                &[&mol],
                &[&rings],
                &ReactionTransformLimits::default(),
            )
            .unwrap();
        assert_eq!(reports.len(), 2);
        assert_eq!(reports[0].variant_index, 0);
        assert!(reports[0].diagnostics.accepted_matches > 0);
        assert_eq!(
            reports
                .iter()
                .map(|report| report.diagnostics.applied_products)
                .sum::<usize>(),
            prepared
                .run_reactants_with_rings_and_limits(
                    &[&mol],
                    &[&rings],
                    &ReactionTransformLimits::default(),
                )
                .unwrap()
                .len()
        );
    }

    /// `find_reaction_matches` and `apply_reaction_match` must propagate
    /// the same `ReactantCountMismatch` error `run_reactants` does when the
    /// number of input molecules doesn't match the SMIRKS's reactant-slot
    /// count.
    #[test]
    fn find_and_apply_match_reactant_count_mismatch_errors() {
        let mol = parse("C").unwrap();
        let smirks = "[N:1].[C:2]>>[N:1][C:2]";

        let find_err = find_reaction_matches(smirks, &[&mol]);
        assert!(matches!(
            find_err,
            Err(TransformError::ReactantCountMismatch {
                expected: 2,
                got: 1
            })
        ));

        let dummy_match = ReactionMatch {
            per_reactant: vec![FxHashMap::default()],
        };
        let apply_err = apply_reaction_match(smirks, &[&mol], &dummy_match, true);
        assert!(matches!(
            apply_err,
            Err(TransformError::ReactantCountMismatch {
                expected: 2,
                got: 1
            })
        ));
    }

    /// `apply_reaction_match` must also reject a `ReactionMatch` whose own
    /// `per_reactant` shape doesn't match the SMIRKS being applied against
    /// (e.g. a match obtained from a different SMIRKS), rather than
    /// panicking on an out-of-bounds index into `template_atom_maps`.
    #[test]
    fn apply_reaction_match_rejects_mismatched_match_shape() {
        let n_mol = parse("N").unwrap();
        let c_mol = parse("C").unwrap();
        let reactants = [&n_mol, &c_mol];
        let smirks = "[N:1].[C:2]>>[N:1][C:2]";

        // A match shaped for a single-reactant SMIRKS, applied against a
        // two-reactant one.
        let mismatched_match = ReactionMatch {
            per_reactant: vec![FxHashMap::default()],
        };
        let err = apply_reaction_match(smirks, &reactants, &mismatched_match, true);
        assert!(matches!(
            err,
            Err(TransformError::ReactantCountMismatch {
                expected: 2,
                got: 1
            })
        ));
    }

    /// `ReactionMatch::atom_map_positions` must likewise reject a
    /// shape mismatch rather than panicking.
    #[test]
    fn atom_map_positions_rejects_mismatched_match_shape() {
        let mismatched_match = ReactionMatch {
            per_reactant: vec![FxHashMap::default()],
        };
        let err = mismatched_match.atom_map_positions("[N:1].[C:2]>>[N:1][C:2]");
        assert!(matches!(
            err,
            Err(TransformError::ReactantCountMismatch {
                expected: 2,
                got: 1
            })
        ));
    }

    #[test]
    fn atom_tag_survives_apply_atom_map_does_not() {
        let mut mol = parse("CCO").unwrap();
        for i in 0..mol.atom_count() {
            mol.set_tag(AtomIdx(i as u32), Some(200 + i as u16));
        }
        let matches = find_reaction_matches("[C:1]>>[C:1]O", &[&mol]).expect("match");
        assert!(!matches.is_empty());
        let products =
            apply_reaction_match("[C:1]>>[C:1]O", &[&mol], &matches[0], true).expect("apply");
        let product = &products.expect("valence")[0];
        let tags: Vec<_> = product
            .atoms()
            .map(|(i, _)| product.atom_tag(i).map(core::num::NonZeroU16::get))
            .collect();
        assert!(
            tags.contains(&Some(200)) || tags.contains(&Some(201)) || tags.contains(&Some(202))
        );
        assert!(product.atoms().all(|(_, a)| a.atom_map.is_none()));
        // Born oxygen has no tag.
        assert!(
            product
                .atoms()
                .any(|(i, a)| product.atom_tag(i).is_none() && a.element.atomic_number() == 8)
        );
    }
}
