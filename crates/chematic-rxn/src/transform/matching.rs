use super::*;

/// Steps 1–2 of the original `run_reactants_impl`: VF2-match every reactant
/// template against its input molecule, take the cartesian product across
/// template slots, and keep only combinations that survive the
/// chirality/E-Z stereo post-checks. Does not build any product.
pub(super) fn find_matches_impl(
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
pub(super) fn find_matches_profile(
    prepared: &PreparedReaction,
    reactants: &[&Molecule],
    limits: &ReactionTransformLimits,
    rings: Option<&[&RingSet]>,
    profile: Profile,
) -> Result<(Vec<ReactionMatch>, Option<ReactionCompatibilityUnsupported>), TransformError> {
    let views = aromatic_views(reactants);
    let perceived = perceived_refs(reactants, &views);
    let reactants: &[&Molecule] = &perceived;
    validate_match_inputs(prepared.rxn.reactants.len(), reactants, rings)?;
    let queries = match (profile, &prepared.rdkit_queries) {
        (Profile::Rdkit, Some(rdkit)) => rdkit,
        _ => &prepared.queries,
    };
    let Some(all_match_sets) = collect_match_sets(queries, reactants, rings, profile, limits)?
    else {
        return Ok((
            vec![],
            Some(ReactionCompatibilityUnsupported::RingModelUnresolved),
        ));
    };
    // No matches when any template has no match.
    if all_match_sets.iter().any(|ms| ms.is_empty()) {
        return Ok((vec![], None));
    }
    validate_combination_count(&all_match_sets, limits)?;

    let mut matches: Vec<ReactionMatch> = Vec::new();
    let mut ez_rejected = false;
    for combo in cartesian_product(&all_match_sets) {
        crate::perf_counters::record_match_combination();
        if !combo_chirality_ok(prepared, reactants, &combo, profile) {
            continue;
        }
        if !combo_ez_ok(prepared, reactants, &combo) {
            ez_rejected = true;
            continue;
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

type AtomMatch = FxHashMap<usize, AtomIdx>;

fn validate_match_inputs(
    expected: usize,
    reactants: &[&Molecule],
    rings: Option<&[&RingSet]>,
) -> Result<(), TransformError> {
    if reactants.len() != expected {
        return Err(TransformError::ReactantCountMismatch {
            expected,
            got: reactants.len(),
        });
    }
    if let Some(rings) = rings
        && rings.len() != expected
    {
        return Err(TransformError::ReactantCountMismatch {
            expected,
            got: rings.len(),
        });
    }
    Ok(())
}

fn reaction_match_config() -> MatchConfig {
    MatchConfig {
        // Reactions cannot deduplicate embeddings by target atom set: swapping
        // template atoms may produce distinct products.
        uniquify: false,
        // Reactant templates are exact molecular specifications.
        use_isotopes: true,
        ..MatchConfig::default()
    }
}

fn collect_match_sets(
    queries: &[QueryMolecule],
    reactants: &[&Molecule],
    rings: Option<&[&RingSet]>,
    profile: Profile,
    limits: &ReactionTransformLimits,
) -> Result<Option<Vec<Vec<AtomMatch>>>, TransformError> {
    let config = reaction_match_config();
    let mut all = Vec::with_capacity(queries.len());
    for (index, (query, mol)) in queries.iter().zip(reactants).enumerate() {
        let Some(matches) = match_query(query, mol, rings.map(|r| r[index]), profile, &config)
        else {
            return Ok(None);
        };
        crate::perf_counters::record_reactant_query_match_call(matches.len());
        if matches.len() > limits.max_matches {
            return Err(TransformError::ResourceLimit {
                resource: "reaction matches",
                actual: matches.len(),
                limit: limits.max_matches,
            });
        }
        all.push(matches);
    }
    Ok(Some(all))
}

fn match_query(
    query: &QueryMolecule,
    mol: &Molecule,
    rings: Option<&RingSet>,
    profile: Profile,
    config: &MatchConfig,
) -> Option<Vec<AtomMatch>> {
    if let Some(rings) = rings {
        return Some(find_matches_with_rings_and_config(
            query, mol, rings, config,
        ));
    }
    if profile != Profile::Rdkit || !query_uses_ring_model(query) || rdkit_ring_model_is_sssr(mol) {
        return Some(find_matches_with_config(query, mol, config));
    }
    let config = chematic_smarts::RdkitParityConfig {
        base: config.clone(),
        ..chematic_smarts::RdkitParityConfig::default()
    };
    chematic_smarts::find_matches_rdkit_parity(query, mol, &config)
        .ok()
        .map(|(matches, _)| matches)
}

fn validate_combination_count(
    match_sets: &[Vec<AtomMatch>],
    limits: &ReactionTransformLimits,
) -> Result<(), TransformError> {
    let actual = match_sets
        .iter()
        .map(Vec::len)
        .try_fold(1usize, usize::checked_mul)
        .unwrap_or(usize::MAX);
    if actual > limits.max_matches {
        return Err(TransformError::ResourceLimit {
            resource: "reaction match combinations",
            actual,
            limit: limits.max_matches,
        });
    }
    Ok(())
}

fn combo_chirality_ok(
    prepared: &PreparedReaction,
    reactants: &[&Molecule],
    combo: &[AtomMatch],
    profile: Profile,
) -> bool {
    !prepared.has_stereo
        || profile != Profile::Native
        || prepared
            .rxn
            .reactants
            .iter()
            .enumerate()
            .all(|(index, template)| smirks_chirality_ok(template, reactants[index], &combo[index]))
}

fn combo_ez_ok(prepared: &PreparedReaction, reactants: &[&Molecule], combo: &[AtomMatch]) -> bool {
    !prepared.has_ez_stereo
        || prepared
            .rxn
            .reactants
            .iter()
            .enumerate()
            .all(|(index, template)| smirks_ez_stereo_ok(template, reactants[index], &combo[index]))
}

/// Whether RDKit's ring models give the same rings as the SSSR for `mol`
/// (no extra symmetrized rings, no organometallic view), so the ordinary
/// matcher answers `R<n>` and ring-size queries as RDKit does. Memoized on
/// the (cached) aromatic view; only cages such as quinuclidine or adamantane
/// need the slower RDKit-parity matcher.
pub(super) fn rdkit_ring_model_is_sssr(mol: &Molecule) -> bool {
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
pub(super) fn query_uses_ring_model(query: &chematic_smarts::QueryMolecule) -> bool {
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
pub(super) type AromaticView =
    Option<std::sync::Arc<Result<Molecule, chematic_perception::AromaticityError>>>;

pub(super) fn aromatic_views(reactants: &[&Molecule]) -> Vec<AromaticView> {
    reactants
        .iter()
        .map(|m| {
            (!chematic_perception::rdkit_parity_view_is_identity(m))
                .then(|| chematic_perception::apply_aromaticity_rdkit_parity_shared(m))
        })
        .collect()
}

pub(super) fn perceived_refs<'a>(
    reactants: &[&'a Molecule],
    views: &'a [AromaticView],
) -> Vec<&'a Molecule> {
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
pub(super) fn apply_match_impl(
    prepared: &PreparedReaction,
    reactants: &[&Molecule],
    m: &ReactionMatch,
    carry_substituents: bool,
) -> Option<Vec<Molecule>> {
    apply_match_traced_impl(prepared, reactants, m, carry_substituents)
        .map(|products| products.into_iter().map(|p| p.molecule).collect())
}

pub(super) fn apply_match_traced_impl(
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

pub(super) fn apply_match_profile(
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
            clear_product_template_unspecified_ez(&mut product, pt, &template_idx_to_new);
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
