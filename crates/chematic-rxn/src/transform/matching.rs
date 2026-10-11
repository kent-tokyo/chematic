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
