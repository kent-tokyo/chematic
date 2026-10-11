use super::*;

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
    pub(super) fn rdkit_stereo(&self) -> Option<&crate::rdkit_stereo::ReactionStereo> {
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
    pub(super) fn new_with_reading(
        smirks: &str,
        rdkit_reading: bool,
    ) -> Result<Self, TransformError> {
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

    pub(super) fn run_reactants_impl(
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

    pub(super) fn run_reactants_with_diagnostics_impl(
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
                    rejection_diagnostics.push(set_diagnostics.clone());
                    if profile == Profile::Rdkit {
                        set.iter_mut().for_each(restore_sources);
                        rejected_products.push(
                            set.into_iter()
                                .zip(set_diagnostics.iter().map(|item| item.accepted))
                                .collect(),
                        );
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

    /// Apply this native template while retaining typed atom-level evidence
    /// for every product set rejected during sanitization.
    pub fn run_reactants_traced_with_rejection_diagnostics(
        &self,
        reactants: &[&Molecule],
        limits: &ReactionTransformLimits,
    ) -> Result<NativeDetailedReactionTransformReport, TransformError> {
        crate::perf_counters::record_run_reactants_call();
        let variants: Vec<&PreparedReaction> = match &self.variants {
            Some(variants) => variants.iter().collect(),
            None => vec![self],
        };
        let mut products = Vec::new();
        let mut rejection_diagnostics = Vec::new();
        let mut diagnostics = ReactionTransformDiagnostics {
            accepted_matches: 0,
            applied_products: 0,
            valence_rejected_matches: 0,
            truncated_matches: false,
        };
        for variant in variants {
            let report = variant
                .run_traced_profile(reactants, true, limits, Profile::Native)?
                .unwrap_or_else(|_| unreachable!("the native profile never declines"));
            diagnostics.accepted_matches += report.report.diagnostics.accepted_matches;
            diagnostics.applied_products += report.report.diagnostics.applied_products;
            diagnostics.valence_rejected_matches +=
                report.report.diagnostics.valence_rejected_matches;
            diagnostics.truncated_matches |= report.report.diagnostics.truncated_matches;
            products.extend(report.report.products);
            rejection_diagnostics.extend(report.rejection_diagnostics);
        }
        Ok(NativeDetailedReactionTransformReport {
            report: TracedReactionTransformReport {
                products,
                diagnostics,
            },
            rejection_diagnostics,
        })
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
