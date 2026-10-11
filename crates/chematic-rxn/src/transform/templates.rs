use super::*;

pub(super) fn aggregate_variant_reports(
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

pub(super) fn variant_diagnostics(
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

pub(super) fn template_atom_maps_of(queries: &[QueryMolecule]) -> Vec<Vec<Option<u16>>> {
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
pub(super) type SmirksTemplates = (
    crate::reaction::Reaction,
    Vec<QueryMolecule>,
    Vec<Vec<ProductAtomSpec>>,
    Option<Vec<QueryMolecule>>,
);

pub(super) fn parse_smirks_templates(
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
pub(super) struct ProductAtomSpec {
    pub(super) charge: bool,
    pub(super) hcount: bool,
}

/// [`ProductAtomSpec`] for each atom of one SMILES product component, in
/// parse order.
pub(super) fn product_atom_specs(component: &str) -> Vec<ProductAtomSpec> {
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
pub(super) fn bracket_atom_spec(inner: &str) -> ProductAtomSpec {
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
pub(super) fn parse_reactant_templates(
    side: &str,
    rdkit_reading: bool,
) -> Result<(Vec<Molecule>, Vec<QueryMolecule>), TransformError> {
    parse_reactant_templates_with_limits(
        side,
        rdkit_reading,
        &crate::reaction::ReactionParseLimits::default(),
    )
}

pub(super) fn parse_reactant_templates_with_limits(
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
pub(super) fn rdkit_reactant_queries(
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
pub(super) fn strip_chirality(q: &AtomQuery) -> Option<AtomQuery> {
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
pub(super) fn simplify_reactant_query(query: &mut QueryMolecule) {
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
pub(super) fn has_stereo(mol: &Molecule) -> bool {
    mol.atoms().any(|(_, a)| a.chirality != Chirality::None)
        || mol
            .bonds()
            .any(|(_, b)| matches!(b.order, BondOrder::Up | BondOrder::Down))
}

/// Dot-separated components of one reaction side, ignoring dots inside
/// bracket atoms and recursive SMARTS parentheses.
pub(super) fn split_components(side: &str) -> Vec<&str> {
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
pub(super) fn shadow_molecule_of(
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
pub(super) fn rdkit_product_bonds(products: &str) -> String {
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
pub(super) fn rdkit_bond_symbol(expr: &str) -> char {
    expr.trim_start_matches('!')
        .chars()
        .next()
        .filter(|c| "-=#:$/\\".contains(*c))
        .unwrap_or('~')
}

/// Whether RDKit gives the product bond `q` no order: its first primitive
/// is a ring primitive or `~` inside a larger expression (a lone `~` is
/// RDKit's "any" bond, which keeps a reactant bond).
pub(super) fn zero_order_product_bond(q: &BondQuery) -> bool {
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
pub(super) fn zero_order_product_bonds(component: &mut Molecule, original: &str) {
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
pub(super) fn grouped_product_components(
    products: &str,
) -> Result<Option<Vec<String>>, TransformError> {
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
pub(super) fn normalize_product_templates(products: &str) -> Result<String, TransformError> {
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
pub(super) fn global_map_of(
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
