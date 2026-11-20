//! Traverse VF2 embeddings without allocating result maps.

use super::*;

/// Visit every embedding [`find_matches_with_config`] returns for `config`
/// with `uniquify: false` and `max_matches: None`, in the same order, without
/// materializing result maps. `visit` receives the mapping as a slice indexed
/// by query atom (`mapping[q]` is the target atom index of query atom `q`).
///
/// `config.uniquify` and `config.max_matches` are ignored; a configured
/// `max_visit_budget` is honored exactly as in the map-returning search, and
/// the returned flag reports whether it was exhausted.
pub fn for_each_embedding(
    query: &QueryMolecule,
    mol: &Molecule,
    config: &MatchConfig,
    mut visit: impl FnMut(&[u32]),
) -> bool {
    visit_embeddings(query, mol, config, |mapping| {
        visit(mapping);
        true
    })
}

/// Whether the default matcher finds enough different target-atom sets.
/// Embeddings with the same target atoms count once, regardless of their order.
pub fn has_at_least_distinct_embeddings(
    query: &QueryMolecule,
    mol: &Molecule,
    count: usize,
) -> bool {
    if count == 0 {
        return true;
    }
    let config = MatchConfig::default();
    if count == 1 {
        return has_match_with_config(query, mol, &config);
    }
    let mut sets = FxHashSet::default();
    visit_embeddings(query, mol, &config, |mapping| {
        let mut atoms = mapping.to_vec();
        atoms.sort_unstable();
        sets.insert(atoms);
        sets.len() < count
    });
    sets.len() == count
}

fn visit_embeddings(
    query: &QueryMolecule,
    mol: &Molecule,
    config: &MatchConfig,
    mut visit: impl FnMut(&[u32]) -> bool,
) -> bool {
    if query.atoms.is_empty() || query.atoms.len() > mol.atom_count() {
        return false;
    }
    if config.max_visit_budget.is_none() && !element_counts_allow_match(query, mol) {
        return false;
    }
    let ctx = EvalCtx {
        mol,
        rings_given: None,
        rings_lazy: std::cell::OnceCell::new(),
        ring_atoms_lazy: std::cell::OnceCell::new(),
        hybridization_lazy: std::cell::OnceCell::new(),
        ez_facts_lazy: std::sync::OnceLock::new(),
        config,
        visit_budget: std::cell::Cell::new(config.max_visit_budget.unwrap_or(u64::MAX)),
        budget_exhausted: std::cell::Cell::new(false),
        min_ring_size_by_atom: std::cell::RefCell::new(None),
        cycle_atoms: std::cell::OnceCell::new(),
        prune_warmup: std::cell::Cell::new(0),
        plans: std::cell::RefCell::new(FxHashMap::default()),
        scratch: std::cell::RefCell::new(Vec::new()),
    };
    let mut st = MapState::new(query, mol.atom_count(), None, &ctx);
    visit_recursive(query, &ctx, &mut st, &mut visit);
    ctx.budget_exhausted.get()
}

/// [`match_recursive`] without result maps: the same candidates, checks and
/// visiting order, stopping when `visit` declines another embedding.
fn visit_recursive(
    query: &QueryMolecule,
    ctx: &EvalCtx<'_>,
    st: &mut MapState,
    visit: &mut impl FnMut(&[u32]) -> bool,
) -> bool {
    let remaining = ctx.visit_budget.get();
    if remaining == 0 {
        ctx.budget_exhausted.set(true);
        return false;
    }
    ctx.visit_budget.set(remaining - 1);

    if st.count == query.atoms.len() {
        return !mapped_chirality_compatible(query, &st.q2t, ctx) || visit(&st.q2t);
    }
    let q_next = st.info.plan[st.count];
    let candidates = Candidates::new(q_next, st, query, ctx);
    for t in candidates.iter() {
        let t_idx = AtomIdx(t);
        if st.used[t as usize] {
            continue;
        }
        if st.info.implied[q_next]
            .is_some_and(|z| effective_atomic_number(ctx.mol.atom(t_idx)) != z)
        {
            continue;
        }
        // Exact pruning (dead branches only); no result-map layout to keep.
        if ctx.prunes_low_degree(st, q_next, t_idx)
            || ctx.prunes_acyclic_target(query, st, q_next, t_idx)
        {
            continue;
        }
        if !eval_atom_query(&query.atoms[q_next].query, t_idx, ctx) {
            continue;
        }
        if !bonds_compatible(q_next, t_idx, st, query, ctx) {
            continue;
        }
        st.set(q_next, t_idx);
        let keep_visiting = visit_recursive(query, ctx, st, visit);
        st.unset(q_next, t_idx);
        if !keep_visiting {
            return false;
        }
    }
    true
}
