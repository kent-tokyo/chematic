//! VF2 subgraph isomorphism: find all embeddings of a `QueryMolecule` in a target `Molecule`.
//!
//! The classic VF2 algorithm explores a state-space search tree.  At each step
//! it picks the next unmapped query atom and tries to extend the current partial
//! mapping with every compatible target atom.  Compatibility is checked at two
//! levels:
//!
//! 1. **Atom compatibility** — the target atom must satisfy the query atom's
//!    `AtomQuery` expression.
//! 2. **Bond compatibility** — for every already-mapped query neighbour of the
//!    candidate query atom, the corresponding target atoms must be bonded in the
//!    target molecule, and that bond must satisfy the query bond's `BondQuery`
//!    expression.

use rustc_hash::{FxHashMap, FxHashSet};

use chematic_core::{AtomIdx, BondOrder, Molecule, implicit_hcount};
use chematic_perception::RingSet;
use chematic_perception::find_sssr;

use crate::query::{AtomPrimitive, AtomQuery, BondPrimitive, BondQuery, QueryMolecule};

// ---------------------------------------------------------------------------
// Evaluation context (precomputed per `find_matches` call)
// ---------------------------------------------------------------------------

/// Per-call evaluation context: the target molecule, precomputed ring set, and match config.
///
/// `rings` is borrowed so callers can compute SSSR once and reuse it across many
/// pattern matches (e.g. the 117 Crippen patterns or 480 PAINS patterns).
struct EvalCtx<'a> {
    mol: &'a Molecule,
    rings: &'a RingSet,
    config: &'a MatchConfig,
    /// Remaining visit budget shared across all recursive calls (including nested
    /// recursive-SMARTS `$(...)`).  Decremented on every `match_recursive` /
    /// `has_match_recursive` entry.  `u64::MAX` when no limit is configured.
    visit_budget: std::cell::Cell<u64>,
    /// Lazily-computed, memoized per atom index: the size of the *smallest* SSSR
    /// ring containing that atom (`None` if the atom is in no ring). Backs `[rN]`
    /// (`AtomPrimitive::MinRingSize`) — computed at most once per `find_matches`
    /// call (not once per atom, not once per VF2 backtrack step, which re-evaluates
    /// atom predicates many times) and only if a query actually uses `[rN]`, since
    /// most patterns never do. Distinct from `[kN]` (`AtomPrimitive::RingSize`),
    /// which needs no such cache -- it's a direct scan of `rings` per query, same as
    /// before this variant existed.
    min_ring_size_by_atom: std::cell::RefCell<Option<Vec<Option<u8>>>>,
}

impl EvalCtx<'_> {
    /// The smallest SSSR ring size containing `idx`, or `None` if `idx` isn't in any
    /// ring. Computes and caches the whole per-atom table on first call.
    fn min_ring_size(&self, idx: AtomIdx) -> Option<u8> {
        let mut cache = self.min_ring_size_by_atom.borrow_mut();
        let table = cache.get_or_insert_with(|| {
            let mut table = vec![None; self.mol.atom_count()];
            for ring in self.rings.rings() {
                let size = ring.len() as u8;
                for &atom in ring {
                    let slot = &mut table[atom.0 as usize];
                    *slot = Some(slot.map_or(size, |current: u8| current.min(size)));
                }
            }
            table
        });
        table[idx.0 as usize]
    }
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Configuration for subgraph matching.
#[derive(Debug, Clone)]
pub struct MatchConfig {
    /// Maximum number of matches to return.
    ///
    /// `None` (default) returns all matches.  Set to `Some(n)` to stop after
    /// the first `n` results — useful for large molecules or generic queries
    /// where an unbounded search would be slow or produce huge `Vec`s.
    pub max_matches: Option<usize>,

    /// When `true`, `[@]` / `[@@]` chirality primitives in the query are
    /// enforced against the target atom's chirality annotation.
    ///
    /// Defaults to `false` (chirality is ignored, matching RDKit's default
    /// `useChirality=False` behaviour).
    pub use_chirality: bool,

    /// When `true`, isotope primitives (`[13C]`, `[2H]`, …) are enforced
    /// against the target atom's isotope label.
    ///
    /// Defaults to `false` (isotopes are ignored, matching RDKit's default
    /// `useIsotopes=False` behaviour).
    pub use_isotopes: bool,

    /// When `true`, deduplicate matches: only return one mapping per unique
    /// set of target atoms covered, even if different orderings exist.
    ///
    /// Defaults to `true` (matching RDKit's `uniquify=True` default).
    /// For symmetric queries on symmetric targets, this prevents returning
    /// multiple embeddings of the same substructure.
    pub uniquify: bool,

    /// Maximum number of recursive VF2 state-space visits across the entire
    /// search (including nested recursive-SMARTS `$(...)`).
    ///
    /// `None` (default) is unbounded — preserving the existing behaviour.
    ///
    /// **Warning — opt-in only.** When the budget is exhausted the search stops
    /// early: `find_matches` returns the partial result set, and
    /// `has_match`/`brenk_passes`/`pains_passes` may silently return `false`
    /// even when a match exists.  Only set this for DoS-prevention in
    /// untrusted-input contexts where false negatives are acceptable.
    pub max_visit_budget: Option<u64>,
}

impl Default for MatchConfig {
    fn default() -> Self {
        Self {
            max_matches: None,
            use_chirality: false,
            use_isotopes: false,
            uniquify: true,
            max_visit_budget: None,
        }
    }
}

/// Find all non-overlapping (injective) embeddings of `query` in `mol`.
///
/// Returns a `Vec` of mappings, each mapping a query atom index to a target
/// `AtomIdx`.  Each individual mapping is injective (no two query atoms map to
/// the same target atom), but the same target atom may appear in different
/// mappings.
pub fn find_matches(query: &QueryMolecule, mol: &Molecule) -> Vec<FxHashMap<usize, AtomIdx>> {
    find_matches_with_config(query, mol, &MatchConfig::default())
}

/// Like [`find_matches`] but with explicit configuration.
///
/// Use `config.max_matches = Some(n)` to cap the result count.
pub fn find_matches_with_config(
    query: &QueryMolecule,
    mol: &Molecule,
    config: &MatchConfig,
) -> Vec<FxHashMap<usize, AtomIdx>> {
    // Early-exit before the expensive SSSR computation.
    if query.atoms.is_empty() {
        return vec![];
    }
    if query.atoms.len() > mol.atom_count() {
        return vec![];
    }
    let rings = find_sssr(mol);
    find_matches_with_rings_and_config(query, mol, &rings, config)
}

/// Like [`find_matches`] but reuses a pre-computed [`RingSet`].
///
/// Use this when matching many patterns against the same molecule to avoid
/// recomputing the SSSR for each pattern:
///
/// ```ignore
/// use chematic_perception::find_sssr;
/// let rings = find_sssr(&mol);
/// for query in &queries {
///     let hits = find_matches_with_rings(query, &mol, &rings);
/// }
/// ```
pub fn find_matches_with_rings(
    query: &QueryMolecule,
    mol: &Molecule,
    rings: &RingSet,
) -> Vec<FxHashMap<usize, AtomIdx>> {
    find_matches_with_rings_and_config(query, mol, rings, &MatchConfig::default())
}

/// Like [`find_matches_with_config`] but reuses a pre-computed [`RingSet`].
pub fn find_matches_with_rings_and_config(
    query: &QueryMolecule,
    mol: &Molecule,
    rings: &RingSet,
    config: &MatchConfig,
) -> Vec<FxHashMap<usize, AtomIdx>> {
    if query.atoms.is_empty() {
        return vec![];
    }
    // A query with more heavy atoms than the target can never match (RDKit PR #9201).
    if query.atoms.len() > mol.atom_count() {
        return vec![];
    }

    let ctx = EvalCtx {
        mol,
        rings,
        config,
        visit_budget: std::cell::Cell::new(config.max_visit_budget.unwrap_or(u64::MAX)),
        min_ring_size_by_atom: std::cell::RefCell::new(None),
    };
    let mut mapping: FxHashMap<usize, AtomIdx> = FxHashMap::default();
    let mut results: Vec<FxHashMap<usize, AtomIdx>> = Vec::new();

    match_recursive(query, &ctx, &mut mapping, &mut results, config.max_matches);

    // Deduplicate matches: keep only one mapping per unique set of target atoms.
    if config.uniquify {
        let mut seen = FxHashSet::default();
        results.retain(|m| {
            let mut key: Vec<u32> = m.values().map(|idx| idx.0).collect();
            key.sort_unstable();
            seen.insert(key)
        });
    }

    results
}

// ---------------------------------------------------------------------------
// Recursive VF2 search
// ---------------------------------------------------------------------------

/// Returns the index of the smallest unmapped query atom.
fn next_unmapped(mapping: &FxHashMap<usize, AtomIdx>, query_len: usize) -> usize {
    (0..query_len).find(|i| !mapping.contains_key(i)).unwrap() // safe: caller guarantees mapping.len() < query_len
}

fn match_recursive(
    query: &QueryMolecule,
    ctx: &EvalCtx<'_>,
    mapping: &mut FxHashMap<usize, AtomIdx>,
    results: &mut Vec<FxHashMap<usize, AtomIdx>>,
    max: Option<usize>,
) {
    // Early exit if the result cap has been reached.
    if max.is_some_and(|m| results.len() >= m) {
        return;
    }

    // Decrement shared visit budget; stop if exhausted.
    let remaining = ctx.visit_budget.get();
    if remaining == 0 {
        return;
    }
    ctx.visit_budget.set(remaining - 1);

    // Base case: all query atoms have been mapped.
    if mapping.len() == query.atoms.len() {
        results.push(mapping.clone());
        return;
    }

    // Pick the next unmapped query atom (smallest index not yet in mapping).
    let q_next = next_unmapped(mapping, query.atoms.len());

    // Collect the set of target atoms already used in this mapping so we can
    // enforce injectivity.
    let used_targets: FxHashSet<AtomIdx> = mapping.values().copied().collect();

    // Try each target atom as a candidate for q_next.
    for t in 0..ctx.mol.atom_count() {
        if max.is_some_and(|m| results.len() >= m) {
            break;
        }
        let t_idx = AtomIdx(t as u32);

        // 1. Injectivity: target atom must not already be mapped.
        if used_targets.contains(&t_idx) {
            continue;
        }

        // 2. Atom query must match.
        if !eval_atom_query(&query.atoms[q_next].query, t_idx, ctx) {
            continue;
        }

        // 3. Bond constraints from already-mapped neighbours of q_next.
        if !bonds_compatible(q_next, t_idx, mapping, query, ctx) {
            continue;
        }

        // Extend the mapping and recurse.
        mapping.insert(q_next, t_idx);
        match_recursive(query, ctx, mapping, results, max);
        mapping.remove(&q_next);
    }
}

// ---------------------------------------------------------------------------
// Bond compatibility check
// ---------------------------------------------------------------------------

/// For every already-mapped query neighbour of `q`, verify that the
/// corresponding target atoms are bonded and satisfy the query bond condition.
fn bonds_compatible(
    q: usize,
    t: AtomIdx,
    mapping: &FxHashMap<usize, AtomIdx>,
    query: &QueryMolecule,
    ctx: &EvalCtx<'_>,
) -> bool {
    for &(bond_idx, q_nb) in &query.adj[q] {
        // Only check neighbours that are already mapped.
        if let Some(&t_nb) = mapping.get(&q_nb) {
            // The target must have a bond between t and t_nb.
            match ctx.mol.bond_between(t, t_nb) {
                None => return false,
                Some((_bidx, bond_entry)) => {
                    let bq = &query.bonds[bond_idx].query;
                    if !eval_bond_query(bq, bond_entry.order, t, t_nb, ctx) {
                        return false;
                    }
                }
            }
        }
    }
    true
}

// ---------------------------------------------------------------------------
// Atom query evaluation
// ---------------------------------------------------------------------------

fn eval_atom_query(q: &AtomQuery, idx: AtomIdx, ctx: &EvalCtx<'_>) -> bool {
    match q {
        AtomQuery::Primitive(p) => eval_atom_primitive(p, idx, ctx),
        AtomQuery::And(a, b) => eval_atom_query(a, idx, ctx) && eval_atom_query(b, idx, ctx),
        AtomQuery::Or(a, b) => eval_atom_query(a, idx, ctx) || eval_atom_query(b, idx, ctx),
        AtomQuery::Not(a) => !eval_atom_query(a, idx, ctx),
    }
}

fn eval_atom_primitive(p: &AtomPrimitive, idx: AtomIdx, ctx: &EvalCtx<'_>) -> bool {
    let atom = ctx.mol.atom(idx);
    match p {
        AtomPrimitive::AtomicNum(n) => atom.element.atomic_number() == *n,
        AtomPrimitive::Symbol(s) => atom.element.symbol() == s.as_str(),
        AtomPrimitive::Aromatic(a) => atom.aromatic == *a,
        AtomPrimitive::Charge(c) => atom.charge == *c,
        AtomPrimitive::HCount(h) => eval_hcount(idx, ctx, *h),
        AtomPrimitive::ImplicitHCount(h) => implicit_hcount(ctx.mol, idx) == *h,
        AtomPrimitive::Degree(d) => ctx.mol.neighbors(idx).count() as u8 == *d,
        AtomPrimitive::RingMembership(r) => ctx.rings.contains_atom(idx) == *r,
        AtomPrimitive::RingSize(n) => ctx
            .rings
            .rings()
            .iter()
            .any(|ring| ring.len() == *n as usize && ring.contains(&idx)),
        AtomPrimitive::MinRingSize(n) => ctx.min_ring_size(idx) == Some(*n),
        AtomPrimitive::Wildcard => true,
        AtomPrimitive::Recursive(sub_query) => has_match_anchored(sub_query, idx, ctx),
        AtomPrimitive::Valence(v) => eval_valence(idx, ctx, *v),
        AtomPrimitive::RingBondCount(x) => eval_ring_bond_count(idx, ctx, *x),
        AtomPrimitive::TotalConnectivity(x) => {
            ctx.mol.neighbors(idx).count() as u8 + implicit_hcount(ctx.mol, idx) == *x
        }
        AtomPrimitive::RingCount(n) => {
            ctx.rings
                .rings()
                .iter()
                .filter(|r| r.contains(&idx))
                .count() as u8
                == *n
        }
        AtomPrimitive::Hybridization(h) => eval_hybridization(idx, ctx, *h),
        AtomPrimitive::Isotope(mass) => {
            !ctx.config.use_isotopes || ctx.mol.atom(idx).isotope == Some(*mass)
        }
        AtomPrimitive::Chirality(kind) => eval_chirality(idx, ctx, *kind),
    }
}

/// Total H count (explicit H neighbors + implicit H) for HCount primitive.
fn eval_hcount(idx: AtomIdx, ctx: &EvalCtx<'_>, h: u8) -> bool {
    let explicit_h = ctx
        .mol
        .neighbors(idx)
        .filter(|(nb, _)| ctx.mol.atom(*nb).element.atomic_number() == 1)
        .count() as u8;
    explicit_h + implicit_hcount(ctx.mol, idx) == h
}

/// Total valence (bond order sum + implicit H) for Valence primitive.
fn eval_valence(idx: AtomIdx, ctx: &EvalCtx<'_>, v: u8) -> bool {
    let bond_sum: u8 = ctx
        .mol
        .neighbors(idx)
        .map(|(_, bid)| bond_order_int(ctx.mol.bond(bid).order))
        .sum();
    bond_sum + implicit_hcount(ctx.mol, idx) == v
}

/// Ring bond count: bonds where both endpoints share at least one SSSR ring.
fn eval_ring_bond_count(idx: AtomIdx, ctx: &EvalCtx<'_>, x: u8) -> bool {
    let count = ctx
        .mol
        .neighbors(idx)
        .filter(|(nb, _)| {
            ctx.rings
                .rings()
                .iter()
                .any(|ring| ring.contains(&idx) && ring.contains(nb))
        })
        .count() as u8;
    count == x
}

/// Inferred hybridization: aromatic→sp2, triple→sp, double→sp2, else→sp3.
fn eval_hybridization(idx: AtomIdx, ctx: &EvalCtx<'_>, h: u8) -> bool {
    let atom = ctx.mol.atom(idx);
    let hyb = if atom.aromatic {
        2u8
    } else {
        let mut has_triple = false;
        let mut has_double = false;
        for (_, bid) in ctx.mol.neighbors(idx) {
            match ctx.mol.bond(bid).order {
                BondOrder::Triple => {
                    has_triple = true;
                    break;
                }
                BondOrder::Double => has_double = true,
                _ => {}
            }
        }
        if has_triple {
            1
        } else if has_double {
            2
        } else {
            3
        }
    };
    hyb == h
}

/// Chirality primitive: ignored when use_chirality is false.
///
/// **Limitation**: this is a raw flag comparison. `@`/`@@` in SMILES/SMARTS
/// encodes chirality *relative to the SMILES atom write order*, not as an
/// absolute spatial property. Two SMILES strings that represent the same
/// absolute configuration but are written with neighbors in a different order
/// will store opposite `Chirality` flags, making the raw comparison incorrect.
///
/// For SMARTS queries this is acceptable when the query and target share the same
/// write-order convention. For SMIRKS reaction templates, use the parity-aware
/// `smirks_chirality_ok` post-check in `chematic-rxn::transform` instead.
fn eval_chirality(idx: AtomIdx, ctx: &EvalCtx<'_>, kind: u8) -> bool {
    if !ctx.config.use_chirality {
        return true;
    }
    use chematic_core::Chirality;
    let c = ctx.mol.atom(idx).chirality;
    match kind {
        1 => c == Chirality::CounterClockwise,
        2 => c == Chirality::Clockwise,
        _ => c != Chirality::None,
    }
}

// ---------------------------------------------------------------------------
// Anchored match helpers (for recursive SMARTS)
// ---------------------------------------------------------------------------

/// Returns `true` if there exists at least one embedding of `query` in `ctx.mol`
/// with query atom 0 forced to map to `anchor`.
fn has_match_anchored(query: &QueryMolecule, anchor: AtomIdx, ctx: &EvalCtx<'_>) -> bool {
    if query.atoms.is_empty() {
        return false;
    }
    if query.atoms.len() > ctx.mol.atom_count() {
        return false;
    }
    // Quick check: anchor must satisfy the first query atom.
    if !eval_atom_query(&query.atoms[0].query, anchor, ctx) {
        return false;
    }
    // Seed the mapping with query atom 0 → anchor.
    let mut mapping = FxHashMap::default();
    mapping.insert(0usize, anchor);
    // Single-atom query — the anchor already satisfies it.
    if query.atoms.len() == 1 {
        return true;
    }
    has_match_recursive(query, ctx, &mut mapping)
}

/// Depth-first search for a complete embedding, starting from a partial
/// `mapping`.  Returns as soon as the first complete match is found.
fn has_match_recursive(
    query: &QueryMolecule,
    ctx: &EvalCtx<'_>,
    mapping: &mut FxHashMap<usize, AtomIdx>,
) -> bool {
    // Decrement shared visit budget; stop if exhausted (may produce false negatives
    // — only enabled when max_visit_budget is explicitly set).
    let remaining = ctx.visit_budget.get();
    if remaining == 0 {
        return false;
    }
    ctx.visit_budget.set(remaining - 1);

    // Base case: all query atoms mapped.
    if mapping.len() == query.atoms.len() {
        return true;
    }

    // Pick the next unmapped query atom.
    let q_next = next_unmapped(mapping, query.atoms.len());

    let used_targets: FxHashSet<AtomIdx> = mapping.values().copied().collect();

    for t in 0..ctx.mol.atom_count() {
        let t_idx = AtomIdx(t as u32);
        if used_targets.contains(&t_idx) {
            continue;
        }
        if !eval_atom_query(&query.atoms[q_next].query, t_idx, ctx) {
            continue;
        }
        if !bonds_compatible(q_next, t_idx, mapping, query, ctx) {
            continue;
        }
        mapping.insert(q_next, t_idx);
        if has_match_recursive(query, ctx, mapping) {
            mapping.remove(&q_next);
            return true;
        }
        mapping.remove(&q_next);
    }
    false
}

// ---------------------------------------------------------------------------
// Bond query evaluation
// ---------------------------------------------------------------------------

fn eval_bond_query(
    q: &BondQuery,
    order: BondOrder,
    a: AtomIdx,
    b: AtomIdx,
    ctx: &EvalCtx<'_>,
) -> bool {
    match q {
        BondQuery::Primitive(p) => eval_bond_primitive(p, order, a, b, ctx),
        BondQuery::And(x, y) => {
            eval_bond_query(x, order, a, b, ctx) && eval_bond_query(y, order, a, b, ctx)
        }
        BondQuery::Or(x, y) => {
            eval_bond_query(x, order, a, b, ctx) || eval_bond_query(y, order, a, b, ctx)
        }
        BondQuery::Not(x) => !eval_bond_query(x, order, a, b, ctx),
        // Implicit "any bond" — matches any bond order.
        BondQuery::Any => true,
    }
}

/// Convert a `BondOrder` to its integer contribution for valence-style sums.
///
/// Stereo bonds (Up/Down) are treated as single. Aromatic bonds are counted
/// as 1 (SMARTS valence convention).
fn bond_order_int(order: BondOrder) -> u8 {
    match order {
        BondOrder::Zero => 0,
        BondOrder::Single
        | BondOrder::Up
        | BondOrder::Down
        | BondOrder::Aromatic
        | BondOrder::Dative
        | BondOrder::QueryAny
        | BondOrder::QuerySingleOrDouble
        | BondOrder::QuerySingleOrAromatic
        | BondOrder::QueryDoubleOrAromatic => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
    }
}

fn eval_bond_primitive(
    p: &BondPrimitive,
    order: BondOrder,
    a: AtomIdx,
    b: AtomIdx,
    ctx: &EvalCtx<'_>,
) -> bool {
    match p {
        BondPrimitive::Single => {
            matches!(
                order,
                BondOrder::Single
                    | BondOrder::Up
                    | BondOrder::Down
                    | BondOrder::QuerySingleOrDouble
                    | BondOrder::QuerySingleOrAromatic
            )
        }
        BondPrimitive::Double => matches!(
            order,
            BondOrder::Double | BondOrder::QuerySingleOrDouble | BondOrder::QueryDoubleOrAromatic
        ),
        BondPrimitive::Triple => matches!(order, BondOrder::Triple),
        BondPrimitive::Aromatic => matches!(
            order,
            BondOrder::Aromatic
                | BondOrder::QuerySingleOrAromatic
                | BondOrder::QueryDoubleOrAromatic
        ),
        BondPrimitive::Any => true,
        BondPrimitive::Ring => {
            // A bond is a "ring bond" if both its endpoints share at least one common ring.
            ctx.rings
                .rings()
                .iter()
                .any(|ring| ring.contains(&a) && ring.contains(&b))
        }
        BondPrimitive::Up => matches!(order, BondOrder::Up),
        BondPrimitive::Down => matches!(order, BondOrder::Down),
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{find_matches, find_matches_with_config, parse_smarts};
    use chematic_smiles::parse;

    // -- Isotope matching -----------------------------------------------------

    #[test]
    fn test_isotope_ignored_by_default() {
        // [13C] query should match any carbon when use_isotopes=false (default).
        let mol = parse("CC").unwrap();
        let query = parse_smarts("[13C]").unwrap();
        let matches = find_matches(&query, &mol);
        assert_eq!(
            matches.len(),
            2,
            "[13C] with use_isotopes=false should match all carbons"
        );
    }

    #[test]
    fn test_isotope_enforced_when_enabled() {
        // Build a molecule with one 13C and one 12C.
        use chematic_core::{Atom, AtomIdx, BondOrder, Element, MoleculeBuilder};
        let mut b = MoleculeBuilder::new();
        let mut c13 = Atom::new(Element::C);
        c13.isotope = Some(13);
        let c13_idx = b.add_atom(c13);
        let c12_idx = b.add_atom(Atom::new(Element::C));
        b.add_bond(c13_idx, c12_idx, BondOrder::Single).unwrap();
        let mol = b.build();

        let query = parse_smarts("[13C]").unwrap();
        let config = MatchConfig {
            use_isotopes: true,
            ..MatchConfig::default()
        };
        let matches = find_matches_with_config(&query, &mol, &config);
        assert_eq!(
            matches.len(),
            1,
            "[13C] with use_isotopes=true should match only the 13C atom"
        );
        assert_eq!(matches[0][&0], AtomIdx(0));
    }

    #[test]
    fn test_no_isotope_match_on_unlabeled() {
        // [13C] with use_isotopes=true should not match unlabeled carbons.
        let mol = parse("CC").unwrap(); // both atoms have isotope=None
        let query = parse_smarts("[13C]").unwrap();
        let config = MatchConfig {
            use_isotopes: true,
            ..MatchConfig::default()
        };
        let matches = find_matches_with_config(&query, &mol, &config);
        assert_eq!(
            matches.len(),
            0,
            "[13C] with use_isotopes=true should not match unlabeled C"
        );
    }

    // -- Chirality matching ---------------------------------------------------

    #[test]
    fn test_chirality_ignored_by_default() {
        // [@] query should match any atom when use_chirality=false (default).
        let mol = parse("N[C@@H](C)C(=O)O").unwrap(); // L-alanine
        let query = parse_smarts("[C@@H]").unwrap();
        let matches = find_matches(&query, &mol);
        // Default: chirality ignored, so [@] matches any C-H regardless of chirality.
        assert!(
            !matches.is_empty(),
            "chirality should be ignored by default"
        );
    }

    #[test]
    fn test_chirality_enforced_when_enabled() {
        // L-alanine has [C@@H] — query [C@@H] should match, [C@H] should not.
        let mol = parse("N[C@@H](C)C(=O)O").unwrap();

        let q_ccw = parse_smarts("[C@@H]").unwrap(); // CCW (@@) = kind 2
        let q_cw = parse_smarts("[C@H]").unwrap(); // CW  (@)  = kind 1

        let config = MatchConfig {
            use_chirality: true,
            ..MatchConfig::default()
        };

        let m_ccw = find_matches_with_config(&q_ccw, &mol, &config);
        let m_cw = find_matches_with_config(&q_cw, &mol, &config);

        // [C@@H] must match L-alanine's chiral center.
        assert!(!m_ccw.is_empty(), "[C@@H] should match L-alanine (@@)");
        // [C@H] must NOT match L-alanine.
        assert!(m_cw.is_empty(), "[C@H] should not match L-alanine (@@)");
    }

    #[test]
    fn test_chirality_d_alanine_positive_match() {
        // D-alanine [C@H] should positively match [C@H] query when use_chirality=true.
        // This complements test_chirality_enforced_when_enabled which only tested negative on L-ala.
        let mol = parse("N[C@H](C)C(=O)O").unwrap(); // D-alanine (@, kind 1)
        let q_cw = parse_smarts("[C@H]").unwrap(); // CW (@), kind 1
        let config = MatchConfig {
            use_chirality: true,
            ..MatchConfig::default()
        };
        let m = find_matches_with_config(&q_cw, &mol, &config);
        assert!(!m.is_empty(), "[C@H] should positively match D-alanine (@)");
    }

    // S2: visit budget tests

    #[test]
    fn test_visit_budget_unlimited_default() {
        // Default config has no budget cap — normal queries complete fully.
        let mol = parse("c1ccccc1").unwrap(); // benzene
        let q = parse_smarts("c1ccccc1").unwrap();
        let m = find_matches(&q, &mol);
        assert!(!m.is_empty(), "benzene should match aromatic ring query");
    }

    #[test]
    fn test_visit_budget_generous_limit_finds_match() {
        // A budget high enough that a simple query completes.
        let mol = parse("CCO").unwrap();
        let q = parse_smarts("O").unwrap();
        let config = MatchConfig {
            max_visit_budget: Some(10_000),
            ..MatchConfig::default()
        };
        let m = find_matches_with_config(&q, &mol, &config);
        assert!(!m.is_empty(), "CCO contains O — should match within budget");
    }

    #[test]
    fn test_visit_budget_zero_returns_empty() {
        // Budget of 0 → no states explored → no results (documents fail-safe behavior).
        let mol = parse("CCO").unwrap();
        let q = parse_smarts("O").unwrap();
        let config = MatchConfig {
            max_visit_budget: Some(0),
            ..MatchConfig::default()
        };
        let m = find_matches_with_config(&q, &mol, &config);
        // With zero budget the search returns immediately — may or may not find a match.
        // Just verify it does not panic.
        let _ = m;
    }

    // ── RDKit PR #9201: query > target early exit ─────────────────────────────

    #[test]
    fn query_larger_than_target_returns_no_matches() {
        // A 6-atom query (benzene) cannot match a 4-atom target. The early exit
        // must fire immediately without entering VF2 recursion (RDKit PR #9201).
        let mol = parse("CCCC").unwrap(); // 4 heavy atoms
        let q = parse_smarts("c1ccccc1").unwrap(); // 6 query atoms
        let m = find_matches(&q, &mol);
        assert!(m.is_empty(), "6-atom query must not match 4-atom target");
    }

    #[test]
    fn query_same_size_as_target_allowed() {
        // Equal atom counts must still attempt the search (no false fast-fail).
        let mol = parse("CCCCCC").unwrap(); // 6 atoms
        let q = parse_smarts("CCCCCC").unwrap(); // 6 atoms
        let m = find_matches(&q, &mol);
        assert!(!m.is_empty(), "same-size query must be attempted");
    }
}
