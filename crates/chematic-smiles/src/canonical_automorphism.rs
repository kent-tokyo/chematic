//! Exact colored-graph automorphism checking.
//!
//! `has_colored_automorphism_mapping` decides, via a real bijection search
//! (never a substructure/subgraph match, never a hash-based shortcut),
//! whether an automorphism of the current search node's colored graph
//! exists mapping `from` to `to`. Molecules handled by this crate are small
//! (even cage/cubane/coronene fixtures are well under 100 atoms), so a
//! pairwise exact backtracking search is used, per the task's own guidance
//! to prefer correctness over unproven speedups.
//!
//! Absolute invariant: a `false` result may cost performance (a missed
//! prune); a `true` result must always be a genuine automorphism. Every
//! candidate mapping is re-verified in full (`verify_full_bijection`) after
//! backtracking finds a complete assignment, independent of the incremental
//! `feasible` pruning used during search.

use chematic_core::AtomIdx;
use smallvec::{SmallVec, smallvec};

use crate::canonical_partition::{CanonicalColoredGraph, Partition};

/// Does a color-preserving, edge-color-preserving bijection of the whole
/// graph exist that maps `from` to `to` and keeps every already-singleton
/// (individualized) cell fixed?
///
/// `coloring` must already reflect the current search node's individualized
/// and refined state (see `canonical_search::search_canonical`). Atoms in a
/// singleton cell can only ever map to themselves, since the candidate
/// search below is restricted to same-cell atoms and a singleton's only
/// same-cell candidate is itself. That is exactly what "keep
/// already-individualized singletons fixed" means in practice.
/// Hard, always-on ceiling on total backtracking steps within one
/// [`has_colored_automorphism_mapping`] call (issue #421). `extend_mapping`'s
/// candidate-vertex backtracking has no depth bound of its own: `feasible`
/// only checks edges to already-assigned neighbors (a purely local,
/// one-hop-away consistency check), so on a molecule with several
/// simultaneously-unresolved large symmetric regions (e.g. 3+ near-identical
/// repeated substituent arms, all still non-singleton cells at the same
/// search node) it can explore a combinatorially large space. Observed to
/// still be running after 100+ seconds (never confirmed to terminate) on a
/// real 94-atom ChEMBL molecule reordered into `canonical_atom_order`'s own
/// output order -- the outer `SearchBudget` in `canonical_search.rs` only
/// counts *calls* to this function, not work done *inside* one call, so it
/// could not catch this.
///
/// Exceeding this ceiling makes [`extend_mapping`] return `false` (via
/// `steps` below), which is always a *safe* fallback per this module's own
/// documented invariant ("a false result may cost performance... a true
/// result must always be a genuine automorphism"): it can only ever cost a
/// missed prune (redundant-but-still-correct exploration one level up in
/// `canonical_search.rs`), never a wrong canonical answer. Chosen generously
/// relative to this crate's actual fixtures/corpus (the entire cage/cubane/
/// coronene test suite in this module needs a tiny fraction of this per
/// call) while still bounding wall-clock to a small fraction of a second even
/// in an unoptimized debug build.
const MAX_EXTEND_MAPPING_STEPS: usize = 200_000;

pub(crate) fn has_colored_automorphism_mapping(
    graph: &CanonicalColoredGraph,
    coloring: &Partition,
    from: AtomIdx,
    to: AtomIdx,
) -> bool {
    if from == to {
        return true;
    }
    if graph.vertex_color(from) != graph.vertex_color(to) {
        return false;
    }
    if coloring.cell_of[from.0 as usize] != coloring.cell_of[to.0 as usize] {
        return false;
    }

    let n = graph.n();
    let mut image: SmallVec<[Option<u32>; 64]> = smallvec![None; n];
    let mut used: SmallVec<[bool; 64]> = smallvec![false; n];
    image[from.0 as usize] = Some(to.0);
    used[to.0 as usize] = true;

    let mut steps = 0usize;
    if !extend_mapping(graph, coloring, &mut image, &mut used, &mut steps) {
        return false;
    }

    verify_full_bijection(graph, &image)
}

fn extend_mapping(
    graph: &CanonicalColoredGraph,
    coloring: &Partition,
    image: &mut [Option<u32>],
    used: &mut [bool],
    steps: &mut usize,
) -> bool {
    *steps += 1;
    if *steps > MAX_EXTEND_MAPPING_STEPS {
        return false;
    }
    let n = image.len();
    let Some(u) = (0..n).find(|&i| image[i].is_none()) else {
        return true;
    };
    let cell = coloring.cell_of[u];
    let u_color = graph.vertex_color(AtomIdx(u as u32));
    // Filter candidates by BOTH partition cell (respects the current search
    // node's individualization state -- e.g. keeps already-individualized
    // singletons fixed) AND raw vertex color directly. In production,
    // `coloring` always comes from `initial_partition` (which bakes vertex
    // color into the composite cell key), so the two filters agree and this
    // is redundant defense-in-depth. But this function must not silently
    // rely on that caller invariant: without the explicit color check here,
    // a coarser partition (e.g. a single-cell "ignore individualization"
    // partition, as used by several of this module's own unit tests) could
    // let the search commit to a color-mismatched candidate, which
    // `verify_full_bijection` would then reject -- and since this function
    // does not backtrack past a `verify_full_bijection` failure, that could
    // spuriously report `false` for a genuinely automorphic pair reachable
    // via a different candidate. Checking color directly here removes that
    // dependency entirely.
    let candidates: SmallVec<[u32; 8]> = (0..n as u32)
        .filter(|&v| {
            !used[v as usize]
                && coloring.cell_of[v as usize] == cell
                && graph.vertex_color(AtomIdx(v)) == u_color
        })
        .collect();

    for v in candidates {
        if !feasible(graph, image, u as u32, v) {
            continue;
        }
        image[u] = Some(v);
        used[v as usize] = true;
        if extend_mapping(graph, coloring, image, used, steps) {
            return true;
        }
        image[u] = None;
        used[v as usize] = false;
    }
    false
}

/// Incremental feasibility check for tentatively mapping `u -> v`: every
/// edge from `u` to an already-assigned neighbor must correspond to an
/// edge of the same color from `v` to that neighbor's image (and,
/// symmetrically, every edge from `v` to an already-assigned vertex's image
/// must correspond to an edge from `u` to that vertex) -- both directions,
/// so a missing/extra edge relative to any already-committed part of the
/// mapping is caught immediately rather than only at final verification.
fn feasible(graph: &CanonicalColoredGraph, image: &[Option<u32>], u: u32, v: u32) -> bool {
    let ua = AtomIdx(u);
    let va = AtomIdx(v);

    for (nb, bidx) in graph.neighbors(ua) {
        if let Some(mapped) = image[nb.0 as usize] {
            let want = graph.edge_color(ua, bidx);
            let has = graph
                .mol()
                .bond_between(va, AtomIdx(mapped))
                .map(|(bidx2, _)| graph.edge_color(va, bidx2));
            if has != Some(want) {
                return false;
            }
        }
    }

    for (nb, bidx) in graph.neighbors(va) {
        if let Some(preimage) = (0..image.len() as u32).find(|&i| image[i as usize] == Some(nb.0)) {
            let want = graph.edge_color(va, bidx);
            let has = graph
                .mol()
                .bond_between(ua, AtomIdx(preimage))
                .map(|(bidx2, _)| graph.edge_color(ua, bidx2));
            if has != Some(want) {
                return false;
            }
        }
    }

    true
}

/// Full, from-scratch re-verification of a complete candidate mapping:
/// bijectivity, vertex-color preservation, and edge existence + edge-color
/// preservation checked from every vertex's own perspective (so both
/// directions of every edge are independently checked). Never accepts a
/// partial or subgraph match -- every vertex of the graph must be mapped.
fn verify_full_bijection(graph: &CanonicalColoredGraph, image: &[Option<u32>]) -> bool {
    let n = image.len();
    let mut seen: SmallVec<[bool; 64]> = smallvec![false; n];
    for (i, &mapped) in image.iter().enumerate() {
        let Some(v) = mapped else {
            return false;
        };
        if seen[v as usize] {
            return false; // not injective
        }
        seen[v as usize] = true;
        if graph.vertex_color(AtomIdx(i as u32)) != graph.vertex_color(AtomIdx(v)) {
            return false;
        }
    }
    if seen.iter().any(|&s| !s) {
        return false; // not surjective
    }

    for i in 0..n {
        let img_i = image[i].expect("fully assigned");
        if graph.mol().degree(AtomIdx(i as u32)) != graph.mol().degree(AtomIdx(img_i)) {
            return false;
        }
        for (nb, bidx) in graph.neighbors(AtomIdx(i as u32)) {
            let img_nb = image[nb.0 as usize].expect("fully assigned");
            let want = graph.edge_color(AtomIdx(i as u32), bidx);
            match graph.mol().bond_between(AtomIdx(img_i), AtomIdx(img_nb)) {
                Some((bidx2, _)) if graph.edge_color(AtomIdx(img_i), bidx2) == want => {}
                _ => return false,
            }
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse;
    use chematic_core::{Chirality, Molecule, MoleculeBuilder};

    fn trivial_partition(graph: &CanonicalColoredGraph) -> Partition {
        // A partition with every atom in the same single cell (used by
        // tests that want the automorphism checker to consider the whole
        // vertex color / edge structure alone, not any current-search-node
        // individualization state).
        Partition {
            cell_of: vec![0; graph.n()],
        }
    }

    // --- Positive controls -------------------------------------------------

    #[test]
    fn benzene_ring_rotation_is_automorphic() {
        let mol = parse("c1ccccc1").unwrap();
        let graph = CanonicalColoredGraph::new(&mol);
        let part = trivial_partition(&graph);
        for i in 1..6 {
            assert!(
                has_colored_automorphism_mapping(&graph, &part, AtomIdx(0), AtomIdx(i)),
                "benzene atom 0 and atom {i} must be automorphic"
            );
        }
    }

    #[test]
    fn cf3_fluorines_are_automorphic() {
        // FC(F)(F)C -- the three fluorines of a CF3 group are a true-twin
        // orbit (identical neighbor set), the exact shape RENKIN's
        // regression corpus hits repeatedly (docs/rfcs/reaction_transform_perf.md).
        let mol = parse("FC(F)(F)C").unwrap();
        let graph = CanonicalColoredGraph::new(&mol);
        let part = trivial_partition(&graph);
        assert!(has_colored_automorphism_mapping(
            &graph,
            &part,
            AtomIdx(0),
            AtomIdx(2)
        ));
        assert!(has_colored_automorphism_mapping(
            &graph,
            &part,
            AtomIdx(0),
            AtomIdx(3)
        ));
    }

    // --- Negative controls (section 19: must be provably breakable) --------

    #[test]
    fn different_element_never_automorphic() {
        let mol = parse("FC(Cl)(Br)C").unwrap();
        let graph = CanonicalColoredGraph::new(&mol);
        let part = trivial_partition(&graph);
        assert!(!has_colored_automorphism_mapping(
            &graph,
            &part,
            AtomIdx(0),
            AtomIdx(2)
        ));
    }

    #[test]
    fn isotope_difference_never_automorphic() {
        // Two otherwise-identical fluorines, one isotope-labeled.
        let mol = parse("[19F]C(F)(F)C").unwrap();
        let graph = CanonicalColoredGraph::new(&mol);
        let part = trivial_partition(&graph);
        assert!(!has_colored_automorphism_mapping(
            &graph,
            &part,
            AtomIdx(0),
            AtomIdx(2)
        ));
    }

    #[test]
    fn charge_difference_never_automorphic() {
        let mut b = MoleculeBuilder::new();
        // Two terminal N atoms off a central C, one charged +1 (as NH3+),
        // one neutral -- otherwise same degree/element.
        let c = b.add_atom(chematic_core::Atom::organic(chematic_core::Element::C));
        let mut n1 = chematic_core::Atom::organic(chematic_core::Element::N);
        n1.charge = 1;
        let n1i = b.add_atom(n1);
        let n2 = b.add_atom(chematic_core::Atom::organic(chematic_core::Element::N));
        b.add_bond(c, n1i, chematic_core::BondOrder::Single)
            .unwrap();
        b.add_bond(c, n2, chematic_core::BondOrder::Single).unwrap();
        let mol: Molecule = b.build();
        let graph = CanonicalColoredGraph::new(&mol);
        let part = trivial_partition(&graph);
        assert!(!has_colored_automorphism_mapping(&graph, &part, n1i, n2));
    }

    #[test]
    fn stereo_difference_never_automorphic() {
        // A tetrahedral stereocenter is always stereo_unique -- never
        // automorphic with anything, even a chemically-identical-looking
        // atom, per this PR's conservative stereo-handling judgment call.
        let mol = parse("C[C@H](N)O").unwrap();
        let graph = CanonicalColoredGraph::new(&mol);
        let part = trivial_partition(&graph);
        // Atom 1 is the stereocenter; no other atom shares its color.
        for other in [0u32, 2, 3] {
            assert!(!has_colored_automorphism_mapping(
                &graph,
                &part,
                AtomIdx(1),
                AtomIdx(other)
            ));
        }
        assert_eq!(mol.atom(AtomIdx(1)).chirality, Chirality::CounterClockwise);
    }

    #[test]
    fn partial_mapping_never_accepted() {
        // Manually construct an `image` that is deliberately incomplete
        // (only maps 0 -> 1) and confirm verify_full_bijection rejects it.
        let mol = parse("c1ccccc1").unwrap();
        let graph = CanonicalColoredGraph::new(&mol);
        let mut image = vec![None; graph.n()];
        image[0] = Some(1);
        assert!(!verify_full_bijection(&graph, &image));
    }

    #[test]
    fn non_bijective_mapping_never_accepted() {
        // Two source atoms mapped to the SAME target atom.
        let mol = parse("c1ccccc1").unwrap();
        let graph = CanonicalColoredGraph::new(&mol);
        let mut image: Vec<Option<u32>> = (0..graph.n() as u32).map(Some).collect();
        image[1] = Some(0); // 0 -> 0 and 1 -> 0: not injective
        assert!(!verify_full_bijection(&graph, &image));
    }

    #[test]
    fn different_cell_never_automorphic_even_if_color_matches() {
        // Force a partition that pins atom 0 into its own singleton cell
        // (simulating "already individualized") -- atom 0 must then never
        // be considered automorphic with any other same-colored atom, even
        // though a trivial (single-cell) partition would allow it.
        let mol = parse("c1ccccc1").unwrap();
        let graph = CanonicalColoredGraph::new(&mol);
        let mut cell_of = vec![1u32; graph.n()];
        cell_of[0] = 0; // atom 0 alone in cell 0, everyone else in cell 1
        let part = Partition { cell_of };
        assert!(!has_colored_automorphism_mapping(
            &graph,
            &part,
            AtomIdx(0),
            AtomIdx(1)
        ));
    }

    // --- Exhaustive small-graph checks (section 13) -------------------------

    /// Build an unlabeled (all-same-color) cycle graph C_n as a Molecule
    /// (single bonds only -- color content is irrelevant to this check, only
    /// topology is exercised).
    fn cycle_molecule(n: usize) -> Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..n)
            .map(|_| b.add_atom(chematic_core::Atom::organic(chematic_core::Element::C)))
            .collect();
        for i in 0..n {
            b.add_bond(
                atoms[i],
                atoms[(i + 1) % n],
                chematic_core::BondOrder::Single,
            )
            .unwrap();
        }
        b.build()
    }

    #[test]
    fn cycles_up_to_8_are_fully_vertex_transitive() {
        for n in 3..=8 {
            let mol = cycle_molecule(n);
            let graph = CanonicalColoredGraph::new(&mol);
            let part = trivial_partition(&graph);
            for i in 0..n as u32 {
                for j in 0..n as u32 {
                    assert!(
                        has_colored_automorphism_mapping(&graph, &part, AtomIdx(i), AtomIdx(j)),
                        "C{n}: {i} and {j} must be automorphic (rotation)"
                    );
                }
            }
        }
    }

    fn complete_graph(n: usize) -> Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..n)
            .map(|_| b.add_atom(chematic_core::Atom::organic(chematic_core::Element::C)))
            .collect();
        for i in 0..n {
            for j in (i + 1)..n {
                b.add_bond(atoms[i], atoms[j], chematic_core::BondOrder::Single)
                    .unwrap();
            }
        }
        b.build()
    }

    #[test]
    fn complete_graphs_up_to_6_are_fully_vertex_transitive() {
        for n in 2..=6 {
            let mol = complete_graph(n);
            let graph = CanonicalColoredGraph::new(&mol);
            let part = trivial_partition(&graph);
            for i in 0..n as u32 {
                for j in 0..n as u32 {
                    assert!(has_colored_automorphism_mapping(
                        &graph,
                        &part,
                        AtomIdx(i),
                        AtomIdx(j)
                    ));
                }
            }
        }
    }

    fn complete_bipartite(a: usize, b: usize) -> Molecule {
        let mut bld = MoleculeBuilder::new();
        let left: Vec<_> = (0..a)
            .map(|_| bld.add_atom(chematic_core::Atom::organic(chematic_core::Element::C)))
            .collect();
        let right: Vec<_> = (0..b)
            .map(|_| bld.add_atom(chematic_core::Atom::organic(chematic_core::Element::N)))
            .collect();
        for &l in &left {
            for &r in &right {
                bld.add_bond(l, r, chematic_core::BondOrder::Single)
                    .unwrap();
            }
        }
        bld.build()
    }

    #[test]
    fn complete_bipartite_respects_the_two_sides() {
        let mol = complete_bipartite(3, 4);
        let graph = CanonicalColoredGraph::new(&mol);
        let part = trivial_partition(&graph);
        // Any two atoms on the same (differently-colored, C vs N) side are automorphic...
        assert!(has_colored_automorphism_mapping(
            &graph,
            &part,
            AtomIdx(0),
            AtomIdx(1)
        ));
        assert!(has_colored_automorphism_mapping(
            &graph,
            &part,
            AtomIdx(3),
            AtomIdx(4)
        ));
        // ...but never across sides (different element => different color).
        assert!(!has_colored_automorphism_mapping(
            &graph,
            &part,
            AtomIdx(0),
            AtomIdx(3)
        ));
    }

    /// The section-13 positive-control witness: a graph where 1-WL /
    /// Morgan-style refinement puts two vertices in the same cell, but they
    /// are in *different* automorphism orbits -- i.e. a graph that is
    /// regular (so plain degree/neighbor-hash refinement never splits it)
    /// but NOT vertex-transitive.
    ///
    /// The disjoint union of a triangle (C3) and a square (C4) is exactly
    /// this: every vertex is degree 2 with degree-2 neighbors, so 1-WL /
    /// Morgan-rank refinement (which only ever looks at local neighbor
    /// colors, never global path length or component identity) can never
    /// split the two components apart -- this is the textbook fact that
    /// 1-WL cannot distinguish regular graphs of different girth/size. A
    /// naive Morgan-rank-only pruning scheme would therefore treat all 7
    /// vertices as one orbit and wrongly merge e.g. a triangle vertex with a
    /// square vertex. But no automorphism can map a component of order 3 to
    /// a component of order 4 (an automorphism must map connected
    /// components to connected components of the same order); the true
    /// orbit structure is `{0,1,2}` (triangle) and `{3,4,5,6}` (square).
    ///
    /// (An earlier draft of this test used two disjoint *triangles*
    /// instead -- that premise was actually wrong: `Aut(C3 |_| C3)` DOES
    /// include a whole-component swap, since the two components are
    /// isomorphic, so it IS transitive on all 6 vertices. C3 |_| C4 has no
    /// such swap since the components differ in size, which is exactly why
    /// it works as a same-WL-cell/different-orbit witness and two
    /// isomorphic components do not.)
    #[test]
    fn triangle_and_square_same_wl_cell_different_orbits() {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..7)
            .map(|_| b.add_atom(chematic_core::Atom::organic(chematic_core::Element::C)))
            .collect();
        // Triangle: 0-1-2-0. Square: 3-4-5-6-3.
        for &(i, j) in &[(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 6), (6, 3)] {
            b.add_bond(atoms[i], atoms[j], chematic_core::BondOrder::Single)
                .unwrap();
        }
        let mol = b.build();

        // Verify the WL/Morgan-rank premise empirically (don't assert what
        // we want -- read it off the actual implementation): all 7 atoms
        // must tie under plain Morgan refinement.
        let ranks = crate::canonical::morgan_ranks(&mol);
        assert_eq!(
            ranks.iter().collect::<std::collections::HashSet<_>>().len(),
            1,
            "triangle and square must tie under plain Morgan refinement (both 2-regular)"
        );

        let graph = CanonicalColoredGraph::new(&mol);
        let part = trivial_partition(&graph);

        // Within-component: automorphic.
        assert!(has_colored_automorphism_mapping(
            &graph,
            &part,
            AtomIdx(0),
            AtomIdx(1)
        ));
        assert!(has_colored_automorphism_mapping(
            &graph,
            &part,
            AtomIdx(3),
            AtomIdx(5)
        ));
        // Across components: NOT automorphic, despite the same WL cell --
        // this is exactly what a rank-hash-only pruning scheme would get
        // wrong (it would wrongly prune away the square's individualization
        // entirely, believing it redundant with the triangle's).
        for i in 0..3u32 {
            for j in 3..7u32 {
                assert!(
                    !has_colored_automorphism_mapping(&graph, &part, AtomIdx(i), AtomIdx(j)),
                    "atom {i} (triangle) and atom {j} (square) must NOT be automorphic"
                );
            }
        }
    }

    // --- Exhaustive/randomized small-graph checks against a brute-force ----
    // --- ground truth (section 13) -----------------------------------------

    /// Ground truth: is there ANY permutation of `0..n` that fixes color
    /// classes, maps `a -> b`, and preserves edge presence + edge color?
    /// O(n!) -- only used in tests, only for tiny `n`.
    fn brute_force_automorphic(
        n: usize,
        edges: &std::collections::HashMap<(usize, usize), u8>,
        colors: &[u32],
        a: usize,
        b: usize,
    ) -> bool {
        if colors[a] != colors[b] {
            return false;
        }
        let mut perm: Vec<usize> = (0..n).collect();

        fn edge_color(
            edges: &std::collections::HashMap<(usize, usize), u8>,
            x: usize,
            y: usize,
        ) -> Option<u8> {
            let k = if x <= y { (x, y) } else { (y, x) };
            edges.get(&k).copied()
        }

        fn permutations<F: FnMut(&[usize]) -> bool>(
            perm: &mut Vec<usize>,
            k: usize,
            f: &mut F,
        ) -> bool {
            if k == perm.len() {
                return f(perm);
            }
            for i in k..perm.len() {
                perm.swap(k, i);
                if permutations(perm, k + 1, f) {
                    perm.swap(k, i);
                    return true;
                }
                perm.swap(k, i);
            }
            false
        }

        permutations(&mut perm, 0, &mut |p: &[usize]| {
            if p[a] != b {
                return false;
            }
            for i in 0..n {
                if colors[i] != colors[p[i]] {
                    return false;
                }
            }
            for x in 0..n {
                for y in (x + 1)..n {
                    if edge_color(edges, x, y) != edge_color(edges, p[x], p[y]) {
                        return false;
                    }
                }
            }
            true
        })
    }

    /// Build a `Molecule` from an explicit edge set (keyed by (min,max) pair
    /// -> an edge-color tag 0/1, mapped to Single/Double) and a per-atom
    /// color tag (mapped to a small distinct-element palette), for exact
    /// cross-checking against `brute_force_automorphic`.
    fn build_colored_graph(
        n: usize,
        edges: &std::collections::HashMap<(usize, usize), u8>,
        colors: &[u32],
    ) -> Molecule {
        let palette = [
            chematic_core::Element::C,
            chematic_core::Element::N,
            chematic_core::Element::O,
            chematic_core::Element::F,
        ];
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..n)
            .map(|i| {
                b.add_atom(chematic_core::Atom::organic(
                    palette[colors[i] as usize % palette.len()],
                ))
            })
            .collect();
        for x in 0..n {
            for y in (x + 1)..n {
                if let Some(&ec) = edges.get(&(x, y)) {
                    let order = if ec == 0 {
                        chematic_core::BondOrder::Single
                    } else {
                        chematic_core::BondOrder::Double
                    };
                    b.add_bond(atoms[x], atoms[y], order).unwrap();
                }
            }
        }
        b.build()
    }

    /// Exhaustive: every simple (uncolored, unweighted) graph on `n <= 5`
    /// vertices (all `2^C(n,2)` edge subsets), every vertex pair, cross-
    /// checked against brute-force ground truth. `n=6` would need
    /// `2^15 * 15` permutation-searches (each up to `6!`), too slow for a
    /// unit test at debug-build speed -- covered instead by the structured
    /// n<=8 cases above (cycles/complete/bipartite) plus the randomized
    /// n<=8 fuzz test below, per this task's "wherever feasible" guidance.
    #[test]
    fn exhaustive_simple_graphs_n_le_5_match_brute_force() {
        for n in 1..=5usize {
            let pairs: Vec<(usize, usize)> = (0..n)
                .flat_map(|x| ((x + 1)..n).map(move |y| (x, y)))
                .collect();
            let m = pairs.len();
            for mask in 0u32..(1u32 << m) {
                let mut edges = std::collections::HashMap::new();
                for (i, &(x, y)) in pairs.iter().enumerate() {
                    if mask & (1 << i) != 0 {
                        edges.insert((x, y), 0u8);
                    }
                }
                let colors = vec![0u32; n]; // uncolored: all same color class
                let mol = build_colored_graph(n, &edges, &colors);
                let graph = CanonicalColoredGraph::new(&mol);
                let part = trivial_partition(&graph);
                for a in 0..n {
                    for bb in 0..n {
                        let expected = brute_force_automorphic(n, &edges, &colors, a, bb);
                        let got = has_colored_automorphism_mapping(
                            &graph,
                            &part,
                            AtomIdx(a as u32),
                            AtomIdx(bb as u32),
                        );
                        assert_eq!(
                            got, expected,
                            "n={n} mask={mask:b} a={a} b={bb}: got {got}, brute force says {expected}"
                        );
                    }
                }
            }
        }
    }

    /// Tiny deterministic xorshift PRNG -- avoids adding a `rand`
    /// dependency (and the `Cargo.lock` churn that would bring) for what a
    /// few lines of arithmetic can do just as well for a seeded fuzz test.
    struct Xorshift64(u64);
    impl Xorshift64 {
        fn next_u64(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x
        }
        fn next_range(&mut self, bound: usize) -> usize {
            (self.next_u64() % bound as u64) as usize
        }
    }

    /// Randomized: colored graphs AND edge-colored graphs, `n` in `6..=8`,
    /// disconnected graphs included (edge density low enough that some
    /// generated graphs are disconnected), cross-checked against brute
    /// force. Seeded PRNG for reproducibility.
    #[test]
    fn random_colored_edge_colored_graphs_n_6_to_8_match_brute_force() {
        let mut rng = Xorshift64(0x9E3779B97F4A7C15);
        for n in 6..=8usize {
            for _trial in 0..40 {
                let pairs: Vec<(usize, usize)> = (0..n)
                    .flat_map(|x| ((x + 1)..n).map(move |y| (x, y)))
                    .collect();
                let mut edges = std::collections::HashMap::new();
                for &(x, y) in &pairs {
                    if rng.next_range(10) < 4 {
                        // ~40% edge density -- some trials will be disconnected.
                        let ec = (rng.next_range(2)) as u8; // edge color 0 or 1 (Single/Double)
                        edges.insert((x, y), ec);
                    }
                }
                let num_colors = 1 + rng.next_range(3); // 1..=3 vertex color classes
                let colors: Vec<u32> = (0..n).map(|_| rng.next_range(num_colors) as u32).collect();

                let mol = build_colored_graph(n, &edges, &colors);
                let graph = CanonicalColoredGraph::new(&mol);
                let part = trivial_partition(&graph);
                for a in 0..n {
                    for bb in 0..n {
                        let expected = brute_force_automorphic(n, &edges, &colors, a, bb);
                        let got = has_colored_automorphism_mapping(
                            &graph,
                            &part,
                            AtomIdx(a as u32),
                            AtomIdx(bb as u32),
                        );
                        assert_eq!(
                            got, expected,
                            "n={n} trial={_trial} edges={edges:?} colors={colors:?} a={a} b={bb}: \
                             got {got}, brute force says {expected}"
                        );
                    }
                }
            }
        }
    }

    /// Explicitly-disconnected graphs (two separate components of unequal
    /// size, mirroring the triangle+square witness above but with random
    /// internal edges) -- brute-force cross-check.
    #[test]
    fn disconnected_graphs_match_brute_force() {
        // Component A: 3 vertices, a path 0-1-2. Component B: 4 vertices, a
        // path 3-4-5-6. No automorphism should ever cross components (of
        // different order), matching the earlier structural witness test.
        let mut edges = std::collections::HashMap::new();
        edges.insert((0, 1), 0u8);
        edges.insert((1, 2), 0u8);
        edges.insert((3, 4), 0u8);
        edges.insert((4, 5), 0u8);
        edges.insert((5, 6), 0u8);
        let colors = vec![0u32; 7];
        let mol = build_colored_graph(7, &edges, &colors);
        let graph = CanonicalColoredGraph::new(&mol);
        let part = trivial_partition(&graph);
        for a in 0..7 {
            for b in 0..7 {
                let expected = brute_force_automorphic(7, &edges, &colors, a, b);
                let got = has_colored_automorphism_mapping(
                    &graph,
                    &part,
                    AtomIdx(a as u32),
                    AtomIdx(b as u32),
                );
                assert_eq!(got, expected, "a={a} b={b}");
            }
        }
    }
}
