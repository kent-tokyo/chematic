//! Structural property tests for [`CipDigraph`] -- Milestone 1's own acceptance bar
//! (see `docs/rfcs/cip_accurate_rfc.md`): does the digraph build correctly, deterministically,
//! and finitely. Whether it *ranks* substituents correctly is Milestone 2's concern and
//! is deliberately not tested here (there is no ranking logic in this crate yet).

use chematic_core::{AtomIdx, Molecule, MoleculeBuilder, STEREO_H_SENTINEL};
use chematic_smiles::parse;

use crate::budget::CipBudget;
use crate::digraph::CipDigraph;
use crate::node::CipNodeKind;

/// Rebuild `mol` with atoms inserted in the order given by `perm` (`perm[new_idx] =
/// old_idx`), remapping bond endpoints accordingly. Same shape as the existing
/// `permute_molecule` test helper in `crates/chematic-perception/src/sssr.rs`, extended
/// to also remap `stereo_neighbor_order` -- Milestone 1's helper didn't need this (ring
/// perception doesn't consult it), but `assign_cip_accurate_experimental` does, so a
/// renumbering-invariance test for R/S needs it carried over too, or every permuted
/// molecule silently loses its stereocenters (no `@`/`@@` interpretation is possible
/// without it). Returns the permuted molecule plus an `old_idx -> new_idx` map, so a
/// caller can locate where a specific atom ended up.
fn permute_molecule(mol: &Molecule, perm: &[usize]) -> (Molecule, Vec<u32>) {
    let mut old_to_new = vec![0u32; perm.len()];
    for (new_idx, &old_idx) in perm.iter().enumerate() {
        old_to_new[old_idx] = new_idx as u32;
    }
    let mut builder = MoleculeBuilder::new();
    for &old_idx in perm {
        builder.add_atom(mol.atom(AtomIdx(old_idx as u32)).clone());
    }
    for (_, bond) in mol.bonds() {
        let a = AtomIdx(old_to_new[bond.atom1.0 as usize]);
        let b = AtomIdx(old_to_new[bond.atom2.0 as usize]);
        let _ = builder.add_bond(a, b, bond.order);
    }
    for (old_idx, _) in mol.atoms() {
        if let Some(order) = mol.stereo_neighbor_order(old_idx) {
            let remapped: Vec<u32> = order
                .iter()
                .map(|&v| {
                    if v == STEREO_H_SENTINEL {
                        v
                    } else {
                        old_to_new[v as usize]
                    }
                })
                .collect();
            builder.set_stereo_neighbor_order(AtomIdx(old_to_new[old_idx.0 as usize]), remapped);
        }
    }
    (builder.build(), old_to_new)
}

fn find_atom_by_map(mol: &Molecule, map_num: u16) -> AtomIdx {
    mol.atoms()
        .find(|(_, a)| a.atom_map == Some(map_num))
        .map(|(idx, _)| idx)
        .expect("atom map tag not found")
}

#[test]
fn test_atom_renumbering_invariance() {
    // A molecule with a ring, a branch, and an ester -- enough shape to be a
    // meaningful check, not so much that hand-verifying the test itself is hard.
    let mol = parse("CC(=O)OC1CCCCC1").unwrap();
    let root_old = AtomIdx(1); // the carbonyl carbon
    let n = mol.atom_count();
    let perm: Vec<usize> = (0..n).rev().collect();
    let (permuted, old_to_new) = permute_molecule(&mol, &perm);
    let root_new = AtomIdx(old_to_new[root_old.0 as usize]);

    let mut g1 = CipDigraph::new(&mol, root_old, CipBudget::default_budget()).unwrap();
    let mut g2 = CipDigraph::new(&permuted, root_new, CipBudget::default_budget()).unwrap();
    g1.expand_all(g1.root()).unwrap();
    g2.expand_all(g2.root()).unwrap();

    // Cheap pre-check: same total node count and same multiset of node kinds, before
    // the real (order-invariant) signature check -- a more readable first failure.
    assert_eq!(g1.nodes().len(), g2.nodes().len());

    let sig1 = g1.branch_signature(g1.root()).unwrap();
    let sig2 = g2.branch_signature(g2.root()).unwrap();
    assert_eq!(
        sig1, sig2,
        "branch signature must not depend on atom numbering"
    );
}

#[test]
fn test_smiles_respelling_invariance() {
    // Same molecule (cyclohexyl methyl ketone), two different SMILES traversal orders;
    // the atom-map tag pins the "same" physical atom as root in both.
    let mol_a = parse("C1CCCCC1[C:1](=O)C").unwrap();
    let mol_b = parse("C[C:1](=O)C1CCCCC1").unwrap();
    let root_a = find_atom_by_map(&mol_a, 1);
    let root_b = find_atom_by_map(&mol_b, 1);

    let mut g1 = CipDigraph::new(&mol_a, root_a, CipBudget::default_budget()).unwrap();
    let mut g2 = CipDigraph::new(&mol_b, root_b, CipBudget::default_budget()).unwrap();
    g1.expand_all(g1.root()).unwrap();
    g2.expand_all(g2.root()).unwrap();

    assert_eq!(g1.nodes().len(), g2.nodes().len());

    let sig1 = g1.branch_signature(g1.root()).unwrap();
    let sig2 = g2.branch_signature(g2.root()).unwrap();
    assert_eq!(
        sig1, sig2,
        "branch signature must not depend on which SMILES respelling was parsed"
    );
}

#[test]
fn test_double_bond_duplication_structure() {
    // C=C, rooted at atom 0: departure-side (atom0's own list, iterating its
    // neighbors) contributes 1 duplicate of atom1; arrival-side (atom1's own list,
    // since its incoming edge was the double bond) contributes 1 duplicate of atom0.
    // 2 duplicates total -- both halves of the symmetric rule, not just one.
    let mol = parse("C=C").unwrap();
    let mut g = CipDigraph::new(&mol, AtomIdx(0), CipBudget::default_budget()).unwrap();
    g.expand_all(g.root()).unwrap();
    let dup_count = g
        .nodes()
        .iter()
        .filter(|n| matches!(n.kind, CipNodeKind::MultipleBondDuplicate { .. }))
        .count();
    assert_eq!(
        dup_count, 2,
        "a double bond must duplicate its partner into BOTH atoms' own lists"
    );
}

#[test]
fn test_triple_bond_duplication_structure() {
    // C#C: 2 duplicates on each side (k-1=2), 4 total.
    let mol = parse("C#C").unwrap();
    let mut g = CipDigraph::new(&mol, AtomIdx(0), CipBudget::default_budget()).unwrap();
    g.expand_all(g.root()).unwrap();
    let dup_count = g
        .nodes()
        .iter()
        .filter(|n| matches!(n.kind, CipNodeKind::MultipleBondDuplicate { .. }))
        .count();
    assert_eq!(
        dup_count, 4,
        "a triple bond must duplicate its partner twice into BOTH atoms' own lists"
    );
}

#[test]
fn test_ring_and_cage_termination() {
    // A simple ring, a fused bicyclic, and a cage (adamantane) -- must all expand
    // fully and finitely with no panic and no budget error under the default budget.
    for smi in [
        "C1CCCCC1",           // simple 6-ring
        "C1CC2CCC1CC2",       // fused bicyclic
        "C1C2CC3CC1CC(C2)C3", // adamantane (cage)
    ] {
        let mol = parse(smi).unwrap();
        let mut g = CipDigraph::new(&mol, AtomIdx(0), CipBudget::default_budget()).unwrap();
        g.expand_all(g.root())
            .expect("ring/cage system must expand fully and finitely");
        let ring_dups = g
            .nodes()
            .iter()
            .filter(|n| matches!(n.kind, CipNodeKind::RingDuplicate { .. }))
            .count();
        assert!(
            ring_dups > 0,
            "a ring system must produce at least one RingDuplicate: {smi}"
        );
    }
}

#[test]
fn test_determinism() {
    let mol = parse("COc1ccc2c3c1OC1C(O)C(CO)(CCCCc4ccccc4)CC4C(C2)N(C)CCC341").unwrap();
    let mut g1 = CipDigraph::new(&mol, AtomIdx(0), CipBudget::default_budget()).unwrap();
    let mut g2 = CipDigraph::new(&mol, AtomIdx(0), CipBudget::default_budget()).unwrap();
    g1.expand_all(g1.root()).unwrap();
    g2.expand_all(g2.root()).unwrap();
    let json1 = crate::debug::to_json(&g1);
    let json2 = crate::debug::to_json(&g2);
    assert_eq!(
        json1, json2,
        "identical input must produce byte-identical debug JSON"
    );
}

// --- Milestone 2: Rule 1a/1b/2 comparator tests -----------------------------------
//
// Expected R/S values below were independently verified against RDKit's modern
// rdCIPLabeler oracle (rdkit.Chem.rdCIPLabeler.AssignCIPLabels), not just
// hand-derived, before being written into these assertions.

use crate::assign::{SkipReason, assign_cip_accurate_experimental};
use crate::compare::{BranchComparison, CompareContext, compare_ligands, rank_children};
use chematic_core::CipCode;

fn assign_one(smi: &str) -> (chematic_core::AtomIdx, CipCode) {
    let mol = parse(smi).unwrap();
    let result = assign_cip_accurate_experimental(&mol, CipBudget::default_budget()).unwrap();
    assert_eq!(
        result.assignments.len(),
        1,
        "{smi}: expected exactly 1 assignment, got assignments={:?} skipped={:?}",
        result.assignments,
        result.skipped
    );
    result.assignments[0]
}

#[test]
fn test_rule_1a_atomic_number_ordering() {
    // Br(35) > Cl(17) > F(9) > H(1).
    let (_, code) = assign_one("[C@@H](F)(Cl)Br");
    assert_eq!(code, CipCode::R);
    let (_, code) = assign_one("[C@H](F)(Cl)Br");
    assert_eq!(code, CipCode::S);
}

#[test]
fn test_rule_1a_renumbering_invariance() {
    let mol = parse("[C@@H](F)(Cl)Br").unwrap();
    let n = mol.atom_count();
    let perm: Vec<usize> = (0..n).rev().collect();
    let (permuted, _) = permute_molecule(&mol, &perm);

    let a = assign_cip_accurate_experimental(&mol, CipBudget::default_budget()).unwrap();
    let b = assign_cip_accurate_experimental(&permuted, CipBudget::default_budget()).unwrap();
    assert_eq!(a.assignments.len(), 1);
    assert_eq!(b.assignments.len(), 1);
    assert_eq!(
        a.assignments[0].1, b.assignments[0].1,
        "R/S must not depend on atom numbering"
    );
}

#[test]
fn test_rule_1b_duplicate_resolves_via_1a_alone() {
    // CHO branch (real-O + duplicate-O at rank 2) vs CH2OH branch (real-O + H at rank
    // 2) -- Rule 1a alone decides this, per the worked trace in compare.rs's module
    // docs. Verified against RDKit: atom 2 is R.
    let (idx, code) = assign_one("OC[C@@H](O)C=O");
    assert_eq!(idx.0, 2);
    assert_eq!(code, CipCode::R);
}

#[test]
fn test_rule_1b_nitrile_duplicate_symmetry() {
    // Triple-bond N: 2 duplicates on both the C and N sides (arrival + departure).
    // Verified against RDKit: atom 2 is R.
    let (idx, code) = assign_one("N#C[C@@H](C)N");
    assert_eq!(idx.0, 2);
    assert_eq!(code, CipCode::R);
}

#[test]
fn test_ethynyl_vs_carboxymethyl_decided_at_shallow_sphere() {
    // The corpus's lone triple_bond_dup case, and the concrete regression that exposed
    // compare_ligands' original depth-first bug (see compare.rs's module docs): ethynyl
    // (-C#CH, own children [C,C,C]) vs carboxymethyl (-CH2COOH, own children [C,H,H])
    // must be decided by the *second* child (C beats H) without ever descending into
    // the carboxyl group's oxygens under the first (tied, C-vs-C) child. A depth-first
    // comparator reaches those oxygens first and gets this backwards. Verified against
    // RDKit's modern rdCIPLabeler: atom 2 is S.
    let (idx, code) = assign_one("C#C[C@H](CC(=O)O)NC(=O)c1cc2n(n1)CCN(CCC1CCNCC1)C2=O");
    assert_eq!(idx.0, 2);
    assert_eq!(code, CipCode::S);
}

#[test]
fn test_rule_2_deuterium_beats_hydrogen() {
    // D > H (isotope Some(2) > None). Verified against RDKit: S.
    let (_, code) = assign_one("[C@H]([2H])(C)O");
    assert_eq!(code, CipCode::S);
}

#[test]
fn test_rule_2_carbon_13_beats_carbon_12() {
    // 13C > 12C (isotope Some(13) > None). Verified against RDKit: R.
    let (_, code) = assign_one("[C@]([13CH3])(C)(F)Cl");
    assert_eq!(code, CipCode::R);
}

#[test]
fn test_rule_2_isotope_at_depth() {
    // Isotope label 3 bonds out from the stereocenter -- proves the recursive
    // comparator carries isotope comparison down, not just at rank-1. Two otherwise
    // structurally-identical branches, one ending in a labeled carbon.
    let mol = parse("[C@H](CC[13CH3])(CCC)F").unwrap();
    let result = assign_cip_accurate_experimental(&mol, CipBudget::default_budget()).unwrap();
    assert_eq!(
        result.assignments.len(),
        1,
        "isotope 3 bonds out must still break the tie between two otherwise-identical \
         propyl-shaped branches: assignments={:?} skipped={:?}",
        result.assignments,
        result.skipped
    );
}

#[test]
fn test_smiles_respelling_invariance_accurate() {
    let mol_a = parse("N[C@@H](C)C(=O)O").unwrap(); // L-alanine
    let mol_b = parse("OC(=O)[C@H](C)N").unwrap(); // same molecule, respelled
    let a = assign_cip_accurate_experimental(&mol_a, CipBudget::default_budget()).unwrap();
    let b = assign_cip_accurate_experimental(&mol_b, CipBudget::default_budget()).unwrap();
    assert_eq!(a.assignments.len(), 1);
    assert_eq!(b.assignments.len(), 1);
    assert_eq!(
        a.assignments[0].1, b.assignments[0].1,
        "R/S must not depend on which SMILES respelling was parsed"
    );
}

#[test]
fn test_compare_ligands_antisymmetry() {
    let mol = parse("[C@@H](F)(Cl)Br").unwrap();
    let mut g = CipDigraph::new(&mol, AtomIdx(0), CipBudget::default_budget()).unwrap();
    let children = g.expand_children(g.root()).unwrap();
    let mut ctx = CompareContext::new();
    let ab = compare_ligands(&mut g, children[0], children[1], &mut ctx).unwrap();
    let mut ctx2 = CompareContext::new();
    let ba = compare_ligands(&mut g, children[1], children[0], &mut ctx2).unwrap();
    let expected_inverse = match ab {
        BranchComparison::Higher => BranchComparison::Lower,
        BranchComparison::Lower => BranchComparison::Higher,
        other => other,
    };
    assert_eq!(
        ba, expected_inverse,
        "compare(a,b) and compare(b,a) must be inverses"
    );
}

#[test]
fn test_compare_ligands_reflexivity() {
    let mol = parse("[C@@H](F)(Cl)Br").unwrap();
    let mut g = CipDigraph::new(&mol, AtomIdx(0), CipBudget::default_budget()).unwrap();
    let children = g.expand_children(g.root()).unwrap();
    let mut ctx = CompareContext::new();
    let aa = compare_ligands(&mut g, children[0], children[0], &mut ctx).unwrap();
    assert_eq!(
        aa,
        BranchComparison::Equal,
        "compare(a,a) must always be Equal"
    );
}

/// Issue #107's pairwise memoization: a repeated `compare_ligands(a, b)` call within
/// the same `CompareContext` must be a cache hit, and that hit's value must be
/// byte-identical to what an independent, uncached (fresh-context) computation of the
/// exact same comparison produces -- proving the cache never returns a stale or
/// incorrect answer, not just that it returns *something* fast.
#[test]
fn test_pairwise_cache_hit_matches_fresh_computation() {
    let mol = parse("COc1ccc2c3c1OC1C(O)C(CO)(CCCCc4ccccc4)CC4C(C2)N(C)CCC341").unwrap();
    let mut g = CipDigraph::new(&mol, AtomIdx(0), CipBudget::default_budget()).unwrap();
    let children = g.expand_children(g.root()).unwrap();
    assert!(children.len() >= 2, "need at least 2 children to compare");
    let (a, b) = (children[0], children[1]);

    let mut ctx = CompareContext::new();
    let first = compare_ligands(&mut g, a, b, &mut ctx).unwrap();
    assert!(
        ctx.cache_misses >= 1,
        "first-ever call must populate at least one cache entry"
    );
    let hits_before = ctx.cache_hits;
    let second = compare_ligands(&mut g, a, b, &mut ctx).unwrap();
    assert_eq!(
        ctx.cache_hits,
        hits_before + 1,
        "identical repeated (a, b) call must be a cache hit"
    );
    assert_eq!(
        first, second,
        "a cache hit must match the original computation"
    );

    // Independent, uncached computation (fresh graph + fresh context) of the exact
    // same comparison must agree.
    let mut g2 = CipDigraph::new(&mol, AtomIdx(0), CipBudget::default_budget()).unwrap();
    let children2 = g2.expand_children(g2.root()).unwrap();
    let mut ctx2 = CompareContext::new();
    let independent = compare_ligands(&mut g2, children2[0], children2[1], &mut ctx2).unwrap();
    assert_eq!(
        first, independent,
        "cached and uncached computations of the same comparison must agree"
    );
}

/// The cache's "outcome-direction normalization": `compare_ligands(a, b)` and
/// `compare_ligands(b, a)` within the SAME context must share one cache entry (the
/// second call is a hit, not a miss) and the second call's returned value must be the
/// correctly-inverted opposite of the first, never the raw (un-inverted) cached value.
#[test]
fn test_pairwise_cache_respects_direction() {
    let mol = parse("[C@@H](F)(Cl)Br").unwrap();
    let mut g = CipDigraph::new(&mol, AtomIdx(0), CipBudget::default_budget()).unwrap();
    let children = g.expand_children(g.root()).unwrap();
    let (a, b) = (children[0], children[1]);

    let mut ctx = CompareContext::new();
    let ab = compare_ligands(&mut g, a, b, &mut ctx).unwrap();
    let hits_before = ctx.cache_hits;
    let ba = compare_ligands(&mut g, b, a, &mut ctx).unwrap();
    assert_eq!(
        ctx.cache_hits,
        hits_before + 1,
        "(b, a) must reuse (a, b)'s cache entry, not miss and recompute"
    );
    let expected_ba = match ab {
        BranchComparison::Higher => BranchComparison::Lower,
        BranchComparison::Lower => BranchComparison::Higher,
        other => other,
    };
    assert_eq!(
        ba, expected_ba,
        "a cached hit in the reversed direction must invert the stored outcome"
    );
}

#[test]
fn test_compare_ligands_determinism() {
    let mol = parse("COc1ccc2c3c1OC1C(O)C(CO)(CCCCc4ccccc4)CC4C(C2)N(C)CCC341").unwrap();
    let mut g1 = CipDigraph::new(&mol, AtomIdx(0), CipBudget::default_budget()).unwrap();
    let mut g2 = CipDigraph::new(&mol, AtomIdx(0), CipBudget::default_budget()).unwrap();
    let c1 = g1.expand_children(g1.root()).unwrap();
    let c2 = g2.expand_children(g2.root()).unwrap();
    let mut ctx1 = CompareContext::new();
    let mut ctx2 = CompareContext::new();
    let r1 = rank_children(&mut g1, &c1, &mut ctx1).unwrap();
    let r2 = rank_children(&mut g2, &c2, &mut ctx2).unwrap();
    let sizes1: Vec<usize> = r1.iter().map(|g| g.len()).collect();
    let sizes2: Vec<usize> = r2.iter().map(|g| g.len()).collect();
    assert_eq!(
        sizes1, sizes2,
        "identical input must rank into identically-shaped groups"
    );
}

#[test]
fn test_budget_increase_does_not_change_result() {
    let mol = parse("N[C@@H](C)C(=O)O").unwrap();
    let small_budget = CipBudget {
        max_nodes: 1_000,
        max_depth: 32,
        max_expansions: 1_000,
    };
    let a = assign_cip_accurate_experimental(&mol, small_budget).unwrap();
    let b = assign_cip_accurate_experimental(&mol, CipBudget::default_budget()).unwrap();
    assert_eq!(
        a.assignments, b.assignments,
        "a larger budget must not change the result for a molecule that fits comfortably in both"
    );
}

#[test]
fn test_no_forced_resolution_on_symmetric_substituents() {
    // A center with two IDENTICAL substituents (two methyls) cannot be resolved by
    // Rules 1a/1b/2 -- must be reported as Tied, never guessed.
    let mol = parse("C[C@](C)(F)Cl").unwrap();
    let result = assign_cip_accurate_experimental(&mol, CipBudget::default_budget()).unwrap();
    assert!(
        result.assignments.is_empty(),
        "identical substituents must not produce an assignment"
    );
    assert!(
        result
            .skipped
            .iter()
            .any(|(_, reason)| *reason == SkipReason::Tied),
        "identical substituents must be reported as Tied, not silently dropped: {:?}",
        result.skipped
    );
}

/// Compare the two real (non-hydrogen) children of a "hub" atom -- `smi`'s atom 0 must
/// have exactly 2 heavy-atom substituents (plus implicit hydrogens, which are filtered
/// out here since these tests care only about the two named branches).
fn compare_hub_branches(smi: &str) -> BranchComparison {
    let mol = parse(smi).unwrap();
    let mut g = CipDigraph::new(&mol, AtomIdx(0), CipBudget::default_budget()).unwrap();
    let root_children = g.expand_children(g.root()).unwrap();
    let real: Vec<_> = root_children
        .iter()
        .copied()
        .filter(|&n| matches!(g.node(n).kind, CipNodeKind::Atom { .. }))
        .collect();
    assert_eq!(
        real.len(),
        2,
        "{smi}: expected exactly 2 real (non-hydrogen) hub substituents, got {}",
        real.len()
    );
    let mut ctx = CompareContext::new();
    compare_ligands(&mut g, real[0], real[1], &mut ctx).unwrap()
}

#[test]
fn test_shallow_dominance_regardless_of_chain_depth() {
    // Generalizes the ethynyl-vs-carboxymethyl regression (see compare.rs's module
    // docs and the triple_bond_dup corpus case): a shallow difference must decide a
    // comparison even when an much deeper, otherwise-decisive difference sits under an
    // earlier-ranked sibling branch. Parameterized over chain length so this isn't
    // pinned to one specific depth.
    //
    // LEFT's first atom carries an extra methyl (ranks it above a plain chain
    // continuation at depth 2, since a branch with a real substituent beats one with
    // just hydrogens); RIGHT's first atom has no such branch. That depth-2 difference
    // must decide LEFT > RIGHT regardless of how long the remaining chain is, even
    // though RIGHT's chain ends in Br (which *would* outrank LEFT's terminal Cl if the
    // comparison ever reached that deep -- it must not).
    for chain_len in [1usize, 3, 5] {
        let chain = "C".repeat(chain_len);
        let left = format!("C(C){chain}Cl");
        let right = format!("C{chain}Br");
        let smi = format!("C({left})({right})");
        assert_eq!(
            compare_hub_branches(&smi),
            BranchComparison::Higher,
            "chain_len={chain_len}: {smi} -- shallow branch-vs-no-branch difference must \
             decide this regardless of the deeper Cl-vs-Br difference at the chain ends"
        );
    }
}

#[test]
fn test_phantom_padding_is_local_to_its_own_position() {
    // Regression test for a hazard the sphere-by-sphere fix introduced alongside its
    // own fix (see compare.rs's LevelSlot::Phantom docs): naively flattening multiple
    // ranked sibling positions into one combined list and comparing raw total lengths
    // corrupts *later* positions once an *earlier* position's substituent count
    // mismatches, because the later position's data silently shifts to fill the gap.
    //
    // X = hub(isopropyl, methyl), Y = hub(phenyl, ethyl). Isopropyl/phenyl both rank
    // above methyl/ethyl on their respective sides (a branch with 2 real substituents
    // beats one with fewer), so isopropyl-vs-phenyl is compared *first*. Phenyl (an
    // aromatic ipso carbon, fully substituted) has only 2 digraph children where
    // isopropyl has 3 (2 branches + H) -- exactly the child-count mismatch that needs
    // local phantom padding, not a global length fallback -- and correctly decides the
    // whole comparison in X's favor (a real hydrogen beats a phantom) *before* the
    // methyl-vs-ethyl position (which would favor Y, ethyl > methyl) is ever reached.
    //
    // A naive global-flatten-then-compare-total-length implementation gets this
    // backwards: it would misalign isopropyl's own trailing H against ethyl's leading
    // carbon (Y's phenyl only contributes 2 slots to the flattened list, so position 2
    // already belongs to Y's *second* branch), reporting C > H and picking Y.
    let smi = "C(C(C(C)C)C)(C(c1ccccc1)CC)";
    assert_eq!(
        compare_hub_branches(smi),
        BranchComparison::Higher,
        "{smi} -- isopropyl-vs-phenyl's local padding must decide this in X's favor \
         without being corrupted by (or itself corrupting) the methyl-vs-ethyl position"
    );
}

#[test]
fn rule5_resolves_the_two_verified_milestone_4a_target_rows() {
    use chematic_core::CipCode;

    let cases: &[(&str, u32)] = &[
        ("N[C@]1(C(=O)O)C[C@H](C(=O)O)[C@H](C(=O)O)C1", 1),
        (
            "CCCCc1cn([C@H]2[C@H](C)CCC[C@@H]2C)c(=O)n1Cc1ccc(-c2ccccc2-c2nn[nH]n2)nc1",
            7,
        ),
    ];

    for (smiles, atom_idx) in cases {
        let mol = chematic_smiles::parse(smiles).expect("valid SMILES");
        let assignment =
            crate::assign_cip_accurate_experimental(&mol, crate::CipBudget::default_budget())
                .expect("assignment succeeds");
        let code = assignment
            .assignments
            .iter()
            .find(|(idx, _)| idx.0 == *atom_idx)
            .map(|(_, code)| *code);
        assert_eq!(
            code,
            Some(CipCode::LowerR),
            "{smiles} atom {atom_idx}: expected lowercase r (RDKit oracle), got {code:?}"
        );
    }
}

/// Pseudoasymmetric r/s labels have a genuinely counter-intuitive property, verified
/// directly against the live RDKit oracle (not assumed from textbook analogy -- this
/// project's own standing lesson from Milestone 4B-1.5's discarded "S precedes R"
/// overfit): under a *global* molecular mirror (every `@`/`@@` in the molecule swapped,
/// including both this center's own tag and the two embedded reference centers that
/// make it pseudoasymmetric), the lowercase label is **invariant**, not covariant like
/// ordinary uppercase R/S. `rdCIPLabeler` confirms this directly on both Milestone-4A
/// target rows (checked, not derived by hand): `orig=r, mirr=r` for both. This is not a
/// corpus quirk -- swapping which of the two constitutionally-identical branches holds
/// the `R`/`S` auxiliary descriptor is itself the local mirror image, and it exactly
/// cancels the center's own tag flip. An earlier version of this test asserted the
/// *opposite* (mirrored expectation = the flipped label) on the assumption that CIP
/// labels always invert under mirroring; that assumption was wrong specifically for
/// pseudoasymmetric centers and was caught by checking the live oracle before trusting
/// it, not by argument.
#[test]
fn rule5_two_row_target_pseudoasymmetric_label_is_mirror_invariant() {
    use chematic_core::CipCode;
    let cases: &[(&str, u32)] = &[
        ("N[C@@]1(C(=O)O)C[C@@H](C(=O)O)[C@@H](C(=O)O)C1", 1),
        (
            "CCCCc1cn([C@@H]2[C@@H](C)CCC[C@H]2C)c(=O)n1Cc1ccc(-c2ccccc2-c2nn[nH]n2)nc1",
            7,
        ),
    ];
    for (smiles, atom_idx) in cases {
        let mol = chematic_smiles::parse(smiles).expect("valid SMILES");
        let assignment =
            crate::assign_cip_accurate_experimental(&mol, crate::CipBudget::default_budget())
                .expect("assignment succeeds");
        let code = assignment
            .assignments
            .iter()
            .find(|(idx, _)| idx.0 == *atom_idx)
            .map(|(_, code)| *code);
        assert_eq!(
            code,
            Some(CipCode::LowerR),
            "{smiles} atom {atom_idx} (mirrored): expected the SAME lowercase r \
             (RDKit oracle: mirroring is invariant for pseudoasymmetric centers, not \
             covariant), got {code:?}"
        );
    }
}

#[test]
fn diagnose_m4a0_quinic_residual_constitutional_identity() {
    // M4A-0: for each quinic/gallic-ester residual molecule, check whether the two
    // tied physical branches at the tied atom are genuinely constitutionally
    // isomorphic (branch_signature equal) -- if so, the tie is a real
    // stereo-dependent case (Rule 4/5 territory, not a comparator bug); if the
    // signatures differ, Rule 1a/1b *should* have broken the tie and didn't, which
    // is a bug, not a missing-rule situation. See advisor's caution in this
    // session: "skip:tied only means the current comparator returned Equal" is not
    // by itself proof of genuine constitutional identity.
    use crate::budget::CipBudget;
    use crate::compare::{CompareContext, rank_children};
    use crate::digraph::CipDigraph;
    use chematic_core::{AtomIdx, Chirality, STEREO_H_SENTINEL};
    use std::collections::HashSet;

    let cases: &[(&str, u32)] = &[
        (
            "O=C(O[C@H]1[C@H](O)C[C@](O)(C(=O)O)C[C@H]1O)c1cc(O)c(O)c(O)c1",
            3,
        ),
        (
            "O=C(O[C@H]1[C@H](O)C[C@](OC(=O)c2cc(O)c(O)c(O)c2)(C(=O)O)C[C@H]1O)c1cc(O)c(O)c(O)c1",
            3,
        ),
        (
            "O=C(O[C@@H]1C[C@](OC(=O)c2cc(O)c(O)c(O)c2)(C(=O)O)C[C@@H](OC(=O)c2cc(O)c(O)c(O)c2)[C@H]1OC(=O)c1cc(O)c(O)c(O)c1)c1cc(O)c(O)c(O)c1",
            5,
        ),
        (
            "O=C(Oc1c(O)cc(C(=O)O[C@H]2[C@H](OC(=O)c3cc(O)c(O)c(O)c3)C[C@](O)(C(=O)O)C[C@H]2OC(=O)c2cc(O)c(O)c(O)c2)cc1O)c1cc(O)c(O)c(O)c1",
            11,
        ),
    ];

    for (smi, atom_idx) in cases {
        let mol = chematic_smiles::parse(smi).expect("valid SMILES");
        let idx = AtomIdx(*atom_idx);
        let atom = mol.atom(idx);
        assert_ne!(atom.chirality, Chirality::None);
        let stereo_order = mol.stereo_neighbor_order(idx).expect("4 substituents");
        assert_eq!(stereo_order.len(), 4);

        let mut graph = CipDigraph::new(&mol, idx, CipBudget::default_budget()).unwrap();
        let root = graph.root();
        let root_children = graph.expand_children(root).unwrap();

        let position_nodes: Vec<_> = stereo_order
            .iter()
            .map(|&pos_val| {
                if pos_val == STEREO_H_SENTINEL {
                    root_children
                        .iter()
                        .copied()
                        .find(|&id| {
                            matches!(
                                graph.node(id).kind,
                                crate::node::CipNodeKind::ImplicitHydrogen
                            )
                        })
                        .unwrap()
                } else {
                    let a = AtomIdx(pos_val);
                    root_children
                        .iter()
                        .copied()
                        .find(|&id| {
                            matches!(graph.node(id).kind, crate::node::CipNodeKind::Atom { atom_idx } if atom_idx == a)
                        })
                        .unwrap()
                }
            })
            .collect();
        let position_set: HashSet<_> = position_nodes.iter().copied().collect();

        let mut ctx = CompareContext::new();
        let groups = rank_children(&mut graph, &root_children, &mut ctx).unwrap();

        let tied: Vec<_> = groups
            .iter()
            .find(|g| g.iter().filter(|n| position_set.contains(n)).count() > 1)
            .map(|g| {
                g.iter()
                    .copied()
                    .filter(|n| position_set.contains(n))
                    .collect::<Vec<_>>()
            })
            .expect("this atom is tied");
        assert_eq!(tied.len(), 2, "expected exactly a 2-way physical tie");

        let sig_a = graph.branch_signature(tied[0]).unwrap();
        let sig_b = graph.branch_signature(tied[1]).unwrap();
        println!(
            "{smi} atom {atom_idx}: branch_signature a={sig_a:#x} b={sig_b:#x} equal={}",
            sig_a == sig_b
        );
    }
}

/// Milestone 4A-2: the three-armed, locally-symmetric adamantane-cage pseudoasymmetric
/// family (`validation/cip_residual_classification_corpus.jsonl`, `engine=accurate`,
/// `bucket=pseudoasymmetric`, 15 rows across 5 distinct molecules -- every one of them,
/// not a sample). Each `(smiles, atom_idx, expected)` pair's `expected` is the RDKit
/// `rdCIPLabeler` oracle label (RDKit 2026.03.3), reproduced via
/// `scripts/cip_pseudoasymmetric_oracle.py`. Before this milestone's fix, every one of
/// these 15 atoms was `SkipReason::Tied` (a declined answer) -- see module docs for why
/// (the embedded reference used for the tiebreak was itself Pass-1/Rule-4b tied, so the
/// old molecular-descriptor `provisional` lookup could never find it).
#[test]
fn rule5_resolves_the_15_row_cage_family() {
    let cases: &[(&str, u32, chematic_core::CipCode)] = &[
        (
            "CCCCCCCC/C=C/CCCCCCCC(=O)OCc1cc(=O)c(OC(=O)[C@]23C[C@H]4C[C@H](C[C@H](C4)C2)C3)co1",
            31,
            chematic_core::CipCode::LowerS,
        ),
        (
            "CCCCCCCC/C=C/CCCCCCCC(=O)OCc1cc(=O)c(OC(=O)[C@]23C[C@H]4C[C@H](C[C@H](C4)C2)C3)co1",
            33,
            chematic_core::CipCode::LowerS,
        ),
        (
            "CCCCCCCC/C=C/CCCCCCCC(=O)OCc1cc(=O)c(OC(=O)[C@]23C[C@H]4C[C@H](C[C@H](C4)C2)C3)co1",
            35,
            chematic_core::CipCode::LowerS,
        ),
        (
            "COc1ccccc1N1CCN(CCCCNC(=O)C23C[C@H]4C[C@@H](C2)C[C@@H](C3)C4)CC1",
            21,
            chematic_core::CipCode::LowerS,
        ),
        (
            "COc1ccccc1N1CCN(CCCCNC(=O)C23C[C@H]4C[C@@H](C2)C[C@@H](C3)C4)CC1",
            23,
            chematic_core::CipCode::LowerS,
        ),
        (
            "COc1ccccc1N1CCN(CCCCNC(=O)C23C[C@H]4C[C@@H](C2)C[C@@H](C3)C4)CC1",
            26,
            chematic_core::CipCode::LowerS,
        ),
        (
            "O=C(NCCCCN1CCCC(/C=C\\c2ccccc2)C1)C12C[C@H]3C[C@@H](C1)C[C@@H](C2)C3",
            23,
            chematic_core::CipCode::LowerS,
        ),
        (
            "O=C(NCCCCN1CCCC(/C=C\\c2ccccc2)C1)C12C[C@H]3C[C@@H](C1)C[C@@H](C2)C3",
            25,
            chematic_core::CipCode::LowerS,
        ),
        (
            "O=C(NCCCCN1CCCC(/C=C\\c2ccccc2)C1)C12C[C@H]3C[C@@H](C1)C[C@@H](C2)C3",
            28,
            chematic_core::CipCode::LowerS,
        ),
        (
            "O=C(NCCCCN1CCN(c2cccc3ccccc23)CC1)C12C[C@H]3C[C@@H](C1)C[C@@H](C2)C3",
            25,
            chematic_core::CipCode::LowerS,
        ),
        (
            "O=C(NCCCCN1CCN(c2cccc3ccccc23)CC1)C12C[C@H]3C[C@@H](C1)C[C@@H](C2)C3",
            27,
            chematic_core::CipCode::LowerS,
        ),
        (
            "O=C(NCCCCN1CCN(c2cccc3ccccc23)CC1)C12C[C@H]3C[C@@H](C1)C[C@@H](C2)C3",
            30,
            chematic_core::CipCode::LowerS,
        ),
        (
            "O=c1cc(COC2CCOCC2)occ1OC(=O)[C@]12C[C@H]3C[C@H](C[C@H](C3)C1)C2",
            20,
            chematic_core::CipCode::LowerS,
        ),
        (
            "O=c1cc(COC2CCOCC2)occ1OC(=O)[C@]12C[C@H]3C[C@H](C[C@H](C3)C1)C2",
            22,
            chematic_core::CipCode::LowerS,
        ),
        (
            "O=c1cc(COC2CCOCC2)occ1OC(=O)[C@]12C[C@H]3C[C@H](C[C@H](C3)C1)C2",
            24,
            chematic_core::CipCode::LowerS,
        ),
    ];

    for (smi, atom_idx, expected) in cases {
        let mol = chematic_smiles::parse(smi).expect("valid SMILES");
        let assignment = crate::assign_cip_accurate_experimental(&mol, CipBudget::default_budget())
            .expect("assignment succeeds");
        let code = assignment
            .assignments
            .iter()
            .find(|(idx, _)| idx.0 == *atom_idx)
            .map(|(_, code)| *code);
        assert_eq!(
            code,
            Some(*expected),
            "{smi} atom {atom_idx}: expected {expected:?} (RDKit oracle), got {code:?}"
        );
    }
}

/// Companion control for the 15-row cage family above (see module docs, "Sign
/// convention: your corpus supplies ~1 bit" concern, and
/// `rule5_two_row_target_pseudoasymmetric_label_is_mirror_invariant`'s doc comment for
/// the underlying property this relies on). Every row in the corpus is oracle-labeled
/// lowercase `s`, which alone can't rule out a same-shape overfit as the one Milestone
/// 4B-1.5 found and discarded ("S precedes R" scored 8/8 on its own one-sided corpus
/// and 0/8 on the mirrored set) -- **but** for pseudoasymmetric centers specifically,
/// checked directly against the live RDKit oracle on all 15 rows (not assumed), a
/// global molecular mirror (every `@`/`@@` swapped) leaves the lowercase label
/// **unchanged** (`s` stays `s`), unlike ordinary R/S. This test's expectation is
/// therefore the *same* label as the original set, not its flip -- an oracle-verified
/// invariance check, not the mirror-antisymmetry check an uppercase R/S case would need.
#[test]
fn rule5_resolves_the_15_row_cage_family_mirrored() {
    let cases: &[(&str, u32, chematic_core::CipCode)] = &[
        (
            "CCCCCCCC/C=C/CCCCCCCC(=O)OCc1cc(=O)c(OC(=O)[C@@]23C[C@@H]4C[C@@H](C[C@@H](C4)C2)C3)co1",
            31,
            chematic_core::CipCode::LowerS,
        ),
        (
            "CCCCCCCC/C=C/CCCCCCCC(=O)OCc1cc(=O)c(OC(=O)[C@@]23C[C@@H]4C[C@@H](C[C@@H](C4)C2)C3)co1",
            33,
            chematic_core::CipCode::LowerS,
        ),
        (
            "CCCCCCCC/C=C/CCCCCCCC(=O)OCc1cc(=O)c(OC(=O)[C@@]23C[C@@H]4C[C@@H](C[C@@H](C4)C2)C3)co1",
            35,
            chematic_core::CipCode::LowerS,
        ),
        (
            "COc1ccccc1N1CCN(CCCCNC(=O)C23C[C@@H]4C[C@H](C2)C[C@H](C3)C4)CC1",
            21,
            chematic_core::CipCode::LowerS,
        ),
        (
            "COc1ccccc1N1CCN(CCCCNC(=O)C23C[C@@H]4C[C@H](C2)C[C@H](C3)C4)CC1",
            23,
            chematic_core::CipCode::LowerS,
        ),
        (
            "COc1ccccc1N1CCN(CCCCNC(=O)C23C[C@@H]4C[C@H](C2)C[C@H](C3)C4)CC1",
            26,
            chematic_core::CipCode::LowerS,
        ),
        (
            "O=C(NCCCCN1CCCC(/C=C\\c2ccccc2)C1)C12C[C@@H]3C[C@H](C1)C[C@H](C2)C3",
            23,
            chematic_core::CipCode::LowerS,
        ),
        (
            "O=C(NCCCCN1CCCC(/C=C\\c2ccccc2)C1)C12C[C@@H]3C[C@H](C1)C[C@H](C2)C3",
            25,
            chematic_core::CipCode::LowerS,
        ),
        (
            "O=C(NCCCCN1CCCC(/C=C\\c2ccccc2)C1)C12C[C@@H]3C[C@H](C1)C[C@H](C2)C3",
            28,
            chematic_core::CipCode::LowerS,
        ),
        (
            "O=C(NCCCCN1CCN(c2cccc3ccccc23)CC1)C12C[C@@H]3C[C@H](C1)C[C@H](C2)C3",
            25,
            chematic_core::CipCode::LowerS,
        ),
        (
            "O=C(NCCCCN1CCN(c2cccc3ccccc23)CC1)C12C[C@@H]3C[C@H](C1)C[C@H](C2)C3",
            27,
            chematic_core::CipCode::LowerS,
        ),
        (
            "O=C(NCCCCN1CCN(c2cccc3ccccc23)CC1)C12C[C@@H]3C[C@H](C1)C[C@H](C2)C3",
            30,
            chematic_core::CipCode::LowerS,
        ),
        (
            "O=c1cc(COC2CCOCC2)occ1OC(=O)[C@@]12C[C@@H]3C[C@@H](C[C@@H](C3)C1)C2",
            20,
            chematic_core::CipCode::LowerS,
        ),
        (
            "O=c1cc(COC2CCOCC2)occ1OC(=O)[C@@]12C[C@@H]3C[C@@H](C[C@@H](C3)C1)C2",
            22,
            chematic_core::CipCode::LowerS,
        ),
        (
            "O=c1cc(COC2CCOCC2)occ1OC(=O)[C@@]12C[C@@H]3C[C@@H](C[C@@H](C3)C1)C2",
            24,
            chematic_core::CipCode::LowerS,
        ),
    ];

    for (smi, atom_idx, expected) in cases {
        let mol = chematic_smiles::parse(smi).expect("valid SMILES");
        let assignment = crate::assign_cip_accurate_experimental(&mol, CipBudget::default_budget())
            .expect("assignment succeeds");
        let code = assignment
            .assignments
            .iter()
            .find(|(idx, _)| idx.0 == *atom_idx)
            .map(|(_, code)| *code);
        assert_eq!(
            code,
            Some(*expected),
            "{smi} atom {atom_idx}: expected {expected:?} (mirrored RDKit oracle), got {code:?}"
        );
    }
}

/// Determinism gate for the Milestone 4A-2 fix: `rank_children`'s equivalence classes
/// come from a `HashMap<usize, Vec<usize>>` (`groups_map` in `compare.rs`), and
/// `assign_one_with_rule5`'s `physical_in_tied[0]`/`[1]` pick is therefore sensitive, in
/// principle, to `mol.neighbors()`'s raw adjacency/bond-creation order -- exactly the
/// class of order-dependence this project has been burned by before (see this repo's
/// standing "never use atom index or HashMap iteration order as a tie-break" policy).
/// The fix's `if is_r_a { (pos_a, pos_b) } else { (pos_b, pos_a) }` choice is a genuine
/// chemical comparison (which branch's auxiliary sign is R), not an index tie-break, so
/// it should be renumbering-invariant by construction -- but that is an argument, not a
/// test, until checked here: worst-of-30 atom-renumbering permutations per molecule,
/// covering all 5 distinct cage molecules (15 target atoms total, 450 checks), using the
/// same `permute_molecule` helper (which also remaps `stereo_neighbor_order`, without
/// which no permuted molecule's `@`/`@@` tags could be reinterpreted at all) the existing
/// Rules-1a/2 renumbering tests use above.
#[test]
fn rule5_15_row_cage_family_is_renumbering_invariant_worst_of_30() {
    let cases: &[(&str, u32)] = &[
        (
            "CCCCCCCC/C=C/CCCCCCCC(=O)OCc1cc(=O)c(OC(=O)[C@]23C[C@H]4C[C@H](C[C@H](C4)C2)C3)co1",
            31,
        ),
        (
            "CCCCCCCC/C=C/CCCCCCCC(=O)OCc1cc(=O)c(OC(=O)[C@]23C[C@H]4C[C@H](C[C@H](C4)C2)C3)co1",
            33,
        ),
        (
            "CCCCCCCC/C=C/CCCCCCCC(=O)OCc1cc(=O)c(OC(=O)[C@]23C[C@H]4C[C@H](C[C@H](C4)C2)C3)co1",
            35,
        ),
        (
            "COc1ccccc1N1CCN(CCCCNC(=O)C23C[C@H]4C[C@@H](C2)C[C@@H](C3)C4)CC1",
            21,
        ),
        (
            "COc1ccccc1N1CCN(CCCCNC(=O)C23C[C@H]4C[C@@H](C2)C[C@@H](C3)C4)CC1",
            23,
        ),
        (
            "COc1ccccc1N1CCN(CCCCNC(=O)C23C[C@H]4C[C@@H](C2)C[C@@H](C3)C4)CC1",
            26,
        ),
        (
            "O=C(NCCCCN1CCCC(/C=C\\c2ccccc2)C1)C12C[C@H]3C[C@@H](C1)C[C@@H](C2)C3",
            23,
        ),
        (
            "O=C(NCCCCN1CCCC(/C=C\\c2ccccc2)C1)C12C[C@H]3C[C@@H](C1)C[C@@H](C2)C3",
            25,
        ),
        (
            "O=C(NCCCCN1CCCC(/C=C\\c2ccccc2)C1)C12C[C@H]3C[C@@H](C1)C[C@@H](C2)C3",
            28,
        ),
        (
            "O=C(NCCCCN1CCN(c2cccc3ccccc23)CC1)C12C[C@H]3C[C@@H](C1)C[C@@H](C2)C3",
            25,
        ),
        (
            "O=C(NCCCCN1CCN(c2cccc3ccccc23)CC1)C12C[C@H]3C[C@@H](C1)C[C@@H](C2)C3",
            27,
        ),
        (
            "O=C(NCCCCN1CCN(c2cccc3ccccc23)CC1)C12C[C@H]3C[C@@H](C1)C[C@@H](C2)C3",
            30,
        ),
        (
            "O=c1cc(COC2CCOCC2)occ1OC(=O)[C@]12C[C@H]3C[C@H](C[C@H](C3)C1)C2",
            20,
        ),
        (
            "O=c1cc(COC2CCOCC2)occ1OC(=O)[C@]12C[C@H]3C[C@H](C[C@H](C3)C1)C2",
            22,
        ),
        (
            "O=c1cc(COC2CCOCC2)occ1OC(=O)[C@]12C[C@H]3C[C@H](C[C@H](C3)C1)C2",
            24,
        ),
    ];

    // Small inline xorshift PRNG -- no `rand` dev-dependency exists in this crate and
    // adding one for a single deterministic-shuffle test isn't worth it (see the crate's
    // own `Cargo.toml`); a fixed seed keeps this test itself reproducible.
    fn next(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }
    fn shuffled(n: usize, state: &mut u64) -> Vec<usize> {
        let mut perm: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            let j = (next(state) as usize) % (i + 1);
            perm.swap(i, j);
        }
        perm
    }

    const PERMUTATIONS_PER_MOLECULE: usize = 30;
    let mut checked = 0usize;
    let mut seed: u64 = 0x9E3779B97F4A7C15;

    for (smi, atom_idx) in cases {
        let mol = chematic_smiles::parse(smi).expect("valid SMILES");
        let n = mol.atom_count();
        for trial in 0..PERMUTATIONS_PER_MOLECULE {
            let perm = shuffled(n, &mut seed);
            let (permuted, old_to_new) = permute_molecule(&mol, &perm);
            let new_idx = old_to_new[*atom_idx as usize];

            let assignment =
                crate::assign_cip_accurate_experimental(&permuted, CipBudget::default_budget())
                    .expect("assignment succeeds");
            let code = assignment
                .assignments
                .iter()
                .find(|(idx, _)| idx.0 == new_idx)
                .map(|(_, code)| *code);
            assert_eq!(
                code,
                Some(chematic_core::CipCode::LowerS),
                "{smi} original atom {atom_idx} (trial {trial}, new idx {new_idx}): \
                 expected LowerS under every renumbering, got {code:?}"
            );
            checked += 1;
        }
    }

    assert_eq!(
        checked,
        cases.len() * PERMUTATIONS_PER_MOLECULE,
        "sanity: every (molecule, permutation) pair must have been checked"
    );
    println!(
        "rule5 cage-family renumbering invariance: {checked}/{checked} identical \
         ({} molecule-atoms x {PERMUTATIONS_PER_MOLECULE} permutations)",
        cases.len()
    );
}

/// Companion determinism gate for the element-level guard in `assign_one_with_rule5`
/// (see `assign.rs` module docs, "Element-level guard: phosphorus stays tied"): the
/// guard is a plain `mol.atom(idx).element == Element::P` check on the *original*
/// atom identity, so it should stay `SkipReason::Tied` under every renumbering by
/// construction -- checked here the same way `rule5_15_row_cage_family_is_renumbering_invariant_worst_of_30`
/// checks the carbon cage family, rather than assumed. Confirms the 2 cyclophosphazene
/// phosphorus stereocenters from Milestone 4C-1 never flap to a resolved label on any
/// of 30 renumbering permutations.
#[test]
fn rule5_phosphorus_ties_stay_tied_across_renumbering_worst_of_30() {
    fn next(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }
    fn shuffled(n: usize, state: &mut u64) -> Vec<usize> {
        let mut perm: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            let j = (next(state) as usize) % (i + 1);
            perm.swap(i, j);
        }
        perm
    }

    const PERMUTATIONS: usize = 30;
    let smi = "CNP1(NC)=N[P@](NC)(N2CC2)=NP(NC)(NC)=N[P@@](NC)(N2CC2)=N1";
    let mol = chematic_smiles::parse(smi).expect("valid SMILES");
    let n = mol.atom_count();
    let mut seed: u64 = 0x9E3779B97F4A7C15;
    let mut checked = 0usize;

    for trial in 0..PERMUTATIONS {
        let perm = shuffled(n, &mut seed);
        let (permuted, old_to_new) = permute_molecule(&mol, &perm);

        let assignment =
            crate::assign_cip_accurate_experimental(&permuted, CipBudget::default_budget())
                .expect("assignment succeeds");

        for atom_idx in [6u32, 19u32] {
            let new_idx = old_to_new[atom_idx as usize];
            let resolved = assignment
                .assignments
                .iter()
                .any(|(idx, _)| idx.0 == new_idx);
            assert!(
                !resolved,
                "original atom {atom_idx} (trial {trial}, new idx {new_idx}): \
                 phosphorus must never resolve to a label"
            );
            let tied = assignment
                .skipped
                .iter()
                .any(|(idx, reason)| idx.0 == new_idx && *reason == SkipReason::Tied);
            assert!(
                tied,
                "original atom {atom_idx} (trial {trial}, new idx {new_idx}): \
                 expected SkipReason::Tied under every renumbering, got {:?}",
                assignment.skipped.iter().find(|(idx, _)| idx.0 == new_idx)
            );
            checked += 1;
        }
    }

    assert_eq!(
        checked,
        PERMUTATIONS * 2,
        "sanity: both atoms x every permutation checked"
    );
    println!(
        "phosphorus tied-stability under renumbering: {checked}/{checked} stably unresolved \
         (2 atoms x {PERMUTATIONS} permutations)"
    );
}

/// Regression for the square-planar-stereo PR's required CIP safety fix:
/// `assign_all`'s per-atom loop gated on `chirality == Chirality::None`
/// (equality), not an exhaustive match -- adding
/// `chematic_core::Chirality::SquarePlanar` did NOT force a compile error
/// there, so an `@SP1`-tagged 4-neighbor Pt center would have silently
/// reached the tetrahedral digraph-comparator machinery and produced a
/// bogus CIP code. Must be skipped (not assigned any code, tetrahedral or
/// otherwise), same as any other non-tetrahedral atom.
#[test]
fn square_planar_center_is_skipped_not_assigned_a_bogus_cip_code() {
    let mol = parse("N->[Pt@SP1](<-N)(Cl)Cl").unwrap();
    let pt = (0..mol.atom_count())
        .map(|i| AtomIdx(i as u32))
        .find(|&i| mol.atom(i).element == chematic_core::Element::PT)
        .expect("fixture has a Pt atom");
    assert!(
        matches!(
            mol.atom(pt).chirality,
            chematic_core::Chirality::SquarePlanar(_)
        ),
        "fixture must actually carry a SquarePlanar tag"
    );

    let result = assign_cip_accurate_experimental(&mol, CipBudget::default_budget()).unwrap();
    assert!(
        result.assignments.iter().all(|&(idx, _)| idx != pt),
        "square-planar center must never appear in assignments: {:?}",
        result.assignments
    );
}
