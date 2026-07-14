//! Structural property tests for [`CipDigraph`] -- Milestone 1's own acceptance bar
//! (see `docs/cip_accurate_rfc.md`): does the digraph build correctly, deterministically,
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
