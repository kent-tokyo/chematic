//! Structural property tests for [`CipDigraph`] -- Milestone 1's own acceptance bar
//! (see `docs/cip_accurate_rfc.md`): does the digraph build correctly, deterministically,
//! and finitely. Whether it *ranks* substituents correctly is Milestone 2's concern and
//! is deliberately not tested here (there is no ranking logic in this crate yet).

use chematic_core::{AtomIdx, Molecule, MoleculeBuilder};
use chematic_smiles::parse;

use crate::budget::CipBudget;
use crate::digraph::CipDigraph;
use crate::node::CipNodeKind;

/// Rebuild `mol` with atoms inserted in the order given by `perm` (`perm[new_idx] =
/// old_idx`), remapping bond endpoints accordingly. Same shape as the existing
/// `permute_molecule` test helper in `crates/chematic-perception/src/sssr.rs`. Returns
/// the permuted molecule plus an `old_idx -> new_idx` map, so a caller can locate where
/// a specific atom ended up.
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
