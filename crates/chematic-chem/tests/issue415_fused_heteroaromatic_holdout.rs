//! New, not-previously-used holdout for issue #415's Phase 1 fix
//! (`chematic-chem/src/tautomer.rs`'s `atom_rank`-based tie-breaking,
//! replacing raw-`AtomIdx`-order candidate/ring selection in
//! `find_matches`/`find_exocyclic_lactam_shift_matches`/
//! `find_direct_aromatic_matches`).
//!
//! Six distinct, real fused/bridged heteroaromatic scaffolds (not variants
//! of the one molecule issue #415 itself was diagnosed on) exercising the
//! same class of defect: multiple structurally-equivalent-looking
//! donor/bridge/acceptor or ring candidates, where the tautomer search must
//! reach the same result regardless of the caller's atom ordering. For each
//! fixture: parse the canonical scaffold, generate several genuine atom-
//! order permutations via `MoleculeBuilder` (never a hand-written alternate
//! SMILES respelling -- that risks accidentally encoding a different
//! molecule; a permutation of an already-parsed, already-verified molecule
//! cannot), and assert `canonical_tautomer` agrees exactly across every
//! permutation and is idempotent on its own output.
//!
//! Deliberately not included: a purine/guanine/cytosine-shaped fixture
//! (carbonyl ring atom flanked by two ring nitrogens) -- that is exactly
//! Phase 2 (RFC section 1.7)'s own known, still-open residual
//! (`tp2_07_09_dual_flank_residual_documented_not_silently_fixed` in
//! `tautomer.rs`), unrelated to and not fixed by this round; including it
//! here would make this holdout red for a reason this PR doesn't address.
//! Also not included: standalone hydrazone/imine-enamine fixtures -- these
//! are single-candidate 1,3-shifts already covered by `test_15_shift_*` in
//! `tautomer.rs` and don't exercise the multi-candidate order-dependence
//! this fix targets.

use std::collections::HashMap;

use chematic_chem::canonical_tautomer;
use chematic_core::{AtomIdx, BondIdx, Molecule, MoleculeBuilder, STEREO_H_SENTINEL, StereoGroup};
use chematic_smiles::{canonical_smiles, parse};

/// Rebuild `mol` with atoms placed in `new_order` (`new_order[i]` is the
/// `AtomIdx` from `mol` that becomes atom `i` in the result). Same
/// value-substitution-in-place pattern as `chematic-chem`'s own
/// `hydrogen.rs::remove_hydrogens` remap, generalized to "atoms permuted,
/// none removed" -- see `tautomer.rs`'s own `#[cfg(test)]` copy of this
/// helper for the full rationale.
fn reorder_atoms(mol: &Molecule, new_order: &[AtomIdx]) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let mut atom_map: HashMap<AtomIdx, AtomIdx> = HashMap::with_capacity(new_order.len());
    for &old_idx in new_order {
        let new_idx = builder.add_atom(mol.atom(old_idx).clone());
        atom_map.insert(old_idx, new_idx);
    }

    let mut bond_map: HashMap<BondIdx, BondIdx> = HashMap::with_capacity(mol.bond_count());
    for i in 0..mol.bond_count() {
        let old_bidx = BondIdx(i as u32);
        let bond = mol.bond(old_bidx);
        let new_a = atom_map[&bond.atom1];
        let new_b = atom_map[&bond.atom2];
        let new_bidx = builder
            .add_bond(new_a, new_b, bond.order)
            .expect("reorder_atoms: bond from a valid molecule must be re-addable");
        bond_map.insert(old_bidx, new_bidx);
    }

    for (&old_bidx, &new_bidx) in &bond_map {
        if let Some(direction) = mol.bond_direction(old_bidx) {
            builder.set_bond_direction(new_bidx, direction);
        }
    }

    for &old_idx in new_order {
        let new_idx = atom_map[&old_idx];
        if let Some(order) = mol.stereo_neighbor_order(old_idx) {
            let remapped: Vec<u32> = order
                .iter()
                .map(|&v| {
                    if v == STEREO_H_SENTINEL {
                        STEREO_H_SENTINEL
                    } else {
                        atom_map[&AtomIdx(v)].0
                    }
                })
                .collect();
            builder.set_stereo_neighbor_order(new_idx, remapped);
        }
    }

    for g in mol.stereo_groups() {
        builder.add_stereo_group(StereoGroup::new(
            g.kind.clone(),
            g.atom_indices.iter().map(|&a| atom_map[&a]).collect(),
        ));
    }

    builder.build()
}

/// Assert `canonical_tautomer` agrees across `n` distinct atom-order
/// permutations of `smi` (identity, full reversal, and evenly-spaced
/// rotations), and that the result is idempotent under repeated
/// application, for every permutation.
fn assert_order_invariant_and_idempotent(name: &str, smi: &str) {
    let base = parse(smi).unwrap_or_else(|e| panic!("{name}: failed to parse '{smi}': {e}"));
    let n = base.atom_count();
    assert!(n >= 4, "{name}: fixture too small to meaningfully permute");

    let mut orderings: Vec<Vec<AtomIdx>> = vec![
        (0..n as u32).map(AtomIdx).collect(),
        (0..n as u32).rev().map(AtomIdx).collect(),
    ];
    for shift in [1usize, (n / 3).max(1), (2 * n / 3).max(2)] {
        let mut order: Vec<AtomIdx> = (0..n as u32).map(AtomIdx).collect();
        order.rotate_left(shift % n);
        orderings.push(order);
    }

    let mut expected: Option<String> = None;
    for (i, order) in orderings.iter().enumerate() {
        let permuted = reorder_atoms(&base, order);
        assert_eq!(
            canonical_smiles(&permuted),
            canonical_smiles(&base),
            "{name}: ordering {i}: reorder_atoms must not change the molecule's own identity"
        );

        let result = canonical_tautomer(&permuted);
        let reapplied = canonical_tautomer(&result);
        let result_smi = canonical_smiles(&result);
        assert_eq!(
            result_smi,
            canonical_smiles(&reapplied),
            "{name}: ordering {i}: canonical_tautomer's own output is not a fixed point of itself"
        );

        match &expected {
            None => expected = Some(result_smi),
            Some(exp) => assert_eq!(
                &result_smi, exp,
                "{name}: ordering {i} reached a different canonical tautomer than ordering 0 \
                 (expected '{exp}', got '{result_smi}')"
            ),
        }
    }
}

#[test]
fn quinolin_2_1h_one_fused_benzo_lactam() {
    // Fused benzo-pyridone lactam (carbostyril); mobile H on the ring
    // nitrogen, exocyclic carbonyl on the adjacent ring carbon.
    assert_order_invariant_and_idempotent("quinolin-2(1H)-one", "O=c1ccc2ccccc2[nH]1");
}

#[test]
fn pyridazin_3_2h_one() {
    // Adjacent ring N-N (not a single ring nitrogen like pyridone) with an
    // exocyclic carbonyl and a mobile H on one of the two nitrogens.
    assert_order_invariant_and_idempotent("pyridazin-3(2H)-one", "O=c1cccn[nH]1");
}

#[test]
fn phthalazin_1_2h_one_fused_benzo_pyridazinone() {
    // Fused benzo-pyridazinone: combines the adjacent-ring-N-N shape above
    // with ring fusion, structurally close to (but a distinct scaffold
    // from) issue #415's own diagnosis molecule.
    assert_order_invariant_and_idempotent("phthalazin-1(2H)-one", "O=c1[nH]ncc2ccccc21");
}

#[test]
fn amino_1h_1_2_4_triazole() {
    // Two non-equivalent ring nitrogens (the amino substituent breaks the
    // symmetry a bare 1,2,4-triazole would have), either of which can carry
    // the mobile H.
    assert_order_invariant_and_idempotent("3-amino-1H-1,2,4-triazole", "Nc1nc[nH]n1");
}

#[test]
fn nitro_1h_benzimidazole_fused_imidazole() {
    // Fused benzimidazole with a symmetry-breaking nitro substituent so the
    // two ring-imidazole nitrogens are chemically distinct, not automorphic.
    assert_order_invariant_and_idempotent(
        "5-nitro-1H-benzimidazole",
        "[O-][N+](=O)c1ccc2[nH]cnc2c1",
    );
}

#[test]
fn quinazolin_4_3h_one_fused_pyrimidinone() {
    // Fused benzo-pyrimidinone: a third distinct fused-lactam ring topology
    // (two ring nitrogens, only one adjacent to the exocyclic carbonyl).
    assert_order_invariant_and_idempotent("quinazolin-4(3H)-one", "O=c1[nH]cnc2ccccc12");
}
