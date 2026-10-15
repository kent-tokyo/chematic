//! Storage-only normalization oracle frozen from the immediate baseline.
//! All chemical inputs come from the existing licensed and pinned 10k corpus.
use super::*;

fn reference_finish(
    mol: &Molecule,
    atom_aromatic: &FxHashSet<AtomIdx>,
    bond_orders: &FxHashMap<BondIdx, BondOrder>,
) -> Molecule {
    use chematic_core::{MoleculeBuilder, implicit_hcount};

    // Implicit-H counts computed BEFORE bond orders are normalized below, for
    // organic-subset atoms without an explicit bracket H count. Needed because
    // normalizing every aromatic-model bond to `BondOrder::Aromatic` (below)
    // discards the Kekule Single/Double pattern that distinguishes a
    // lone-pair-donating "pyrrole-type" heteroatom (2 ring single bonds pre-
    // normalization, needs 1 implicit H) from a "pyridine-type" one (1 ring
    // single + 1 ring double, needs 0) -- post-normalization both look
    // identical (aromatic, 2 aromatic-order ring bonds, no substituent), so
    // `implicit_hcount`'s aromatic-path heuristic (correct for SMILES that
    // was aromatic-written from the start, per OpenSMILES convention: bare
    // aromatic `n` is pyridine-type, pyrrole-type is always `[nH]`) silently
    // returns the wrong value for atoms that reach this function via
    // Kekule-then-perceive instead. This under-counts molecular weight and
    // formula, not just fingerprints/canonical SMILES.
    let pre_h: Vec<Option<u8>> = mol
        .atoms()
        .map(|(idx, atom)| {
            if atom.hydrogen_count.is_some() {
                None // already explicit; nothing to preserve
            } else {
                Some(implicit_hcount(mol, idx))
            }
        })
        .collect();

    let mut builder = MoleculeBuilder::new();
    for (idx, atom) in mol.atoms() {
        let mut a = atom.clone();
        a.aromatic = atom_aromatic.contains(&idx);
        builder.add_atom(a);
    }
    let mut stash = Vec::new();
    for (bidx, bond) in mol.bonds() {
        let order = bond_orders[&bidx];
        if let Ok(new_bidx) = builder.add_bond(bond.atom1, bond.atom2, order)
            && order == BondOrder::Aromatic
            && matches!(bond.order, BondOrder::Up | BondOrder::Down)
        {
            // Kekule input promoted to Aromatic here loses its E/Z direction
            // the same way the SMILES parser's aromatic-aromatic coercion
            // does — stash it so an exocyclic double bond anchored on this
            // ring bond still round-trips through the canonical writer.
            stash.push((new_bidx, bond.order, bond.atom1));
        }
    }
    // Atoms/bonds above are re-added in `mol`'s own enumeration order with
    // none skipped, so indices line up 1:1 — safe to copy side-channel
    // metadata wholesale. (This rebuild previously dropped stereo_groups and
    // stereo_neighbor_order silently; closing that here too.)
    builder.copy_r_groups_from(mol);
    builder.copy_atom_tags_from(mol);
    builder.copy_stereo_groups_from(mol);
    builder.copy_stereo_from(mol);
    builder.copy_bond_directions_from(mol);
    // After the wholesale copy, which would otherwise erase them.
    for (bidx, dir, anchor) in stash {
        builder.set_bond_direction(bidx, dir);
        builder.set_bond_direction_anchor(bidx, anchor);
    }
    let normalized = builder.build();

    // Compare the pre-normalization implicit H against what the same
    // (already-tested, unmodified) `implicit_hcount` computes on the
    // normalized bonds; only atoms where normalization actually changed the
    // answer get an explicit H frozen in. Benzene CH and pyridine-type N
    // (heuristic already agrees) are left untouched -- no spurious bracket
    // notation for atoms that didn't need it.
    let needs_patch: Vec<(chematic_core::AtomIdx, u8)> = normalized
        .atoms()
        .filter_map(|(idx, _)| {
            let pre = pre_h[idx.0 as usize]?;
            let post = implicit_hcount(&normalized, idx);
            (pre != post).then_some((idx, pre))
        })
        .collect();
    if needs_patch.is_empty() {
        return normalized;
    }

    let mut patched = MoleculeBuilder::new();
    for (idx, atom) in normalized.atoms() {
        let mut a = atom.clone();
        if let Some(&(_, h)) = needs_patch.iter().find(|(pidx, _)| *pidx == idx) {
            a.hydrogen_count = Some(h);
        }
        patched.add_atom(a);
    }
    for (_bond_idx, bond) in normalized.bonds() {
        let _ = patched.add_bond(bond.atom1, bond.atom2, bond.order);
    }
    patched.copy_stereo_groups_from(&normalized);
    patched.copy_stereo_from(&normalized);
    patched.copy_atom_tags_from(&normalized);
    patched.copy_bond_directions_from(&normalized);
    patched.build()
}

fn compare_graph(expected: &Molecule, actual: &Molecule) {
    assert_eq!(actual.atom_count(), expected.atom_count());
    assert_eq!(actual.bond_count(), expected.bond_count());
    assert_eq!(actual.stereo_groups(), expected.stereo_groups());
    assert_eq!(actual.rdkit_bond_order(), expected.rdkit_bond_order());
    for group in 0..expected.stereo_groups().len() {
        assert_eq!(
            actual.stereo_group_bond_indices(group),
            expected.stereo_group_bond_indices(group)
        );
    }
    for (idx, atom) in expected.atoms() {
        assert_eq!(actual.atom(idx), atom);
        assert_eq!(actual.atom_tag(idx), expected.atom_tag(idx));
        assert_eq!(actual.r_group_label(idx), expected.r_group_label(idx));
        assert_eq!(
            actual.neighbors(idx).collect::<Vec<_>>(),
            expected.neighbors(idx).collect::<Vec<_>>()
        );
        assert_eq!(
            actual.stereo_neighbor_order(idx),
            expected.stereo_neighbor_order(idx)
        );
        assert_eq!(
            chematic_core::implicit_hcount(actual, idx),
            chematic_core::implicit_hcount(expected, idx)
        );
    }
    for (idx, bond) in expected.bonds() {
        assert_eq!(actual.bond(idx), bond);
        assert_eq!(actual.bond_direction(idx), expected.bond_direction(idx));
        assert_eq!(
            actual.bond_direction_anchor(idx),
            expected.bond_direction_anchor(idx)
        );
        assert_eq!(
            actual.smiles_ring_closure_metadata(idx),
            expected.smiles_ring_closure_metadata(idx)
        );
    }
}

fn compare_normalization(input: &Molecule, model: &AromaticityModel) -> usize {
    let orders = compute_bond_orders(input, model);
    let flags = input
        .atoms()
        .filter_map(|(idx, atom)| (model.is_atom_aromatic(idx) || atom.aromatic).then_some(idx))
        .collect();
    let expected = reference_finish(input, &flags, &orders);
    let actual = finish_molecule_with_flags(input, &flags, &orders);
    compare_graph(&expected, &actual);
    input
        .atoms()
        .filter(|(idx, atom)| atom.hydrogen_count != actual.atom(*idx).hydrogen_count)
        .count()
}

#[test]
fn fixed_corpus_normalization_matches_original_storage() {
    let corpus = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"
    ));
    let mut rows = 0;
    let mut patches = 0;
    for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
        let input =
            chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
        for algorithm in [
            AromaticityAlgorithm::Huckel,
            AromaticityAlgorithm::RdkitLike,
        ] {
            patches += compare_normalization(&input, &assign_aromaticity_ex(&input, algorithm));
            if let Ok(kekule) = chematic_core::kekulize(&input) {
                let kekulized = chematic_core::apply_kekule(&input, &kekule);
                patches += compare_normalization(
                    &kekulized,
                    &assign_aromaticity_ex(&kekulized, algorithm),
                );
            }
        }
        rows += 1;
    }
    assert_eq!(rows, 10000);
    assert!(patches > 0, "exercise the indexed hydrogen patch map");
}
