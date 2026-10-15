//! Compare normalized descriptor graphs with the pre-change copy path.
use super::*;

fn frozen_descriptor_view(mol: &Molecule) -> DescriptorView<'_> {
    // The RDKit-parity view is often an exact copy of `mol`; use `mol` (and
    // its caches) directly then.
    // The view also carries RDKit's sanitization rewrites (perchlorate,
    // organometallic dative bonds), so it is never `mol` when those apply.
    if chematic_perception::rdkit_parity_view_is_identity(mol) {
        return DescriptorView::Same(mol);
    }
    // Memoized on `mol`: TPSA, Crippen LogP/MR, HBA and friends all start from
    // the same perceived copy, and that copy carries its own SSSR cache.
    DescriptorView::Perceived(
        mol.derived(chematic_core::DerivedSlot::DescriptorAromatic, || {
            let mut view = chematic_perception::apply_aromaticity_rdkit_parity_experimental(mol)
                .unwrap_or_else(|_| {
                    let cleaned = chematic_perception::rdkit_sanitize_cleanup(mol);
                    chematic_perception::apply_aromaticity(cleaned.as_ref().unwrap_or(mol))
                });
            // A ring triple bond RDKit flags aromatic (`c1ccccc#1`) is kept
            // as a triple bond in the view; RDKit's descriptors read its
            // aromatic flag.
            let triples: Vec<chematic_core::BondIdx> = view
                .bonds()
                .filter(|(_, b)| {
                    b.order == chematic_core::BondOrder::Triple
                        && view.atom(b.atom1).aromatic
                        && view.atom(b.atom2).aromatic
                })
                .map(|(i, _)| i)
                .collect();
            for b in triples {
                view.set_bond_order(b, chematic_core::BondOrder::Aromatic);
            }
            view
        }),
    )
}

fn compare(mol: &Molecule) {
    let reference_input = mol.clone();
    let expected = frozen_descriptor_view(&reference_input);
    let actual = descriptor_aromaticity(mol);
    assert_eq!(
        expected.atoms().collect::<Vec<_>>(),
        actual.atoms().collect::<Vec<_>>()
    );
    assert_eq!(
        expected.bonds().collect::<Vec<_>>(),
        actual.bonds().collect::<Vec<_>>()
    );
    assert_eq!(expected.rdkit_bond_order(), actual.rdkit_bond_order());
    assert_eq!(expected.stereo_groups(), actual.stereo_groups());
    for (i, _) in expected.atoms() {
        assert_eq!(expected.neighbor_slice(i), actual.neighbor_slice(i));
        assert_eq!(
            expected.stereo_neighbor_order(i),
            actual.stereo_neighbor_order(i)
        );
        assert_eq!(expected.atom_tag(i), actual.atom_tag(i));
        assert_eq!(expected.r_group_label(i), actual.r_group_label(i));
        assert_eq!(
            chematic_core::implicit_hcount(&expected, i),
            chematic_core::implicit_hcount(&actual, i)
        );
    }
    for (i, _) in expected.bonds() {
        assert_eq!(expected.bond_direction(i), actual.bond_direction(i));
        assert_eq!(
            expected.bond_direction_anchor(i),
            actual.bond_direction_anchor(i)
        );
        assert_eq!(
            expected.smiles_ring_closure_metadata(i),
            actual.smiles_ring_closure_metadata(i)
        );
    }
}

#[test]
fn fixed_ten_thousand_inputs_keep_normalized_descriptor_graphs() {
    let corpus = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"
    ));
    let mut rows = 0;
    for line in corpus.lines().filter(|l| !l.trim().is_empty()) {
        let mol = chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
        compare(&mol);
        compare(&mol);
        rows += 1;
    }
    assert_eq!(rows, 10000);
}

#[test]
fn shared_and_normalizing_views_keep_stereo_and_metadata() {
    for smiles in [
        "C1=CC=CC=C1",
        "c1ccccc#1",
        "C1=C=CC=CC1",
        "c1cccc1",
        "n1nnnn1",
        "C/C=C/C1=CC=CC=C1",
        "[C@H](F)(Cl)C1=CC=CC=C1",
        "[*:1]C1=CC=CC=C1",
        "CO(C)[Cu]",
    ] {
        let mut mol = chematic_smiles::parse_template(smiles).unwrap();
        mol.set_tag(chematic_core::AtomIdx(0), Some(23));
        compare(&mol);
        compare(&mol);
    }
}

#[test]
fn externally_seeded_view_keeps_dative_marker_copy() {
    let mol = chematic_smiles::parse_template("C1=CC=CC=C1").unwrap();
    let view = chematic_perception::apply_aromaticity_rdkit_parity_experimental(&mol).unwrap();
    let input = mol.clone();
    input.seed_derived(
        chematic_core::DerivedSlot::RdkitParityAromatic,
        std::sync::Arc::new(Ok::<_, chematic_perception::AromaticityError>(view.clone())),
    );
    let actual = descriptor_aromaticity(&input);
    assert!(
        actual
            .derived_if_computed::<Vec<bool>>(chematic_core::DerivedSlot::RdkitDativeBonds)
            .is_some_and(|flags| !flags.iter().any(|&d| d))
    );
    compare(&input);
}
