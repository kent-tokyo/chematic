//! Exact comparison against the pre-change aromatic view path.
use super::*;

fn frozen_unseeded_reference(mol: &Molecule) -> Result<Molecule, AromaticityError> {
    match parity_shortcut(mol) {
        ParityShortcut::Identity => return Ok(mol.clone()),
        ParityShortcut::ClearFlags => return Ok(clear_aromatic_flags(mol)),
        ParityShortcut::Full => {}
    }
    let (kekulized, verdict) = kekulize_for_rdkit_parity_with_verdict(mol)?;
    if kekulized
        .bonds()
        .any(|(_, bond)| bond.order == BondOrder::Aromatic)
    {
        return Ok(kekulized);
    }
    let model = match verdict {
        Some(v) => assign_from_verdict(&kekulized, v)?,
        None => assign_from_kekulized(&kekulized)?,
    };
    Ok(crate::aromaticity::build_molecule_from_model(
        &kekulized, &model,
    ))
}


fn compare(input: &Molecule) {
    let expected = frozen_unseeded_reference(input);
    let actual = apply_aromaticity_rdkit_parity_unseeded(input);
    match (expected, actual) {
        (Err(a), Err(b)) => assert_eq!(a, b),
        (Ok(a), Ok(b)) => {
            assert_eq!(a.atoms().collect::<Vec<_>>(), b.atoms().collect::<Vec<_>>());
            assert_eq!(a.bonds().collect::<Vec<_>>(), b.bonds().collect::<Vec<_>>());
            assert_eq!(a.rdkit_bond_order(), b.rdkit_bond_order());
            assert_eq!(a.stereo_groups(), b.stereo_groups());
            for (i, _) in a.atoms() {
                assert_eq!(a.neighbor_slice(i), b.neighbor_slice(i));
                assert_eq!(a.stereo_neighbor_order(i), b.stereo_neighbor_order(i));
                assert_eq!(a.r_group_label(i), b.r_group_label(i));
                assert_eq!(a.atom_tag(i), b.atom_tag(i));
                assert_eq!(a.smiles_ring_closure_count(i), b.smiles_ring_closure_count(i));
            }
            for (i, _) in a.bonds() {
                assert_eq!(a.bond_direction(i), b.bond_direction(i));
                assert_eq!(a.bond_direction_anchor(i), b.bond_direction_anchor(i));
                assert_eq!(a.smiles_ring_closure_metadata(i), b.smiles_ring_closure_metadata(i));
            }
            for i in 0..a.stereo_groups().len() {
                assert_eq!(a.stereo_group_bond_indices(i), b.stereo_group_bond_indices(i));
            }
        }
        _ => panic!("aromatic view success/error boundary changed"),
    }
}

#[test]
fn fixed_ten_thousand_inputs_match_frozen_view() {
    let corpus = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"));
    let mut rows = 0;
    let mut borrowed = 0;
    for line in corpus.lines().filter(|l| !l.trim().is_empty()) {
        let mol = chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
        borrowed += usize::from(parity_shortcut(&mol) == ParityShortcut::Full && already_kekulized_in_bond_order(&mol));
        compare(&mol);
        rows += 1;
    }
    assert_eq!(rows, 10000);
    assert!(borrowed > 1000, "exercise the new full-perception path");
}

#[test]
fn flags_stereo_tags_and_r_groups_match_reference() {
    for smiles in ["C1=CC=CC=C1", "c1ccccc1", "C1=Cc2ccccc2C1", "O=C1NNC(=O)N1", "[C@H](F)(Cl)C1=CC=CC=C1", "C/C=C/C1=CC=CC=C1", "[*:1]C1=CC=CC=C1", "C1=C=CC=CC1", "[O-][Cl](=O)(=O)=O", "CO(C)[Cu]"] {
        let mut mol = chematic_smiles::parse_template(smiles).unwrap();
        mol.set_tag(AtomIdx(0), Some(23));
        mol.add_stereo_group_with_bonds(
            chematic_core::StereoGroup::new(chematic_core::StereoGroupKind::Or(9), vec![AtomIdx(0)]),
            vec![BondIdx(0)],
        );
        if mol.atom(AtomIdx(0)).wildcard {mol.set_r_group(AtomIdx(0), chematic_core::RGroupLabel::unnumbered());}
        compare(&mol);
        mol.set_atom_aromatic(AtomIdx(0), true);
        assert!(!already_kekulized_in_bond_order(&mol));
        compare(&mol);
    }
}

#[test]
fn removed_and_readded_bond_matches_reference() {
    let mut mol = chematic_smiles::parse_template("C1=CC=CC=C1").unwrap();
    let first = mol.bond(BondIdx(0)).clone();
    mol.remove_bond(BondIdx(0));
    mol.add_bond(first.atom1, first.atom2, first.order).unwrap();
    compare(&mol);
}
