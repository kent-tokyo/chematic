use super::*;
use chematic_smiles::{canonical_smiles, parse};

#[test]
fn core_removal_records_shifted_atoms_and_attachment_labels() {
    let mol = parse("CCOC").unwrap();
    let core = parse("CO").unwrap();
    let result = replace_core_checked(&mol, &core, &[AtomIdx(1), AtomIdx(2)]).unwrap();
    assert_eq!(
        result.original_to_output,
        [Some(AtomIdx(0)), None, None, Some(AtomIdx(1))]
    );
    assert_eq!(result.molecule.atom_count(), 4);
    assert_eq!(
        result
            .attachments
            .iter()
            .map(|a| (a.core_atom.0, a.source_anchor.0, a.anchor.0, a.dummy.0))
            .collect::<Vec<_>>(),
        [(0, 0, 0, 2), (1, 3, 1, 3)]
    );
    assert!(result.molecule.atom(AtomIdx(2)).wildcard);
    assert_eq!(result.molecule.atom(AtomIdx(3)).isotope, Some(1));
    assert_eq!(mol.atom_count(), 4);
}

#[test]
fn core_removal_keeps_whole_unrelated_components_and_remaining_attributes() {
    let mut mol = parse("[13CH3:7]CO.C[NH3+]").unwrap();
    mol.set_tag(AtomIdx(4), Some(91));
    let core = parse("O").unwrap();
    let result = replace_core_checked(&mol, &core, &[AtomIdx(2)]).unwrap();
    for (index, atom) in mol.atoms() {
        if let Some(new) = result.original_to_output[index.0 as usize] {
            assert_eq!(atom, result.molecule.atom(new));
            assert_eq!(mol.atom_tag(index), result.molecule.atom_tag(new));
            assert_eq!(
                mol.implicit_hydrogen_count(index),
                result.molecule.implicit_hydrogen_count(new)
            );
        }
    }
    assert!(canonical_smiles(&result.molecule).contains("."));
}

#[test]
fn core_removal_whole_ring_and_whole_source_are_supported() {
    let mol = parse("Cc1ccccc1").unwrap();
    let core = parse("c1ccccc1").unwrap();
    let matching = (1..7).map(AtomIdx).collect::<Vec<_>>();
    assert_eq!(
        replace_core_checked(&mol, &core, &matching)
            .unwrap()
            .molecule
            .atom_count(),
        2
    );
    let all = (0..mol.atom_count() as u32)
        .map(AtomIdx)
        .collect::<Vec<_>>();
    let result = replace_core_checked(&mol, &mol, &all).unwrap();
    assert_eq!(result.molecule.atom_count(), 0);
    assert!(result.attachments.is_empty());
}

#[test]
fn core_removal_refuses_bad_maps_and_unsupported_chemistry() {
    for (source, core, matching) in [
        ("CC", "C", vec![]),
        ("CC", "C", vec![AtomIdx(2)]),
        ("CC", "CC", vec![AtomIdx(0), AtomIdx(0)]),
        ("CC", "O", vec![AtomIdx(0)]),
        ("C=C", "C", vec![AtomIdx(0)]),
        ("C1CC1", "C", vec![AtomIdx(0)]),
        ("*CC", "C", vec![AtomIdx(1)]),
        ("N[C@H](C)O", "N", vec![AtomIdx(0)]),
        ("CC", "", vec![]),
        ("CC", "C=C", vec![AtomIdx(0), AtomIdx(1)]),
    ] {
        let source = parse(source).unwrap();
        let core = if core.is_empty() {
            chematic_core::MoleculeBuilder::new().build()
        } else {
            parse(core).unwrap()
        };
        assert!(replace_core_checked(&source, &core, &matching).is_err());
    }
}
