use chematic_core::{AtomIdx, STEREO_H_SENTINEL};
use chematic_smiles::{canonical_smiles, parse};

#[test]
fn ring_digits_take_their_written_place_at_an_allene_terminal() {
    let opening_first = parse("CC=[C@]=C1CC[C@H](CC1)CC").unwrap();
    let opening_last = parse("CC=[C@@]=C(CC[C@H](CC1)CC)1").unwrap();
    assert_eq!(
        opening_first.stereo_neighbor_order(AtomIdx(2)).unwrap(),
        &[1, 0, STEREO_H_SENTINEL, 3, 8, 4]
    );
    assert_eq!(
        opening_last.stereo_neighbor_order(AtomIdx(2)).unwrap(),
        &[1, 0, STEREO_H_SENTINEL, 3, 4, 8]
    );
    let expected = canonical_smiles(&opening_first);
    assert_eq!(canonical_smiles(&opening_last), expected);
    assert_eq!(canonical_smiles(&parse(&expected).unwrap()), expected);
}

#[test]
fn disconnected_syntax_preserves_the_initial_hydrogen_position() {
    let initial = parse("FC1.C1=[C@]=C(C)Cl").unwrap();
    let connected = parse("FCC=[C@@]=C(C)Cl").unwrap();
    assert_eq!(canonical_smiles(&initial), canonical_smiles(&connected));
}

#[test]
fn a_directional_single_bond_is_a_valid_allene_ligand() {
    let molecule = parse("C/C=C\\C=C=[C@@]=C=CO").unwrap();
    assert!(molecule.is_allene_stereo(AtomIdx(5)));
    let canonical = canonical_smiles(&molecule);
    let reparsed = parse(&canonical).unwrap();
    assert!(
        reparsed
            .atoms()
            .any(|(atom, _)| reparsed.is_allene_stereo(atom))
    );
    assert_eq!(canonical_smiles(&reparsed), canonical);
}
