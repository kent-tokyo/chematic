use chematic_core::{AtomIdx, BondOrder};
use chematic_mol::{parse_ket_3d, write_ket_3d};

#[test]
fn ket_roundtrip_preserves_coordinates_charge_isotope_and_stereo_bonds() {
    for (bond, expected) in [
        (r#"{"type":1,"stereo":1,"atoms":[0,1]}"#, BondOrder::Up),
        (r#"{"type":1,"stereo":6,"atoms":[0,1]}"#, BondOrder::Down),
        (r#"{"type":2,"atoms":[0,1]}"#, BondOrder::Double),
        (r#"{"type":3,"atoms":[0,1]}"#, BondOrder::Triple),
        (r#"{"type":4,"atoms":[0,1]}"#, BondOrder::Aromatic),
    ] {
        let input = format!(
            r#"{{"atoms":[{{"label":"C","isotope":13,"location":[1,2,3]}},{{"label":"N","charge":1,"location":[4,5,6]}}],"bonds":[{bond}]}}"#
        );
        let (molecule, coords) = parse_ket_3d(&input).unwrap();
        let (reread, roundtrip_coords) = parse_ket_3d(&write_ket_3d(&molecule, &coords)).unwrap();
        assert_eq!(roundtrip_coords, [(1.0, 2.0, 3.0), (4.0, 5.0, 6.0)]);
        assert_eq!(reread.atom(AtomIdx(0)).isotope, Some(13));
        assert_eq!(reread.atom(AtomIdx(1)).charge, 1);
        assert_eq!(reread.bonds().next().unwrap().1.order, expected);
    }
}
