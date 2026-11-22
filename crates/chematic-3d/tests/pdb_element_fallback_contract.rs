use chematic_3d::pdb::{PdbAtom, pdb_to_molecule};
use chematic_core::{AtomIdx, Element};

fn atom(name: &str, element: &str, x: f64) -> PdbAtom {
    PdbAtom {
        serial: 1,
        name: name.into(),
        res_name: "LIG".into(),
        chain_id: 'A',
        res_seq: 1,
        x,
        y: 0.,
        z: 0.,
        occupancy: 1.,
        temp_factor: 0.,
        element: element.into(),
    }
}

#[test]
fn pdb_conversion_uses_element_then_name_then_carbon_fallback() {
    let input = [
        atom("C1", "Cl", 0.),
        atom(" N1 ", "", 4.),
        atom("X1", "", 8.),
        atom(" ", "", 12.),
    ];
    let (molecule, coords) = pdb_to_molecule(&input);
    assert_eq!(
        molecule.atoms().map(|(_, a)| a.element).collect::<Vec<_>>(),
        [Element::CL, Element::N, Element::C, Element::C]
    );
    assert_eq!(molecule.bond_count(), 0);
    for (i, entry) in input.iter().enumerate() {
        assert_eq!(coords.get(AtomIdx(i as u32)).x, entry.x);
    }
}
