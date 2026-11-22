use chematic_3d::coords::{Coords3D, Point3};
use chematic_3d::minimize::{ForceField, MinimizeConfig, minimize_with_config};
use chematic_core::{Atom, AtomIdx, BondOrder, Element, MoleculeBuilder};

fn config() -> MinimizeConfig {
    MinimizeConfig {
        max_steps: 5,
        step_size: 0.01,
        convergence: 1e-6,
        force_field: ForceField::UFF,
    }
}

#[test]
fn stretched_bonds_contract_for_varied_elements_and_orders() {
    let pairs = [
        ("C", "H", BondOrder::Single),
        ("C", "S", BondOrder::Aromatic),
        ("H", "N", BondOrder::Single),
        ("H", "O", BondOrder::Single),
        ("H", "P", BondOrder::Single),
        ("N", "N", BondOrder::Aromatic),
        ("N", "O", BondOrder::Single),
        ("N", "O", BondOrder::Aromatic),
        ("Cl", "Cl", BondOrder::Double),
        ("Cl", "Cl", BondOrder::Triple),
        ("Cl", "Cl", BondOrder::Aromatic),
        ("C", "C", BondOrder::Triple),
        ("C", "N", BondOrder::Double),
        ("C", "O", BondOrder::Aromatic),
        ("C", "S", BondOrder::Double),
        ("C", "Br", BondOrder::Single),
        ("C", "I", BondOrder::Single),
        ("C", "Si", BondOrder::Single),
        ("H", "H", BondOrder::Single),
        ("H", "S", BondOrder::Single),
        ("N", "N", BondOrder::Triple),
        ("N", "O", BondOrder::Double),
        ("O", "O", BondOrder::Double),
        ("S", "S", BondOrder::Double),
        ("P", "P", BondOrder::Single),
        ("Cl", "Cl", BondOrder::Quadruple),
    ];
    for (a, b, order) in pairs {
        let mut builder = MoleculeBuilder::new();
        let i = builder.add_atom(Atom::new(Element::from_symbol(a).unwrap()));
        let j = builder.add_atom(Atom::new(Element::from_symbol(b).unwrap()));
        builder.add_bond(i, j, order).unwrap();
        let mol = builder.build();
        let mut coords = Coords3D::new_zeroed(2);
        coords.set(AtomIdx(1), Point3::new(3., 0., 0.));
        let result = minimize_with_config(&mol, coords, &config());
        let d = (result.get(i).x - result.get(j).x).abs();
        assert!(d.is_finite() && d < 3., "{a}-{b} {order:?}: {d}");
    }
}

#[test]
fn nonbonded_and_heteroatom_angles_remain_finite() {
    for center in ["O", "N", "S", "P", "Se"] {
        let mut builder = MoleculeBuilder::new();
        let a = builder.add_atom(Atom::new(Element::C));
        let c = builder.add_atom(Atom::new(Element::from_symbol(center).unwrap()));
        let b = builder.add_atom(Atom::new(Element::C));
        builder.add_bond(a, c, BondOrder::Single).unwrap();
        builder.add_bond(c, b, BondOrder::Single).unwrap();
        let mol = builder.build();
        let mut coords = Coords3D::new_zeroed(3);
        coords.set(a, Point3::new(1.5, 0., 0.));
        coords.set(b, Point3::new(0., 1.5, 0.));
        let result = minimize_with_config(&mol, coords, &config());
        for idx in [a, b, c] {
            let p = result.get(idx);
            assert!(
                p.x.is_finite() && p.y.is_finite() && p.z.is_finite(),
                "{center}"
            );
        }
    }
}

#[test]
fn close_unbonded_atoms_repel_for_supported_and_fallback_elements() {
    for symbol in [
        "H", "C", "N", "O", "F", "Si", "P", "S", "Cl", "Br", "I", "Se", "Te", "Ne",
    ] {
        let mut builder = MoleculeBuilder::new();
        let atom = Atom::new(Element::from_symbol(symbol).unwrap());
        let a = builder.add_atom(atom.clone());
        let b = builder.add_atom(atom);
        let mut coords = Coords3D::new_zeroed(2);
        coords.set(b, Point3::new(2., 0., 0.));
        let result = minimize_with_config(&builder.build(), coords, &config());
        let distance = result.get(a).distance(&result.get(b));
        assert!(
            distance.is_finite() && distance > 2.,
            "{symbol}: {distance}"
        );
    }
}
