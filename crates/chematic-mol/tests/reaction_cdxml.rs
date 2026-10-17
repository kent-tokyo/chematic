use chematic_mol::{parse_cdxml_all, write_reaction_cdxml_checked as write};
use chematic_rxn::{Reaction, parse_reaction};
const ARROW: [[f64; 2]; 2] = [[100., 0.], [200., 0.]];
fn coords(r: &Reaction) -> Vec<Vec<[f64; 2]>> {
    r.reactants
        .iter()
        .chain(&r.agents)
        .chain(&r.products)
        .enumerate()
        .map(|(i, m)| {
            (0..m.atom_count())
                .map(|j| [i as f64 * 100. + j as f64 * 30., (j % 2) as f64 * 15.])
                .collect()
        })
        .collect()
}
fn check_roundtrip(r: &Reaction) -> String {
    let c = coords(r);
    let text = write(r, &c, ARROW).unwrap();
    let parsed = parse_cdxml_all(&text).unwrap();
    assert_eq!(parsed.len(), c.len());
    for (((before, points), (after, read)), i) in r
        .reactants
        .iter()
        .chain(&r.agents)
        .chain(&r.products)
        .zip(&c)
        .zip(parsed)
        .zip(0..)
    {
        assert_eq!(before.atom_count(), after.atom_count(), "component {i}");
        assert_eq!(before.bond_count(), after.bond_count());
        for ((_, a), (_, b)) in before.atoms().zip(after.atoms()) {
            assert_eq!(a, b);
        }
        for ((_, a), (_, b)) in before.bonds().zip(after.bonds()) {
            assert_eq!(a, b);
        }
        assert_eq!(
            points.iter().map(|p| (p[0], p[1])).collect::<Vec<_>>(),
            read
        );
    }
    text
}
#[test]
fn explicit_roles_arrow_and_global_ids() {
    let r = parse_reaction("CCO.O>Cl>CC=O").unwrap();
    let text = check_roundtrip(&r);
    assert!(text.contains("ReactionStepReactants=\"2 3\""));
    assert!(text.contains("ReactionStepObjectsAboveArrow=\"4\""));
    assert!(text.contains("ReactionStepProducts=\"5\""));
    assert_eq!(text.matches("<arrow ").count(), 1);
    let mut ids = std::collections::HashSet::new();
    for line in text.lines().filter(|line| line.contains(" id=\"")) {
        let id = line
            .split("id=\"")
            .nth(1)
            .unwrap()
            .split('"')
            .next()
            .unwrap();
        assert!(ids.insert(id), "duplicate object ID {id}");
    }
    let arrow = text
        .lines()
        .find(|line| line.starts_with("<arrow "))
        .unwrap()
        .split('"')
        .nth(1)
        .unwrap();
    assert!(text.contains(&format!("ReactionStepArrows=\"{arrow}\"")));
}
#[test]
fn duplicate_components_and_empty_agents_are_retained() {
    let r = parse_reaction("O.O>>O.O").unwrap();
    let text = check_roundtrip(&r);
    assert_eq!(text.matches("<fragment ").count(), 4);
    assert!(text.contains("ReactionStepProducts=\"4 5\""));
    assert!(text.contains("ReactionStepObjectsAboveArrow=\"\""));
}
#[test]
fn charge_isotope_explicit_hydrogen_and_triple_bonds() {
    check_roundtrip(&parse_reaction("[13CH3]O.[Na+]>[Cl-]>C#N").unwrap());
}
#[test]
fn missing_sides_and_empty_components_are_refused() {
    for smiles in [">>O", "O>>", ">>"] {
        let r = parse_reaction(smiles).unwrap();
        assert!(write(&r, &coords(&r), ARROW).is_err());
    }
    let mut r = parse_reaction("C>>O").unwrap();
    r.agents.push(chematic_core::MoleculeBuilder::new().build());
    assert!(write(&r, &coords(&r), ARROW).is_err());
}
#[test]
fn exact_coordinate_shapes_and_finite_bounds() {
    let r = parse_reaction("CC>>O").unwrap();
    let c = coords(&r);
    assert!(write(&r, &c[..1], ARROW).is_err());
    let mut bad = c.clone();
    bad[0].pop();
    assert!(write(&r, &bad, ARROW).is_err());
    let mut bad = c.clone();
    bad[1].push([0., 0.]);
    assert!(write(&r, &bad, ARROW).is_err());
    for value in [
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
        1_000_001.,
        -1_000_001.,
    ] {
        let mut bad = c.clone();
        bad[0][0][0] = value;
        assert!(write(&r, &bad, ARROW).is_err());
        assert!(write(&r, &c, [[0., 0.], [value, 1.]]).is_err());
    }
    assert!(write(&r, &c, [[0., 0.]; 2]).is_err());
    assert!(write(&r, &c, [[-1e6, 1e6], [1e6, -1e6]]).is_ok());
}
#[test]
fn unsupported_annotations_and_nonterminal_double_bonds_fail() {
    for smiles in [
        "c1ccccc1>>C",
        "[C@H](F)(Cl)Br>>C",
        "[CH3:1]>>C",
        "*>>C",
        "C/C=C/C>>C",
        "CC=CC>>C",
        "[O]>>O",
        "[N]>>O",
        "[CH3]>>O",
        "[H]>>O",
        "[SiH4]>>O",
    ] {
        let r = parse_reaction(smiles).unwrap();
        assert!(write(&r, &coords(&r), ARROW).is_err(), "{smiles}");
    }
}
#[test]
fn terminal_double_bonds_do_not_invent_geometric_stereo() {
    check_roundtrip(&parse_reaction("CC=O>>C=C").unwrap());
}
#[test]
fn component_cap_exact_and_next_item() {
    let mut r = parse_reaction("C>>O").unwrap();
    r.agents = vec![chematic_smiles::parse("Cl").unwrap(); 30];
    assert!(write(&r, &coords(&r), ARROW).is_ok());
    r.agents.push(chematic_smiles::parse("Cl").unwrap());
    assert!(write(&r, &coords(&r), ARROW).is_err());
}
#[test]
fn atom_caps_exact_and_next_atom() {
    let large = chematic_smiles::parse(&"C".repeat(256)).unwrap();
    let mut r = Reaction {
        reactants: vec![large.clone(); 3],
        agents: vec![],
        products: vec![large],
    };
    assert!(write(&r, &coords(&r), ARROW).is_ok());
    r.agents.push(chematic_smiles::parse("O").unwrap());
    assert!(write(&r, &coords(&r), ARROW).is_err());
    let r = parse_reaction(&format!("{}>>O", "C".repeat(257))).unwrap();
    assert!(write(&r, &coords(&r), ARROW).is_err());
}

#[test]
fn unsupported_bond_types_and_stereo_groups_never_become_plain_bonds() {
    use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder, StereoGroup, StereoGroupKind};
    for order in [
        BondOrder::Quadruple,
        BondOrder::Zero,
        BondOrder::Dative,
        BondOrder::QueryAny,
        BondOrder::QuerySingleOrDouble,
        BondOrder::QuerySingleOrAromatic,
        BondOrder::QueryDoubleOrAromatic,
        BondOrder::Aromatic,
        BondOrder::Up,
        BondOrder::Down,
    ] {
        let mut builder = MoleculeBuilder::new();
        let a = builder.add_atom(Atom::new(Element::C));
        let b = builder.add_atom(Atom::new(Element::N));
        builder.add_bond(a, b, order).unwrap();
        let r = Reaction {
            reactants: vec![builder.build()],
            agents: vec![],
            products: vec![chematic_smiles::parse("O").unwrap()],
        };
        assert!(write(&r, &coords(&r), ARROW).is_err());
    }
    let mut r = parse_reaction("CC>>O").unwrap();
    r.reactants[0].set_stereo_groups(vec![StereoGroup::new(
        StereoGroupKind::Absolute,
        vec![chematic_core::AtomIdx(0)],
    )]);
    assert!(write(&r, &coords(&r), ARROW).is_err());
}
