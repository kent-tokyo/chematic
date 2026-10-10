use crate::{Layout, Point, RenderOptions, render_eps_opts, render_svg_opts};
use chematic_core::{Atom, AtomIdx, BondIdx, BondOrder, Element, MoleculeBuilder};
fn bond_molecule(order: BondOrder) -> chematic_core::Molecule {
    let mut b = MoleculeBuilder::new();
    let a = b.add_atom(Atom::new(Element::C));
    let c = b.add_atom(Atom::new(Element::C));
    b.add_bond(a, c, order).unwrap();
    b.build()
}
#[test]
fn eps_bond_primitives_have_distinct_paths_and_finite_degenerate_output() {
    for (order, strokes, filled) in [
        (BondOrder::Single, 1, 0),
        (BondOrder::Double, 2, 0),
        (BondOrder::Triple, 3, 0),
        (BondOrder::Up, 0, 1),
        (BondOrder::Down, 7, 0),
        (BondOrder::Dative, 1, 1),
        (BondOrder::Quadruple, 1, 0),
        (BondOrder::Zero, 1, 0),
    ] {
        let mol = bond_molecule(order);
        let layout = Layout {
            coords: vec![Point::new(0., 0.), Point::new(40., 0.)],
        };
        let opts = RenderOptions {
            background: "transparent".into(),
            ..Default::default()
        };
        let output = render_eps_opts(&mol, &layout, &opts);
        assert_eq!(
            output.matches("lineto stroke").count(),
            strokes,
            "{order:?}"
        );
        assert_eq!(
            output.matches("closepath fill").count(),
            filled,
            "{order:?}"
        );
        assert!(!output.contains("rectfill"));
        let collapsed = Layout {
            coords: vec![Point::new(0., 0.); 2],
        };
        let output = render_eps_opts(&mol, &collapsed, &opts);
        assert!(
            !output.contains("NaN") && !output.contains("inf"),
            "{order:?}: {output}"
        );
    }
    for order in [
        BondOrder::QueryAny,
        BondOrder::QuerySingleOrDouble,
        BondOrder::QuerySingleOrAromatic,
        BondOrder::QueryDoubleOrAromatic,
    ] {
        let mol = bond_molecule(order);
        let layout = Layout {
            coords: vec![Point::new(0., 0.), Point::new(40., 0.)],
        };
        let output = render_eps_opts(&mol, &layout, &Default::default());
        assert!(output.contains("setdash"));
        assert!(output.contains("lineto stroke"));
    }
}
#[test]
fn eps_explicit_size_highlights_and_named_colors_survive_export() {
    let mol = bond_molecule(BondOrder::Single);
    let layout = Layout {
        coords: vec![Point::new(-10., 5.), Point::new(30., 5.)],
    };
    for (color, rgb) in [
        ("black", "0.0000 0.0000 0.0000"),
        ("white", "1.0000 1.0000 1.0000"),
        ("red", "1.0000 0.0000 0.0000"),
        ("green", "0.0000 0.5020 0.0000"),
        ("blue", "0.0000 0.0000 1.0000"),
        ("#123456", "0.0706 0.2039 0.3373"),
        ("invalid", "1.0000 1.0000 1.0000"),
    ] {
        let mut opts = RenderOptions {
            width: Some(320),
            height: Some(240),
            background: color.into(),
            dark: true,
            ..Default::default()
        };
        opts.highlight_atoms
            .extend([AtomIdx(0), AtomIdx(1), AtomIdx(99)]);
        opts.highlight_bonds.insert(BondIdx(0));
        opts.atom_color_map.insert(AtomIdx(0), "red".into());
        let output = render_eps_opts(&mol, &layout, &opts);
        assert!(output.contains("%%BoundingBox: 0 0 320 240"));
        assert!(
            output.contains(&format!("{rgb} setrgbcolor")),
            "{color}: {output}"
        );
        assert_eq!(output.matches("arc fill").count(), 2);
        assert!(output.contains("1.0000 0.5490 0.0000 setrgbcolor"));
        assert!(output.contains("4.0000 setlinewidth"));
    }
}
#[test]
fn eps_labels_convert_all_hydrogen_subscripts_to_ascii() {
    for hydrogens in 0..=19 {
        let mut b = MoleculeBuilder::new();
        let mut atom = Atom::new(Element::C);
        atom.hydrogen_count = Some(hydrogens);
        b.add_atom(atom);
        let mol = b.build();
        let layout = Layout {
            coords: vec![Point::new(0., 0.)],
        };
        let output = render_eps_opts(
            &mol,
            &layout,
            &RenderOptions {
                dark: true,
                ..Default::default()
            },
        );
        let expected = match hydrogens {
            0 => "C".into(),
            1 => "CH".into(),
            n => format!("CH{n}"),
        };
        assert!(output.contains(&format!("({expected}) show")), "{output}");
        assert!(output.is_ascii());
    }
}
#[test]
fn svg_query_and_degenerate_bonds_keep_valid_numeric_output() {
    for order in [
        BondOrder::Triple,
        BondOrder::Up,
        BondOrder::Down,
        BondOrder::Dative,
        BondOrder::Zero,
        BondOrder::Quadruple,
        BondOrder::QueryAny,
        BondOrder::QuerySingleOrDouble,
        BondOrder::QuerySingleOrAromatic,
        BondOrder::QueryDoubleOrAromatic,
    ] {
        let mol = bond_molecule(order);
        for separation in [0., 40.] {
            let layout = Layout {
                coords: vec![Point::new(0., 0.), Point::new(separation, 0.)],
            };
            let opts = RenderOptions {
                background: "transparent".into(),
                dark: true,
                atom_ids: true,
                show_atom_indices: true,
                ..Default::default()
            };
            let svg = render_svg_opts(&mol, &layout, &opts);
            assert!(svg.starts_with("<svg") && svg.ends_with("</svg>"));
            assert!(
                !svg.contains("NaN") && !svg.contains("inf"),
                "{order:?}: {svg}"
            );
        }
    }
}
