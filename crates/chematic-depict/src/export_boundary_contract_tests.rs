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

#[test]
fn metadata_export_keeps_graph_styles_and_escapes_authored_text() {
    let mol = chematic_smiles::parse("c1ccncc1.[NH4+]").unwrap();
    let layout = crate::compute_layout(&mol);
    for kekulize in [false, true] {
        let mut opts = RenderOptions {
            kekulize,
            atom_ids: true,
            show_atom_indices: true,
            background: "transparent".into(),
            ..Default::default()
        };
        opts.highlight_atoms.extend([AtomIdx(0), AtomIdx(999)]);
        opts.highlight_bonds.insert(BondIdx(0));
        opts.atom_color_map.insert(AtomIdx(0), "red".into());
        opts.atom_color_map.insert(AtomIdx(1), "blue".into());
        opts.atom_color_map.insert(AtomIdx(999), "invalid".into());
        for text in ["", "C<&\"' >"] {
            let svg = crate::render_svg_with_metadata(&mol, &layout, &opts, text);
            if text.is_empty() {
                let metadata = svg
                    .split("<smiles>")
                    .nth(1)
                    .unwrap()
                    .split("</smiles>")
                    .next()
                    .unwrap();
                let recovered = chematic_smiles::parse(metadata).unwrap();
                assert_eq!(
                    chematic_smiles::rdkit_canonical_smiles(&recovered).unwrap(),
                    chematic_smiles::rdkit_canonical_smiles(&mol).unwrap()
                );
            } else {
                assert!(svg.contains("<smiles>C&lt;&amp;&quot;' &gt;</smiles>"));
            }
            assert_eq!(svg.matches("<circle ").count(), 2);
            assert!(svg.contains("#FF8C00"));
            assert!(svg.contains("data-atom-idx"));
            assert!(!svg.contains("NaN"));
        }
    }
}

#[test]
fn styled_grid_defaults_match_plain_grid_and_keep_per_cell_styles() {
    let a = chematic_smiles::parse("CO").unwrap();
    let b = chematic_smiles::parse("N#N").unwrap();
    let layouts = vec![crate::compute_layout(&a), crate::compute_layout(&b)];
    let default = RenderOptions::default();
    assert_eq!(
        crate::depict_svg_grid_with_opts_and_layouts(
            &[(&a, None), (&b, Some(&default))],
            &layouts,
            1
        ),
        crate::depict_svg_grid_with_layouts(&[&a, &b], &layouts, 1)
    );
    let dark = RenderOptions {
        dark: true,
        ..Default::default()
    };
    let grid = crate::depict_svg_grid_with_opts_and_layouts(
        &[(&a, Some(&dark)), (&b, None)],
        &layouts,
        99,
    );
    assert_eq!(grid.matches("<g id=\"mol-").count(), 2);
    assert!(grid.contains("mol-0") && grid.contains("mol-1"));
    assert!(!grid.contains("NaN"));
    for (mols, cols) in [(vec![], 2), (vec![(&a, Some(&dark))], 0)] {
        assert!(
            crate::depict_svg_grid_with_opts_and_layouts(&mols, &layouts, cols)
                .contains("width=\"0\"")
        );
    }
}

#[test]
fn preflight_reports_graph_limits_overlapping_labels_and_crossings_at_their_paths() {
    use crate::{Layout, Point, RenderOptions, preflight::*};
    let mol = chematic_smiles::parse("NN.NN").unwrap();
    let layout = Layout {
        coords: vec![
            Point::new(0., 0.),
            Point::new(20., 20.),
            Point::new(0., 20.),
            Point::new(20., 0.),
        ],
    };
    let options = RenderOptions {
        padding: 0.,
        ..Default::default()
    };
    let report = preflight_svg(
        &mol,
        &layout,
        &options,
        &PreflightLimits {
            max_atoms: 3,
            max_bonds: 1,
            max_labels: 3,
        },
    );
    assert!(!report.ready);
    for path in ["molecule.atoms", "molecule.bonds", "molecule.labels"] {
        assert!(
            report
                .diagnostics
                .iter()
                .any(|d| d.code == PreflightCode::ResourceLimit && d.path == path),
            "{path}: {:?}",
            report.diagnostics
        );
    }
    assert!(
        report
            .diagnostics
            .iter()
            .any(|d| d.code == PreflightCode::CrossingBond && d.path == "bonds[0]/bonds[1]")
    );
    assert!(
        report
            .diagnostics
            .iter()
            .any(|d| d.code == PreflightCode::ClippedAtom)
    );
    assert!(
        report
            .diagnostics
            .iter()
            .any(|d| d.code == PreflightCode::ClippedLabel)
    );
    let crowded = Layout {
        coords: vec![
            Point::new(0., 0.),
            Point::new(1., 0.),
            Point::new(2., 0.),
            Point::new(3., 0.),
        ],
    };
    let report = preflight_svg(
        &mol,
        &crowded,
        &RenderOptions::default(),
        &PreflightLimits::default(),
    );
    assert!(
        report
            .diagnostics
            .iter()
            .any(|d| d.code == PreflightCode::LabelOverlap)
    );
    assert!(
        report
            .diagnostics
            .iter()
            .any(|d| d.code == PreflightCode::AtomOverlap)
    );
    let short = Layout {
        coords: vec![Point::new(0., 0.)],
    };
    let report = preflight_svg(
        &mol,
        &short,
        &RenderOptions::default(),
        &PreflightLimits::default(),
    );
    assert!(
        report
            .diagnostics
            .iter()
            .any(|d| d.code == PreflightCode::ResourceLimit && d.path == "layout.coords")
    );
}
