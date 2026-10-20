use super::*;
use chematic_smiles::parse;

fn mol(s: &str) -> chematic_core::Molecule {
    parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"))
}

#[test]
fn test_perp_unit_horizontal() {
    let (px, py) = perp_unit(Point::new(0.0, 0.0), Point::new(1.0, 0.0));
    assert!((px.abs() - 0.0).abs() < 1e-9);
    assert!((py.abs() - 1.0).abs() < 1e-9);
}

#[test]
fn isolated_water_label_h2o() {
    let m = mol("O");
    let layout = crate::layout::compute_layout(&m);
    let svg = render_svg(&m, &layout);
    assert!(
        svg.contains("H2O"),
        "water 'O' should show H2O label, got: {}",
        svg
    );
}

#[test]
fn isolated_methane_label_ch4() {
    let m = mol("C");
    let layout = crate::layout::compute_layout(&m);
    let svg = render_svg(&m, &layout);
    assert!(
        svg.contains("CH4"),
        "methane 'C' should show CH4 label, got: {}",
        svg
    );
    assert!(svg.contains("<text"), "single C must have a text label now");
}

#[test]
fn isolated_ammonia_label_h3n() {
    let m = mol("N");
    let layout = crate::layout::compute_layout(&m);
    let svg = render_svg(&m, &layout);
    assert!(svg.contains("H3N"), "ammonia 'N' should show H3N label");
}

#[test]
fn multi_atom_carbon_no_label() {
    let m = mol("CC");
    let layout = crate::layout::compute_layout(&m);
    let svg = render_svg(&m, &layout);
    assert!(!svg.contains("<text"), "ethane should have no atom labels");
}

#[test]
fn render_opts_transparent_no_bg_rect() {
    let m = mol("c1ccccc1");
    let layout = crate::layout::compute_layout(&m);
    let opts = RenderOptions {
        background: "transparent".into(),
        ..Default::default()
    };
    let svg = render_svg_opts(&m, &layout, &opts);
    assert!(svg.contains("<svg"), "must be valid SVG");
    // No background rect when transparent.
    assert!(
        !svg.contains("fill=\"transparent\""),
        "no bg rect fill for transparent"
    );
}

#[test]
fn render_opts_custom_size() {
    let m = mol("CCO");
    let layout = crate::layout::compute_layout(&m);
    let opts = RenderOptions {
        width: Some(300),
        height: Some(200),
        ..Default::default()
    };
    let svg = render_svg_opts(&m, &layout, &opts);
    assert!(svg.contains("width=\"300\""), "SVG width should be 300");
    assert!(svg.contains("height=\"200\""), "SVG height should be 200");
}

#[test]
fn render_opts_dark_theme_white_bonds() {
    let m = mol("CC");
    let layout = crate::layout::compute_layout(&m);
    let opts = RenderOptions {
        dark: true,
        background: "#0f172a".into(),
        ..Default::default()
    };
    let svg = render_svg_opts(&m, &layout, &opts);
    assert!(
        svg.contains("stroke=\"white\""),
        "dark theme bonds should be white"
    );
}

#[test]
fn render_opts_highlight_atoms() {
    let m = mol("c1ccncc1");
    let layout = crate::layout::compute_layout(&m);
    let n_idx = m
        .atoms()
        .find(|(_, a)| a.element.atomic_number() == 7)
        .map(|(idx, _)| idx)
        .expect("pyridine has N");
    let mut hl = std::collections::HashSet::new();
    hl.insert(n_idx);
    let opts = RenderOptions {
        highlight_atoms: hl,
        ..Default::default()
    };
    let svg = render_svg_opts(&m, &layout, &opts);
    assert!(svg.contains("<circle"), "highlight must produce a circle");
}

#[test]
fn test_atom_color_rgb_matches_hex() {
    // Verify atom_color_rgb values match atom_color hex strings.
    for atomic_number in [7u8, 8, 9, 15, 16, 17, 35, 53, 6] {
        let hex = atom_color(atomic_number);
        let rgb = atom_color_rgb(atomic_number);
        let r = u8::from_str_radix(&hex[1..3], 16).unwrap();
        let g = u8::from_str_radix(&hex[3..5], 16).unwrap();
        let b = u8::from_str_radix(&hex[5..7], 16).unwrap();
        assert_eq!(
            rgb,
            [r, g, b],
            "mismatch for atomic_number {atomic_number}: hex={hex}"
        );
    }
}

#[test]
fn test_atom_color_rgb_nitrogen_blue() {
    assert_eq!(atom_color_rgb(7), [0x30, 0x50, 0xF8]);
}

#[test]
fn test_atom_display_label_terminal_carbon() {
    let m = mol("CC");
    // Terminal carbon has 3 implicit H
    let label = atom_display_label(&m, chematic_core::AtomIdx(0));
    assert_eq!(label, "CH₃", "terminal carbon should show CH₃");
}

#[test]
fn test_atom_display_label_interior_carbon() {
    let m = mol("CCC");
    // Interior carbons have 2 implicit H
    let label = atom_display_label(&m, chematic_core::AtomIdx(1));
    assert_eq!(label, "", "interior carbon should have empty label");
}

#[test]
fn test_atom_display_label_isolated_oxygen() {
    let m = mol("O");
    let label = atom_display_label(&m, chematic_core::AtomIdx(0));
    assert_eq!(label, "H₂O", "isolated oxygen should show H₂O");
}

#[test]
fn test_atom_label_with_h_nitrogen() {
    let m = mol("CCN");
    let label = atom_label_with_h(&m, chematic_core::AtomIdx(2));
    assert_eq!(label.symbol, "N");
    assert_eq!(label.h_count, 2);
    assert_eq!(label.h_position, HPosition::Right);
}

#[test]
fn test_atom_color_rgb_unknown_black() {
    assert_eq!(atom_color_rgb(0), [0, 0, 0]);
    assert_eq!(atom_color_rgb(118), [0, 0, 0]);
}

// C5 Tests: Dative and query bond rendering

#[test]
fn test_render_dative_bond() {
    // Test that dative bonds render with arrow notation
    let p1 = Point::new(0.0, 0.0);
    let p2 = Point::new(10.0, 0.0);
    let svg = render_dative_bond(p1, p2, "black");

    // Should contain line SVG element
    assert!(
        svg.contains("<line"),
        "dative bond should include line element"
    );
    // Should contain polygon for arrowhead
    assert!(
        svg.contains("<polygon"),
        "dative bond should include arrowhead polygon"
    );
}

#[test]
fn test_render_query_any_bond() {
    // Test that query_any bonds render with dotted line and asterisk
    let p1 = Point::new(0.0, 0.0);
    let p2 = Point::new(10.0, 0.0);
    let svg = render_query_any_bond(p1, p2, "black");

    // Should contain dashed line
    assert!(
        svg.contains("stroke-dasharray"),
        "query_any should have dashed line"
    );
    // Should contain asterisk marker
    assert!(
        svg.contains("*"),
        "query_any should include asterisk marker"
    );
    assert!(
        svg.contains("<text"),
        "query_any should include text for asterisk"
    );
}

#[test]
fn test_render_query_single_or_double() {
    // Test query single-or-double bond with short dashes
    let p1 = Point::new(0.0, 0.0);
    let p2 = Point::new(10.0, 0.0);
    let svg = render_query_dashed_bond(p1, p2, "black", "2,2");

    assert!(
        svg.contains("stroke-dasharray=\"2,2\""),
        "should have correct dash pattern"
    );
}

#[test]
fn test_render_query_single_or_aromatic() {
    // Test query single-or-aromatic bond with medium dashes
    let p1 = Point::new(0.0, 0.0);
    let p2 = Point::new(10.0, 0.0);
    let svg = render_query_dashed_bond(p1, p2, "black", "4,2");

    assert!(
        svg.contains("stroke-dasharray=\"4,2\""),
        "should have correct dash pattern"
    );
}

#[test]
fn test_render_query_double_or_aromatic() {
    // Test query double-or-aromatic bond with long dashes
    let p1 = Point::new(0.0, 0.0);
    let p2 = Point::new(10.0, 0.0);
    let svg = render_query_dashed_bond(p1, p2, "black", "3,3");

    assert!(
        svg.contains("stroke-dasharray=\"3,3\""),
        "should have correct dash pattern"
    );
}

#[test]
fn test_bond_rendering_all_types() {
    // Comprehensive test ensuring all bond types render without errors
    let positions = vec![
        (BondOrder::Single, "single"),
        (BondOrder::Double, "double"),
        (BondOrder::Triple, "triple"),
        (BondOrder::Aromatic, "aromatic"),
        (BondOrder::Up, "wedge up"),
        (BondOrder::Down, "down"),
        (BondOrder::Dative, "dative"),
        (BondOrder::QueryAny, "query any"),
        (BondOrder::QuerySingleOrDouble, "query single-or-double"),
        (BondOrder::QuerySingleOrAromatic, "query single-or-aromatic"),
        (BondOrder::QueryDoubleOrAromatic, "query double-or-aromatic"),
    ];

    let p1 = Point::new(0.0, 0.0);
    let p2 = Point::new(10.0, 0.0);

    for (bond_order, name) in positions {
        let svg = render_bond_c(bond_order, p1, p2, "black");
        assert!(
            !svg.is_empty(),
            "{} bond should produce non-empty SVG",
            name
        );
        assert!(
            svg.contains("<"),
            "{} bond should contain SVG elements",
            name
        );
    }
}

#[test]
fn test_dative_bond_color() {
    // Test that dative bond respects color parameter
    let p1 = Point::new(0.0, 0.0);
    let p2 = Point::new(10.0, 0.0);
    let svg = render_dative_bond(p1, p2, "#FF0000");

    assert!(
        svg.contains("#FF0000"),
        "dative bond should use specified color"
    );
}

#[test]
fn test_query_any_bond_midpoint_marker() {
    // Test that asterisk appears approximately at bond midpoint
    let p1 = Point::new(0.0, 0.0);
    let p2 = Point::new(20.0, 0.0);
    let svg = render_query_any_bond(p1, p2, "black");

    // Should have text element with coordinates near midpoint (10.0, 0.0)
    assert!(
        svg.contains("x=\"10.00\""),
        "asterisk should be near midpoint"
    );
}

#[test]
fn test_bond_rendering_coverage() {
    // Ensure all BondOrder variants are handled
    let orders = vec![
        BondOrder::Zero,
        BondOrder::Single,
        BondOrder::Double,
        BondOrder::Triple,
        BondOrder::Aromatic,
        BondOrder::Up,
        BondOrder::Down,
        BondOrder::Quadruple,
        BondOrder::Dative,
        BondOrder::QueryAny,
        BondOrder::QuerySingleOrDouble,
        BondOrder::QuerySingleOrAromatic,
        BondOrder::QueryDoubleOrAromatic,
    ];

    let p1 = Point::new(0.0, 0.0);
    let p2 = Point::new(10.0, 0.0);

    for order in orders {
        // Each should render without panicking
        let _ = render_bond_c(order, p1, p2, "black");
    }
}
