//! PostScript EPS renderer for molecular 2D layouts.
//!
//! Generates a self-contained EPS string directly from a `Layout` and `Molecule`
//! without any external dependencies.  Coordinate convention: SVG uses Y-down;
//! EPS uses Y-up.  Transformation: `eps_x = x − view_x`, `eps_y = view_h − (y − view_y)`.

use chematic_core::{BondOrder, Molecule};

use crate::layout::{BOND_LEN, Layout, Point};
use crate::svg::{RenderOptions, atom_color_rgb, atom_display_label};

// ---------------------------------------------------------------------------
// Public entry points
// ---------------------------------------------------------------------------

/// Compute 2D layout and render as a self-contained EPS string.
pub fn render_eps(mol: &Molecule, layout: &Layout) -> String {
    render_eps_opts(mol, layout, &RenderOptions::default())
}

struct EpsViewport {
    view_x: f64,
    view_y: f64,
    view_h: f64,
    width: f64,
    height: f64,
    bond_color: (f64, f64, f64),
    background: Option<(f64, f64, f64)>,
}

impl EpsViewport {
    fn new(layout: &Layout, opts: &RenderOptions) -> Self {
        let (min_x, min_y, max_x, max_y) = layout.bounding_box();
        let raw_w = (max_x - min_x).max(BOND_LEN);
        let raw_h = (max_y - min_y).max(BOND_LEN);
        let view_w = raw_w + 2.0 * opts.padding;
        let view_h = raw_h + 2.0 * opts.padding;
        Self {
            view_x: min_x - opts.padding,
            view_y: min_y - opts.padding,
            view_h,
            width: opts.width.unwrap_or(view_w.round() as u32) as f64,
            height: opts.height.unwrap_or(view_h.round() as u32) as f64,
            bond_color: if opts.dark {
                (1.0, 1.0, 1.0)
            } else {
                (0.0, 0.0, 0.0)
            },
            background: (opts.background != "transparent")
                .then(|| parse_color_or(&opts.background, (1.0, 1.0, 1.0))),
        }
    }

    fn point(&self, point: Point) -> Point {
        Point::new(point.x - self.view_x, self.view_h - (point.y - self.view_y))
    }
}

fn write_eps_header(out: &mut String, viewport: &EpsViewport) {
    out.push_str("%!PS-Adobe-3.0 EPSF-3.0\n");
    out.push_str(&format!(
        "%%BoundingBox: 0 0 {} {}\n",
        viewport.width.ceil() as u32,
        viewport.height.ceil() as u32
    ));
    out.push_str("%%EndComments\n");
    out.push_str("/Helvetica findfont 12 scalefont setfont\n");
    out.push_str("1 setlinecap\n1 setlinejoin\n");
    if let Some((r, g, b)) = viewport.background {
        out.push_str(&format!("{r:.4} {g:.4} {b:.4} setrgbcolor\n"));
        out.push_str(&format!(
            "0 0 {:.2} {:.2} rectfill\n",
            viewport.width, viewport.height
        ));
    }
}

fn write_eps_highlights(
    out: &mut String,
    mol: &Molecule,
    layout: &Layout,
    opts: &RenderOptions,
    viewport: &EpsViewport,
) {
    let default_color = parse_color_or(&opts.highlight_color, (1.0, 1.0, 0.0));
    for atom in &opts.highlight_atoms {
        if atom.0 as usize >= mol.atom_count() {
            continue;
        }
        let point = viewport.point(layout.get(*atom));
        let (r, g, b) = opts
            .atom_color_map
            .get(atom)
            .map(|color| parse_color_or(color, default_color))
            .unwrap_or(default_color);
        out.push_str(&format!("{r:.4} {g:.4} {b:.4} setrgbcolor\n"));
        out.push_str(&format!(
            "{:.2} {:.2} 16 0 360 arc fill\n",
            point.x, point.y
        ));
    }
}

fn write_eps_bonds(
    out: &mut String,
    mol: &Molecule,
    layout: &Layout,
    opts: &RenderOptions,
    viewport: &EpsViewport,
) {
    set_color(out, viewport.bond_color);
    for (bond_index, bond) in mol.bonds() {
        let (start, end) =
            crate::svg::linear_carbon::bond_endpoints(mol, layout, bond.atom1, bond.atom2);
        let color = if opts.highlight_bonds.contains(&bond_index) {
            (1.0_f64, 0.549, 0.0)
        } else {
            viewport.bond_color
        };
        render_bond_eps(
            out,
            bond.order,
            viewport.point(start),
            viewport.point(end),
            color,
            if opts.highlight_bonds.contains(&bond_index) {
                4.0
            } else {
                1.5
            },
        );
    }
}

fn write_eps_labels(
    out: &mut String,
    mol: &Molecule,
    layout: &Layout,
    opts: &RenderOptions,
    viewport: &EpsViewport,
) {
    for (atom_index, atom) in mol.atoms() {
        let raw_label = atom_display_label(mol, atom_index);
        if raw_label.is_empty() {
            continue;
        }
        let label = unicode_subscripts_to_ascii(&raw_label);
        let point = viewport.point(layout.get(atom_index));
        if let Some((r, g, b)) = viewport.background {
            out.push_str(&format!("{r:.4} {g:.4} {b:.4} setrgbcolor\n"));
            out.push_str(&format!(
                "{:.2} {:.2} {:.2} {:.2} rectfill\n",
                point.x - 8.0,
                point.y - 7.0,
                16.0,
                14.0,
            ));
        }
        let text_color = if opts.dark && atom.element.atomic_number() == 6 {
            (1.0_f64, 1.0, 1.0)
        } else {
            let [r, g, b] = atom_color_rgb(atom.element.atomic_number());
            (r as f64 / 255.0, g as f64 / 255.0, b as f64 / 255.0)
        };
        set_color(out, text_color);
        let approximate_width = label.chars().count() as f64 * 3.5;
        out.push_str(&format!(
            "{:.2} {:.2} moveto ({}) show\n",
            point.x - approximate_width,
            point.y - 4.0,
            escape_ps(&label),
        ));
    }
}

/// Compute 2D layout and render as EPS with full style control.
pub fn render_eps_opts(mol: &Molecule, layout: &Layout, opts: &RenderOptions) -> String {
    let viewport = EpsViewport::new(layout, opts);
    let mut out = String::with_capacity(8192);
    write_eps_header(&mut out, &viewport);
    write_eps_highlights(&mut out, mol, layout, opts, &viewport);
    write_eps_bonds(&mut out, mol, layout, opts, &viewport);
    write_eps_labels(&mut out, mol, layout, opts, &viewport);
    out.push_str("showpage\n%%EOF\n");
    out
}

// ---------------------------------------------------------------------------
// Bond rendering helpers
// ---------------------------------------------------------------------------

fn render_bond_eps(
    out: &mut String,
    order: BondOrder,
    p1: Point,
    p2: Point,
    color: (f64, f64, f64),
    base_w: f64,
) {
    match order {
        BondOrder::Single => eps_line(out, p1, p2, base_w, color),
        BondOrder::Up => eps_wedge(out, p1, p2, color),
        BondOrder::Down => eps_dash_bond(out, p1, p2, color),
        BondOrder::Double => eps_double_bond(out, p1, p2, color),
        BondOrder::Triple => eps_triple_bond(out, p1, p2, color),
        BondOrder::Aromatic => eps_aromatic_bond(out, p1, p2, color),
        BondOrder::Quadruple => eps_line(out, p1, p2, 3.0, color),
        BondOrder::Zero => eps_line(out, p1, p2, 1.0, color),
        BondOrder::Dative => eps_dative_bond(out, p1, p2, color),
        // Query bonds: dashed lines
        BondOrder::QueryAny
        | BondOrder::QuerySingleOrDouble
        | BondOrder::QuerySingleOrAromatic
        | BondOrder::QueryDoubleOrAromatic => eps_dashed_line(out, p1, p2, base_w, color),
    }
}

fn perp_unit(p1: Point, p2: Point) -> (f64, f64) {
    let dx = p2.x - p1.x;
    let dy = p2.y - p1.y;
    let len = (dx * dx + dy * dy).sqrt();
    if len < 1e-10 {
        return (0.0, 1.0);
    }
    (-dy / len, dx / len)
}

fn eps_line(out: &mut String, p1: Point, p2: Point, w: f64, color: (f64, f64, f64)) {
    set_color(out, color);
    out.push_str(&format!(
        "{:.4} setlinewidth\n{:.2} {:.2} moveto {:.2} {:.2} lineto stroke\n",
        w, p1.x, p1.y, p2.x, p2.y
    ));
}

fn eps_dashed_line(out: &mut String, p1: Point, p2: Point, w: f64, color: (f64, f64, f64)) {
    set_color(out, color);
    out.push_str(&format!("{:.4} setlinewidth\n", w));
    out.push_str("[4 3] 0 setdash\n");
    out.push_str(&format!(
        "{:.2} {:.2} moveto {:.2} {:.2} lineto stroke\n",
        p1.x, p1.y, p2.x, p2.y
    ));
    out.push_str("[] 0 setdash\n");
}

fn eps_double_bond(out: &mut String, p1: Point, p2: Point, color: (f64, f64, f64)) {
    let offset = 2.0;
    let (px, py) = perp_unit(p1, p2);
    for sign in [-1.0_f64, 1.0] {
        eps_line(
            out,
            Point::new(p1.x + px * offset * sign, p1.y + py * offset * sign),
            Point::new(p2.x + px * offset * sign, p2.y + py * offset * sign),
            1.5,
            color,
        );
    }
}

fn eps_triple_bond(out: &mut String, p1: Point, p2: Point, color: (f64, f64, f64)) {
    let (px, py) = perp_unit(p1, p2);
    for &offset in &[0.0_f64, -3.0, 3.0] {
        eps_line(
            out,
            Point::new(p1.x + px * offset, p1.y + py * offset),
            Point::new(p2.x + px * offset, p2.y + py * offset),
            1.5,
            color,
        );
    }
}

fn eps_aromatic_bond(out: &mut String, p1: Point, p2: Point, color: (f64, f64, f64)) {
    let offset = 2.0;
    let (px, py) = perp_unit(p1, p2);
    eps_line(
        out,
        Point::new(p1.x - px * offset, p1.y - py * offset),
        Point::new(p2.x - px * offset, p2.y - py * offset),
        1.5,
        color,
    );
    eps_dashed_line(
        out,
        Point::new(p1.x + px * offset, p1.y + py * offset),
        Point::new(p2.x + px * offset, p2.y + py * offset),
        1.5,
        color,
    );
}

fn eps_wedge(out: &mut String, p1: Point, p2: Point, color: (f64, f64, f64)) {
    let (px, py) = perp_unit(p1, p2);
    let half_w = 3.0;
    let x2a = p2.x - px * half_w;
    let y2a = p2.y - py * half_w;
    let x2b = p2.x + px * half_w;
    let y2b = p2.y + py * half_w;
    set_color(out, color);
    out.push_str(&format!(
        "{:.2} {:.2} moveto {:.2} {:.2} lineto {:.2} {:.2} lineto closepath fill\n",
        p1.x, p1.y, x2a, y2a, x2b, y2b,
    ));
}

fn eps_dash_bond(out: &mut String, p1: Point, p2: Point, color: (f64, f64, f64)) {
    let dx = p2.x - p1.x;
    let dy = p2.y - p1.y;
    let len = (dx * dx + dy * dy).sqrt();
    if len < 1e-10 {
        return;
    }
    let (px, py) = perp_unit(p1, p2);
    let steps = 6usize;
    set_color(out, color);
    out.push_str("1.0 setlinewidth\n");
    for i in 0..=steps {
        let t = i as f64 / steps as f64;
        let cx = p1.x + t * dx;
        let cy = p1.y + t * dy;
        let hw = t * 3.0 + 0.5;
        out.push_str(&format!(
            "{:.2} {:.2} moveto {:.2} {:.2} lineto stroke\n",
            cx - px * hw,
            cy - py * hw,
            cx + px * hw,
            cy + py * hw,
        ));
    }
}

fn eps_dative_bond(out: &mut String, p1: Point, p2: Point, color: (f64, f64, f64)) {
    eps_line(out, p1, p2, 1.5, color);
    // Arrowhead
    let dx = p2.x - p1.x;
    let dy = p2.y - p1.y;
    let len = (dx * dx + dy * dy).sqrt();
    if len < 1e-10 {
        return;
    }
    let (dxn, dyn_) = (dx / len, dy / len);
    let arrow_len = 8.0;
    let arrow_w = 6.0;
    let tip_x = p2.x - dxn * 2.0;
    let tip_y = p2.y - dyn_ * 2.0;
    let base_x = tip_x - dxn * arrow_len;
    let base_y = tip_y - dyn_ * arrow_len;
    let (px, py) = (-dyn_, dxn);
    let lx = base_x + px * arrow_w / 2.0;
    let ly = base_y + py * arrow_w / 2.0;
    let rx = base_x - px * arrow_w / 2.0;
    let ry = base_y - py * arrow_w / 2.0;
    set_color(out, color);
    out.push_str(&format!(
        "{:.2} {:.2} moveto {:.2} {:.2} lineto {:.2} {:.2} lineto closepath fill\n",
        p2.x, p2.y, lx, ly, rx, ry,
    ));
}

// ---------------------------------------------------------------------------
// Utilities
// ---------------------------------------------------------------------------

fn set_color(out: &mut String, (r, g, b): (f64, f64, f64)) {
    out.push_str(&format!("{:.4} {:.4} {:.4} setrgbcolor\n", r, g, b));
}

fn parse_color_or(s: &str, fallback: (f64, f64, f64)) -> (f64, f64, f64) {
    let s = s.trim().trim_start_matches('#');
    if s.len() == 6
        && let (Ok(r), Ok(g), Ok(b)) = (
            u8::from_str_radix(&s[0..2], 16),
            u8::from_str_radix(&s[2..4], 16),
            u8::from_str_radix(&s[4..6], 16),
        )
    {
        return (r as f64 / 255.0, g as f64 / 255.0, b as f64 / 255.0);
    }
    // Named colors
    match s.to_ascii_lowercase().as_str() {
        "white" => (1.0, 1.0, 1.0),
        "black" => (0.0, 0.0, 0.0),
        "red" => (1.0, 0.0, 0.0),
        "green" => (0.0, 0.502, 0.0),
        "blue" => (0.0, 0.0, 1.0),
        _ => fallback,
    }
}

/// Convert Unicode subscripts (₀–₉) to ASCII digits for PostScript compatibility.
fn unicode_subscripts_to_ascii(s: &str) -> String {
    s.chars()
        .map(|c| match c {
            '₀' => '0',
            '₁' => '1',
            '₂' => '2',
            '₃' => '3',
            '₄' => '4',
            '₅' => '5',
            '₆' => '6',
            '₇' => '7',
            '₈' => '8',
            '₉' => '9',
            other => other,
        })
        .collect()
}

/// Escape PostScript special characters in a label string.
fn escape_ps(s: &str) -> String {
    s.replace('\\', "\\\\")
        .replace('(', "\\(")
        .replace(')', "\\)")
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::layout::compute_layout;
    use chematic_smiles::parse;

    fn mol(s: &str) -> Molecule {
        parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"))
    }

    #[test]
    fn eps_benzene_has_header() {
        let m = mol("c1ccccc1");
        let layout = compute_layout(&m);
        let eps = render_eps(&m, &layout);
        assert!(
            eps.starts_with("%!PS-Adobe-3.0 EPSF-3.0"),
            "missing EPS header"
        );
        assert!(eps.contains("%%BoundingBox:"), "missing BoundingBox");
        assert!(eps.contains("lineto"), "missing bond lines");
        assert!(eps.ends_with("showpage\n%%EOF\n"), "missing EPS footer");
    }

    #[test]
    fn eps_pyridine_contains_label() {
        let m = mol("c1ccncc1");
        let layout = compute_layout(&m);
        let eps = render_eps(&m, &layout);
        assert!(eps.contains("(N)"), "pyridine EPS must contain N label");
    }

    #[test]
    fn eps_ethanol_contains_oh_label() {
        let m = mol("CCO");
        let layout = compute_layout(&m);
        let eps = render_eps(&m, &layout);
        assert!(eps.contains("(OH)"), "ethanol EPS must contain OH label");
    }

    #[test]
    fn eps_double_bond_has_two_lines() {
        let m = mol("C=C");
        let layout = compute_layout(&m);
        let eps = render_eps(&m, &layout);
        let lineto_count = eps.matches("lineto").count();
        assert!(
            lineto_count >= 2,
            "C=C should have >= 2 lineto, got {lineto_count}"
        );
    }

    #[test]
    fn eps_unicode_subscripts_stripped() {
        assert_eq!(unicode_subscripts_to_ascii("H₂O"), "H2O");
        assert_eq!(unicode_subscripts_to_ascii("CH₄"), "CH4");
        assert_eq!(unicode_subscripts_to_ascii("NH₃"), "NH3");
    }

    #[test]
    fn eps_single_carbon_shows_ch4() {
        let m = mol("C");
        let layout = compute_layout(&m);
        let eps = render_eps(&m, &layout);
        assert!(eps.contains("(CH4)"), "single C EPS should show CH4");
    }

    #[test]
    fn eps_co2_shows_middle_carbon_with_a_bond_gap() {
        let m = mol("O=C=O");
        let layout = compute_layout(&m);
        let opts = RenderOptions {
            background: "transparent".into(),
            ..Default::default()
        };
        let eps = render_eps_opts(&m, &layout, &opts);
        assert!(eps.contains("(C) show"));
        assert!(eps.contains("20.00 58.00 moveto 51.00 58.00 lineto"));
        assert!(eps.contains("69.00 58.00 moveto 100.00 58.00 lineto"));
    }
}
