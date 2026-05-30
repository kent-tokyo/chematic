//! SVG serializer for molecular 2D layouts.
//!
//! Converts a `Layout` (atom coordinates) plus a `Molecule` (atoms, bonds)
//! into a self-contained SVG string suitable for embedding in HTML or saving
//! as a `.svg` file.

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule};

use crate::layout::{Layout, Point, BOND_LEN};

/// Font size used for atom labels, in SVG pixels.
const FONT_SIZE: f64 = 12.0;

/// Approximate half-width of a label background rectangle.
const LABEL_HALF_W: f64 = 8.0;

/// Approximate half-height of a label background rectangle.
const LABEL_HALF_H: f64 = 7.0;

// ---------------------------------------------------------------------------
// Public options
// ---------------------------------------------------------------------------

/// Options controlling SVG rendering style.
#[derive(Clone, Debug)]
pub struct RenderOptions {
    /// Override SVG `width` attribute (px). `None` = auto from bounding box.
    pub width: Option<u32>,
    /// Override SVG `height` attribute (px). `None` = auto from bounding box.
    pub height: Option<u32>,
    /// Padding around the molecule bounding box (SVG user units, default 20.0).
    pub padding: f64,
    /// Background fill color.  `"transparent"` suppresses the background rect
    /// and also removes per-label background rectangles so bonds show through.
    /// Default: `"white"`.
    pub background: String,
    /// Dark-theme mode: bonds rendered in white, carbon text in white.
    /// Default: `false`.
    pub dark: bool,
    /// Atom indices to highlight (yellow circle background by default).
    pub highlight_atoms: std::collections::HashSet<AtomIdx>,
    /// Bond indices to highlight (orange thick line by default).
    pub highlight_bonds: std::collections::HashSet<BondIdx>,
    /// Highlight color (CSS hex, default `"#FFFF00"`).
    pub highlight_color: String,
}

impl Default for RenderOptions {
    fn default() -> Self {
        Self {
            width: None,
            height: None,
            padding: 20.0,
            background: "white".into(),
            dark: false,
            highlight_atoms: std::collections::HashSet::new(),
            highlight_bonds: std::collections::HashSet::new(),
            highlight_color: "#FFFF00".into(),
        }
    }
}

// ---------------------------------------------------------------------------
// Internal draw context
// ---------------------------------------------------------------------------

struct DrawCtx<'a> {
    bond_color: &'a str,
    label_rect_fill: Option<&'a str>, // None = skip background rect
    dark: bool,
}

impl<'a> DrawCtx<'a> {
    fn from_opts(opts: &'a RenderOptions) -> Self {
        let bond_color = if opts.dark { "white" } else { "black" };
        let label_rect_fill = if opts.background == "transparent" {
            None
        } else {
            Some(opts.background.as_str())
        };
        DrawCtx { bond_color, label_rect_fill, dark: opts.dark }
    }

    fn text_color(&self, atomic_number: u8) -> &str {
        if self.dark && atomic_number == 6 {
            "white"
        } else {
            atom_color(atomic_number)
        }
    }
}

// ---------------------------------------------------------------------------
// Entry points
// ---------------------------------------------------------------------------

/// Render just the bonds and atom labels for `mol` without an SVG wrapper.
///
/// Used by the grid renderer to compose multiple molecules into one SVG.
pub(crate) fn render_mol_body(mol: &Molecule, layout: &Layout) -> String {
    let ctx = DrawCtx { bond_color: "black", label_rect_fill: Some("white"), dark: false };
    let mut body = String::new();
    for (_, bond) in mol.bonds() {
        let p1 = layout.get(bond.atom1);
        let p2 = layout.get(bond.atom2);
        body.push_str(&render_bond_c(bond.order, p1, p2, ctx.bond_color));
    }
    write_atom_labels_ctx(mol, layout, &ctx, &mut body);
    body
}

/// Render `mol` with the given `layout` as a self-contained SVG string.
pub fn render_svg(mol: &Molecule, layout: &Layout) -> String {
    render_svg_opts(mol, layout, &RenderOptions::default())
}

/// Render `mol` with highlighted atoms and bonds.
///
/// An empty `highlight_atoms`/`highlight_bonds` produces the same output
/// as [`render_svg`].
pub fn render_svg_highlighted(
    mol: &Molecule,
    layout: &Layout,
    highlight_atoms: &std::collections::HashSet<AtomIdx>,
    highlight_bonds: &std::collections::HashSet<BondIdx>,
) -> String {
    let opts = RenderOptions {
        highlight_atoms: highlight_atoms.clone(),
        highlight_bonds: highlight_bonds.clone(),
        ..RenderOptions::default()
    };
    render_svg_opts(mol, layout, &opts)
}

/// Render `mol` with full control over style via [`RenderOptions`].
pub fn render_svg_opts(mol: &Molecule, layout: &Layout, opts: &RenderOptions) -> String {
    let ctx = DrawCtx::from_opts(opts);
    let mut svg = String::new();

    write_svg_header_opts(layout, opts, &mut svg);

    // Highlight atom circles (beneath bonds).
    for idx in &opts.highlight_atoms {
        let p = layout.get(*idx);
        svg.push_str(&format!(
            "  <circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"16\" fill=\"{}\" opacity=\"0.5\"/>\n",
            p.x, p.y, opts.highlight_color
        ));
    }

    // Bonds.
    for (bond_idx, bond) in mol.bonds() {
        let p1 = layout.get(bond.atom1);
        let p2 = layout.get(bond.atom2);
        if opts.highlight_bonds.contains(&bond_idx) {
            svg.push_str(&format!(
                "  <line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" \
                 stroke=\"#FF8C00\" stroke-width=\"4.0\" fill=\"none\"/>\n",
                p1.x, p1.y, p2.x, p2.y
            ));
        } else {
            svg.push_str(&render_bond_c(bond.order, p1, p2, ctx.bond_color));
        }
    }

    write_atom_labels_ctx(mol, layout, &ctx, &mut svg);

    svg.push_str("</svg>");
    svg
}

// ---------------------------------------------------------------------------
// SVG header
// ---------------------------------------------------------------------------

fn write_svg_header_opts(layout: &Layout, opts: &RenderOptions, svg: &mut String) {
    let padding = opts.padding;
    let (min_x, min_y, max_x, max_y) = layout.bounding_box();

    let raw_w = (max_x - min_x).max(BOND_LEN);
    let raw_h = (max_y - min_y).max(BOND_LEN);

    let view_x = min_x - padding;
    let view_y = min_y - padding;
    let view_w = raw_w + 2.0 * padding;
    let view_h = raw_h + 2.0 * padding;

    let display_w = opts.width.unwrap_or_else(|| view_w.round() as u32);
    let display_h = opts.height.unwrap_or_else(|| view_h.round() as u32);

    svg.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" \
         width=\"{}\" height=\"{}\" \
         viewBox=\"{:.2} {:.2} {:.2} {:.2}\">\n",
        display_w, display_h, view_x, view_y, view_w, view_h
    ));

    if opts.background != "transparent" {
        svg.push_str(&format!(
            "  <rect x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" fill=\"{}\"/>\n",
            view_x, view_y, view_w, view_h, opts.background
        ));
    }
}

// ---------------------------------------------------------------------------
// Atom labels
// ---------------------------------------------------------------------------

fn write_atom_labels_ctx(mol: &Molecule, layout: &Layout, ctx: &DrawCtx, svg: &mut String) {
    for (idx, _atom) in mol.atoms() {
        let label = atom_label(mol, idx);
        if label.is_empty() {
            continue;
        }
        let p = layout.get(idx);

        if let Some(fill) = ctx.label_rect_fill {
            svg.push_str(&format!(
                "  <rect x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" fill=\"{}\"/>\n",
                p.x - LABEL_HALF_W,
                p.y - LABEL_HALF_H,
                LABEL_HALF_W * 2.0,
                LABEL_HALF_H * 2.0,
                fill,
            ));
        }

        svg.push_str(&format!(
            "  <text x=\"{:.2}\" y=\"{:.2}\" \
             font-family=\"sans-serif\" font-size=\"{}\" \
             text-anchor=\"middle\" dominant-baseline=\"central\" \
             fill=\"{}\">{}</text>\n",
            p.x,
            p.y,
            FONT_SIZE as u32,
            ctx.text_color(mol.atom(idx).element.atomic_number()),
            escape_xml(&label)
        ));
    }
}

// ---------------------------------------------------------------------------
// Bond rendering
// ---------------------------------------------------------------------------

fn render_bond_c(order: BondOrder, p1: Point, p2: Point, color: &str) -> String {
    match order {
        BondOrder::Single   => render_line(p1, p2, "1.5", color),
        BondOrder::Up       => render_wedge_up(p1, p2, color),
        BondOrder::Down     => render_dash_bond(p1, p2, color),
        BondOrder::Double   => render_double_bond(p1, p2, color),
        BondOrder::Triple   => render_triple_bond(p1, p2, color),
        BondOrder::Aromatic => render_aromatic_bond(p1, p2, color),
        BondOrder::Quadruple => render_line(p1, p2, "3.0", color),
    }
}

fn render_line(p1: Point, p2: Point, stroke_width: &str, color: &str) -> String {
    format!(
        "  <line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" \
         stroke=\"{}\" stroke-width=\"{}\" fill=\"none\"/>\n",
        p1.x, p1.y, p2.x, p2.y, color, stroke_width
    )
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

fn render_double_bond(p1: Point, p2: Point, color: &str) -> String {
    let offset = 2.0;
    let (px, py) = perp_unit(p1, p2);
    let mut s = String::new();
    for sign in [-1.0_f64, 1.0] {
        let ox = px * offset * sign;
        let oy = py * offset * sign;
        s.push_str(&format!(
            "  <line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" \
             stroke=\"{}\" stroke-width=\"1.5\" fill=\"none\"/>\n",
            p1.x + ox, p1.y + oy, p2.x + ox, p2.y + oy, color
        ));
    }
    s
}

fn render_triple_bond(p1: Point, p2: Point, color: &str) -> String {
    let (px, py) = perp_unit(p1, p2);
    let mut s = String::new();
    for &offset in &[0.0_f64, -3.0, 3.0] {
        let ox = px * offset;
        let oy = py * offset;
        s.push_str(&format!(
            "  <line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" \
             stroke=\"{}\" stroke-width=\"1.5\" fill=\"none\"/>\n",
            p1.x + ox, p1.y + oy, p2.x + ox, p2.y + oy, color
        ));
    }
    s
}

fn render_aromatic_bond(p1: Point, p2: Point, color: &str) -> String {
    let offset = 2.0;
    let (px, py) = perp_unit(p1, p2);
    let mut s = String::new();
    s.push_str(&format!(
        "  <line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" \
         stroke=\"{}\" stroke-width=\"1.5\" fill=\"none\"/>\n",
        p1.x - px * offset, p1.y - py * offset,
        p2.x - px * offset, p2.y - py * offset,
        color
    ));
    s.push_str(&format!(
        "  <line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" \
         stroke=\"{}\" stroke-width=\"1.5\" fill=\"none\" stroke-dasharray=\"4,3\"/>\n",
        p1.x + px * offset, p1.y + py * offset,
        p2.x + px * offset, p2.y + py * offset,
        color
    ));
    s
}

fn render_wedge_up(p1: Point, p2: Point, color: &str) -> String {
    let (px, py) = perp_unit(p1, p2);
    let half_w = 3.0;
    let x2a = p2.x - px * half_w;
    let y2a = p2.y - py * half_w;
    let x2b = p2.x + px * half_w;
    let y2b = p2.y + py * half_w;
    format!(
        "  <polygon points=\"{:.2},{:.2} {:.2},{:.2} {:.2},{:.2}\" \
         fill=\"{}\" stroke=\"{}\" stroke-width=\"0.5\"/>\n",
        p1.x, p1.y, x2a, y2a, x2b, y2b, color, color
    )
}

fn render_dash_bond(p1: Point, p2: Point, color: &str) -> String {
    let dx = p2.x - p1.x;
    let dy = p2.y - p1.y;
    let len = (dx * dx + dy * dy).sqrt();
    if len < 1e-10 {
        return String::new();
    }
    let (px, py) = perp_unit(p1, p2);
    let steps = 6usize;
    let mut s = String::new();
    for i in 0..=steps {
        let t = i as f64 / steps as f64;
        let cx = p1.x + t * dx;
        let cy = p1.y + t * dy;
        let hw = t * 3.0 + 0.5;
        s.push_str(&format!(
            "  <line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" \
             stroke=\"{}\" stroke-width=\"1.0\" fill=\"none\"/>\n",
            cx - px * hw, cy - py * hw,
            cx + px * hw, cy + py * hw,
            color
        ));
    }
    s
}

// ---------------------------------------------------------------------------
// Atom coloring (CPK palette)
// ---------------------------------------------------------------------------

fn atom_color(atomic_number: u8) -> &'static str {
    match atomic_number {
        7  => "#3050F8", // N  blue
        8  => "#FF0D0D", // O  red
        16 => "#FFFF30", // S  yellow
        17 => "#1FF01F", // Cl green
        9  => "#90E050", // F  light-green
        35 => "#A62929", // Br dark-red/brown
        53 => "#940094", // I  purple
        15 => "#FF8000", // P  orange
        _  => "#000000", // default black
    }
}

// ---------------------------------------------------------------------------
// Atom labels
// ---------------------------------------------------------------------------

/// Compute the display label for an atom.
///
/// Isolated single-atom molecules use molecular-formula style (H2O, CH4, NH3).
/// In multi-atom molecules, plain carbons (no charge, no isotope) return "".
fn atom_label(mol: &Molecule, idx: AtomIdx) -> String {
    let atom = mol.atom(idx);
    let is_carbon = atom.element.atomic_number() == 6;
    let has_charge = atom.charge != 0;
    let has_isotope = atom.isotope.is_some();

    // Single-atom molecule: display as molecular formula (H2O, CH4, NH3...).
    if mol.atom_count() == 1 {
        let h = chematic_core::implicit_hcount(mol, idx);
        return build_isolated_label(
            atom.element.symbol(),
            atom.element.atomic_number(),
            h,
            atom.charge,
        );
    }

    // Plain carbon in a multi-atom molecule: no label in skeletal structure.
    if is_carbon && !has_charge && !has_isotope {
        return String::new();
    }

    let mut label = atom.element.symbol().to_string();

    // Implicit H count for non-carbon atoms.
    if !is_carbon {
        let h = chematic_core::implicit_hcount(mol, idx);
        if h == 1 {
            label.push('H');
        } else if h > 1 {
            label.push('H');
            label.push_str(&h.to_string());
        }
    }

    // Charge.
    if has_charge {
        let c = atom.charge;
        if c == 1 {
            label.push('+');
        } else if c == -1 {
            label.push('-');
        } else if c > 1 {
            label.push_str(&format!("{c}+"));
        } else {
            label.push_str(&format!("{}−", -c));
        }
    }

    label
}

/// Build a molecular-formula-style label for an isolated single atom.
///
/// Uses Hill notation convention: C first, then H, then alphabetical.
/// Examples: C→"CH4", O→"H2O", N→"H3N", S→"H2S", noble gas→"He"
fn build_isolated_label(symbol: &str, atomic_number: u8, h: u8, charge: i8) -> String {
    let base = match atomic_number {
        6 => match h {
            0 => symbol.to_string(),
            1 => format!("{}H", symbol),
            n => format!("{}H{}", symbol, n),
        },
        _ => match h {
            0 => symbol.to_string(),
            1 => format!("H{}", symbol),
            n => format!("H{}{}", n, symbol),
        },
    };

    if charge == 0 {
        return base;
    }
    let charge_str = match charge {
        1  => "+".to_string(),
        -1 => "-".to_string(),
        c if c > 1 => format!("{c}+"),
        c => format!("{}−", -c),
    };
    format!("{}{}", base, charge_str)
}

/// Escape XML special characters in a label string.
fn escape_xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
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
        assert!(svg.contains("H2O"), "water 'O' should show H2O label, got: {}", svg);
    }

    #[test]
    fn isolated_methane_label_ch4() {
        let m = mol("C");
        let layout = crate::layout::compute_layout(&m);
        let svg = render_svg(&m, &layout);
        assert!(svg.contains("CH4"), "methane 'C' should show CH4 label, got: {}", svg);
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
        let opts = RenderOptions { background: "transparent".into(), ..Default::default() };
        let svg = render_svg_opts(&m, &layout, &opts);
        assert!(svg.contains("<svg"), "must be valid SVG");
        // No background rect when transparent.
        assert!(!svg.contains("fill=\"transparent\""), "no bg rect fill for transparent");
    }

    #[test]
    fn render_opts_custom_size() {
        let m = mol("CCO");
        let layout = crate::layout::compute_layout(&m);
        let opts = RenderOptions { width: Some(300), height: Some(200), ..Default::default() };
        let svg = render_svg_opts(&m, &layout, &opts);
        assert!(svg.contains("width=\"300\""), "SVG width should be 300");
        assert!(svg.contains("height=\"200\""), "SVG height should be 200");
    }

    #[test]
    fn render_opts_dark_theme_white_bonds() {
        let m = mol("CC");
        let layout = crate::layout::compute_layout(&m);
        let opts = RenderOptions { dark: true, background: "#0f172a".into(), ..Default::default() };
        let svg = render_svg_opts(&m, &layout, &opts);
        assert!(svg.contains("stroke=\"white\""), "dark theme bonds should be white");
    }

    #[test]
    fn render_opts_highlight_atoms() {
        let m = mol("c1ccncc1");
        let layout = crate::layout::compute_layout(&m);
        let n_idx = m.atoms()
            .find(|(_, a)| a.element.atomic_number() == 7)
            .map(|(idx, _)| idx)
            .expect("pyridine has N");
        let mut hl = std::collections::HashSet::new();
        hl.insert(n_idx);
        let opts = RenderOptions { highlight_atoms: hl, ..Default::default() };
        let svg = render_svg_opts(&m, &layout, &opts);
        assert!(svg.contains("<circle"), "highlight must produce a circle");
    }
}
