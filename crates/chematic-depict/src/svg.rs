//! SVG serializer for molecular 2D layouts.
//!
//! Converts a `Layout` (atom coordinates) plus a `Molecule` (atoms, bonds)
//! into a self-contained SVG string suitable for embedding in HTML or saving
//! as a `.svg` file.

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule, apply_kekule, kekulize};

use crate::layout::{BOND_LEN, Layout, Point};

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
    /// Per-atom color overrides.  Keys are atom indices; values are CSS colors.
    /// Takes precedence over [`highlight_color`] when both apply to the same atom.
    /// Atoms in this map are highlighted even if absent from [`highlight_atoms`].
    pub atom_color_map: std::collections::HashMap<AtomIdx, String>,
    /// Attach `data-atom-idx`, `data-element`, `data-charge` attributes to
    /// atom label `<text>` elements, enabling JS hover/click tooltips.
    /// Default: `false`.
    pub atom_ids: bool,
    /// Overlay atom index numbers (0-based) as small grey text near each atom.
    /// Default: `false`.
    pub show_atom_indices: bool,
    /// Render aromatic bonds as Kekulé (alternating single/double) instead of
    /// the default dashed aromatic style.
    /// Default: `false`.
    pub kekulize: bool,
}

impl RenderOptions {
    /// Return `RenderOptions` with per-atom CPK colors pre-filled for every
    /// non-carbon atom in `mol`.  Carbon is left black (default label color).
    ///
    /// This is a convenience constructor for quick coloured depictions; all
    /// other options retain their defaults.
    pub fn with_cpk_colors_for(mol: &chematic_core::Molecule) -> Self {
        let mut opts = Self::default();
        for (idx, atom) in mol.atoms() {
            let color = atom_color(atom.element.atomic_number());
            if color != "#000000" {
                opts.atom_color_map.insert(idx, color.to_string());
            }
        }
        opts
    }
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
            atom_color_map: std::collections::HashMap::new(),
            atom_ids: false,
            show_atom_indices: false,
            kekulize: false,
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
    atom_ids: bool,
    show_atom_indices: bool,
}

impl<'a> DrawCtx<'a> {
    fn from_opts(opts: &'a RenderOptions) -> Self {
        let bond_color = if opts.dark { "white" } else { "black" };
        let label_rect_fill = if opts.background == "transparent" {
            None
        } else {
            Some(opts.background.as_str())
        };
        DrawCtx {
            bond_color,
            label_rect_fill,
            dark: opts.dark,
            atom_ids: opts.atom_ids,
            show_atom_indices: opts.show_atom_indices,
        }
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
    let ctx = DrawCtx {
        bond_color: "black",
        label_rect_fill: Some("white"),
        dark: false,
        atom_ids: false,
        show_atom_indices: false,
    };
    let mut body = String::new();
    for (_, bond) in mol.bonds() {
        let p1 = layout.get(bond.atom1);
        let p2 = layout.get(bond.atom2);
        body.push_str(&render_bond_c(bond.order, p1, p2, ctx.bond_color));
    }
    write_atom_labels_ctx(mol, layout, &ctx, &mut body);
    body
}

/// Render bonds, highlight circles, and atom labels without an SVG wrapper.
///
/// Used by the highlighted grid renderer.
pub(crate) fn render_mol_body_opts(
    mol: &Molecule,
    layout: &Layout,
    opts: &RenderOptions,
) -> String {
    let ctx = DrawCtx::from_opts(opts);
    let mut body = String::new();

    // Highlight atom circles (drawn beneath bonds).
    let default_hc = escape_xml(&opts.highlight_color);
    let atom_count = mol.atom_count();
    let mut circles: Vec<(AtomIdx, String)> = Vec::new();
    for idx in &opts.highlight_atoms {
        if idx.0 as usize >= atom_count {
            continue;
        }
        let color = opts
            .atom_color_map
            .get(idx)
            .map(|c| escape_xml(c))
            .unwrap_or_else(|| default_hc.clone());
        circles.push((*idx, color));
    }
    circles.sort_unstable_by_key(|(idx, _)| *idx);
    for (idx, color) in &circles {
        let p = layout.get(*idx);
        body.push_str(&format!(
            "  <circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"16\" fill=\"{}\" opacity=\"0.5\"/>\n",
            p.x, p.y, color,
        ));
    }

    // Bonds.
    for (bond_idx, bond) in mol.bonds() {
        let p1 = layout.get(bond.atom1);
        let p2 = layout.get(bond.atom2);
        if opts.highlight_bonds.contains(&bond_idx) {
            body.push_str(&render_line(p1, p2, "4.0", "#FF8C00"));
        } else {
            body.push_str(&render_bond_c(bond.order, p1, p2, ctx.bond_color));
        }
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
    // R4: cleaner than the uninit-borrow dance; kekulize failure silently falls back to aromatic.
    let maybe_kekule: Option<Molecule> = if opts.kekulize {
        kekulize(mol).ok().map(|kr| apply_kekule(mol, &kr))
    } else {
        None
    };
    let mol: &Molecule = maybe_kekule.as_ref().unwrap_or(mol);

    let ctx = DrawCtx::from_opts(opts);
    let mut svg = String::new();

    write_svg_header_opts(layout, opts, &mut svg);

    // Highlight atom circles (beneath bonds).
    // Union of highlight_atoms and atom_color_map keys; per-atom color overrides highlight_color.
    let atom_count = mol.atom_count();
    let default_hc = escape_xml(&opts.highlight_color);
    let mut atom_circles: Vec<(AtomIdx, String)> = Vec::new();
    for idx in &opts.highlight_atoms {
        if idx.0 as usize >= atom_count {
            continue;
        }
        let color = opts
            .atom_color_map
            .get(idx)
            .map(|c| escape_xml(c))
            .unwrap_or_else(|| default_hc.clone());
        atom_circles.push((*idx, color));
    }
    for (idx, color) in &opts.atom_color_map {
        if idx.0 as usize >= atom_count {
            continue;
        }
        if !opts.highlight_atoms.contains(idx) {
            atom_circles.push((*idx, escape_xml(color)));
        }
    }
    atom_circles.sort_unstable_by_key(|(idx, _)| *idx);
    for (idx, color) in &atom_circles {
        let p = layout.get(*idx);
        svg.push_str(&format!(
            "  <circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"16\" fill=\"{}\" opacity=\"0.5\"/>\n",
            p.x, p.y, color
        ));
    }

    // Bonds.
    for (bond_idx, bond) in mol.bonds() {
        let p1 = layout.get(bond.atom1);
        let p2 = layout.get(bond.atom2);
        if opts.highlight_bonds.contains(&bond_idx) {
            svg.push_str(&render_line(p1, p2, "4.0", "#FF8C00")); // R1
        } else {
            svg.push_str(&render_bond_c(bond.order, p1, p2, ctx.bond_color));
        }
    }

    write_atom_labels_ctx(mol, layout, &ctx, &mut svg);

    svg.push_str("</svg>");
    svg
}

/// Render `mol` with embedded SMILES metadata in the SVG.
///
/// The canonical SMILES is embedded in a `<metadata><smiles>...</smiles></metadata>`
/// element immediately after the opening `<svg>` tag. This allows the structure to
/// be recovered from the image without external data.
pub fn render_svg_with_metadata(
    mol: &Molecule,
    layout: &Layout,
    opts: &RenderOptions,
    smiles: &str,
) -> String {
    use chematic_smiles::canonical_smiles;

    let maybe_kekule: Option<Molecule> = if opts.kekulize {
        kekulize(mol).ok().map(|kr| apply_kekule(mol, &kr))
    } else {
        None
    };
    let mol: &Molecule = maybe_kekule.as_ref().unwrap_or(mol);

    let ctx = DrawCtx::from_opts(opts);
    let mut svg = String::new();

    write_svg_header_opts(layout, opts, &mut svg);

    // Embed metadata with SMILES
    let display_smiles = if smiles.is_empty() {
        canonical_smiles(mol)
    } else {
        smiles.to_string()
    };
    svg.push_str(&format!(
        "  <metadata><smiles>{}</smiles></metadata>\n",
        escape_xml(&display_smiles)
    ));

    // Rest of rendering is identical to render_svg_opts
    let atom_count = mol.atom_count();
    let default_hc = escape_xml(&opts.highlight_color);
    let mut atom_circles: Vec<(AtomIdx, String)> = Vec::new();
    for idx in &opts.highlight_atoms {
        if idx.0 as usize >= atom_count {
            continue;
        }
        let color = opts
            .atom_color_map
            .get(idx)
            .map(|c| escape_xml(c))
            .unwrap_or_else(|| default_hc.clone());
        atom_circles.push((*idx, color));
    }
    for (idx, color) in &opts.atom_color_map {
        if idx.0 as usize >= atom_count {
            continue;
        }
        if !opts.highlight_atoms.contains(idx) {
            atom_circles.push((*idx, escape_xml(color)));
        }
    }
    atom_circles.sort_unstable_by_key(|(idx, _)| *idx);
    for (idx, color) in &atom_circles {
        let p = layout.get(*idx);
        svg.push_str(&format!(
            "  <circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"16\" fill=\"{}\" opacity=\"0.5\"/>\n",
            p.x, p.y, color
        ));
    }

    for (bond_idx, bond) in mol.bonds() {
        let p1 = layout.get(bond.atom1);
        let p2 = layout.get(bond.atom2);
        if opts.highlight_bonds.contains(&bond_idx) {
            svg.push_str(&render_line(p1, p2, "4.0", "#FF8C00"));
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
            view_x,
            view_y,
            view_w,
            view_h,
            escape_xml(&opts.background) // S1: escape user input
        ));
    }
}

// ---------------------------------------------------------------------------
// Atom labels
// ---------------------------------------------------------------------------

fn write_atom_labels_ctx(mol: &Molecule, layout: &Layout, ctx: &DrawCtx, svg: &mut String) {
    for (idx, atom) in mol.atoms() {
        let label = atom_label(mol, idx);
        let p = layout.get(idx);

        if !label.is_empty() {
            if let Some(fill) = ctx.label_rect_fill {
                svg.push_str(&format!(
                    "  <rect x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" fill=\"{}\"/>\n",
                    p.x - LABEL_HALF_W,
                    p.y - LABEL_HALF_H,
                    LABEL_HALF_W * 2.0,
                    LABEL_HALF_H * 2.0,
                    escape_xml(fill), // S1: escape background colour
                ));
            }

            // R2: build optional data attrs, then single format!
            let data_attrs = if ctx.atom_ids {
                format!(
                    " data-atom-idx=\"{}\" data-element=\"{}\" data-charge=\"{}\"",
                    idx.0,
                    atom.element.symbol(),
                    atom.charge
                )
            } else {
                String::new()
            };
            svg.push_str(&format!(
                "  <text x=\"{:.2}\" y=\"{:.2}\" \
                 font-family=\"sans-serif\" font-size=\"{}\" \
                 text-anchor=\"middle\" dominant-baseline=\"central\" \
                 fill=\"{}\"{}>{}</text>\n",
                p.x,
                p.y,
                FONT_SIZE as u32,
                ctx.text_color(atom.element.atomic_number()),
                data_attrs,
                escape_xml(&label)
            ));
        } else if ctx.atom_ids {
            // B3: unlabelled atoms (plain carbons) get an invisible anchor so JS can address them.
            svg.push_str(&format!(
                "  <text x=\"{:.2}\" y=\"{:.2}\" font-size=\"0\" \
                 data-atom-idx=\"{}\" data-element=\"{}\" data-charge=\"{}\"/>\n",
                p.x,
                p.y,
                idx.0,
                atom.element.symbol(),
                atom.charge
            ));
        }

        // Atom index overlay (all atoms, including unlabelled carbons).
        if ctx.show_atom_indices {
            svg.push_str(&format!(
                "  <text x=\"{:.2}\" y=\"{:.2}\" \
                 font-family=\"sans-serif\" font-size=\"8\" \
                 text-anchor=\"start\" dominant-baseline=\"auto\" \
                 fill=\"#8b92a9\">{}</text>\n",
                p.x + LABEL_HALF_W,
                p.y - LABEL_HALF_H,
                idx.0
            ));
        }
    }
}

// ---------------------------------------------------------------------------
// Bond rendering
// ---------------------------------------------------------------------------

fn render_bond_c(order: BondOrder, p1: Point, p2: Point, color: &str) -> String {
    match order {
        BondOrder::Single => render_line(p1, p2, "1.5", color),
        BondOrder::Up => render_wedge_up(p1, p2, color),
        BondOrder::Down => render_dash_bond(p1, p2, color),
        BondOrder::Double => render_double_bond(p1, p2, color),
        BondOrder::Triple => render_triple_bond(p1, p2, color),
        BondOrder::Aromatic => render_aromatic_bond(p1, p2, color),
        BondOrder::Quadruple => render_line(p1, p2, "3.0", color),
        BondOrder::Zero => render_line(p1, p2, "1.0", color), // Zero-order as thin line
        BondOrder::Dative => render_dative_bond(p1, p2, color),
        BondOrder::QueryAny => render_query_any_bond(p1, p2, color),
        BondOrder::QuerySingleOrDouble => render_query_dashed_bond(p1, p2, color, "2,2"),
        BondOrder::QuerySingleOrAromatic => render_query_dashed_bond(p1, p2, color, "4,2"),
        BondOrder::QueryDoubleOrAromatic => render_query_dashed_bond(p1, p2, color, "3,3"),
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
        s.push_str(&render_line(
            Point::new(p1.x + px * offset * sign, p1.y + py * offset * sign),
            Point::new(p2.x + px * offset * sign, p2.y + py * offset * sign),
            "1.5",
            color,
        ));
    }
    s
}

fn render_triple_bond(p1: Point, p2: Point, color: &str) -> String {
    let (px, py) = perp_unit(p1, p2);
    let mut s = String::new();
    for &offset in &[0.0_f64, -3.0, 3.0] {
        s.push_str(&render_line(
            Point::new(p1.x + px * offset, p1.y + py * offset),
            Point::new(p2.x + px * offset, p2.y + py * offset),
            "1.5",
            color,
        ));
    }
    s
}

fn render_line_dashed(
    p1: Point,
    p2: Point,
    stroke_width: &str,
    color: &str,
    dasharray: &str,
) -> String {
    format!(
        "  <line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" \
         stroke=\"{}\" stroke-width=\"{}\" fill=\"none\" stroke-dasharray=\"{}\"/>\n",
        p1.x, p1.y, p2.x, p2.y, color, stroke_width, dasharray
    )
}

fn render_aromatic_bond(p1: Point, p2: Point, color: &str) -> String {
    let offset = 2.0;
    let (px, py) = perp_unit(p1, p2);
    let mut s = String::new();
    s.push_str(&render_line(
        Point::new(p1.x - px * offset, p1.y - py * offset),
        Point::new(p2.x - px * offset, p2.y - py * offset),
        "1.5",
        color,
    ));
    s.push_str(&render_line_dashed(
        Point::new(p1.x + px * offset, p1.y + py * offset),
        Point::new(p2.x + px * offset, p2.y + py * offset),
        "1.5",
        color,
        "4,3",
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
        s.push_str(&render_line(
            Point::new(cx - px * hw, cy - py * hw),
            Point::new(cx + px * hw, cy + py * hw),
            "1.0",
            color,
        ));
    }
    s
}

/// Render a dative bond (coordinate covalent) as an arrow from donor to acceptor.
/// Typically indicates a bond where both electrons come from one atom.
fn render_dative_bond(p1: Point, p2: Point, color: &str) -> String {
    let dx = p2.x - p1.x;
    let dy = p2.y - p1.y;
    let len = (dx * dx + dy * dy).sqrt();
    if len < 1e-10 {
        return String::new();
    }

    let mut s = String::new();

    // Draw the line from p1 to p2
    s.push_str(&render_line(p1, p2, "1.5", color));

    // Draw arrowhead at p2 pointing in the direction of the bond
    let arrow_len = 8.0;
    let arrow_width = 6.0;

    // Normalize direction
    let dx_norm = dx / len;
    let dy_norm = dy / len;

    // Arrowhead point is slightly before p2
    let arrow_tip_x = p2.x - dx_norm * 2.0;
    let arrow_tip_y = p2.y - dy_norm * 2.0;

    // Perpendicular direction
    let perp_x = -dy_norm;
    let perp_y = dx_norm;

    // Arrow base
    let base_x = arrow_tip_x - dx_norm * arrow_len;
    let base_y = arrow_tip_y - dy_norm * arrow_len;

    // Arrowhead polygon (filled triangle)
    let left_x = base_x + perp_x * arrow_width / 2.0;
    let left_y = base_y + perp_y * arrow_width / 2.0;
    let right_x = base_x - perp_x * arrow_width / 2.0;
    let right_y = base_y - perp_y * arrow_width / 2.0;

    s.push_str(&format!(
        "  <polygon points=\"{:.2},{:.2} {:.2},{:.2} {:.2},{:.2}\" \
         fill=\"{}\" stroke=\"{}\" stroke-width=\"0.5\"/>\n",
        p2.x, p2.y, left_x, left_y, right_x, right_y, color, color
    ));

    s
}

/// Render a query bond for "any bond type" as a dotted line with asterisks.
fn render_query_any_bond(p1: Point, p2: Point, color: &str) -> String {
    let dx = p2.x - p1.x;
    let dy = p2.y - p1.y;
    let len = (dx * dx + dy * dy).sqrt();
    if len < 1e-10 {
        return String::new();
    }

    let mut s = String::new();

    // Draw main dotted line
    s.push_str(&render_line_dashed(p1, p2, "1.5", color, "3,3"));

    // Add asterisk marker in the middle
    let mid_x = (p1.x + p2.x) / 2.0;
    let mid_y = (p1.y + p2.y) / 2.0;

    s.push_str(&format!(
        "  <text x=\"{:.2}\" y=\"{:.2}\" font-family=\"serif\" font-size=\"14\" \
         text-anchor=\"middle\" dominant-baseline=\"central\" fill=\"{}\" \
         font-weight=\"bold\">*</text>\n",
        mid_x, mid_y, color
    ));

    s
}

/// Render query bonds with different dash patterns to distinguish bond types.
/// - QuerySingleOrDouble: short dashes (2,2)
/// - QuerySingleOrAromatic: medium dashes (4,2)
/// - QueryDoubleOrAromatic: long dashes (3,3)
fn render_query_dashed_bond(p1: Point, p2: Point, color: &str, dasharray: &str) -> String {
    render_line_dashed(p1, p2, "1.5", color, dasharray)
}

// ---------------------------------------------------------------------------
// Atom coloring (CPK palette)
// ---------------------------------------------------------------------------

/// CPK standard color for an element (black for carbon/hydrogen/default).
pub fn atom_color(atomic_number: u8) -> &'static str {
    match atomic_number {
        7 => "#3050F8",  // N  blue
        8 => "#FF0D0D",  // O  red
        16 => "#FFFF30", // S  yellow
        17 => "#1FF01F", // Cl green
        9 => "#90E050",  // F  light-green
        35 => "#A62929", // Br dark-red/brown
        53 => "#940094", // I  purple
        15 => "#FF8000", // P  orange
        _ => "#000000",  // default black
    }
}

/// CPK color for `atomic_number` as an `[r, g, b]` byte triple.
///
/// Values match [`atom_color`]; elements not in the CPK table return `[0, 0, 0]`.
/// Use this instead of parsing the hex string when integrating with a GUI
/// framework (e.g. `egui::Color32::from_rgb(r, g, b)`).
pub fn atom_color_rgb(atomic_number: u8) -> [u8; 3] {
    match atomic_number {
        7 => [0x30, 0x50, 0xF8],  // N  blue
        8 => [0xFF, 0x0D, 0x0D],  // O  red
        16 => [0xFF, 0xFF, 0x30], // S  yellow
        17 => [0x1F, 0xF0, 0x1F], // Cl green
        9 => [0x90, 0xE0, 0x50],  // F  light-green
        35 => [0xA6, 0x29, 0x29], // Br brown
        53 => [0x94, 0x00, 0x94], // I  purple
        15 => [0xFF, 0x80, 0x00], // P  orange
        _ => [0x00, 0x00, 0x00],  // default black
    }
}

// ---------------------------------------------------------------------------
// Atom labels
// ---------------------------------------------------------------------------

/// Format a formal charge as a suffix string.
/// Centralised here so the character used for minus is consistent everywhere (B5, R3).
fn format_charge(c: i8) -> String {
    match c {
        0 => String::new(),
        1 => "+".to_string(),
        -1 => "-".to_string(),
        c if c > 1 => format!("{c}+"),
        c => format!("{}-", -c), // e.g. "2-", "3-"
    }
}

/// Compute the display label for an atom.
///
/// Isolated atoms (single-atom molecule or degree-0 in a multi-atom molecule)
/// use Hill-notation molecular-formula style: H2O, CH4, NH3.
/// Plain carbons in bonded multi-atom molecules return "".
fn atom_label(mol: &Molecule, idx: AtomIdx) -> String {
    let atom = mol.atom(idx);
    let is_carbon = atom.element.atomic_number() == 6;
    let has_charge = atom.charge != 0;
    let has_isotope = atom.isotope.is_some();

    // Isolated atom: use molecular-formula style (H2O, CH4, etc.).
    // This covers both single-atom molecules and degree-0 atoms in disconnected SMILES (B4).
    if mol.atom_count() == 1 || mol.degree(idx) == 0 {
        let h = chematic_core::implicit_hcount(mol, idx);
        return build_isolated_label(
            atom.element.symbol(),
            atom.element.atomic_number(),
            h,
            atom.charge,
        );
    }

    // Plain carbon in a bonded multi-atom molecule: no label in skeletal structure.
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

    // Charge (R3: use shared helper; B5: consistent ASCII minus).
    if has_charge {
        label.push_str(&format_charge(atom.charge));
    }

    label
}

/// Build a molecular-formula-style label for an isolated atom.
///
/// Hill notation: C first, then H, then alphabetical.
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
    format!("{}{}", base, format_charge(charge)) // R3, B5
}

/// Position where H count should be displayed relative to the atom symbol.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HPosition {
    Left,
    Right,
    Up,
    Down,
}

/// Hydrogen count and position for condensed notation (e.g., "CH₃", "NH₂").
#[derive(Debug, Clone)]
pub struct AtomLabel {
    pub symbol: String,
    pub h_count: u8,
    pub h_position: HPosition,
}

/// Return the structured display label for an atom with implicit H position.
///
/// Terminal atoms (degree 1) and isolated atoms show condensed notation:
/// - "CH₃" for terminal carbon
/// - "NH₂" for terminal nitrogen
/// - "OH" for terminal oxygen
///
/// Interior carbons (degree ≥ 2) with no charge/isotope return empty label.
pub fn atom_label_with_h(mol: &Molecule, idx: AtomIdx) -> AtomLabel {
    let atom = mol.atom(idx);
    let is_carbon = atom.element.atomic_number() == 6;
    let degree = mol.degree(idx);

    // Single-atom molecule: always show as molecular formula
    if mol.atom_count() == 1 {
        let h = chematic_core::implicit_hcount(mol, idx);
        return AtomLabel {
            symbol: atom.element.symbol().to_string(),
            h_count: h,
            h_position: HPosition::Right,
        };
    }

    // Isolated degree-0 atom: show label
    if degree == 0 {
        let h = chematic_core::implicit_hcount(mol, idx);
        return AtomLabel {
            symbol: atom.element.symbol().to_string(),
            h_count: h,
            h_position: HPosition::Right,
        };
    }

    // Interior carbon with no charge/isotope: no label (skeletal structure)
    if is_carbon && atom.charge == 0 && atom.isotope.is_none() && degree >= 2 {
        return AtomLabel {
            symbol: String::new(),
            h_count: 0,
            h_position: HPosition::Right,
        };
    }

    // Terminal atom (degree 1): show with implicit H count
    if degree == 1 {
        let h = chematic_core::implicit_hcount(mol, idx);
        return AtomLabel {
            symbol: atom.element.symbol().to_string(),
            h_count: h,
            h_position: HPosition::Right,
        };
    }

    // Non-terminal non-carbon with implicit H: show element symbol and H count
    if !is_carbon {
        let h = chematic_core::implicit_hcount(mol, idx);
        return AtomLabel {
            symbol: atom.element.symbol().to_string(),
            h_count: h,
            h_position: HPosition::Right,
        };
    }

    // Fallback: carbon with charge or isotope
    AtomLabel {
        symbol: atom.element.symbol().to_string(),
        h_count: 0,
        h_position: HPosition::Right,
    }
}

/// Return the display label string for an atom in condensed notation.
///
/// Examples:
/// - Terminal carbon with 3 H: "CH₃"
/// - Non-terminal nitrogen with 2 H: "NH₂"
/// - Hydroxyl (O with 1 H): "OH"
/// - Isolated oxygen with 2 H: "H₂O" (Hill notation)
/// - Interior carbon in skeleton: "" (empty)
pub fn atom_display_label(mol: &Molecule, idx: AtomIdx) -> String {
    let label = atom_label_with_h(mol, idx);

    if label.symbol.is_empty() {
        return String::new();
    }

    if label.h_count == 0 {
        return label.symbol;
    }

    // For single-atom molecules, use Hill notation (H first for non-carbon)
    let is_isolated = mol.atom_count() == 1;
    let is_carbon = mol.atom(idx).element.atomic_number() == 6;

    if is_isolated && !is_carbon {
        // Hill notation: H first for non-carbon atoms (e.g., H₂O not OH₂)
        let h_subscript = match label.h_count {
            1 => "H".to_string(),
            2 => "H₂".to_string(),
            3 => "H₃".to_string(),
            4 => "H₄".to_string(),
            5 => "H₅".to_string(),
            n => format!("H{}", n),
        };
        return format!("{}{}", h_subscript, label.symbol);
    }

    // For other cases, element symbol followed by H subscript
    let h_subscript = match label.h_count {
        1 => "H".to_string(),
        2 => "H₂".to_string(),
        3 => "H₃".to_string(),
        4 => "H₄".to_string(),
        5 => "H₅".to_string(),
        n => format!("H{}", n),
    };

    format!("{}{}", label.symbol, h_subscript)
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
}
