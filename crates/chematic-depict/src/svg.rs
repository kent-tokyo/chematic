//! SVG serializer for molecular 2D layouts.
//!
//! Converts a `Layout` (atom coordinates) plus a `Molecule` (atoms, bonds)
//! into a self-contained SVG string suitable for embedding in HTML or saving
//! as a `.svg` file.

use chematic_core::{AtomIdx, BondOrder, Molecule};
use chematic_perception::find_sssr;

use crate::layout::{Layout, Point, BOND_LEN};

/// Padding around the bounding box, in SVG pixels.
const PADDING: f64 = 20.0;

/// Font size used for atom labels, in SVG pixels.
const FONT_SIZE: f64 = 12.0;

/// Approximate half-width of a label background rectangle.
const LABEL_HALF_W: f64 = 8.0;

/// Approximate half-height of a label background rectangle.
const LABEL_HALF_H: f64 = 7.0;

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Render `mol` with the given `layout` as a self-contained SVG string.
///
/// Rendering order:
/// 1. Bond lines (drawn first, beneath labels).
/// 2. White background rectangles for atom labels (drawn over bond lines).
/// 3. Atom label text (drawn on top).
pub fn render_svg(mol: &Molecule, layout: &Layout) -> String {
    // Compute bounding box with padding.
    let (min_x, min_y, max_x, max_y) = layout.bounding_box();

    // For a single atom the bounding box collapses; ensure a minimum viewport.
    let raw_w = (max_x - min_x).max(BOND_LEN);
    let raw_h = (max_y - min_y).max(BOND_LEN);

    let view_x = min_x - PADDING;
    let view_y = min_y - PADDING;
    let view_w = raw_w + 2.0 * PADDING;
    let view_h = raw_h + 2.0 * PADDING;

    let width = view_w.round() as u32;
    let height = view_h.round() as u32;

    // Detect aromatic bonds using SSSR + the Aromatic bond order flag.
    // Bonds with BondOrder::Aromatic are treated as aromatic regardless of perception.
    let ring_set = find_sssr(mol);
    let _ = ring_set; // Available if needed for ring membership tests.

    let mut svg = String::new();

    // SVG header.
    svg.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" \
         width=\"{}\" height=\"{}\" \
         viewBox=\"{:.2} {:.2} {:.2} {:.2}\">\n",
        width, height, view_x, view_y, view_w, view_h
    ));

    // 1. Draw all bonds.
    for (_, bond) in mol.bonds() {
        let p1 = layout.get(bond.atom1);
        let p2 = layout.get(bond.atom2);
        let bond_svg = render_bond(bond.order, p1, p2);
        svg.push_str(&bond_svg);
    }

    // 2. Draw label backgrounds, then labels.
    for (idx, _atom) in mol.atoms() {
        let label = atom_label(mol, idx);
        if label.is_empty() {
            continue;
        }
        let p = layout.get(idx);
        // White background rect.
        svg.push_str(&format!(
            "  <rect x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" fill=\"white\"/>\n",
            p.x - LABEL_HALF_W,
            p.y - LABEL_HALF_H,
            LABEL_HALF_W * 2.0,
            LABEL_HALF_H * 2.0,
        ));
        // Label text.
        svg.push_str(&format!(
            "  <text x=\"{:.2}\" y=\"{:.2}\" \
             font-family=\"sans-serif\" font-size=\"{}\" \
             text-anchor=\"middle\" dominant-baseline=\"central\" \
             fill=\"black\">{}</text>\n",
            p.x, p.y, FONT_SIZE as u32, escape_xml(&label)
        ));
    }

    svg.push_str("</svg>");
    svg
}

// ---------------------------------------------------------------------------
// Bond rendering
// ---------------------------------------------------------------------------

/// Render a single bond as SVG elements, returned as a string fragment.
fn render_bond(order: BondOrder, p1: Point, p2: Point) -> String {
    match order {
        BondOrder::Single | BondOrder::Up | BondOrder::Down => {
            // For Up/Down we use simple lines here.
            // Full wedge/dash stereo rendering.
            match order {
                BondOrder::Up => render_wedge_up(p1, p2),
                BondOrder::Down => render_dash_bond(p1, p2),
                _ => render_single_line(p1, p2, "1.5"),
            }
        }
        BondOrder::Double => render_double_bond(p1, p2),
        BondOrder::Triple => render_triple_bond(p1, p2),
        BondOrder::Aromatic => render_aromatic_bond(p1, p2),
        BondOrder::Quadruple => render_single_line(p1, p2, "3.0"),
    }
}

/// A simple solid line between two points.
fn render_single_line(p1: Point, p2: Point, stroke_width: &str) -> String {
    format!(
        "  <line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" \
         stroke=\"black\" stroke-width=\"{}\" fill=\"none\"/>\n",
        p1.x, p1.y, p2.x, p2.y, stroke_width
    )
}

/// Perpendicular unit vector to the direction (p2-p1).
fn perp_unit(p1: Point, p2: Point) -> (f64, f64) {
    let dx = p2.x - p1.x;
    let dy = p2.y - p1.y;
    let len = (dx * dx + dy * dy).sqrt();
    if len < 1e-10 {
        return (0.0, 1.0);
    }
    (-dy / len, dx / len)
}

/// Double bond: two parallel lines, offset ±2px perpendicular to bond direction.
fn render_double_bond(p1: Point, p2: Point) -> String {
    let offset = 2.0;
    let (px, py) = perp_unit(p1, p2);
    let mut s = String::new();
    for sign in [-1.0_f64, 1.0] {
        let ox = px * offset * sign;
        let oy = py * offset * sign;
        s.push_str(&format!(
            "  <line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" \
             stroke=\"black\" stroke-width=\"1.5\" fill=\"none\"/>\n",
            p1.x + ox,
            p1.y + oy,
            p2.x + ox,
            p2.y + oy,
        ));
    }
    s
}

/// Triple bond: three parallel lines, center + offset ±3px.
fn render_triple_bond(p1: Point, p2: Point) -> String {
    let (px, py) = perp_unit(p1, p2);
    let mut s = String::new();
    for &offset in &[0.0_f64, -3.0, 3.0] {
        let ox = px * offset;
        let oy = py * offset;
        s.push_str(&format!(
            "  <line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" \
             stroke=\"black\" stroke-width=\"1.5\" fill=\"none\"/>\n",
            p1.x + ox,
            p1.y + oy,
            p2.x + ox,
            p2.y + oy,
        ));
    }
    s
}

/// Aromatic bond: one solid line + one dashed line, offset ±2px.
fn render_aromatic_bond(p1: Point, p2: Point) -> String {
    let offset = 2.0;
    let (px, py) = perp_unit(p1, p2);
    let mut s = String::new();

    // Solid line.
    s.push_str(&format!(
        "  <line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" \
         stroke=\"black\" stroke-width=\"1.5\" fill=\"none\"/>\n",
        p1.x - px * offset,
        p1.y - py * offset,
        p2.x - px * offset,
        p2.y - py * offset,
    ));
    // Dashed line.
    s.push_str(&format!(
        "  <line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" \
         stroke=\"black\" stroke-width=\"1.5\" fill=\"none\" stroke-dasharray=\"4,3\"/>\n",
        p1.x + px * offset,
        p1.y + py * offset,
        p2.x + px * offset,
        p2.y + py * offset,
    ));
    s
}

/// Up (wedge) bond: filled triangle from p1 (narrow) to p2 (wide, ±3px).
fn render_wedge_up(p1: Point, p2: Point) -> String {
    let (px, py) = perp_unit(p1, p2);
    let half_w = 3.0;
    let x1 = p1.x;
    let y1 = p1.y;
    let x2a = p2.x - px * half_w;
    let y2a = p2.y - py * half_w;
    let x2b = p2.x + px * half_w;
    let y2b = p2.y + py * half_w;
    format!(
        "  <polygon points=\"{:.2},{:.2} {:.2},{:.2} {:.2},{:.2}\" \
         fill=\"black\" stroke=\"black\" stroke-width=\"0.5\"/>\n",
        x1, y1, x2a, y2a, x2b, y2b
    )
}

/// Down (dashed) bond: series of short bars across the bond direction.
fn render_dash_bond(p1: Point, p2: Point) -> String {
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
        let hw = t * 3.0 + 0.5; // grows from narrow to wide.
        s.push_str(&format!(
            "  <line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" \
             stroke=\"black\" stroke-width=\"1.0\" fill=\"none\"/>\n",
            cx - px * hw,
            cy - py * hw,
            cx + px * hw,
            cy + py * hw,
        ));
    }
    s
}

// ---------------------------------------------------------------------------
// Atom labels
// ---------------------------------------------------------------------------

/// Compute the text label for an atom.
///
/// Returns an empty string for plain carbon atoms (no charge, no isotope).
/// Returns the element symbol (with optional H count and charge) for other atoms.
fn atom_label(mol: &Molecule, idx: AtomIdx) -> String {
    let atom = mol.atom(idx);

    let is_carbon = atom.element.atomic_number() == 6;
    let has_charge = atom.charge != 0;
    let has_isotope = atom.isotope.is_some();

    // Plain carbon with no charge or isotope: no label.
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

    #[test]
    fn test_perp_unit_horizontal() {
        // Horizontal bond: perpendicular should be (0, 1) or (0, -1).
        let (px, py) = perp_unit(Point::new(0.0, 0.0), Point::new(1.0, 0.0));
        assert!((px.abs() - 0.0).abs() < 1e-9);
        assert!((py.abs() - 1.0).abs() < 1e-9);
    }
}
