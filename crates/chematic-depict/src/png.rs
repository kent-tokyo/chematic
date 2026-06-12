//! PNG depiction via tiny-skia rasterization.
//!
//! Renders molecules directly to PNG using tiny-skia 2D graphics.
//! Draws bonds, atom positions, and element labels from Layout data.

use crate::layout::Layout;
use crate::svg::RenderOptions;
use chematic_core::Molecule;
use tiny_skia::{Color, Paint, Pixmap, Stroke};

const PIXELS_PER_UNIT: f64 = 10.0;
const BOND_WIDTH: f32 = 1.5;
const ATOM_RADIUS: f32 = 3.0;

/// CPK color table: atomic number → RGB tuple
fn cpk_color(atomic_number: u8) -> (u8, u8, u8) {
    match atomic_number {
        1 => (255, 255, 255), // H: white
        6 => (80, 80, 80),    // C: dark grey
        7 => (48, 80, 248),   // N: blue
        8 => (255, 13, 13),   // O: red
        9 => (144, 224, 80),  // F: green
        15 => (255, 128, 0),  // P: orange
        16 => (255, 200, 50), // S: yellow
        17 => (31, 240, 31),  // Cl: light green
        35 => (166, 41, 41),  // Br: dark red
        53 => (148, 0, 148),  // I: purple
        _ => (255, 20, 147),  // Other: pink
    }
}

/// Bitmap patterns for element labels (5×4 pixels each).
fn get_char_bitmap(ch: char) -> Option<[u8; 5]> {
    match ch.to_ascii_uppercase() {
        'C' => Some([0b0110, 0b1001, 0b1000, 0b1001, 0b0110]),
        'N' => Some([0b1001, 0b1011, 0b1101, 0b1001, 0b1001]),
        'O' => Some([0b0110, 0b1001, 0b1001, 0b1001, 0b0110]),
        'H' => Some([0b1001, 0b1001, 0b1111, 0b1001, 0b1001]),
        'F' => Some([0b1111, 0b1000, 0b1110, 0b1000, 0b1000]),
        'P' => Some([0b1110, 0b1001, 0b1110, 0b1000, 0b1000]),
        'S' => Some([0b0110, 0b1000, 0b0110, 0b0001, 0b1110]),
        'B' => Some([0b1110, 0b1001, 0b1110, 0b1001, 0b1110]),
        'I' => Some([0b1111, 0b0100, 0b0100, 0b0100, 0b1111]),
        'K' => Some([0b1001, 0b1010, 0b1100, 0b1010, 0b1001]),
        'X' => Some([0b1001, 0b0110, 0b0110, 0b0110, 0b1001]),
        _ => None,
    }
}

/// Render molecule as PNG bytes using tiny-skia.
pub fn render_png(mol: &Molecule, layout: &Layout) -> Vec<u8> {
    render_png_opts(mol, layout, &RenderOptions::default())
}

/// Render PNG with custom options.
pub fn render_png_opts(mol: &Molecule, layout: &Layout, _opts: &RenderOptions) -> Vec<u8> {
    if mol.atom_count() == 0 {
        return empty_png();
    }

    let bounds = layout.bounding_box();
    let width_units = bounds.2 - bounds.0 + 2.0;
    let height_units = bounds.3 - bounds.1 + 2.0;
    let width = (width_units * PIXELS_PER_UNIT).max(100.0) as u32;
    let height = (height_units * PIXELS_PER_UNIT).max(100.0) as u32;

    let mut pixmap = match Pixmap::new(width, height) {
        Some(pm) => pm,
        None => return empty_png(),
    };

    pixmap.fill(Color::WHITE);
    let offset_x = -bounds.0 + 1.0;
    let offset_y = -bounds.1 + 1.0;

    // Draw bonds (simplified: use pixmap as-is, small dots for atoms)
    let mut paint_bond = Paint::default();
    paint_bond.set_color(Color::BLACK);

    // Draw atoms as markers with CPK colors and element labels
    for (idx, _) in mol.atoms() {
        let p = layout.get(idx);
        let x = ((p.x + offset_x) * PIXELS_PER_UNIT) as f32;
        let y = ((p.y + offset_y) * PIXELS_PER_UNIT) as f32;
        let r = ATOM_RADIUS;

        // Draw filled rectangle with CPK color
        if let Some(rect) = tiny_skia::Rect::from_xywh(x - r, y - r, 2.0 * r, 2.0 * r) {
            let atom = mol.atom(idx);
            let (r_val, g_val, b_val) = cpk_color(atom.element.atomic_number());
            let mut paint = Paint::default();
            paint.set_color(Color::from_rgba8(r_val, g_val, b_val, 255));
            pixmap.fill_rect(rect, &paint, tiny_skia::Transform::default(), None);
        }

        // Draw element symbol as text label
        let atom = mol.atom(idx);
        let symbol = atom.element.symbol();
        for (i, ch) in symbol.chars().enumerate() {
            if let Some(bitmap) = get_char_bitmap(ch) {
                let label_x = x - 2.0 + (i as f32 * 3.0);
                let label_y = y - 2.5;

                let mut text_paint = Paint::default();
                text_paint.set_color(Color::BLACK);

                // Draw 5 rows × 4 columns
                for (row, bits) in bitmap.iter().enumerate() {
                    for col in 0..4 {
                        if (bits >> (3 - col)) & 1 == 1 {
                            let px = label_x + (col as f32 * 0.5);
                            let py = label_y + (row as f32 * 0.5);
                            if let Some(rect) = tiny_skia::Rect::from_xywh(px, py, 0.5, 0.5) {
                                pixmap.fill_rect(
                                    rect,
                                    &text_paint,
                                    tiny_skia::Transform::default(),
                                    None,
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    // Draw bonds as lines between atoms (with wedge/dash styling)
    use chematic_core::BondOrder;

    for (_, bond) in mol.bonds() {
        let p1 = layout.get(bond.atom1);
        let p2 = layout.get(bond.atom2);
        let x1 = ((p1.x + offset_x) * PIXELS_PER_UNIT) as f32;
        let y1 = ((p1.y + offset_y) * PIXELS_PER_UNIT) as f32;
        let x2 = ((p2.x + offset_x) * PIXELS_PER_UNIT) as f32;
        let y2 = ((p2.y + offset_y) * PIXELS_PER_UNIT) as f32;

        // Determine stroke width and style based on bond order
        let stroke_width = match bond.order {
            BondOrder::Up => BOND_WIDTH * 2.5,   // Wedge: thicker
            BondOrder::Down => BOND_WIDTH * 0.8, // Dash: thinner
            _ => BOND_WIDTH,                     // Normal: regular
        };

        let stroke = Stroke {
            width: stroke_width,
            ..Default::default()
        };

        let mut pb = tiny_skia::PathBuilder::new();
        pb.move_to(x1, y1);
        pb.line_to(x2, y2);
        if let Some(path) = pb.finish() {
            pixmap.stroke_path(
                &path,
                &paint_bond,
                &stroke,
                tiny_skia::Transform::default(),
                None,
            );
        }

        // Draw dashed lines for down bonds (approximation)
        if bond.order == BondOrder::Down {
            let dx = x2 - x1;
            let dy = y2 - y1;
            let len = (dx * dx + dy * dy).sqrt();
            if len > 0.0 {
                let ux = dx / len;
                let uy = dy / len;
                let dash_len = 2.0;
                let gap_len = 2.0;

                let mut pos = 0.0;
                while pos < len {
                    let start_pos = (pos).min(len);
                    let end_pos = (pos + dash_len).min(len);

                    let sx = x1 + ux * start_pos;
                    let sy = y1 + uy * start_pos;
                    let ex = x1 + ux * end_pos;
                    let ey = y1 + uy * end_pos;

                    let mut pb = tiny_skia::PathBuilder::new();
                    pb.move_to(sx, sy);
                    pb.line_to(ex, ey);
                    if let Some(path) = pb.finish() {
                        pixmap.stroke_path(
                            &path,
                            &paint_bond,
                            &stroke,
                            tiny_skia::Transform::default(),
                            None,
                        );
                    }

                    pos += dash_len + gap_len;
                }
            }
        }
    }

    pixmap.encode_png().unwrap_or_else(|_| empty_png())
}

fn empty_png() -> Vec<u8> {
    vec![
        0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A, 0x00, 0x00, 0x00, 0x0D, 0x49, 0x48, 0x44,
        0x52, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x08, 0x02, 0x00, 0x00, 0x00, 0x90,
        0x77, 0x53, 0xDE, 0x00, 0x00, 0x00, 0x0C, 0x49, 0x44, 0x41, 0x54, 0x08, 0x99, 0x01, 0x01,
        0x00, 0x00, 0xFE, 0xFF, 0x00, 0x00, 0x00, 0x02, 0x00, 0x01, 0x9A, 0x7E, 0x0B, 0xBB, 0x00,
        0x00, 0x00, 0x00, 0x49, 0x45, 0x4E, 0x44, 0xAE, 0x42, 0x60, 0x82,
    ]
}
