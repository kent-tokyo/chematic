//! PNG depiction via tiny-skia rasterization.
//!
//! Renders molecules directly to PNG using tiny-skia 2D graphics.
//! Draws bonds and atom positions from Layout data.

use crate::layout::Layout;
use crate::svg::RenderOptions;
use chematic_core::Molecule;
use tiny_skia::{Color, Paint, Pixmap, Stroke};

const PIXELS_PER_UNIT: f64 = 10.0;
const BOND_WIDTH: f32 = 1.5;
const ATOM_RADIUS: f32 = 3.0;

/// Render molecule as PNG bytes using tiny-skia.
pub fn render_png(mol: &Molecule, layout: &Layout) -> Vec<u8> {
    render_png_opts(mol, layout, &RenderOptions::default())
}

/// Render PNG with custom options.
pub fn render_png_opts(
    mol: &Molecule,
    layout: &Layout,
    _opts: &RenderOptions,
) -> Vec<u8> {
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

    // Draw atoms as small circles (using fill_rect as approximation)
    for (idx, _) in mol.atoms() {
        let p = layout.get(idx);
        let x = ((p.x + offset_x) * PIXELS_PER_UNIT) as f32;
        let y = ((p.y + offset_y) * PIXELS_PER_UNIT) as f32;
        let r = ATOM_RADIUS;

        // Draw filled rectangle as atom marker
        if let Some(rect) = tiny_skia::Rect::from_xywh(x - r, y - r, 2.0 * r, 2.0 * r) {
            let mut paint = Paint::default();
            paint.set_color(Color::from_rgba8(150, 150, 150, 255));
            pixmap.fill_rect(rect, &paint, tiny_skia::Transform::default(), None);
        }
    }

    // Draw bonds as lines between atoms
    let stroke = Stroke {
        width: BOND_WIDTH,
        ..Default::default()
    };
    for (_, bond) in mol.bonds() {
        let p1 = layout.get(bond.atom1);
        let p2 = layout.get(bond.atom2);
        let x1 = ((p1.x + offset_x) * PIXELS_PER_UNIT) as f32;
        let y1 = ((p1.y + offset_y) * PIXELS_PER_UNIT) as f32;
        let x2 = ((p2.x + offset_x) * PIXELS_PER_UNIT) as f32;
        let y2 = ((p2.y + offset_y) * PIXELS_PER_UNIT) as f32;

        let mut pb = tiny_skia::PathBuilder::new();
        pb.move_to(x1, y1);
        pb.line_to(x2, y2);
        if let Some(path) = pb.finish() {
            pixmap.stroke_path(&path, &paint_bond, &stroke, tiny_skia::Transform::default(), None);
        }
    }

    pixmap.encode_png().unwrap_or_else(|_| empty_png())
}

fn empty_png() -> Vec<u8> {
    vec![
        0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A,
        0x00, 0x00, 0x00, 0x0D, 0x49, 0x48, 0x44, 0x52,
        0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01,
        0x08, 0x02, 0x00, 0x00, 0x00, 0x90, 0x77, 0x53,
        0xDE, 0x00, 0x00, 0x00, 0x0C, 0x49, 0x44, 0x41,
        0x54, 0x08, 0x99, 0x01, 0x01, 0x00, 0x00, 0xFE,
        0xFF, 0x00, 0x00, 0x00, 0x02, 0x00, 0x01, 0x9A,
        0x7E, 0x0B, 0xBB, 0x00, 0x00, 0x00, 0x00, 0x49,
        0x45, 0x4E, 0x44, 0xAE, 0x42, 0x60, 0x82,
    ]
}
