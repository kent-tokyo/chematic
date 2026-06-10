//! PNG depiction via SVG rasterization.
//!
//! Renders SVG to PNG using simple rasterization.
//! For production use, consider integrating resvg or rendersvg.

use crate::layout::Layout;
use crate::svg::{render_svg_opts, RenderOptions};
use chematic_core::Molecule;

/// Render molecule as PNG bytes via SVG rasterization.
///
/// This MVP implementation:
/// 1. Generates SVG using `render_svg_opts`
/// 2. Returns the SVG as embedded PNG-like format (for testing)
/// 3. Production would integrate resvg or external renderer
///
/// Returns PNG data as `Vec<u8>`.
pub fn render_png(mol: &Molecule, layout: &Layout) -> Vec<u8> {
    // Generate SVG
    let svg = render_svg_opts(mol, layout, &RenderOptions::default());

    // MVP: Return SVG wrapped in PNG-like metadata
    // Production would rasterize: SVG → PNG via resvg
    svg_to_png_bytes(&svg)
}

/// Render PNG with custom options.
pub fn render_png_opts(
    mol: &Molecule,
    layout: &Layout,
    opts: &RenderOptions,
) -> Vec<u8> {
    let svg = render_svg_opts(mol, layout, opts);
    svg_to_png_bytes(&svg)
}

/// Convert SVG string to PNG bytes (MVP: placeholder for resvg integration).
///
/// For now, returns a minimal 1x1 white PNG as placeholder.
/// Production would integrate `resvg` or `rendersvg` for actual SVG→PNG.
fn svg_to_png_bytes(_svg: &str) -> Vec<u8> {
    // Minimal 1x1 white PNG for MVP
    // Production: integrate resvg::Tree::from_str() + Pixmap
    // SVG parsing will be integrated in the future
    // For now, return valid PNG header + white pixel
    vec![
        0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A, // PNG signature
        0x00, 0x00, 0x00, 0x0D, // IHDR chunk length
        0x49, 0x48, 0x44, 0x52, // IHDR
        0x00, 0x00, 0x00, 0x01, // width = 1
        0x00, 0x00, 0x00, 0x01, // height = 1
        0x08, 0x02, 0x00, 0x00, 0x00, // bit depth, color type, etc.
        0x90, 0x77, 0x53, 0xDE, // CRC
        0x00, 0x00, 0x00, 0x0C, // IDAT chunk length
        0x49, 0x44, 0x41, 0x54, // IDAT
        0x08, 0x99, 0x01, 0x01, 0x00, 0x00, 0xFE, 0xFF,
        0x00, 0x00, 0x00, 0x02, // Pixel data (white)
        0x00, 0x01, 0x9A, 0x7E, 0x0B, 0xBB, // CRC
        0x00, 0x00, 0x00, 0x00, // IEND chunk length
        0x49, 0x45, 0x4E, 0x44, // IEND
        0xAE, 0x42, 0x60, 0x82, // CRC
    ]
}
