//! SVG grid depiction for multiple molecules.
//!
//! [`depict_svg_grid`] tiles any number of molecules into a single SVG image
//! with uniform cell sizes and centered molecule layout.

use chematic_core::Molecule;

use crate::layout::{BOND_LEN, Layout, Point, compute_layout};
use crate::svg::{RenderOptions, render_mol_body, render_mol_body_opts};

/// Padding around each molecule within its grid cell (SVG pixels).
const PADDING: f64 = 20.0;

/// Render multiple molecules as a grid SVG with `cols` columns.
///
/// Each molecule is centred inside a uniformly-sized cell.  The cell
/// dimensions are determined by the largest molecule bounding box so that
/// every structure fits without clipping.
///
/// # Arguments
/// * `mols` — slice of molecule references to display.
/// * `cols` — number of columns in the grid (clamped to `mols.len()`).
///
/// Returns an empty `<svg>` element when `mols` is empty or `cols` is 0.
pub fn depict_svg_grid(mols: &[&Molecule], cols: usize) -> String {
    if mols.is_empty() || cols == 0 {
        return "<svg xmlns=\"http://www.w3.org/2000/svg\" \
                width=\"0\" height=\"0\"></svg>"
            .to_string();
    }

    let cols = cols.min(mols.len());
    let rows = mols.len().div_ceil(cols);

    // Compute 2D layouts for all molecules.
    let layouts: Vec<Layout> = mols.iter().map(|m| compute_layout(m)).collect();

    // Determine uniform cell size from the largest molecule bounding box.
    let (cell_w, cell_h) = layouts.iter().fold(
        (
            BOND_LEN * 2.0 + 2.0 * PADDING,
            BOND_LEN * 2.0 + 2.0 * PADDING,
        ),
        |(cw, ch), l| {
            let (min_x, min_y, max_x, max_y) = l.bounding_box();
            let w = (max_x - min_x).max(BOND_LEN) + 2.0 * PADDING;
            let h = (max_y - min_y).max(BOND_LEN) + 2.0 * PADDING;
            (cw.max(w), ch.max(h))
        },
    );

    let total_w = cols as f64 * cell_w;
    let total_h = rows as f64 * cell_h;

    let mut svg = format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" \
         width=\"{}\" height=\"{}\" \
         viewBox=\"0 0 {:.2} {:.2}\">\n",
        total_w.round() as u32,
        total_h.round() as u32,
        total_w,
        total_h,
    );

    for (i, (mol, layout)) in mols.iter().zip(layouts.iter()).enumerate() {
        let col = (i % cols) as f64;
        let row = (i / cols) as f64;

        // Centre of this grid cell.
        let cx = col * cell_w + cell_w / 2.0;
        let cy = row * cell_h + cell_h / 2.0;

        // Centre of molecule bounding box.
        let (min_x, min_y, max_x, max_y) = layout.bounding_box();
        let mol_cx = (min_x + max_x) / 2.0;
        let mol_cy = (min_y + max_y) / 2.0;

        // Translate layout so the molecule is centred in its cell.
        let dx = cx - mol_cx;
        let dy = cy - mol_cy;
        let translated = Layout {
            coords: layout
                .coords
                .iter()
                .map(|p| Point {
                    x: p.x + dx,
                    y: p.y + dy,
                })
                .collect(),
        };

        let body = render_mol_body(mol, &translated);
        svg.push_str(&format!("  <g id=\"mol-{i}\">\n{body}  </g>\n"));
    }

    svg.push_str("</svg>");
    svg
}

/// Render multiple molecules as a grid SVG, highlighting per-cell atom sets.
///
/// # Arguments
/// * `mols` — slice of `(molecule, option<render_options>)` pairs.
///   Pass `None` for the options to use the default (no highlighting).
/// * `cols` — number of columns.
///
/// This is the primitive used by [`depict_svg_grid_highlighted_smarts`] and
/// can be called directly when atom indices are already known.
pub fn depict_svg_grid_with_opts(
    mols: &[(&Molecule, Option<&RenderOptions>)],
    cols: usize,
) -> String {
    if mols.is_empty() || cols == 0 {
        return "<svg xmlns=\"http://www.w3.org/2000/svg\" \
                width=\"0\" height=\"0\"></svg>"
            .to_string();
    }

    let default_opts = RenderOptions::default();
    let cols = cols.min(mols.len());
    let rows = mols.len().div_ceil(cols);
    let layouts: Vec<Layout> = mols.iter().map(|(m, _)| compute_layout(m)).collect();

    let (cell_w, cell_h) = layouts.iter().fold(
        (
            BOND_LEN * 2.0 + 2.0 * PADDING,
            BOND_LEN * 2.0 + 2.0 * PADDING,
        ),
        |(cw, ch), l| {
            let (min_x, min_y, max_x, max_y) = l.bounding_box();
            let w = (max_x - min_x).max(BOND_LEN) + 2.0 * PADDING;
            let h = (max_y - min_y).max(BOND_LEN) + 2.0 * PADDING;
            (cw.max(w), ch.max(h))
        },
    );

    let total_w = cols as f64 * cell_w;
    let total_h = rows as f64 * cell_h;

    let mut svg = format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" \
         width=\"{}\" height=\"{}\" \
         viewBox=\"0 0 {:.2} {:.2}\">\n",
        total_w.round() as u32,
        total_h.round() as u32,
        total_w,
        total_h,
    );

    for (i, ((mol, opts), layout)) in mols.iter().zip(layouts.iter()).enumerate() {
        let col = (i % cols) as f64;
        let row = (i / cols) as f64;
        let cx = col * cell_w + cell_w / 2.0;
        let cy = row * cell_h + cell_h / 2.0;
        let (min_x, min_y, max_x, max_y) = layout.bounding_box();
        let mol_cx = (min_x + max_x) / 2.0;
        let mol_cy = (min_y + max_y) / 2.0;
        let dx = cx - mol_cx;
        let dy = cy - mol_cy;
        let translated = Layout {
            coords: layout
                .coords
                .iter()
                .map(|p| Point {
                    x: p.x + dx,
                    y: p.y + dy,
                })
                .collect(),
        };

        let body = render_mol_body_opts(mol, &translated, opts.unwrap_or(&default_opts));
        svg.push_str(&format!("  <g id=\"mol-{i}\">\n{body}  </g>\n"));
    }

    svg.push_str("</svg>");
    svg
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> chematic_core::Molecule {
        parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"))
    }

    #[test]
    fn grid_empty_returns_empty_svg() {
        let svg = depict_svg_grid(&[], 2);
        assert!(svg.contains("<svg"), "empty grid returns SVG element");
        assert!(!svg.contains("<line"), "empty grid has no lines");
    }

    #[test]
    fn grid_zero_cols_returns_empty_svg() {
        let benzene = mol("c1ccccc1");
        let mols: Vec<&chematic_core::Molecule> = vec![&benzene];
        let svg = depict_svg_grid(&mols, 0);
        assert!(svg.contains("<svg"), "zero-col grid returns SVG element");
        assert!(!svg.contains("<line"), "zero-col grid has no lines");
    }

    #[test]
    fn grid_single_molecule_has_lines() {
        let benzene = mol("c1ccccc1");
        let mols = vec![&benzene];
        let svg = depict_svg_grid(&mols, 1);
        assert!(svg.starts_with("<svg"), "starts with <svg");
        assert!(svg.ends_with("</svg>"), "ends with </svg>");
        assert!(svg.contains("<line"), "benzene should have bond lines");
    }

    #[test]
    fn grid_two_molecules_two_cols() {
        let aspirin = mol("CC(=O)Oc1ccccc1C(=O)O");
        let caffeine = mol("Cn1cnc2c1c(=O)n(c(=O)n2C)C");
        let mols = vec![&aspirin, &caffeine];
        let svg = depict_svg_grid(&mols, 2);
        assert!(svg.contains("mol-0"), "contains mol-0 group");
        assert!(svg.contains("mol-1"), "contains mol-1 group");
        // Both molecules contain N labels (caffeine) and O labels (both)
        assert!(svg.contains('N'), "grid should have N label");
        assert!(svg.contains('O'), "grid should have O label");
    }

    #[test]
    fn grid_four_molecules_two_cols_two_rows() {
        let benzene = mol("c1ccccc1");
        let ethanol = mol("CCO");
        let pyridine = mol("c1ccncc1");
        let aspirin = mol("CC(=O)Oc1ccccc1C(=O)O");
        let mols = vec![&benzene, &ethanol, &pyridine, &aspirin];
        let svg = depict_svg_grid(&mols, 2);
        assert!(svg.contains("mol-3"), "4 molecules, last id=mol-3");
    }

    #[test]
    fn grid_cols_clamped_to_mol_count() {
        // Requesting more columns than molecules should not panic.
        let benzene = mol("c1ccccc1");
        let mols = vec![&benzene];
        let svg = depict_svg_grid(&mols, 10);
        assert!(svg.contains("<svg"), "grid with excess cols still renders");
    }
}
