//! `chematic-depict` — 2D SVG depiction engine for chematic.
//!
//! Entry point: `depict_svg(mol)` returns an SVG string.

#![forbid(unsafe_code)]

pub mod eps;
pub mod grid;
pub mod layout;
#[cfg(feature = "png")]
pub mod png;
pub mod similarity_map;
pub mod svg;

use chematic_core::{AtomIdx, BondIdx, BondOrder, Element, Molecule};

pub use eps::{render_eps, render_eps_opts};
pub use grid::{depict_svg_grid, depict_svg_grid_with_opts};
pub use layout::{
    BOND_LEN, Layout, Point, compute_layout, detect_crossings, suggest_bond_direction,
};
#[cfg(feature = "png")]
pub use png::{render_png, render_png_opts};
pub use reaction_svg::{depict_reaction_svg, depict_reaction_svg_opts};
pub use svg::{
    AtomLabel, HPosition, RenderOptions, atom_color, atom_color_rgb, atom_display_label,
    atom_label_with_h, render_svg, render_svg_highlighted, render_svg_opts,
    render_svg_with_metadata,
};

pub mod reaction_svg;

pub use similarity_map::{similarity_map_svg, similarity_map_svg_opts};

// ---------------------------------------------------------------------------
// DepictData — structured drawing data for egui/canvas renderers
// ---------------------------------------------------------------------------

/// Visual bond type for rendering.
#[derive(Debug, Clone, PartialEq)]
pub enum DepictBondKind {
    Single,
    Double,
    Triple,
    Aromatic,
    /// Wedge bond pointing up (solid triangle).
    Up,
    /// Wedge bond pointing down (dashed lines).
    Down,
}

/// A single atom's drawing data.
#[derive(Debug, Clone)]
pub struct DepictAtom {
    /// Index in the parent molecule.
    pub idx: AtomIdx,
    /// Element.
    pub element: Element,
    /// 2D position in layout units (~40 units per Å).
    pub pos: Point,
    /// Formal charge.
    pub charge: i8,
    /// Display label.  `None` = suppress (carbon in skeletal structure).
    pub label: Option<String>,
    /// CPK color as CSS hex string (e.g. `"#3050F8"` for N).
    pub color: String,
}

/// A single bond's drawing data.
#[derive(Debug, Clone)]
pub struct DepictBond {
    pub idx: BondIdx,
    pub atom1: AtomIdx,
    pub atom2: AtomIdx,
    pub kind: DepictBondKind,
}

/// Structured depiction data — use this to drive egui/canvas renderers instead
/// of parsing SVG output.
#[derive(Debug, Clone)]
pub struct DepictData {
    pub atoms: Vec<DepictAtom>,
    pub bonds: Vec<DepictBond>,
}

/// Compute structured depiction data for `mol` using auto-generated 2D layout.
pub fn compute_depict_data(mol: &Molecule) -> DepictData {
    let layout = compute_layout(mol);
    depict_data_from_layout(mol, &layout)
}

/// Compute structured depiction data for `mol` using caller-supplied coordinates.
///
/// `coords[i]` is the `(x, y)` position for atom `i` in the same units as
/// [`compute_layout`] (~40 units per Å).  Atoms beyond `coords.len()` receive
/// position `(0.0, 0.0)`.
pub fn depict_data_with_coords(mol: &Molecule, coords: &[(f64, f64)]) -> DepictData {
    let layout = Layout {
        coords: coords.iter().map(|&(x, y)| Point { x, y }).collect(),
    };
    depict_data_from_layout(mol, &layout)
}

/// Internal helper: build `DepictData` from a pre-computed `Layout`.
fn depict_data_from_layout(mol: &Molecule, layout: &Layout) -> DepictData {
    let atoms: Vec<DepictAtom> = mol
        .atoms()
        .map(|(idx, atom)| {
            let pos = layout.get(idx);
            let color = atom_color(atom.element.atomic_number()).to_string();
            let label = if atom.element.atomic_number() == 6
                && atom.charge == 0
                && atom.isotope.is_none()
                && mol.degree(idx) > 0
            {
                None
            } else {
                Some(atom.element.symbol().to_string())
            };
            DepictAtom {
                idx,
                element: atom.element,
                pos,
                charge: atom.charge,
                label,
                color,
            }
        })
        .collect();

    let bonds: Vec<DepictBond> = mol
        .bonds()
        .map(|(bidx, bond)| {
            let kind = match bond.order {
                BondOrder::Single => DepictBondKind::Single,
                BondOrder::Double => DepictBondKind::Double,
                BondOrder::Triple => DepictBondKind::Triple,
                BondOrder::Aromatic => DepictBondKind::Aromatic,
                BondOrder::Up => DepictBondKind::Up,
                BondOrder::Down => DepictBondKind::Down,
                BondOrder::Quadruple => DepictBondKind::Triple,
                BondOrder::Zero
                | BondOrder::Dative
                | BondOrder::QueryAny
                | BondOrder::QuerySingleOrDouble
                | BondOrder::QuerySingleOrAromatic
                | BondOrder::QueryDoubleOrAromatic => DepictBondKind::Single,
            };
            DepictBond {
                idx: bidx,
                atom1: bond.atom1,
                atom2: bond.atom2,
                kind,
            }
        })
        .collect();

    DepictData { atoms, bonds }
}

// ---------------------------------------------------------------------------
// PDF output (optional feature)
// ---------------------------------------------------------------------------

/// Render a molecule as a PDF document (bytes).
///
/// Requires the `pdf` Cargo feature.
/// Equivalent to `depict_svg` but converts the SVG to a PDF byte vector.
#[cfg(feature = "pdf")]
pub fn depict_pdf(mol: &Molecule) -> Vec<u8> {
    let svg = depict_svg(mol);
    svg_bytes_to_pdf(&svg)
}

/// Render a molecule as a PDF document (bytes) with full style control.
///
/// Requires the `pdf` Cargo feature.
#[cfg(feature = "pdf")]
pub fn depict_pdf_opts(mol: &Molecule, opts: &RenderOptions) -> Vec<u8> {
    let svg = depict_svg_opts(mol, opts);
    svg_bytes_to_pdf(&svg)
}

#[cfg(feature = "pdf")]
fn svg_bytes_to_pdf(svg: &str) -> Vec<u8> {
    let mut options = svg2pdf::usvg::Options::default();
    options.fontdb_mut().load_system_fonts();
    let tree = svg2pdf::usvg::Tree::from_str(svg, &options).expect("generated SVG should parse");
    svg2pdf::to_pdf(
        &tree,
        svg2pdf::ConversionOptions::default(),
        svg2pdf::PageOptions::default(),
    )
    .expect("PDF conversion should not fail for valid SVG")
}

/// Compute a 2D layout and render it as an EPS string.
///
/// Generates a self-contained PostScript EPS document from the molecule's
/// 2D structure without any external dependencies.
/// Suitable for publication-quality vector graphics in LaTeX, Adobe Illustrator, etc.
pub fn depict_eps(mol: &Molecule) -> String {
    let layout = compute_layout(mol);
    render_eps(mol, &layout)
}

/// Compute a 2D layout and render it as EPS with full style control.
pub fn depict_eps_opts(mol: &Molecule, opts: &RenderOptions) -> String {
    let layout = compute_layout(mol);
    render_eps_opts(mol, &layout, opts)
}

/// Compute a 2D layout and render it as an SVG string.
pub fn depict_svg(mol: &Molecule) -> String {
    let layout = compute_layout(mol);
    render_svg(mol, &layout)
}

/// Compute a 2D layout and render it as an SVG with full style control.
pub fn depict_svg_opts(mol: &Molecule, opts: &RenderOptions) -> String {
    let layout = compute_layout(mol);
    render_svg_opts(mol, &layout, opts)
}

/// Compute a 2D layout and render it as an SVG with highlighted atoms/bonds.
pub fn depict_svg_highlighted(
    mol: &Molecule,
    highlight_atoms: &std::collections::HashSet<AtomIdx>,
    highlight_bonds: &std::collections::HashSet<BondIdx>,
) -> String {
    let layout = compute_layout(mol);
    render_svg_highlighted(mol, &layout, highlight_atoms, highlight_bonds)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::layout::BOND_LEN;

    use chematic_smiles::parse;

    // Helper: parse SMILES, panic with a helpful message on failure.
    fn mol(smiles: &str) -> Molecule {
        parse(smiles).unwrap_or_else(|e| panic!("Failed to parse '{}': {:?}", smiles, e))
    }

    // -------------------------------------------------------------------
    // 1. compute_layout — benzene: 6 atoms, all coords distinct,
    //    no two atoms closer than BOND_LEN/2.
    // -------------------------------------------------------------------
    #[test]
    fn test_layout_benzene_six_distinct_coords() {
        let m = mol("c1ccccc1");
        assert_eq!(m.atom_count(), 6);
        let layout = compute_layout(&m);
        assert_eq!(layout.coords.len(), 6);

        // All coordinates must be distinct.
        for i in 0..6 {
            for j in (i + 1)..6 {
                let d = layout.coords[i].dist(&layout.coords[j]);
                assert!(
                    d > BOND_LEN / 2.0,
                    "Atoms {} and {} are too close: {:.2} < {:.2}",
                    i,
                    j,
                    d,
                    BOND_LEN / 2.0
                );
            }
        }
    }

    // -------------------------------------------------------------------
    // 2. compute_layout — single atom: one coord near origin.
    // -------------------------------------------------------------------
    #[test]
    fn test_layout_single_atom() {
        let m = mol("[C]");
        assert_eq!(m.atom_count(), 1);
        let layout = compute_layout(&m);
        assert_eq!(layout.coords.len(), 1);
        // Single atom should be placed at (0, 0).
        let p = layout.coords[0];
        assert!(
            (p.x).abs() < 1.0 && (p.y).abs() < 1.0,
            "Single atom not near origin: {:?}",
            p
        );
    }

    // -------------------------------------------------------------------
    // 3. compute_layout — ethane (CC): 2 atoms, distance ~= BOND_LEN (±1%).
    // -------------------------------------------------------------------
    #[test]
    fn test_layout_ethane_bond_length() {
        let m = mol("CC");
        assert_eq!(m.atom_count(), 2);
        let layout = compute_layout(&m);
        assert_eq!(layout.coords.len(), 2);
        let d = layout.coords[0].dist(&layout.coords[1]);
        let tolerance = BOND_LEN * 0.01;
        assert!(
            (d - BOND_LEN).abs() < tolerance,
            "Ethane bond distance {:.4} != BOND_LEN {:.4} (±{:.4})",
            d,
            BOND_LEN,
            tolerance
        );
    }

    // -------------------------------------------------------------------
    // 4. compute_layout — naphthalene: 10 distinct coords, reasonable bbox.
    // -------------------------------------------------------------------
    #[test]
    fn test_layout_naphthalene_ten_coords() {
        let m = mol("c1ccc2ccccc2c1");
        assert_eq!(m.atom_count(), 10);
        let layout = compute_layout(&m);
        assert_eq!(layout.coords.len(), 10);

        // All coords distinct.
        for i in 0..10 {
            for j in (i + 1)..10 {
                let d = layout.coords[i].dist(&layout.coords[j]);
                assert!(
                    d > BOND_LEN / 2.0,
                    "Naphthalene atoms {} and {} too close: {:.2}",
                    i,
                    j,
                    d
                );
            }
        }

        // Reasonable bounding box: no wider/taller than 6 * BOND_LEN.
        let (min_x, min_y, max_x, max_y) = layout.bounding_box();
        assert!(max_x - min_x < 6.0 * BOND_LEN, "Naphthalene too wide");
        assert!(max_y - min_y < 6.0 * BOND_LEN, "Naphthalene too tall");
    }

    // -------------------------------------------------------------------
    // 5. compute_layout — disconnected mol ("CC.CC"): atoms from different
    //    fragments are farther than BOND_LEN apart.
    // -------------------------------------------------------------------
    #[test]
    fn test_layout_disconnected_fragments_no_overlap() {
        let m = mol("CC.CC");
        assert_eq!(m.atom_count(), 4);
        let layout = compute_layout(&m);

        // Fragment 0 = atoms 0,1; fragment 1 = atoms 2,3.
        // No atom from fragment 0 should be within BOND_LEN of any atom in fragment 1.
        for i in 0..2 {
            for j in 2..4 {
                let d = layout.coords[i].dist(&layout.coords[j]);
                assert!(
                    d >= BOND_LEN,
                    "Atoms from different fragments too close: atoms {} and {}, dist {:.2}",
                    i,
                    j,
                    d
                );
            }
        }
    }

    // -------------------------------------------------------------------
    // 6. render_svg — benzene: SVG contains <svg and <line but no <text.
    // -------------------------------------------------------------------
    #[test]
    fn test_svg_benzene_no_text() {
        let m = mol("c1ccccc1");
        let layout = compute_layout(&m);
        let svg = render_svg(&m, &layout);
        assert!(svg.contains("<svg"), "SVG must start with <svg");
        assert!(svg.contains("<line"), "SVG must have bond lines");
        assert!(
            !svg.contains("<text"),
            "Benzene SVG must have no atom text labels"
        );
    }

    // -------------------------------------------------------------------
    // 7. render_svg — pyridine: SVG contains <text with "N".
    // -------------------------------------------------------------------
    #[test]
    fn test_svg_pyridine_contains_nitrogen_label() {
        let m = mol("c1ccncc1");
        let layout = compute_layout(&m);
        let svg = render_svg(&m, &layout);
        assert!(svg.contains("<text"), "Pyridine SVG must have a text label");
        assert!(svg.contains('N'), "Pyridine SVG must contain 'N' label");
        // N should be colored blue (#3050F8).
        assert!(
            svg.contains("#3050F8"),
            "Pyridine N label should be blue (#3050F8)"
        );
    }

    // -------------------------------------------------------------------
    // 8. render_svg — aspirin: non-empty, contains <line.
    // -------------------------------------------------------------------
    #[test]
    fn test_svg_aspirin_non_empty() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        let layout = compute_layout(&m);
        let svg = render_svg(&m, &layout);
        assert!(!svg.is_empty(), "Aspirin SVG must be non-empty");
        assert!(
            svg.contains("<line"),
            "Aspirin SVG must contain line elements"
        );
    }

    // -------------------------------------------------------------------
    // 9. render_svg — double bond (C=C): SVG contains two <line elements.
    // -------------------------------------------------------------------
    #[test]
    fn test_svg_double_bond_two_lines() {
        let m = mol("C=C");
        let layout = compute_layout(&m);
        let svg = render_svg(&m, &layout);
        let count = svg.matches("<line").count();
        assert!(
            count >= 2,
            "C=C SVG should have >= 2 <line elements, got {}",
            count
        );
    }

    // -------------------------------------------------------------------
    // 10. depict_svg — caffeine: produces valid SVG.
    // -------------------------------------------------------------------
    #[test]
    fn test_depict_svg_caffeine_valid() {
        let m = mol("Cn1cnc2c1c(=O)n(c(=O)n2C)C");
        let svg = depict_svg(&m);
        assert!(svg.starts_with("<svg"), "Caffeine SVG must start with <svg");
        assert!(svg.ends_with("</svg>"), "Caffeine SVG must end with </svg>");
    }

    // -------------------------------------------------------------------
    // 11. depict_svg — water ([OH2]): SVG contains "O" label.
    // -------------------------------------------------------------------
    #[test]
    fn test_depict_svg_water_contains_o() {
        let m = mol("[OH2]");
        let svg = depict_svg(&m);
        assert!(svg.contains('O'), "Water SVG must contain 'O'");
    }

    // -------------------------------------------------------------------
    // 12. depict_svg — single carbon (C): shows molecular formula "CH4".
    // -------------------------------------------------------------------
    #[test]
    fn test_depict_svg_single_carbon_shows_ch4() {
        let m = mol("C");
        let svg = depict_svg(&m);
        assert!(svg.starts_with("<svg"), "Single C SVG must start with <svg");
        assert!(svg.ends_with("</svg>"), "Single C SVG must end with </svg>");
        // Isolated carbon shows molecular formula label.
        assert!(svg.contains("CH4"), "Single C SVG must show CH4 label");
    }

    // -------------------------------------------------------------------
    // 13. depict_svg — single oxygen (O): shows "H2O" label.
    // -------------------------------------------------------------------
    #[test]
    fn test_depict_svg_single_oxygen_shows_h2o() {
        let m = mol("O");
        let svg = depict_svg(&m);
        assert!(svg.contains("H2O"), "Single O SVG must show H2O label");
    }

    // -------------------------------------------------------------------
    // render_svg_highlighted — pyridine with N highlighted.
    // -------------------------------------------------------------------
    #[test]
    fn test_svg_highlighted_pyridine() {
        use std::collections::HashSet;
        let m = mol("c1ccncc1");
        // Find the N atom index.
        let n_idx = m
            .atoms()
            .find(|(_, a)| a.element.atomic_number() == 7)
            .map(|(idx, _)| idx)
            .expect("pyridine has no N");
        let mut hl_atoms = HashSet::new();
        hl_atoms.insert(n_idx);
        let svg = depict_svg_highlighted(&m, &hl_atoms, &HashSet::new());
        assert!(
            svg.contains("circle"),
            "highlighted SVG must contain a circle"
        );
        assert!(svg.contains("FFFF00"), "highlight circle should be yellow");
    }
}
