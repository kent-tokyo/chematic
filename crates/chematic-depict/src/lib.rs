//! `chematic-depict` — 2D SVG depiction engine for chematic.
//!
//! Entry point: `depict_svg(mol)` returns an SVG string.

#![forbid(unsafe_code)]

pub mod layout;
pub mod svg;

use chematic_core::{AtomIdx, BondIdx, Molecule};

pub use layout::{Layout, Point, compute_layout};
pub use svg::{render_svg, render_svg_highlighted};

/// Compute a 2D layout and render it as an SVG string.
pub fn depict_svg(mol: &Molecule) -> String {
    let layout = compute_layout(mol);
    render_svg(mol, &layout)
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
        assert!((p.x).abs() < 1.0 && (p.y).abs() < 1.0, "Single atom not near origin: {:?}", p);
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
        assert!(!svg.contains("<text"), "Benzene SVG must have no atom text labels");
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
        assert!(svg.contains("#3050F8"), "Pyridine N label should be blue (#3050F8)");
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
        assert!(svg.contains("<line"), "Aspirin SVG must contain line elements");
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
        assert!(count >= 2, "C=C SVG should have >= 2 <line elements, got {}", count);
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
    // 12. depict_svg — single carbon (C): SVG produced without error,
    //     no label for plain C.
    // -------------------------------------------------------------------
    #[test]
    fn test_depict_svg_single_carbon_no_label() {
        let m = mol("C");
        let svg = depict_svg(&m);
        assert!(svg.starts_with("<svg"), "Single C SVG must start with <svg");
        assert!(svg.ends_with("</svg>"), "Single C SVG must end with </svg>");
        // Plain carbon has no text label.
        assert!(!svg.contains("<text"), "Single C SVG should have no text label");
    }

    // -------------------------------------------------------------------
    // render_svg_highlighted — pyridine with N highlighted.
    // -------------------------------------------------------------------
    #[test]
    fn test_svg_highlighted_pyridine() {
        use std::collections::HashSet;
        let m = mol("c1ccncc1");
        // Find the N atom index.
        let n_idx = m.atoms()
            .find(|(_, a)| a.element.atomic_number() == 7)
            .map(|(idx, _)| idx)
            .expect("pyridine has no N");
        let mut hl_atoms = HashSet::new();
        hl_atoms.insert(n_idx);
        let svg = depict_svg_highlighted(&m, &hl_atoms, &HashSet::new());
        assert!(svg.contains("circle"), "highlighted SVG must contain a circle");
        assert!(svg.contains("FFFF00"), "highlight circle should be yellow");
    }
}
