//! Keep a degree-two carbon visible when its bonds form a straight axis.

use chematic_core::{AtomIdx, BondOrder, Molecule};

use crate::layout::{Layout, Point};

/// An sp carbon is linear even when a caller supplies its own layout.
pub(super) fn has_sp_axis(mol: &Molecule, idx: AtomIdx) -> bool {
    if mol.degree(idx) != 2 {
        return false;
    }
    let mut orders = mol.neighbors(idx).map(|(_, bond)| mol.bond(bond).order);
    let (Some(first), Some(second)) = (orders.next(), orders.next()) else {
        return false;
    };
    matches!(first, BondOrder::Triple)
        || matches!(second, BondOrder::Triple)
        || (first == BondOrder::Double && second == BondOrder::Double)
}

pub(super) fn needs_label(mol: &Molecule, layout: &Layout, idx: AtomIdx) -> bool {
    let atom = mol.atom(idx);
    if atom.element.atomic_number() != 6
        || atom.charge != 0
        || atom.isotope.is_some()
        || mol.degree(idx) != 2
    {
        return false;
    }
    if has_sp_axis(mol, idx) {
        return true;
    }

    let center = layout.get(idx);
    let mut neighbors = mol.neighbors(idx).map(|(other, _)| layout.get(other));
    let (Some(first), Some(second)) = (neighbors.next(), neighbors.next()) else {
        return false;
    };
    let (ux, uy) = (first.x - center.x, first.y - center.y);
    let (vx, vy) = (second.x - center.x, second.y - center.y);
    let (dot, cross) = (ux * vx + uy * vy, ux * vy - uy * vx);
    let length_product = (ux * ux + uy * uy) * (vx * vx + vy * vy);
    dot < 0.0 && length_product > 0.0 && cross * cross <= 0.0025 * length_product
}

/// Leave room for a visible C, including when the SVG background is transparent.
pub(crate) fn bond_endpoints(
    mol: &Molecule,
    layout: &Layout,
    a: AtomIdx,
    b: AtomIdx,
) -> (Point, Point) {
    let mut start = layout.get(a);
    let mut end = layout.get(b);
    let clip_start = needs_label(mol, layout, a);
    let clip_end = needs_label(mol, layout, b);
    if !clip_start && !clip_end {
        return (start, end);
    }
    let (dx, dy) = (end.x - start.x, end.y - start.y);
    let length = dx.hypot(dy);
    if length == 0.0 {
        return (start, end);
    }
    let fraction = ((super::LABEL_HALF_W + 1.0) / length).min(0.25);
    if clip_start {
        start = Point::new(start.x + dx * fraction, start.y + dy * fraction);
    }
    if clip_end {
        end = Point::new(end.x - dx * fraction, end.y - dy * fraction);
    }
    (start, end)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn co2_is_labeled_and_bonds_stop_before_the_text() {
        let mol = parse("O=C=O").unwrap();
        let layout = crate::compute_layout(&mol);
        let svg = crate::render_svg(&mol, &layout);
        assert_eq!(svg.matches(">C</text>").count(), 1);
        assert_eq!(svg.matches("<line ").count(), 4);
        assert!(svg.contains("x1=\"0.00\" y1=\"-2.00\" x2=\"31.00\""));
        assert!(svg.contains("x1=\"49.00\" y1=\"-2.00\" x2=\"80.00\""));
        assert!(svg.contains("x=\"32.00\" y=\"-7.00\" width=\"16.00\""));

        let transparent = crate::RenderOptions {
            background: "transparent".into(),
            ..Default::default()
        };
        let svg = crate::render_svg_opts(&mol, &layout, &transparent);
        assert_eq!(svg.matches(">C</text>").count(), 1);
        assert!(svg.contains("x2=\"31.00\""));
        assert!(!svg.contains("x=\"32.00\" y=\"-7.00\" width=\"16.00\""));
    }

    #[test]
    fn shared_render_paths_and_structured_label_show_linear_carbons() {
        let mol = parse("O=C=O").unwrap();
        let layout = crate::compute_layout(&mol);
        let options = crate::RenderOptions::default();
        for svg in [
            crate::render_svg_opts(&mol, &layout, &options),
            crate::render_svg_with_metadata(&mol, &layout, &options, "O=C=O"),
            crate::depict_svg_grid(&[&mol], 1),
            crate::depict_svg_grid_with_opts(&[(&mol, Some(&options))], 1),
        ] {
            assert_eq!(svg.matches(">C</text>").count(), 1);
            assert_eq!(svg.matches("<line ").count(), 4);
        }
        assert_eq!(crate::atom_display_label(&mol, AtomIdx(1)), "C");
    }

    #[test]
    fn cumulene_and_nitrile_centers_are_visible_but_bent_chains_stay_skeletal() {
        for (smiles, expected) in [("C=C=C", 1), ("CC#N", 1), ("CC#CC", 2)] {
            let mol = parse(smiles).unwrap();
            let svg = crate::depict_svg(&mol);
            assert_eq!(svg.matches(">C</text>").count(), expected, "{smiles}");
        }
        for smiles in ["CC", "CCC", "CCO"] {
            let mol = parse(smiles).unwrap();
            assert!(!crate::depict_svg(&mol).contains(">C</text>"), "{smiles}");
        }
        let propane = parse("CCC").unwrap();
        assert_eq!(crate::atom_display_label(&propane, AtomIdx(1)), "");

        let straight = Layout {
            coords: vec![
                Point::new(0.0, 0.0),
                Point::new(40.0, 0.0),
                Point::new(80.0, 0.0),
            ],
        };
        assert_eq!(
            crate::render_svg(&propane, &straight)
                .matches(">C</text>")
                .count(),
            1
        );
    }
}
