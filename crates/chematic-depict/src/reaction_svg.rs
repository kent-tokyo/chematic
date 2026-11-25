//! SVG rendering of chemical reactions.
//!
//! Produces a side-by-side SVG layout:
//!
//! ```text
//!  [mol] + [mol] → [mol] + [mol]
//! ```

use chematic_rxn::Reaction;

use crate::layout::{BOND_LEN, compute_layout};
use crate::svg::{RenderOptions, render_svg_opts};

/// Render a `Reaction` as a self-contained SVG string using default options.
pub fn depict_reaction_svg(rxn: &Reaction) -> String {
    depict_reaction_svg_opts(rxn, &RenderOptions::default())
}

/// Render a `Reaction` as a self-contained SVG string with custom options.
///
/// The output places reactants to the left of a "→" arrow and products to
/// the right.  Each component is laid out independently at the same vertical
/// centre line.  A "+" separator is inserted between components on each side.
pub fn depict_reaction_svg_opts(rxn: &Reaction, opts: &RenderOptions) -> String {
    depict_reaction_svg_prepared(rxn, opts, &|mol| (None, compute_layout(mol)))
}

/// [`depict_reaction_svg_opts`] with each component drawn as `prepare`
/// gives it: an optional replacement molecule (e.g. with stereo wedges) and
/// its layout, in [`compute_layout`] units.
pub fn depict_reaction_svg_prepared(
    rxn: &Reaction,
    opts: &RenderOptions,
    prepare: &dyn Fn(&chematic_core::Molecule) -> (Option<chematic_core::Molecule>, crate::Layout),
) -> String {
    let pad = BOND_LEN; // padding between components
    let arrow_w = BOND_LEN * 2.0; // width of the reaction arrow
    let sep_w = BOND_LEN; // width of the "+" separator

    // -----------------------------------------------------------------------
    // Step 1: render each molecule to SVG and measure its bounding box.
    // -----------------------------------------------------------------------

    struct Component {
        svg_inner: String, // SVG content without outer <svg> wrapper
        w: f64,
        h: f64,
        // The content is drawn in layout coordinates: left edge and vertical
        // centre of the layout's bounding box.
        min_x: f64,
        mid_y: f64,
    }

    let render_mol = |mol: &chematic_core::Molecule| -> Component {
        let (copy, layout) = prepare(mol);
        let mol = copy.as_ref().unwrap_or(mol);
        let (min_x, min_y, max_x, max_y) = layout.bounding_box();
        let (mw, mh) = ((max_x - min_x).max(1.0), (max_y - min_y).max(1.0));
        let mol_svg = render_svg_opts(mol, &layout, opts);
        // Strip the outer <svg …> wrapper so we can re-embed with a translate.
        let inner = strip_svg_wrapper(&mol_svg);
        Component {
            svg_inner: inner,
            w: mw + pad,
            h: mh + pad,
            min_x,
            mid_y: (min_y + max_y) / 2.0,
        }
    };
    // Arrow and "+" follow the theme, like the bonds.
    let ink = if opts.dark { "white" } else { "#333" };

    let reactants: Vec<Component> = rxn.reactants.iter().map(render_mol).collect();
    let products: Vec<Component> = rxn.products.iter().map(render_mol).collect();

    let max_h = reactants
        .iter()
        .chain(products.iter())
        .map(|c| c.h)
        .fold(BOND_LEN * 4.0, f64::max);
    let total_h = max_h + pad;
    // Components, "+" and arrow share this horizontal centre line.
    let mid_y = total_h / 2.0;

    // -----------------------------------------------------------------------
    // Step 2: lay out horizontally.
    // -----------------------------------------------------------------------

    let mut svg_parts: Vec<String> = Vec::new();
    let mut cursor_x = pad / 2.0;

    // Each component takes a cell of width `w`, its bounding box centred in
    // it, so every gap between two boxes, or a box and the edge, is `pad` wide.
    let emit_components = |parts: &mut Vec<String>, comps: &[Component], cursor: &mut f64| {
        for (i, comp) in comps.iter().enumerate() {
            if i > 0 {
                // "+" separator
                let cx = *cursor + sep_w / 2.0;
                parts.push(format!(
                    "<text x=\"{cx:.1}\" y=\"{mid_y:.1}\" text-anchor=\"middle\" dominant-baseline=\"middle\" font-family=\"sans-serif\" font-size=\"18\" fill=\"{ink}\">+</text>"
                ));
                *cursor += sep_w;
            }
            let tx = *cursor + pad / 2.0 - comp.min_x;
            let ty = mid_y - comp.mid_y;
            parts.push(format!(
                r#"<g transform="translate({tx:.1},{ty:.1})">{}</g>"#,
                comp.svg_inner
            ));
            *cursor += comp.w;
        }
    };

    emit_components(&mut svg_parts, &reactants, &mut cursor_x);

    // Reaction arrow
    let arrow_x1 = cursor_x + pad / 2.0;
    let arrow_x2 = arrow_x1 + arrow_w;
    svg_parts.push(format!(
        "<line x1=\"{:.1}\" y1=\"{mid_y:.1}\" x2=\"{:.1}\" y2=\"{mid_y:.1}\" stroke=\"{ink}\" stroke-width=\"2\" marker-end=\"url(#arrowhead)\"/>",
        arrow_x1, arrow_x2
    ));
    cursor_x += pad / 2.0 + arrow_w + pad / 2.0;

    emit_components(&mut svg_parts, &products, &mut cursor_x);
    cursor_x += pad / 2.0;

    let total_w = cursor_x;

    // -----------------------------------------------------------------------
    // Step 3: assemble the final SVG.
    // -----------------------------------------------------------------------

    let defs = format!(
        "<defs><marker id=\"arrowhead\" markerWidth=\"10\" markerHeight=\"7\" refX=\"9\" refY=\"3.5\" orient=\"auto\"><polygon points=\"0 0, 10 3.5, 0 7\" fill=\"{ink}\"/></marker></defs>"
    );

    format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{total_w:.0}" height="{total_h:.0}" viewBox="0 0 {total_w:.1} {total_h:.1}">{defs}{}</svg>"#,
        svg_parts.join("")
    )
}

// ---------------------------------------------------------------------------
// Helper: strip the outer <svg …> … </svg> tags
// ---------------------------------------------------------------------------

fn strip_svg_wrapper(svg: &str) -> String {
    // Find the first `>` that closes the opening <svg …> tag.
    let start = svg.find('>').map(|i| i + 1).unwrap_or(0);
    let end = svg.rfind("</svg>").unwrap_or(svg.len());
    svg[start..end].to_string()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_rxn::parse_reaction;

    #[test]
    fn test_depict_reaction_svg_produces_svg() {
        let rxn = parse_reaction("CC>>CCO").unwrap();
        let svg = depict_reaction_svg(&rxn);
        assert!(svg.starts_with("<svg"), "must start with <svg");
        assert!(svg.ends_with("</svg>"), "must end with </svg>");
        assert!(svg.contains("arrowhead"), "must contain arrow marker");
    }

    #[test]
    fn test_depict_reaction_svg_two_reactants() {
        let rxn = parse_reaction("CC.N>>CCN").unwrap();
        let svg = depict_reaction_svg(&rxn);
        assert!(svg.contains('+'), "must contain '+' separator");
    }

    /// Value of the first `name="…"` attribute after `from`.
    fn attr(svg: &str, from: &str, name: &str) -> f64 {
        let rest = &svg[svg.find(from).unwrap()..];
        let key = format!(" {name}=\"");
        let rest = &rest[rest.find(&key).unwrap() + key.len()..];
        rest[..rest.find('"').unwrap()].parse().unwrap()
    }

    /// Attribute of the reaction arrow (bonds are `<line>`s too).
    fn arrow(svg: &str, name: &str) -> f64 {
        let end = svg.find("marker-end").unwrap();
        attr(&svg[svg[..end].rfind("<line").unwrap()..], "<line", name)
    }

    /// Each component's layout bounding box on the canvas, in drawing order.
    fn placed_boxes(rxn: &Reaction, svg: &str) -> Vec<(f64, f64, f64, f64)> {
        let mols: Vec<_> = rxn.reactants.iter().chain(&rxn.products).collect();
        let translations: Vec<(f64, f64)> = svg
            .match_indices("translate(")
            .map(|(i, _)| {
                let rest = &svg[i + "translate(".len()..];
                let (x, y) = rest[..rest.find(')').unwrap()].split_once(',').unwrap();
                (x.parse().unwrap(), y.parse().unwrap())
            })
            .collect();
        assert_eq!(translations.len(), mols.len());
        mols.iter()
            .zip(translations)
            .map(|(mol, (tx, ty))| {
                let (min_x, min_y, max_x, max_y) = compute_layout(mol).bounding_box();
                (tx + min_x, ty + min_y, tx + max_x, ty + max_y)
            })
            .collect()
    }

    /// Esterification: each component has atoms at negative layout
    /// coordinates and labels on its edges.
    const ESTERIFICATION: &str = "CC(=O)O.OCC>>CC(=O)OCC.O";

    #[test]
    fn components_are_inside_the_canvas_and_centred_on_the_arrow() {
        let rxn = parse_reaction(ESTERIFICATION).unwrap();
        let svg = depict_reaction_svg(&rxn);
        let (w, h) = (attr(&svg, "<svg", "width"), attr(&svg, "<svg", "height"));
        let arrow_y = arrow(&svg, "y1");
        for (left, top, right, bottom) in placed_boxes(&rxn, &svg) {
            // Labels reach about FONT_SIZE / 2 past the atom centres.
            assert!(
                left >= 10.0 && right <= w - 10.0,
                "x {left:.1}..{right:.1}, width {w}"
            );
            assert!(
                top >= 10.0 && bottom <= h - 10.0,
                "y {top:.1}..{bottom:.1}, height {h}"
            );
            let centre = (top + bottom) / 2.0;
            assert!(
                (centre - arrow_y).abs() < 0.2,
                "centre {centre:.1}, arrow at {arrow_y}"
            );
        }
    }

    #[test]
    fn separators_sit_midway_between_components() {
        let rxn = parse_reaction(ESTERIFICATION).unwrap();
        let svg = depict_reaction_svg(&rxn);
        let boxes = placed_boxes(&rxn, &svg);
        let arrow_mid = (arrow(&svg, "x1") + arrow(&svg, "x2")) / 2.0;
        let plus: Vec<f64> = svg
            .match_indices(">+</text>")
            .map(|(i, _)| attr(&svg[svg[..i].rfind("<text").unwrap()..], "<text", "x"))
            .collect();
        assert_eq!(plus.len(), 2);
        for (left_edge, mid, right_edge) in [
            (boxes[0].2, plus[0], boxes[1].0),
            (boxes[1].2, arrow_mid, boxes[2].0),
            (boxes[2].2, plus[1], boxes[3].0),
        ] {
            let (before, after) = (mid - left_edge, right_edge - mid);
            assert!(
                (before - after).abs() < 0.2,
                "{before:.1} before, {after:.1} after"
            );
        }
    }

    #[test]
    fn plus_names_a_font_and_dark_theme_recolours_arrow_and_plus() {
        let rxn = parse_reaction(ESTERIFICATION).unwrap();
        let svg = depict_reaction_svg(&rxn);
        // As for atom labels: a renderer without a default font draws nothing.
        assert!(
            svg.contains("font-family=\"sans-serif\" font-size=\"18\""),
            "{svg}"
        );
        let dark = RenderOptions {
            dark: true,
            background: "transparent".into(),
            ..RenderOptions::default()
        };
        let svg = depict_reaction_svg_opts(&rxn, &dark);
        assert!(!svg.contains("#333"), "{svg}");
        assert!(svg.contains("stroke=\"white\" stroke-width=\"2\""), "{svg}");
    }
}
