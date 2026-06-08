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
    }

    let render_mol = |mol: &chematic_core::Molecule| -> Component {
        let layout = compute_layout(mol);
        let (min_x, min_y, max_x, max_y) = layout.bounding_box();
        let (mw, mh) = ((max_x - min_x).max(1.0), (max_y - min_y).max(1.0));
        let mol_svg = render_svg_opts(mol, &layout, opts);
        // Strip the outer <svg …> wrapper so we can re-embed with a translate.
        let inner = strip_svg_wrapper(&mol_svg);
        Component {
            svg_inner: inner,
            w: mw + pad,
            h: mh + pad,
        }
    };

    let reactants: Vec<Component> = rxn.reactants.iter().map(render_mol).collect();
    let products: Vec<Component> = rxn.products.iter().map(render_mol).collect();

    let max_h = reactants
        .iter()
        .chain(products.iter())
        .map(|c| c.h)
        .fold(BOND_LEN * 4.0, f64::max);

    // -----------------------------------------------------------------------
    // Step 2: lay out horizontally.
    // -----------------------------------------------------------------------

    let mut svg_parts: Vec<String> = Vec::new();
    let mut cursor_x = pad;

    let emit_components = |parts: &mut Vec<String>,
                           comps: &[Component],
                           cursor: &mut f64,
                           max_h: f64| {
        for (i, comp) in comps.iter().enumerate() {
            if i > 0 {
                // "+" separator
                let cx = *cursor + sep_w / 2.0;
                let cy = max_h / 2.0;
                parts.push(format!(
                    "<text x=\"{cx:.1}\" y=\"{cy:.1}\" text-anchor=\"middle\" dominant-baseline=\"middle\" font-size=\"18\" fill=\"#333\">+</text>"
                ));
                *cursor += sep_w;
            }
            let tx = *cursor;
            let ty = (max_h - comp.h) / 2.0;
            parts.push(format!(
                r#"<g transform="translate({tx:.1},{ty:.1})">{}</g>"#,
                comp.svg_inner
            ));
            *cursor += comp.w;
        }
    };

    emit_components(&mut svg_parts, &reactants, &mut cursor_x, max_h);

    // Reaction arrow
    let arrow_x1 = cursor_x + pad / 2.0;
    let arrow_x2 = arrow_x1 + arrow_w;
    let arrow_y = max_h / 2.0;
    svg_parts.push(format!(
        "<line x1=\"{:.1}\" y1=\"{arrow_y:.1}\" x2=\"{:.1}\" y2=\"{arrow_y:.1}\" stroke=\"#333\" stroke-width=\"2\" marker-end=\"url(#arrowhead)\"/>",
        arrow_x1, arrow_x2
    ));
    cursor_x += pad / 2.0 + arrow_w + pad / 2.0;

    emit_components(&mut svg_parts, &products, &mut cursor_x, max_h);
    cursor_x += pad;

    let total_w = cursor_x;
    let total_h = max_h + pad;

    // -----------------------------------------------------------------------
    // Step 3: assemble the final SVG.
    // -----------------------------------------------------------------------

    let defs = "<defs><marker id=\"arrowhead\" markerWidth=\"10\" markerHeight=\"7\" refX=\"9\" refY=\"3.5\" orient=\"auto\"><polygon points=\"0 0, 10 3.5, 0 7\" fill=\"#333\"/></marker></defs>";

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
}
