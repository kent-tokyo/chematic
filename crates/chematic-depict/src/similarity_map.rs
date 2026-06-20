//! SimilarityMap: atom-colored SVG where each atom is tinted by a scalar weight.
//!
//! Positive weights → blue tint, negative weights → red tint, zero → white (no tint).
//! Useful for visualising per-atom property contributions such as LogP, TPSA, or
//! fingerprint similarity.
//!
//! # Usage
//!
//! ```
//! # use chematic_smiles::parse;
//! # use chematic_depict::similarity_map_svg;
//! let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
//! // Pretend we have per-atom LogP weights (index = AtomIdx).
//! let weights: Vec<f64> = vec![0.1, -0.2, 0.5, 0.3, 0.0, 0.1, 0.2, 0.3, -0.1, 0.4, 0.6, -0.3, 0.1];
//! let svg = similarity_map_svg(&mol, &weights);
//! assert!(svg.contains("<svg"));
//! ```

use chematic_core::{AtomIdx, Molecule};

use crate::svg::RenderOptions;

/// Render `mol` as an SVG with each atom tinted by the corresponding `weights` value.
///
/// `weights` must have the same length as the number of heavy atoms in `mol`
/// (`mol.atom_count()`). Out-of-range weights are silently ignored.
///
/// - Positive weights are coloured blue (`#7070FF`).
/// - Negative weights are coloured red (`#FF7070`).
/// - Zero weights produce white circles (invisible on the default white background).
///
/// All circles are rendered at **50% opacity** (matching the existing highlight style).
/// Atom labels and bond styles are unaffected.
///
/// If `weights` is shorter than `mol.atom_count()`, missing atoms get no tint.
pub fn similarity_map_svg(mol: &Molecule, weights: &[f64]) -> String {
    similarity_map_svg_opts(mol, weights, &RenderOptions::default())
}

/// Like [`similarity_map_svg`] but accepts custom [`RenderOptions`].
///
/// Any `atom_color_map` already present in `opts` will be **overridden** for atoms
/// covered by `weights`.  Other rendering options (size, padding, dark mode, …)
/// are respected as-is.
pub fn similarity_map_svg_opts(mol: &Molecule, weights: &[f64], opts: &RenderOptions) -> String {
    let mut opts = opts.clone();

    // Find the maximum absolute weight for normalisation.
    let max_abs = weights
        .iter()
        .copied()
        .map(f64::abs)
        .fold(0.0_f64, f64::max);

    if max_abs > 1e-9 {
        for (i, &w) in weights.iter().enumerate() {
            if i >= mol.atom_count() {
                break;
            }
            let t = (w / max_abs).clamp(-1.0, 1.0);
            let color = weight_to_hex(t);
            opts.atom_color_map.insert(AtomIdx(i as u32), color);
        }
    }

    crate::depict_svg_opts(mol, &opts)
}

/// Map a normalised weight `t ∈ [-1, 1]` to a CSS hex colour.
///
/// - `t > 0` → blue tint lerp(#FFFFFF, #7070FF, t)
/// - `t < 0` → red tint lerp(#FFFFFF, #FF7070, |t|)
/// - `t = 0` → `#FFFFFF` (white)
fn weight_to_hex(t: f64) -> String {
    let t = t.clamp(-1.0, 1.0);
    let abs_t = t.abs();

    // Both colour endpoints share the same "mid" channel value (0x70 = 112).
    // The full channel (0xFF = 255) varies by sign.
    let mid = 112.0_f64;
    let full = 255.0_f64;

    // lerp from white (255) toward the accent channel.
    let accent = (full - (full - mid) * abs_t).round() as u8;

    let (r, g, b) = if t >= 0.0 {
        // Blue: R and G dim together; B stays full.
        (accent, accent, 255u8)
    } else {
        // Red: G and B dim together; R stays full.
        (255u8, accent, accent)
    };

    format!("#{r:02X}{g:02X}{b:02X}")
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> Molecule {
        parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"))
    }

    #[test]
    fn test_weight_to_hex_extremes() {
        assert_eq!(weight_to_hex(0.0), "#FFFFFF");
        assert_eq!(weight_to_hex(1.0), "#7070FF");
        assert_eq!(weight_to_hex(-1.0), "#FF7070");
    }

    #[test]
    fn test_weight_to_hex_half() {
        // t = 0.5: accent = round(255 - 143 * 0.5) = round(255 - 71.5) = 184 = 0xB8
        let color = weight_to_hex(0.5);
        assert_eq!(color, "#B8B8FF");
        let color_neg = weight_to_hex(-0.5);
        assert_eq!(color_neg, "#FFB8B8");
    }

    #[test]
    fn test_similarity_map_svg_returns_valid_svg() {
        let m = mol("CC(=O)O");
        let weights = vec![0.1, -0.2, 0.5, 0.3];
        let svg = similarity_map_svg(&m, &weights);
        assert!(
            svg.starts_with("<svg") || svg.contains("<svg"),
            "not SVG: {}",
            &svg[..50.min(svg.len())]
        );
        assert!(svg.contains("fill"), "no fill attributes in SVG");
    }

    #[test]
    fn test_similarity_map_zero_weights_no_circles() {
        let m = mol("c1ccccc1");
        let weights = vec![0.0; 6];
        let svg = similarity_map_svg(&m, &weights);
        // All-zero weights → no atom_color_map entries → no colored circles
        assert!(
            !svg.contains("7070FF"),
            "should be no blue circles for zero weights"
        );
        assert!(
            !svg.contains("FF7070"),
            "should be no red circles for zero weights"
        );
    }

    #[test]
    fn test_similarity_map_short_weights() {
        // Fewer weights than atoms: should not panic.
        let m = mol("c1ccccc1");
        let weights = vec![0.5, -0.5]; // only 2 of 6 atoms
        let svg = similarity_map_svg(&m, &weights);
        assert!(svg.contains("<svg") || svg.contains("svg"));
    }

    #[test]
    fn test_similarity_map_contains_atom_colors() {
        let m = mol("CC");
        let weights = vec![1.0, -1.0];
        let svg = similarity_map_svg(&m, &weights);
        assert!(svg.contains("7070FF"), "blue atom circle missing");
        assert!(svg.contains("FF7070"), "red atom circle missing");
    }
}
