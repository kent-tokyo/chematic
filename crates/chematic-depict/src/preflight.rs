//! Deterministic, conservative publication checks for SVG consumers.
//!
//! This catches malformed geometry, obvious clipping/overlap, and crossings before
//! publication. It is not a substitute for pixel inspection by the target renderer.

use crate::{Layout, Point, RenderOptions, atom_display_label, detect_crossings};
use chematic_core::Molecule;
use serde::{Deserialize, Serialize};

const LABEL_HALF_W: f64 = 8.0;
const LABEL_HALF_H: f64 = 7.0;
const ATOM_RADIUS: f64 = 6.0;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PreflightCode {
    NonFiniteGeometry,
    NegativeStyle,
    ClippedAtom,
    ClippedBond,
    ClippedLabel,
    AtomOverlap,
    LabelOverlap,
    CrossingBond,
    DegenerateBond,
    ResourceLimit,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct PreflightDiagnostic {
    pub code: PreflightCode,
    pub path: String,
    pub message: String,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PreflightStyle {
    pub font_family: String,
    pub font_size: f64,
    pub bond_width: f64,
    pub padding: f64,
    pub viewport_width: u32,
    pub viewport_height: u32,
    pub font_metrics: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct PreflightLimits {
    pub max_atoms: usize,
    pub max_bonds: usize,
    pub max_labels: usize,
}

impl Default for PreflightLimits {
    fn default() -> Self {
        Self {
            max_atoms: 10_000,
            max_bonds: 20_000,
            max_labels: 10_000,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PreflightReport {
    pub schema_version: u32,
    pub ready: bool,
    pub diagnostics: Vec<PreflightDiagnostic>,
    pub style: PreflightStyle,
    pub fingerprint: String,
}

/// Run deterministic conservative checks over a molecule and its layout.
pub fn preflight_svg(
    mol: &Molecule,
    layout: &Layout,
    opts: &RenderOptions,
    limits: &PreflightLimits,
) -> PreflightReport {
    let style = effective_style(layout, opts);
    let mut diagnostics = Vec::new();
    let labels = display_labels(mol, layout);
    check_preflight_limits(mol, layout, &style, limits, labels.len(), &mut diagnostics);
    check_finite_coordinates(layout, &mut diagnostics);
    if geometry_is_complete(mol, layout) {
        check_complete_geometry(mol, layout, &style, &labels, &mut diagnostics);
    }
    let fingerprint = fingerprint(mol, layout, &style);
    PreflightReport {
        schema_version: 1,
        ready: diagnostics.is_empty(),
        diagnostics,
        style,
        fingerprint,
    }
}

type DisplayLabel = (usize, Point, String);

fn diagnostic(
    diagnostics: &mut Vec<PreflightDiagnostic>,
    code: PreflightCode,
    path: String,
    message: &str,
) {
    diagnostics.push(PreflightDiagnostic {
        code,
        path,
        message: message.to_string(),
    });
}

fn display_labels(mol: &Molecule, layout: &Layout) -> Vec<DisplayLabel> {
    mol.atoms()
        .filter_map(|(idx, _)| {
            layout
                .coords
                .get(idx.0 as usize)
                .copied()
                .map(|point| (idx.0 as usize, point, atom_display_label(mol, idx)))
        })
        .filter(|(_, _, label)| !label.is_empty())
        .collect()
}

fn check_preflight_limits(
    mol: &Molecule,
    layout: &Layout,
    style: &PreflightStyle,
    limits: &PreflightLimits,
    label_count: usize,
    diagnostics: &mut Vec<PreflightDiagnostic>,
) {
    if mol.atom_count() > limits.max_atoms {
        diagnostic(
            diagnostics,
            PreflightCode::ResourceLimit,
            "molecule.atoms".into(),
            "atom count exceeds limit",
        );
    }
    if mol.bond_count() > limits.max_bonds {
        diagnostic(
            diagnostics,
            PreflightCode::ResourceLimit,
            "molecule.bonds".into(),
            "bond count exceeds limit",
        );
    }
    if layout.coords.len() != mol.atom_count() {
        diagnostic(
            diagnostics,
            PreflightCode::ResourceLimit,
            "layout.coords".into(),
            "layout coordinate count does not match atom count",
        );
    }
    if !style.padding.is_finite()
        || style.padding < 0.0
        || style.font_size <= 0.0
        || style.bond_width <= 0.0
    {
        diagnostic(
            diagnostics,
            PreflightCode::NegativeStyle,
            "style".into(),
            "style values must be finite and positive",
        );
    }
    if label_count > limits.max_labels {
        diagnostic(
            diagnostics,
            PreflightCode::ResourceLimit,
            "molecule.labels".into(),
            "label count exceeds limit",
        );
    }
}

fn check_finite_coordinates(layout: &Layout, diagnostics: &mut Vec<PreflightDiagnostic>) {
    for (index, point) in layout.coords.iter().enumerate() {
        if !point.x.is_finite() || !point.y.is_finite() {
            diagnostic(
                diagnostics,
                PreflightCode::NonFiniteGeometry,
                format!("atoms[{index}].geometry"),
                "atom coordinates must be finite",
            );
        }
    }
}

fn geometry_is_complete(mol: &Molecule, layout: &Layout) -> bool {
    layout.coords.len() == mol.atom_count()
        && layout
            .coords
            .iter()
            .all(|point| point.x.is_finite() && point.y.is_finite())
}

fn check_complete_geometry(
    mol: &Molecule,
    layout: &Layout,
    style: &PreflightStyle,
    labels: &[DisplayLabel],
    diagnostics: &mut Vec<PreflightDiagnostic>,
) {
    let (min_x, min_y, max_x, max_y) = layout.bounding_box();
    let view = (
        min_x - style.padding,
        min_y - style.padding,
        (max_x - min_x).max(crate::BOND_LEN) + 2.0 * style.padding,
        (max_y - min_y).max(crate::BOND_LEN) + 2.0 * style.padding,
    );
    for (index, point) in layout.coords.iter().enumerate() {
        if outside(view, *point, ATOM_RADIUS) {
            diagnostic(
                diagnostics,
                PreflightCode::ClippedAtom,
                format!("atoms[{index}].geometry"),
                "atom extends outside the SVG viewBox",
            );
        }
    }
    for (index, bond) in mol.bonds() {
        let Some(a) = layout.coords.get(bond.atom1.0 as usize).copied() else {
            continue;
        };
        let Some(b) = layout.coords.get(bond.atom2.0 as usize).copied() else {
            continue;
        };
        if a.dist(&b) <= f64::EPSILON {
            diagnostic(
                diagnostics,
                PreflightCode::DegenerateBond,
                format!("bonds[{}].geometry", index.0),
                "bond endpoints are coincident",
            );
        } else if segment_outside(view, a, b) {
            diagnostic(
                diagnostics,
                PreflightCode::ClippedBond,
                format!("bonds[{}].geometry", index.0),
                "bond extends outside the SVG viewBox",
            );
        }
    }
    for (index, (atom_i, point, _)) in labels.iter().enumerate() {
        if outside(view, *point, LABEL_HALF_W.max(LABEL_HALF_H)) {
            diagnostic(
                diagnostics,
                PreflightCode::ClippedLabel,
                format!("atoms[{atom_i}].label"),
                "label may be clipped by the SVG viewBox",
            );
        }
        for (atom_j, other, _) in labels.iter().skip(index + 1) {
            if rect_overlap(
                (*point, LABEL_HALF_W, LABEL_HALF_H),
                (*other, LABEL_HALF_W, LABEL_HALF_H),
            ) {
                diagnostic(
                    diagnostics,
                    PreflightCode::LabelOverlap,
                    format!("atoms[{atom_i}].label/atoms[{atom_j}].label"),
                    "label boxes overlap",
                );
            }
        }
        if layout
            .coords
            .iter()
            .enumerate()
            .any(|(other_index, other)| {
                other_index != *atom_i && point.dist(other) < ATOM_RADIUS + LABEL_HALF_H
            })
        {
            diagnostic(
                diagnostics,
                PreflightCode::AtomOverlap,
                format!("atoms[{atom_i}].label"),
                "label overlaps a nearby atom",
            );
        }
    }
    for (a, b) in detect_crossings(layout, mol) {
        diagnostic(
            diagnostics,
            PreflightCode::CrossingBond,
            format!("bonds[{}]/bonds[{}]", a.0, b.0),
            "bonds cross",
        );
    }
}

pub fn preflight_svg_json(
    mol: &Molecule,
    layout: &Layout,
    opts: &RenderOptions,
    limits: &PreflightLimits,
) -> String {
    serde_json::to_string(&preflight_svg(mol, layout, opts, limits))
        .expect("preflight report is serializable")
}

fn effective_style(layout: &Layout, opts: &RenderOptions) -> PreflightStyle {
    let (min_x, min_y, max_x, max_y) = layout.bounding_box();
    let width = opts.width.unwrap_or_else(|| {
        ((max_x - min_x).max(crate::BOND_LEN) + 2.0 * opts.padding)
            .round()
            .max(0.0) as u32
    });
    let height = opts.height.unwrap_or_else(|| {
        ((max_y - min_y).max(crate::BOND_LEN) + 2.0 * opts.padding)
            .round()
            .max(0.0) as u32
    });
    PreflightStyle {
        font_family: "sans-serif".into(),
        font_size: 12.0,
        bond_width: 1.5,
        padding: opts.padding,
        viewport_width: width,
        viewport_height: height,
        font_metrics: "approximate; target renderer authoritative".into(),
    }
}

fn outside(view: (f64, f64, f64, f64), p: Point, radius: f64) -> bool {
    p.x - radius < view.0
        || p.y - radius < view.1
        || p.x + radius > view.0 + view.2
        || p.y + radius > view.1 + view.3
}

fn segment_outside(view: (f64, f64, f64, f64), a: Point, b: Point) -> bool {
    outside(view, a, 0.75) || outside(view, b, 0.75)
}

fn rect_overlap(a: (Point, f64, f64), b: (Point, f64, f64)) -> bool {
    a.0.x - a.1 < b.0.x + b.1
        && b.0.x - b.1 < a.0.x + a.1
        && a.0.y - a.2 < b.0.y + b.2
        && b.0.y - b.2 < a.0.y + a.2
}

fn fingerprint(mol: &Molecule, layout: &Layout, style: &PreflightStyle) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    let mut add = |s: String| {
        for byte in s.as_bytes() {
            h ^= u64::from(*byte);
            h = h.wrapping_mul(0x100000001b3);
        }
    };
    add(format!("style:{style:?};"));
    for (i, p) in layout.coords.iter().enumerate() {
        add(format!("a{i}:{:.9},{:.9};", p.x, p.y));
    }
    for (idx, bond) in mol.bonds() {
        add(format!(
            "b{}:{}-{}:{:?};",
            idx.0, bond.atom1.0, bond.atom2.0, bond.order
        ));
    }
    format!("fnv1a64-{h:016x}")
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn clean_report_is_ready_and_stable() {
        let mol = parse("CCO").unwrap();
        let layout = crate::compute_layout(&mol);
        let opts = RenderOptions::default();
        let a = preflight_svg_json(&mol, &layout, &opts, &PreflightLimits::default());
        let b = preflight_svg_json(&mol, &layout, &opts, &PreflightLimits::default());
        assert_eq!(a, b);
        assert!(preflight_svg(&mol, &layout, &opts, &PreflightLimits::default()).ready);
    }

    #[test]
    fn reports_invalid_overlap_and_clipping() {
        let mol = parse("NN").unwrap();
        let layout = Layout {
            coords: vec![Point::new(f64::NAN, 0.0), Point::new(0.0, 0.0)],
        };
        let report = preflight_svg(
            &mol,
            &layout,
            &RenderOptions::default(),
            &PreflightLimits::default(),
        );
        assert!(
            report
                .diagnostics
                .iter()
                .any(|d| d.code == PreflightCode::NonFiniteGeometry)
        );
        let layout = Layout {
            coords: vec![Point::new(0.0, 0.0), Point::new(0.0, 0.0)],
        };
        let opts = RenderOptions {
            padding: -1.0,
            ..RenderOptions::default()
        };
        let report = preflight_svg(&mol, &layout, &opts, &PreflightLimits::default());
        assert!(
            report
                .diagnostics
                .iter()
                .any(|d| d.code == PreflightCode::NegativeStyle)
        );
        assert!(
            report
                .diagnostics
                .iter()
                .any(|d| d.code == PreflightCode::DegenerateBond)
        );
    }
}
