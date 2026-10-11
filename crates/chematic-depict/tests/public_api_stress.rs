//! Depiction public-surface robustness on special chemistry and style modes.

use std::collections::HashSet;

use chematic_core::{AtomIdx, BondIdx};
use chematic_depict::*;
use chematic_rxn::parse_reaction;
use chematic_smiles::parse;

const STRESS: &str = include_str!("../../../validation/cosmolkit_comparison/corpora/stress-v1.smi");

#[test]
fn vector_depiction_surface_handles_special_chemistry_and_styles() {
    let mut molecules = Vec::new();
    for line in STRESS.lines() {
        if molecules.len() == 256 {
            break;
        }
        let smiles = line.split_whitespace().next().unwrap_or("");
        let Ok(mol) = parse(smiles) else { continue };
        if mol.atom_count() == 0 || mol.atom_count() > 100 {
            continue;
        }
        let row = molecules.len();
        let mut opts = if row % 3 == 0 {
            RenderOptions::with_cpk_colors_for(&mol)
        } else {
            RenderOptions::default()
        };
        opts.width = (row % 5 == 0).then_some(320);
        opts.height = (row % 7 == 0).then_some(240);
        opts.padding = if row % 11 == 0 { 0.0 } else { 20.0 };
        opts.background = if row % 13 == 0 {
            "transparent".into()
        } else if row % 17 == 0 {
            "#202124".into()
        } else {
            "white".into()
        };
        opts.dark = row % 17 == 0;
        opts.atom_ids = row % 19 == 0;
        opts.show_atom_indices = row % 23 == 0;
        opts.kekulize = row % 2 == 0;
        if mol.atom_count() > 0 {
            opts.highlight_atoms.insert(AtomIdx(0));
            opts.atom_color_map.insert(AtomIdx(0), "#ff00ff".into());
        }
        if mol.bond_count() > 0 {
            opts.highlight_bonds.insert(BondIdx(0));
        }

        let layout = compute_layout(&mol);
        let _ = layout.bounding_box();
        let _ = detect_crossings(&layout, &mol);
        if mol.atom_count() > 0 {
            let _ = suggest_bond_direction(&mol, AtomIdx(0), &layout);
            let _ = atom_label_with_h(&mol, AtomIdx(0));
            let _ = atom_display_label(&mol, AtomIdx(0));
        }
        let _ = compute_depict_data(&mol);
        let coords = layout.coords.iter().map(|p| (p.x, p.y)).collect::<Vec<_>>();
        let _ = depict_data_with_coords(&mol, &coords);
        let svg = depict_svg_opts(&mol, &opts);
        assert!(svg.starts_with("<svg"));
        let _ = depict_svg(&mol);
        let _ = depict_svg_highlighted(&mol, &opts.highlight_atoms, &opts.highlight_bonds);
        let _ = render_svg(&mol, &layout);
        let _ = render_svg_opts(&mol, &layout, &opts);
        let _ = render_svg_with_metadata(&mol, &layout, &opts, smiles);
        let _ = depict_eps(&mol);
        let _ = depict_eps_opts(&mol, &opts);
        let weights = (0..mol.atom_count())
            .map(|i| i as f64 / mol.atom_count().max(1) as f64 - 0.5)
            .collect::<Vec<_>>();
        let _ = similarity_map_svg(&mol, &weights);
        let _ = similarity_map_svg_opts(&mol, &weights, &opts);
        let _ = similarity_map_options(&mol, &weights, &opts);
        let _ = preflight_svg(&mol, &layout, &opts, &PreflightLimits::default());
        let _ = preflight_svg_json(&mol, &layout, &opts, &PreflightLimits::default());
        for atomic_number in [0, 1, 6, 7, 8, 9, 15, 16, 17, 35, 53, 78, 118] {
            let _ = atom_color(atomic_number);
            let _ = atom_color_rgb(atomic_number);
        }
        molecules.push((mol, layout, opts));
    }
    assert_eq!(molecules.len(), 256);

    for chunk in molecules.chunks(8) {
        let refs = chunk.iter().map(|(mol, _, _)| mol).collect::<Vec<_>>();
        let layouts = chunk
            .iter()
            .map(|(_, layout, _)| Layout {
                coords: layout.coords.clone(),
            })
            .collect::<Vec<_>>();
        let with_opts = chunk
            .iter()
            .map(|(mol, _, opts)| (mol, Some(opts)))
            .collect::<Vec<_>>();
        let _ = depict_svg_grid(&refs, 3);
        let _ = depict_svg_grid_with_layouts(&refs, &layouts, 3);
        let _ = depict_svg_grid_with_opts(&with_opts, 3);
        let _ = depict_svg_grid_with_opts_and_layouts(&with_opts, &layouts, 3);
    }
    assert!(depict_svg_grid(&[], 0).contains("width=\"0\""));
}

#[test]
fn reaction_and_preflight_boundary_modes_are_reported() {
    let reaction = parse_reaction("CCO.O=C=O>>CCOC=O").unwrap();
    let opts = RenderOptions {
        dark: true,
        background: "transparent".into(),
        show_atom_indices: true,
        ..RenderOptions::default()
    };
    assert!(depict_reaction_svg(&reaction).starts_with("<svg"));
    assert!(depict_reaction_svg_opts(&reaction, &opts).starts_with("<svg"));

    let mol = parse("N#N").unwrap();
    let malformed = Layout {
        coords: vec![Point::new(f64::NAN, 0.0)],
    };
    let limits = PreflightLimits {
        max_atoms: 0,
        max_bonds: 0,
        max_labels: 0,
    };
    let report = preflight_svg(&mol, &malformed, &opts, &limits);
    assert!(!report.ready);
    assert!(!report.diagnostics.is_empty());

    let coincident = Layout {
        coords: vec![Point::new(0.0, 0.0), Point::new(0.0, 0.0)],
    };
    let _ = preflight_svg(
        &mol,
        &coincident,
        &RenderOptions::default(),
        &Default::default(),
    );
    let mut crossing = coincident;
    relieve_layout_clashes(&mol, &mut crossing);
    let _ = render_svg_highlighted(&mol, &crossing, &HashSet::new(), &HashSet::new());
}

#[test]
fn layout_engine_handles_the_full_special_chemistry_corpus() {
    let mut parsed = 0usize;
    for line in STRESS.lines() {
        let smiles = line.split_whitespace().next().unwrap_or("");
        let Ok(mol) = parse(smiles) else { continue };
        if mol.atom_count() == 0 || mol.atom_count() > 160 {
            continue;
        }
        parsed += 1;
        let mut layout = compute_layout(&mol);
        assert_eq!(layout.coords.len(), mol.atom_count());
        let _ = layout.bounding_box();
        let _ = detect_crossings(&layout, &mol);
        relieve_layout_clashes(&mol, &mut layout);
        if mol.atom_count() > 0 {
            let _ = suggest_bond_direction(&mol, AtomIdx(0), &layout);
        }
    }
    assert!(parsed > 3_500);
}
