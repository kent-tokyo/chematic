use crate::*;
use chematic_core::{AtomIdx, MoleculeBuilder};

#[test]
fn public_depiction_entrypoints_handle_empty_graphs_and_retain_atom_metadata() {
    let empty = MoleculeBuilder::new().build();
    assert!(compute_layout(&empty).coords.is_empty());
    assert_eq!(compute_layout(&empty).bounding_box(), (0.0, 0.0, 0.0, 0.0));
    assert!(compute_depict_data(&empty).atoms.is_empty());
    assert!(depict_svg(&empty).ends_with("</svg>"));
    assert!(depict_eps(&empty).starts_with("%!PS-Adobe-3.0 EPSF-3.0"));
    for source in [
        "[13CH3][NH2+]C(=O)[O-]",
        "FC(Cl)(Br)I",
        "c1ccncc1",
        "C#N",
        "C=C",
    ] {
        let mol = chematic_smiles::parse(source).unwrap();
        let data = compute_depict_data(&mol);
        assert_eq!(data.atoms.len(), mol.atom_count());
        assert_eq!(data.bonds.len(), mol.bond_count());
        for atom in data.atoms {
            assert_eq!(atom.element, mol.atom(atom.idx).element);
            assert_eq!(atom.charge, mol.atom(atom.idx).charge);
            assert!(atom.pos.x.is_finite() && atom.pos.y.is_finite());
        }
        let mut opts = RenderOptions::with_cpk_colors_for(&mol);
        opts.width = Some(320);
        opts.height = Some(240);
        opts.dark = true;
        opts.background = "#123456".into();
        let eps = depict_eps_opts(&mol, &opts);
        assert!(eps.contains("%%BoundingBox: 0 0 320 240"));
        assert!(!eps.contains("NaN") && !eps.contains("Infinity"));
        let grid = grid::depict_svg_grid_with_opts(&[(&mol, None), (&mol, Some(&opts))], 2);
        assert!(grid.contains("mol-0") && grid.contains("mol-1"));
        assert!(grid::depict_svg_grid_with_opts(&[(&mol, None)], 0).contains("width=\"0\""));
        let partial = depict_data_with_coords(&mol, &[(7.0, -3.0)]);
        assert_eq!(
            (partial.atoms[0].pos.x, partial.atoms[0].pos.y),
            (7.0, -3.0)
        );
        assert_eq!((partial.atoms[1].pos.x, partial.atoms[1].pos.y), (0.0, 0.0));
        for (idx, atom) in mol.atoms() {
            if atom.element.atomic_number() == 6 {
                assert!(!opts.atom_color_map.contains_key(&idx));
            } else {
                assert!(opts.atom_color_map.contains_key(&idx));
            }
        }
        assert_eq!(partial.atoms[0].idx, AtomIdx(0));
    }
}
