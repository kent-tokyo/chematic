use chematic_core::{AtomIdx, BondOrder};
use chematic_depict::{DepictBondKind, RenderOptions, compute_depict_data, depict_svg_opts};
use chematic_smiles::parse;

#[test]
fn alkene_directions_are_drawn_as_geometry() {
    for (smiles, cis) in [("C/C=C/C", false), (r"C/C=C\C", true)] {
        let mol = parse(smiles).unwrap();
        let data = compute_depict_data(&mol);
        let [a, b, c, d] = std::array::from_fn(|i| data.atoms[i].pos);
        let cross =
            |p: chematic_depict::Point| (c.x - b.x) * (p.y - b.y) - (c.y - b.y) * (p.x - b.x);
        assert!(cross(a).abs() > 1.0 && cross(d).abs() > 1.0);
        assert_eq!(cross(a) * cross(d) > 0.0, cis, "{smiles}");
        assert_eq!(
            data.bonds.iter().map(|b| &b.kind).collect::<Vec<_>>(),
            [
                &DepictBondKind::Single,
                &DepictBondKind::Double,
                &DepictBondKind::Single
            ]
        );
        for kekulize in [false, true] {
            let svg = depict_svg_opts(
                &mol,
                &RenderOptions {
                    kekulize,
                    ..Default::default()
                },
            );
            assert!(!svg.contains("<polygon"), "{svg}");
            assert_eq!(svg.matches("<line ").count(), 4, "{svg}");
        }
    }
}

#[test]
fn explicit_wedge_is_retained() {
    let mut mol = parse("F[C@](Cl)(Br)I").unwrap();
    let (idx, _) = mol.bond_between(AtomIdx(0), AtomIdx(1)).unwrap();
    mol.set_bond_order(idx, BondOrder::Up);
    let svg = depict_svg_opts(&mol, &RenderOptions::default());
    assert_eq!(svg.matches("<polygon ").count(), 1, "{svg}");
}

#[test]
fn isolated_labels_show_hydrogen_counts_isotopes_and_charge() {
    for (smiles, expected) in [
        ("[NH4+]", "NH₄⁺"),
        ("O", "H₂O"),
        ("[15NH4+]", "¹⁵NH₄⁺"),
        ("[Ca+2]", "Ca²⁺"),
        ("[Cl-]", "Cl⁻"),
    ] {
        let mol = parse(smiles).unwrap();
        for kekulize in [false, true] {
            let svg = depict_svg_opts(
                &mol,
                &RenderOptions {
                    kekulize,
                    ..Default::default()
                },
            );
            assert!(
                svg.contains(&format!(">{expected}</text>")),
                "{smiles}: {svg}"
            );
        }
    }
}
