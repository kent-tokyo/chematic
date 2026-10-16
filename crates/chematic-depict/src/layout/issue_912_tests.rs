use super::*;

/// The angle at an atom in the structured drawing data, in degrees.
fn depict_angle(data: &crate::DepictData, a: u32, center: u32, b: u32) -> f64 {
    let point = |idx: u32| data.atoms[idx as usize].pos;
    let (a, c, b) = (point(a), point(center), point(b));
    let (ux, uy) = (a.x - c.x, a.y - c.y);
    let (vx, vy) = (b.x - c.x, b.y - c.y);
    ((ux * vx + uy * vy) / ((ux * ux + uy * uy).sqrt() * (vx * vx + vy * vy).sqrt()))
        .clamp(-1.0, 1.0)
        .acos()
        .to_degrees()
}

#[test]
fn degree_two_sp_centers_are_linear_in_structured_depiction() {
    use chematic_smiles::parse;
    // The second half varies atom order and makes a linear centre the
    // first atom, so the rule cannot depend on a particular SMILES walk.
    for (smiles, axes) in [
        ("O=C=O", &[(0, 1, 2)][..]),
        ("CC#N", &[(0, 1, 2)]),
        ("CC#CC", &[(0, 1, 2), (1, 2, 3)]),
        ("C=C=C", &[(0, 1, 2)]),
        ("C(=O)=O", &[(1, 0, 2)]),
        ("N#CC", &[(0, 1, 2)]),
        ("C(#CC)C", &[(1, 0, 3), (0, 1, 2)]),
        ("C(=C)=C", &[(1, 0, 2)]),
    ] {
        let mol = parse(smiles).unwrap();
        let data = crate::compute_depict_data(&mol);
        for &(a, center, b) in axes {
            let angle = depict_angle(&data, a, center, b);
            assert!(
                (angle - 180.0).abs() < 1e-6,
                "{smiles}: center {center} has angle {angle}°"
            );
        }
        let layout = compute_layout(&mol);
        for (_, bond) in mol.bonds() {
            let length = layout.get(bond.atom1).dist(&layout.get(bond.atom2));
            assert!(
                (length - BOND_LEN).abs() < 1e-6,
                "{smiles}: bond length {length}"
            );
        }
    }
}

#[test]
fn ordinary_degree_two_chains_keep_their_bend() {
    use chematic_smiles::parse;
    for smiles in ["CCC", "COC"] {
        let mol = parse(smiles).unwrap();
        let data = crate::compute_depict_data(&mol);
        let angle = depict_angle(&data, 0, 1, 2);
        assert!((angle - 150.0).abs() < 1e-6, "{smiles}: angle {angle}°");
    }
}

#[test]
fn svg_bonds_follow_the_linear_sp_axis() {
    use chematic_smiles::parse;
    // The SVG renderer writes double bonds as parallel lines and a triple
    // bond as three lines. These coordinates check the actual output path.
    let co2 = crate::depict_svg(&parse("O=C=O").unwrap());
    assert!(co2.contains("<line x1=\"0.00\" y1=\"-2.00\" x2=\"31.00\" y2=\"-2.00\""));
    assert!(co2.contains("<line x1=\"49.00\" y1=\"-2.00\" x2=\"80.00\" y2=\"-2.00\""));
    assert!(co2.contains(">C</text>"));

    let nitrile = crate::depict_svg(&parse("CC#N").unwrap());
    assert!(nitrile.contains("<line x1=\"0.00\" y1=\"0.00\" x2=\"31.00\" y2=\"0.00\""));
    assert!(nitrile.contains("<line x1=\"49.00\" y1=\"0.00\" x2=\"80.00\" y2=\"0.00\""));
}
