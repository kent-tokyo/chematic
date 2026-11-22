use chematic_smiles::rdkit_mol_from_mol2_block;

fn alkene(points: [[f64; 3]; 4]) -> String {
    let mut text =
        "@<TRIPOS>MOLECULE\nalkene\n4 3 0 0 0\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n".to_owned();
    for (i, ((name, kind), [x, y, z])) in [("F", "F"), ("C", "C.2"), ("C", "C.2"), ("Cl", "Cl")]
        .into_iter()
        .zip(points)
        .enumerate()
    {
        text.push_str(&format!(
            "{} {name} {x:.3} {y:.3} {z:.3} {kind} 1 RES 0\n",
            i + 1
        ));
    }
    text.push_str("@<TRIPOS>BOND\n1 1 2 1\n2 2 3 2\n3 3 4 1\n");
    text
}

#[test]
fn mol2_alkene_geometry_preserves_or_omits_stereo_as_supported() {
    // Independent reference: RDKit 2025.09.3 MolFromMol2Block with sanitization.
    let cases = [
        (
            [[-1., 1., 0.], [0., 0., 0.], [1., 0., 0.], [2., -1., 0.]],
            r"F/[C]=[C]/Cl",
        ),
        (
            [[-1., 1., 0.], [0., 0., 0.], [1., 0., 0.], [2., 1., 0.]],
            r"F/[C]=[C]\Cl",
        ),
        (
            [[-1., 0., 0.], [0., 0., 0.], [1., 0., 0.], [2., 1., 0.]],
            "F[C]=[C]Cl",
        ),
        (
            [[-1., 1., 0.], [0., 0., 0.], [1., 0., 0.], [2., 0., 0.]],
            "F[C]=[C]Cl",
        ),
    ];
    for (points, expected) in cases {
        let parsed = rdkit_mol_from_mol2_block(&alkene(points), true, true, true).unwrap();
        assert_eq!(parsed.smiles, expected, "{points:?}");
    }
}
