use chematic_smiles::{rdkit_mol_from_pdb_block, rdkit_smiles};

fn atom(serial: usize, element: &str, v: [f64; 3]) -> String {
    format!(
        "HETATM{serial:>5} {element:<4} LIG A   1    {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {element:>2}  ",
        v[0], v[1], v[2], 1., 0.
    )
}

#[test]
fn pdb_non_tetrahedral_stereo_tracks_each_opposite_ligand_pair() {
    // Expected canonical strings were independently generated with RDKit 2025.09.3
    // using MolFromPDBBlock(sanitize=False, removeHs=False, proximityBonding=False).
    let cases = [
        (
            (0, 1),
            "[F][Pt@OH6]([Cl])([Br])[I]",
            "[F][Pt@TB5]([Cl])([Br])[I]",
        ),
        (
            (0, 2),
            "[F][Pt@OH3]([Cl])([Br])[I]",
            "[F][Pt@TB3]([Cl])([Br])[I]",
        ),
        (
            (0, 3),
            "[F][Pt@OH1]([Cl])([Br])[I]",
            "[F][Pt@TB1]([Cl])([Br])[I]",
        ),
        (
            (1, 2),
            "[F][Pt@OH20]([Cl])([Br])[I]",
            "[F][Pt@TB19]([Cl])([Br])[I]",
        ),
        (
            (1, 3),
            "[F][Pt@OH19]([Cl])([Br])[I]",
            "[F][Pt@TB20]([Cl])([Br])[I]",
        ),
        (
            (2, 3),
            "[F][Pt@OH21]([Cl])([Br])[I]",
            "[F][Pt@TB17]([Cl])([Br])[I]",
        ),
    ];
    for ((first, second), octahedral, bipyramidal) in cases {
        for (last, expected) in [([0., 0., 1.], octahedral), ([0., -0.5, 0.866], bipyramidal)] {
            let vectors = [[1., 0., 0.], [-1., 0., 0.], [0., 1., 0.], last];
            let mut points = [[0.; 3]; 4];
            points[first] = vectors[0];
            points[second] = vectors[1];
            let mut next = 2;
            for (i, point) in points.iter_mut().enumerate() {
                if i != first && i != second {
                    *point = vectors[next];
                    next += 1;
                }
            }
            let mut lines = vec![atom(1, "Pt", [0.; 3])];
            for (i, element) in ["F", "Cl", "Br", "I"].iter().enumerate() {
                lines.push(atom(i + 2, element, points[i]));
            }
            lines.push("CONECT    1    2    3    4    5".into());
            lines.push("END".into());
            let parsed = rdkit_mol_from_pdb_block(&lines.join("\n"), false, false, 0, false)
                .unwrap()
                .unwrap();
            let smiles = rdkit_smiles(&parsed.molecule, &Default::default()).unwrap();
            assert_eq!(smiles, expected, "opposite pair {first} {second}");
        }
    }
}
