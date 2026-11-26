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

#[test]
fn pdb_five_coordinate_stereo_tracks_each_axial_ligand_pair() {
    // Independent reference: RDKit 2025.09.3, unsanitized PDB import.
    let cases = [
        ((0, 1), "[NH2][Pt@TB13]([F])([Cl])([Br])[I]"),
        ((0, 2), "[NH2][Pt@TB10]([F])([Cl])([Br])[I]"),
        ((0, 3), "[NH2][Pt@TB9]([F])([Cl])([Br])[I]"),
        ((0, 4), "[NH2][Pt@TB8]([F])([Cl])([Br])[I]"),
        ((1, 2), "[NH2][Pt@TB16]([F])([Cl])([Br])[I]"),
        ((1, 3), "[NH2][Pt@TB15]([F])([Cl])([Br])[I]"),
        ((1, 4), "[NH2][Pt@TB6]([F])([Cl])([Br])[I]"),
        ((2, 3), "[NH2][Pt@TB17]([F])([Cl])([Br])[I]"),
        ((2, 4), "[NH2][Pt@TB4]([F])([Cl])([Br])[I]"),
        ((3, 4), "[NH2][Pt@TB2]([F])([Cl])([Br])[I]"),
    ];
    let vectors = [
        [1., 0., 0.],
        [-1., 0., 0.],
        [0., 1., 0.],
        [0., -0.5, 0.866],
        [0., -0.5, -0.866],
    ];
    for ((first, second), expected) in cases {
        let mut points = [[0.; 3]; 5];
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
        for (i, element) in ["F", "Cl", "Br", "I", "N"].iter().enumerate() {
            lines.push(atom(i + 2, element, points[i]));
        }
        lines.extend([
            "CONECT    1    2    3    4    5".into(),
            "CONECT    1    6".into(),
            "END".into(),
        ]);
        let parsed = rdkit_mol_from_pdb_block(&lines.join("\n"), false, false, 0, false)
            .unwrap()
            .unwrap();
        assert_eq!(
            rdkit_smiles(&parsed.molecule, &Default::default()).unwrap(),
            expected,
            "axial pair {first} {second}"
        );
    }
}

#[test]
fn overlapping_ligand_and_center_are_rejected_during_3d_stereo_assignment() {
    let mut lines = vec![atom(1, "Pt", [0.; 3])];
    for (i, (element, point)) in [
        ("F", [0., 0., 0.]),
        ("Cl", [1., 0., 0.]),
        ("Br", [0., 1., 0.]),
        ("I", [0., 0., 1.]),
    ]
    .into_iter()
    .enumerate()
    {
        lines.push(atom(i + 2, element, point));
    }
    lines.extend(["CONECT    1    2    3    4    5".into(), "END".into()]);
    let error = rdkit_mol_from_pdb_block(&lines.join("\n"), false, false, 0, false)
        .err()
        .unwrap();
    assert!(error.to_string().contains("zero length vector"), "{error}");
}

#[test]
fn seven_coordinate_center_does_not_receive_an_unsupported_stereo_tag() {
    // RDKit 2025.09.3 unsanitized import also leaves the seven-coordinate Pt untagged.
    let mut lines = vec![atom(1, "Pt", [0., 0., 0.])];
    for (i, (element, point)) in [
        ("F", [1., 0., 0.]),
        ("Cl", [-1., 0., 0.]),
        ("Br", [0., 1., 0.]),
        ("I", [0., -1., 0.]),
        ("N", [0., 0., 1.]),
        ("O", [0., 0., -1.]),
        ("S", [1., 1., 1.]),
    ]
    .into_iter()
    .enumerate()
    {
        lines.push(atom(i + 2, element, point));
    }
    lines.extend([
        "CONECT    1    2    3    4    5".into(),
        "CONECT    1    6    7    8".into(),
        "END".into(),
    ]);
    let parsed = rdkit_mol_from_pdb_block(&lines.join("\n"), false, false, 0, false)
        .unwrap()
        .unwrap();
    let center = chematic_core::AtomIdx(0);
    assert_eq!(parsed.molecule.neighbors(center).count(), 7);
    assert_eq!(
        parsed.molecule.atom(center).chirality,
        chematic_core::Chirality::None
    );
}

#[test]
fn ambiguous_opposite_ligand_pairs_are_left_untagged() {
    // RDKit 2025.09.3 unsanitized import returns CHI_UNSPECIFIED for this geometry.
    let mut lines = vec![atom(1, "Pt", [0.; 3])];
    for (i, (element, point)) in [
        ("F", [1., 0., 0.]),
        ("Cl", [-1., 0., 0.]),
        ("Br", [-1., 0., 0.]),
        ("I", [0., 1., 0.]),
    ]
    .into_iter()
    .enumerate()
    {
        lines.push(atom(i + 2, element, point));
    }
    lines.extend(["CONECT    1    2    3    4    5".into(), "END".into()]);
    let parsed = rdkit_mol_from_pdb_block(&lines.join("\n"), false, false, 0, false)
        .unwrap()
        .unwrap();
    assert_eq!(
        parsed.molecule.atom(chematic_core::AtomIdx(0)).chirality,
        chematic_core::Chirality::None
    );
}
