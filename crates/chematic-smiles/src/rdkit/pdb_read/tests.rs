use super::*;

fn pdb_atom(
    record: &str,
    serial: i32,
    name: &str,
    residue: (&str, i32),
    xyz: [f64; 3],
    element_charge: (&str, &str),
) -> String {
    let (residue, residue_no) = residue;
    let (element, charge) = element_charge;
    format!(
        "{record:<6}{serial:>5} {name:4} {residue:>3} A{residue_no:>4}    {x:>8.3}{y:>8.3}{z:>8.3}{occupancy:>6.2}{temperature:>6.2}          {element:>2}{charge:>2}\n",
        x = xyz[0],
        y = xyz[1],
        z = xyz[2],
        occupancy = 1.0,
        temperature = 0.0,
    )
}

fn residue(name: &[u8; 4], res_name: &[u8], res_no: i32, hetero: bool) -> ResidueInfo {
    ResidueInfo {
        name: *name,
        res_name: res_name.to_vec(),
        chain: b'A',
        insertion_code: b' ',
        res_no,
        hetero,
    }
}

#[test]
fn numeric_fields_and_symbols_follow_rdkit_edge_cases() {
    assert_eq!(to_int(b"  -42"), Ok(-42));
    assert_eq!(to_int(b" +12 "), Ok(0));
    assert_eq!(to_int(b"     "), Ok(0));
    assert_eq!(to_int(b"12-3 "), Ok(12));
    assert_eq!(to_int(b"999999999999999999"), Ok(0));
    assert!(to_int(b"12x").is_err());

    assert_eq!(to_double(b" -12.50 "), Ok(-12.5));
    assert_eq!(to_double(b"+.25"), Ok(0.25));
    assert_eq!(to_double(b"12,34"), Ok(12.0));
    assert_eq!(to_double(b" . "), Ok(0.0));
    assert!(to_double(b"1e3").is_err());

    assert_eq!(atom_from_symbol("D").unwrap(), Some((1, 2)));
    assert_eq!(atom_from_symbol("T").unwrap(), Some((1, 3)));
    assert_eq!(atom_from_symbol("Cl").unwrap(), Some((17, 0)));
    assert!(atom_from_symbol("Qq").is_err());
}

#[test]
fn reads_isotopes_charges_and_residue_metadata() {
    let block = [
        pdb_atom(
            "HETATM",
            1,
            " D1 ",
            ("DOD", 7),
            [0.0, 0.0, 0.1],
            ("D", "+2"),
        ),
        pdb_atom(
            "ATOM",
            2,
            " CL ",
            ("ALA", -2),
            [1.2, 0.0, 0.0],
            ("CL", "2-"),
        ),
        "END\n".to_string(),
    ]
    .concat();
    let pdb = mol_from_pdb_block(&block, false, false, 0, false)
        .unwrap()
        .unwrap();
    assert_eq!(pdb.mol.atoms.len(), 2);
    assert_eq!((pdb.mol.atoms[0].anum, pdb.mol.atoms[0].isotope), (1, 2));
    assert_eq!(pdb.mol.atoms[0].charge, 2);
    assert_eq!(pdb.mol.atoms[1].charge, -2);
    assert_eq!(pdb.info[0].res_name, b"DOD");
    assert_eq!(pdb.info[1].res_no, -2);
    assert!(pdb.info[0].hetero);
    assert!(!pdb.info[1].hetero);
    assert!(pdb.is_3d);
}

#[test]
fn pdb_filters_and_flavor_override_are_explicit() {
    let mut alternate = pdb_atom("HETATM", 1, " C1 ", ("UNL", 1), [0.0; 3], ("C", ""));
    alternate.replace_range(16..17, "B");
    assert_eq!(
        mol_from_pdb_block(&alternate, false, false, 0, false)
            .unwrap()
            .unwrap()
            .mol
            .atoms
            .len(),
        0
    );
    assert_eq!(
        mol_from_pdb_block(&alternate, false, false, 1, false)
            .unwrap()
            .unwrap()
            .mol
            .atoms
            .len(),
        1
    );

    let mut dummy = pdb_atom("HETATM", 1, " C1 ", ("DUM", 1), [0.0; 3], ("C", ""));
    dummy.replace_range(18..21, "DUM");
    assert_eq!(
        mol_from_pdb_block(&dummy, false, false, 0, false)
            .unwrap()
            .unwrap()
            .mol
            .atoms
            .len(),
        0
    );
    assert_eq!(
        mol_from_pdb_block(&dummy, false, false, 1, false)
            .unwrap()
            .unwrap()
            .mol
            .atoms
            .len(),
        1
    );
}

#[test]
fn malformed_pdb_fields_return_typed_errors() {
    assert!(
        mol_from_pdb_block("REMARK no atoms\n", true, true, 0, true)
            .unwrap()
            .is_none()
    );
    assert!(
        mol_from_pdb_block("HEADER test\n", false, false, 0, false)
            .unwrap()
            .is_some()
    );

    let mut bad_serial = pdb_atom("HETATM", 1, " C1 ", ("UNL", 1), [0.0; 3], ("C", ""));
    bad_serial.replace_range(6..11, "  1x ");
    assert!(mol_from_pdb_block(&bad_serial, false, false, 0, false).is_err());

    let mut bad_coord = pdb_atom("HETATM", 1, " C1 ", ("UNL", 1), [0.0; 3], ("C", ""));
    bad_coord.replace_range(30..38, "  1e3  ");
    assert!(mol_from_pdb_block(&bad_coord, false, false, 0, false).is_err());

    let mut bad_element = pdb_atom("HETATM", 1, " Z1 ", ("UNL", 1), [0.0; 3], ("QQ", ""));
    bad_element.replace_range(76..78, "QQ");
    assert!(mol_from_pdb_block(&bad_element, false, false, 0, false).is_err());
}

#[test]
fn proximity_bonding_respects_hydrogen_and_residue_rules() {
    let mut mol = Mol::default();
    for anum in [6, 8, 1, 1, 6] {
        mol.add_atom(Atom::new(anum));
    }
    let mut pdb = PdbMol {
        mol,
        coords: vec![
            [0.0, 0.0, 0.0],
            [1.25, 0.0, 0.0],
            [0.8, 0.0, 0.0],
            [0.8, 0.1, 0.0],
            [1.25, 0.0, 0.0],
        ],
        info: vec![
            residue(b" C  ", b"ALA", 1, false),
            residue(b" O  ", b"ALA", 1, false),
            residue(b" H1 ", b"ALA", 1, false),
            residue(b" H2 ", b"ALA", 1, false),
            residue(b" C  ", b"GLY", 2, false),
        ],
        is_3d: false,
    };
    connect_the_dots(&mut pdb);
    assert!(pdb.mol.bond_between(0, 1).is_some());
    assert!(pdb.mol.bond_between(2, 3).is_none());
    assert!(pdb.mol.degree(2) <= 1);
    assert!(pdb.mol.bond_between(1, 4).is_none());
}

#[test]
fn standard_residue_bond_orders_cover_protein_and_nucleic_acid_tables() {
    for (res_name, left, right) in [
        (b"ALA".as_slice(), b" C  ", b" O  "),
        (b"ARG".as_slice(), b" CZ ", b" NH2"),
        (b"HIS".as_slice(), b" CE1", b" ND1"),
        (b" DA".as_slice(), b" C8 ", b" N7 "),
        (b" DT".as_slice(), b" C4 ", b" O4 "),
    ] {
        let mut mol = Mol::default();
        mol.add_atom(Atom::new(6));
        mol.add_atom(Atom::new(8));
        mol.add_bond(Bond::new(0, 1, BondType::Single));
        let mut pdb = PdbMol {
            mol,
            coords: vec![[0.0; 3]; 2],
            info: vec![
                residue(left, res_name, 1, false),
                residue(right, res_name, 1, false),
            ],
            is_3d: false,
        };
        standard_residue_bond_orders(&mut pdb);
        assert_eq!(pdb.mol.bonds[0].bt, BondType::Double, "{res_name:?}");
    }

    assert!(!standard_chiral_atom(b"GLY", *b" CA "));
    assert!(standard_chiral_atom(b"ILE", *b" CB "));
    assert!(standard_chiral_atom(b"ALA", *b" CA "));
    assert!(!standard_chiral_atom(b"ALA", *b" CB "));
}

#[test]
fn three_dimensional_stereo_helpers_cover_tetrahedral_and_square_planar() {
    let mut tetra = Mol::default();
    for anum in [6, 9, 17, 35, 53] {
        tetra.add_atom(Atom::new(anum));
    }
    for nb in 1..5 {
        tetra.add_bond(Bond::new(0, nb, BondType::Single));
    }
    let coords = [
        [0.0, 0.0, 0.0],
        [1.0, 1.0, 1.0],
        [-1.0, -1.0, 1.0],
        [-1.0, 1.0, -1.0],
        [1.0, -1.0, -1.0],
    ];
    assign_chiral_types_from_3d(&mut tetra, &coords).unwrap();
    assert!(tetra.atoms[0].chiral.is_tetrahedral());

    let mut square = Mol::default();
    for anum in [78, 7, 7, 7, 7] {
        square.add_atom(Atom::new(anum));
    }
    for nb in 1..5 {
        square.add_bond(Bond::new(0, nb, BondType::Single));
    }
    let coords = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
    ];
    let tag = nontetrahedral_from_3d(&square, &coords, 0)
        .unwrap()
        .unwrap();
    assert_eq!(tag.0, ChiralTag::SquarePlanar);
    assert!((1..=3).contains(&tag.1));
}

fn coordination_mol(vectors: &[[f64; 3]]) -> (Mol, Vec<[f64; 3]>) {
    let mut mol = Mol::default();
    mol.add_atom(Atom::new(15));
    for _ in vectors {
        let nb = mol.add_atom(Atom::new(6));
        mol.add_bond(Bond::new(0, nb, BondType::Single));
    }
    let mut coords = vec![[0.0, 0.0, 0.0]];
    coords.extend_from_slice(vectors);
    (mol, coords)
}

#[test]
fn non_tetrahedral_geometry_classes_and_fail_closed_cases() {
    let cases: &[(&[[f64; 3]], ChiralTag)] = &[
        (
            &[[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            ChiralTag::SquarePlanar,
        ),
        (
            &[
                [1.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            ChiralTag::Octahedral,
        ),
        (
            &[
                [1.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, -0.5, 0.866],
            ],
            ChiralTag::TrigonalBipyramidal,
        ),
        (
            &[
                [1.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, -1.0, 0.0],
            ],
            ChiralTag::SquarePlanar,
        ),
        (
            &[
                [1.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, -1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            ChiralTag::Octahedral,
        ),
        (
            &[
                [1.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, -1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, -1.0],
            ],
            ChiralTag::Octahedral,
        ),
    ];
    for (vectors, expected) in cases {
        let (mol, coords) = coordination_mol(vectors);
        let (tag, permutation) = nontetrahedral_from_3d(&mol, &coords, 0).unwrap().unwrap();
        assert_eq!(tag, *expected, "{vectors:?}");
        assert!(permutation > 0);
    }

    let (mut carbon, coords) =
        coordination_mol(&[[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
    carbon.atoms[0].anum = 6;
    assert_eq!(nontetrahedral_from_3d(&carbon, &coords, 0).unwrap(), None);

    let (mol, mut coords) = coordination_mol(&[[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
    coords[1] = [0.0; 3];
    assert!(nontetrahedral_from_3d(&mol, &coords, 0).is_err());
}

#[test]
fn conect_repetitions_reach_triple_and_quadruple_orders() {
    for (repeat, expected) in [
        ("    2    2    2", BondType::Triple),
        ("    2    2    2    2", BondType::Quadruple),
    ] {
        let block = [
            pdb_atom("HETATM", 1, " C1 ", ("UNL", 1), [0.0; 3], ("C", "")),
            pdb_atom("HETATM", 2, " N1 ", ("UNL", 1), [1.2, 0.0, 0.0], ("N", "")),
            format!("CONECT    1{repeat}\nEND\n"),
        ]
        .concat();
        let pdb = mol_from_pdb_block(&block, false, false, 0, false)
            .unwrap()
            .unwrap();
        assert_eq!(pdb.mol.bonds[0].bt, expected);
    }
}

#[test]
fn proximity_predicate_covers_distance_and_hydrogen_boundaries() {
    let entry = |x: f32, elem: u32| ProximityEntry {
        x,
        y: 0.0,
        z: 0.0,
        r: 0.7,
        hash: 0,
        next: 0,
        elem,
    };
    assert!(!is_bonded(&entry(0.0, 1), &entry(1.0, 1)));
    assert!(!is_bonded(&entry(0.0, 6), &entry(0.1, 8)));
    assert!(is_bonded(&entry(0.0, 6), &entry(1.3, 8)));
    assert!(!is_bonded(&entry(0.0, 6), &entry(6.0, 8)));

    assert!(is_blacklisted_atom(1));
    assert!(!is_blacklisted_atom(6));
    assert!(!is_blacklisted_atom(16));
    assert!(is_blacklisted_atom(17));
}

#[test]
fn reads_rdkit_written_ethanol() {
    // Chem.MolToPDBBlock of an embedded ethanol (heavy atoms).
    let block = "HETATM    1  C1  UNL     1      -0.889  -0.214   0.011  1.00  0.00           C  \n\
HETATM    2  C2  UNL     1       0.457   0.442  -0.252  1.00  0.00           C  \n\
HETATM    3  O1  UNL     1       1.440  -0.506   0.010  1.00  0.00           O  \n\
CONECT    1    2\nCONECT    2    3\nEND\n";
    let res = crate::rdkit_mol_from_pdb_block(block, true, true, 0, true)
        .unwrap()
        .unwrap();
    assert_eq!(res.smiles, "CCO");
    assert_eq!(res.coords[2], [1.440, -0.506, 0.010]);
    assert_eq!(crate::rdkit_canonical_smiles(&res.molecule).unwrap(), "CCO");
}

#[test]
fn conect_duplicates_give_bond_orders() {
    let block = "HETATM    1  C1  UNL     1       0.000   0.000   0.000  1.00  0.00           C  \n\
HETATM    2  O1  UNL     1       1.200   0.000   0.000  1.00  0.00           O  \n\
CONECT    1    2    2\nCONECT    2    1    1\nEND\n";
    let res = crate::rdkit_mol_from_pdb_block(block, true, true, 0, true)
        .unwrap()
        .unwrap();
    assert_eq!(res.smiles, "C=O");
}

#[test]
fn xyz_block_reads_atoms_without_bonds() {
    let (mol, coords) =
        crate::rdkit_mol_from_xyz_block("2\ncomment\nC 0 0 0\nCL 1.0 2.0 3.5\n").unwrap();
    assert_eq!(mol.atom_count(), 2);
    assert_eq!(mol.bond_count(), 0);
    assert_eq!(coords[1], [1.0, 2.0, 3.5]);
    assert!(crate::rdkit_mol_from_xyz_block("1\nc\nC 1.0 2.0\n").is_err());
}

#[test]
fn numeric_symbol_geometry_and_residue_helpers_cover_boundaries() {
    assert_eq!(to_int(b"12\0bad"), Ok(12));
    assert_eq!(to_int(b"    "), Ok(0));
    assert_eq!(to_int(b"  +"), Ok(0));
    assert_eq!(to_int(b"99999999999999999999"), Ok(0));
    assert!(to_int(b"12x").is_err());

    assert_eq!(to_double(b"  -1.25\0bad"), Ok(-1.25));
    assert_eq!(to_double(b"+.5"), Ok(0.5));
    assert_eq!(to_double(b"."), Ok(0.0));
    assert!(to_double(b"1e3").is_err());

    assert_eq!(atom_from_symbol("D").unwrap(), Some((1, 2)));
    assert_eq!(atom_from_symbol("T").unwrap(), Some((1, 3)));
    assert_eq!(atom_from_symbol("*").unwrap(), None);
    assert!(atom_from_symbol("NotAnElement").is_err());

    assert!(!standard_chiral_atom(b"GLY", *b" CA "));
    assert!(standard_chiral_atom(b"ILE", *b" CB "));
    assert!(!standard_chiral_atom(b"ILE", *b" CG1"));
    assert!(standard_chiral_atom(b"ALA", *b" CA "));
    assert!(!standard_chiral_atom(b"UNK", *b" CA "));

    assert_eq!(
        angle_to([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]),
        std::f64::consts::PI
    );
    assert_eq!(angle_to([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]), 0.0);

    let entry = |x: f32, y: f32| ProximityEntry {
        x,
        y,
        z: 0.0,
        r: 0.7,
        hash: 0,
        next: 0,
        elem: 6,
    };
    assert!(!is_bonded(&entry(0.0, 0.0), &entry(0.0, 6.0)));
}
