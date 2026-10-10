//! Expected strings are RDKit 2026.03.1's
//! `Chem.MolToSmiles(Chem.MolFromSmiles(input))`.

use super::{RdkitSmilesError, rdkit_canonical_smiles};

fn rd(s: &str) -> Result<String, RdkitSmilesError> {
    let mol = crate::parse(s).expect("parses");
    rdkit_canonical_smiles(&mol)
}

fn check(cases: &[(&str, &str)]) {
    for &(input, want) in cases {
        assert_eq!(rd(input).as_deref(), Ok(want), "input {input}");
    }
}

#[test]
fn tetrahedral_centers() {
    check(&[
        ("N[C@@H](C)C(=O)O", "C[C@H](N)C(=O)O"),
        ("OC(=O)[C@@H]1CCCN1", "O=C(O)[C@@H]1CCCN1"),
        ("[C@@H](F)(Cl)Br", "F[C@H](Cl)Br"),
        ("[H][C@](F)(Cl)Br", "F[C@@H](Cl)Br"),
        ("C[C@H](N)[C@@H](C)O", "C[C@H](N)[C@@H](C)O"),
        (
            "C[C@]12CC[C@H]3[C@@H](CC=C4C[C@@H](O)CC[C@@]43C)[C@@H]1CC[C@@H]2O",
            "C[C@]12CC[C@H]3[C@@H](CC=C4C[C@@H](O)CC[C@@]43C)[C@@H]1CC[C@@H]2O",
        ),
        ("C[S@](=O)c1ccccc1", "C[S@](=O)c1ccccc1"),
        // Not a stereocentre: the tag and the bracket H go away.
        ("C[C@H](C)C", "CC(C)C"),
        ("[NH3+][C@@H](C)C(=O)[O-]", "C[C@H]([NH3+])C(=O)[O-]"),
    ]);
}

#[test]
fn ring_stereo() {
    check(&[
        ("C[C@H]1CC[C@@H](C)CC1", "C[C@H]1CC[C@@H](C)CC1"),
        ("C[C@H]1CC[C@H](C)CC1", "C[C@H]1CC[C@H](C)CC1"),
        ("O[C@H]1CC[C@@H](CC1)O", "O[C@H]1CC[C@H](O)CC1"),
    ]);
}

#[test]
fn double_bond_stereo() {
    check(&[
        ("C/C=C/C", "C/C=C/C"),
        ("C/C=C\\C", "C/C=C\\C"),
        ("C(=C/C)\\C", "C/C=C/C"),
        ("F/C(Cl)=C(\\F)Cl", "F/C(Cl)=C(\\F)Cl"),
        ("O/N=C/c1ccccc1", "O/N=C/c1ccccc1"),
        ("[H]/C(C)=C/C", "C/C=C\\C"),
        // Ring double bonds keep stereo from eight members on.
        ("C1CCCC/C=C/CCC1", "C1=C/CCCCCCCC/1"),
        ("C1CC/C=C/CC1", "C1=CCCCCC1"),
    ]);
}

#[test]
fn conjugated_double_bond_stereo() {
    check(&[
        ("F/C=C/C=C/F", "F/C=C/C=C/F"),
        ("F/C=C\\C=C/C=C\\F", "F\\C=C/C=C\\C=C/F"),
        ("C/C=C/C=C/C=C/C", "C/C=C/C=C/C=C/C"),
        ("c1ccccc1/C=C/C=C/C=O", "O=C/C=C/C=C/c1ccccc1"),
    ]);
}

#[test]
fn rings_and_ring_closure_digits() {
    check(&[
        ("C%10CCCCC%10", "C1CCCCC1"),
        ("C%(100)CC%(100)", "C1CC1"),
        ("C12C3C4C1C5C2C3C45", "C12C3C4C1C1C2C3C41"),
        ("C1CC2CCC1C2", "C1CC2CCC1C2"),
        (
            "c1cc2ccc3ccc4ccc1c1c2c3c4cc1",
            "c1cc2ccc3ccc4ccc1c1ccc2c3c41",
        ),
        (
            "C1CCC2(CC1)CCC1(CC2)CCC2(CC1)CCC1(CC2)CCC2(CC1)CCC1(CC2)CCC2(CC1)CCC1(CC2)CCC2(CC1)CCC1(CC2)CCCCC1",
            "C1CCC2(CC1)CCC1(CC2)CCC2(CC1)CCC1(CC2)CCC2(CC1)CCC1(CCC3(CCC4(CCC5(CCC6(CCCCC6)CC5)CC4)CC3)CC1)CC2",
        ),
    ]);
}

#[test]
fn ring_closures_beyond_nine() {
    // Eleven rings open at once: two-digit ring-closure labels.
    let input = "C12CC3CC4CC5CC6CC7CC8CC9CC%10CC%11CC1C%11C%10C9C8C7C6C5C4C3C2";
    let want = "C1C2CC3CC4CC5CC6CC7CC8CC9CC%10CC%11C1CC2C3C4C5C6C7C8C9C%10%11";
    assert_eq!(rd(input).as_deref(), Ok(want));
}

#[test]
fn charges_and_cleanup() {
    check(&[
        ("CN(=O)=O", "C[N+](=O)[O-]"),
        ("O=n1ccccc1", "[O-][n+]1ccccc1"),
        ("C=[N+]=[N-]", "C=[N+]=[N-]"),
        ("OCl(=O)(=O)=O", "[O-][Cl+3]([O-])([O-])O"),
        ("F[P-](F)(F)(F)(F)F", "F[P-](F)(F)(F)(F)F"),
        ("C[N+](C)(C)C", "C[N+](C)(C)C"),
    ]);
}

#[test]
fn isotopes_and_explicit_hydrogens() {
    check(&[
        ("[H]C([H])([H])C", "CC"),
        ("[2H]C([2H])([2H])O", "[2H]C([2H])([2H])O"),
        ("[13CH4]", "[13CH4]"),
        ("[13C@@H](F)(Cl)Br", "F[13C@H](Cl)Br"),
        ("[2H][C@](C)(O)c1ccccc1", "[2H][C@](C)(O)c1ccccc1"),
        ("[H][H]", "[H][H]"),
    ]);
}

#[test]
fn aromatic_and_kekule_inputs() {
    check(&[
        ("C1=CC=CC=C1", "c1ccccc1"),
        ("N1C=CC=C1", "c1cc[nH]c1"),
        ("O=C1C=CC=C[NH]1", "O=c1cccc[nH]1"),
        ("C1=CC=C2C=CC=CC2=C1", "c1ccc2ccccc2c1"),
        ("c1cc2ccc1cc2", "c1cc2ccc1cc2"),
        ("C1=CC=CC=CC=C1", "C1=CC=CC=CC=C1"),
        ("B1=CC=CC=C1", "b1ccccc1"),
        ("Cn1cccc1", "Cn1cccc1"),
    ]);
}

#[test]
fn multiple_fragments() {
    check(&[
        ("[Na+].[Cl-]", "[Cl-].[Na+]"),
        (
            "CC(=O)[O-].CC(=O)[O-].[Ca+2]",
            "CC(=O)[O-].CC(=O)[O-].[Ca+2]",
        ),
        (
            "O.O.O.[Na+].[Na+].[O-]C(=O)C(=O)[O-]",
            "O.O.O.O=C([O-])C(=O)[O-].[Na+].[Na+]",
        ),
        ("c1ccccc1.c1ccccc1.C.C.O", "C.C.O.c1ccccc1.c1ccccc1"),
        (
            "[Fe+2].c1cc[cH-]c1.c1cc[cH-]c1",
            "[Fe+2].c1cc[cH-]c1.c1cc[cH-]c1",
        ),
    ]);
}

#[test]
fn radicals_maps_dummies_and_metals() {
    check(&[
        ("[CH2]", "[CH2]"),
        ("[C]", "[C]"),
        ("[OH-]", "[OH-]"),
        ("[CH3:1][OH:2]", "[CH3:1][OH:2]"),
        ("*c1ccccc1", "*c1ccccc1"),
        ("[1*]CC[2*]", "[1*]CC[2*]"),
        ("N->[Pt](<-N)(Cl)Cl", "[NH3]->[Pt](<-[NH3])([Cl])[Cl]"),
        ("CCO[Ti](OCC)(OCC)OCC", "CC[O][Ti]([O]CC)([O]CC)[O]CC"),
        ("C=C=C", "C=C=C"),
    ]);
}

#[test]
fn rdkit_rejected_inputs_are_errors() {
    // RDKit's `MolFromSmiles` returns None for these.
    for s in ["CN(C)(C)(C)C", "c1cccc1"] {
        assert!(
            matches!(rd(s), Err(RdkitSmilesError::Sanitization(_))),
            "input {s}"
        );
    }
}

#[test]
fn nontetrahedral_stereo_is_written_like_rdkit() {
    // Chem.MolToSmiles(Chem.MolFromSmiles(s)), RDKit 2026.03.1.
    for (s, want) in [
        ("F[Pt@SP1](Cl)(Br)I", "[F][Pt@SP1]([Cl])([Br])[I]"),
        ("F[Pt@SP1](C)(O)Cl", "[CH3][Pt@SP3]([OH])([F])[Cl]"),
        ("S[As@TB1](F)(Cl)(Br)N", "N[As@TB6](F)(S)(Cl)Br"),
        (
            "C[Pt@OH1](F)(O)(N)(Br)Cl",
            "[CH3][Pt@OH16]([NH2])([OH])([F])([Cl])[Br]",
        ),
        ("CC[Pt@SP](C)(O)F", "C[CH2][Pt@SP]([CH3])([OH])[F]"),
        ("F[C@TH2H](C)O", "C[C@@H](O)F"),
    ] {
        assert_eq!(rd(s).unwrap(), want, "input {s}");
    }
}

#[test]
fn native_canonical_smiles_is_unchanged() {
    let mol = crate::parse("OC(=O)[C@@H]1CCCN1").expect("parses");
    let native = crate::canonical_smiles(&mol);
    let _ = rdkit_canonical_smiles(&mol).expect("rdkit smiles");
    assert_eq!(crate::canonical_smiles(&mol), native);
}

#[test]
fn long_chain() {
    let chain = "C".repeat(2000);
    assert_eq!(rd(&chain).as_deref(), Ok(chain.as_str()));
}

#[test]
fn writer_options_match_rdkit() {
    use super::{RdkitSmilesParams, rdkit_smiles};
    let d = RdkitSmilesParams::default();
    let cases: [(&str, RdkitSmilesParams, &str); 7] = [
        (
            "c1ccccc1O",
            RdkitSmilesParams { kekule: true, ..d },
            "OC1=CC=CC=C1",
        ),
        (
            "F/C=C/[C@H](N)C",
            RdkitSmilesParams {
                isomeric: false,
                ..d
            },
            "CC(N)C=CF",
        ),
        (
            "c1ccccc1[13CH2]O",
            RdkitSmilesParams {
                all_bonds_explicit: true,
                all_hs_explicit: true,
                ..d
            },
            "[OH]-[13CH2]-[c]1:[cH]:[cH]:[cH]:[cH]:[cH]:1",
        ),
        (
            "OCC.[Na+]",
            RdkitSmilesParams {
                canonical: false,
                ..d
            },
            "OCC.[Na+]",
        ),
        (
            "OCC(=O)N",
            RdkitSmilesParams {
                rooted_at_atom: Some(3),
                ..d
            },
            "O=C(N)CO",
        ),
        (
            "F/C=C/[C@H](N)C",
            RdkitSmilesParams {
                isomeric: false,
                all_bonds_explicit: true,
                ..d
            },
            "C-C(-N)/C=C/F",
        ),
        (
            "c1ccc[nH]1",
            RdkitSmilesParams {
                kekule: true,
                all_hs_explicit: true,
                ..d
            },
            "[CH]1=[CH][NH][CH]=[CH]1",
        ),
    ];
    for (input, params, want) in cases {
        let mol = crate::parse(input).expect("parses");
        assert_eq!(
            rdkit_smiles(&mol, &params).as_deref(),
            Ok(want),
            "{input} {params:?}"
        );
    }
}

/// Every implicit H as a new graph atom, heavy atom by heavy atom (what
/// chematic-chem's `add_hydrogens` and RDKit's `AddHs` do).
fn add_hs(mol: &chematic_core::Molecule) -> chematic_core::Molecule {
    use chematic_core::{Atom, AtomIdx, BondIdx, BondOrder, Element, MoleculeBuilder};
    let mut b = MoleculeBuilder::new();
    for i in 0..mol.atom_count() {
        let mut a = mol.atom(AtomIdx(i as u32)).clone();
        a.hydrogen_count = Some(0);
        b.add_atom(a);
    }
    for i in 0..mol.bond_count() {
        let bond = mol.bond(BondIdx(i as u32));
        b.add_bond(bond.atom1, bond.atom2, bond.order)
            .expect("bond");
    }
    b.copy_stereo_from(mol);
    b.copy_bond_directions_from(mol);
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        for _ in 0..chematic_core::implicit_hcount(mol, idx) {
            let h = b.add_atom(Atom::new(Element::H));
            b.add_bond(idx, h, BondOrder::Single).expect("bond");
            if let Some(order) = mol.stereo_neighbor_order(idx) {
                let o = order
                    .iter()
                    .map(|&v| {
                        if v == chematic_core::STEREO_H_SENTINEL {
                            h.0
                        } else {
                            v
                        }
                    })
                    .collect();
                b.set_stereo_neighbor_order(idx, o);
            }
        }
    }
    b.build()
}

#[test]
fn explicit_hydrogen_atoms_are_written_like_rdkit_add_hs() {
    // Chem.MolToSmiles(Chem.AddHs(Chem.MolFromSmiles(input)))
    for (input, want) in [
        (
            "N[C@@H](C)C(=O)O",
            "[H]OC(=O)[C@@]([H])(N([H])[H])C([H])([H])[H]",
        ),
        ("c1cc[nH]c1", "[H]c1c([H])c([H])n([H])c1[H]"),
        ("F/C=C/C", "[H]/C(F)=C(/[H])C([H])([H])[H]"),
        ("[2H]C", "[H]C([H])([H])[2H]"),
    ] {
        let mol = add_hs(&crate::parse(input).expect("parses"));
        assert_eq!(rdkit_canonical_smiles(&mol).as_deref(), Ok(want), "{input}");
    }
}

/// RDKit 2026.03.1 `MolToSmarts`, `MolToCXSmarts`,
/// `MurckoScaffold.GetScaffoldForMol`, `GetStereoisomerCount`,
/// `FindMolChiralCenters(includeUnassigned=True)` and `MolToPDBBlock`.
#[test]
fn rdkit_writers_match_rdkit() {
    let p = |s: &str| crate::parse(s).expect("parses");
    let cases: &[(&str, &str, &str, &str)] = &[
        (
            "C[C@H](N)C(=O)O",
            "[#6]-[#6@H](-[#7])-[#6](=[#8])-[#8]",
            "[#6]-[#6@H](-[#7])-[#6](=[#8])-[#8]",
            "",
        ),
        ("N->[Cu]", "[#7]->[Cu]", "[#7]-[Cu] |C:0.0|", ""),
        (
            "CCCCOC1=CC=C(NC[S](=O)=O)C=N1",
            "[#6]-[#6]-[#6]-[#6]-[#8]-[#6]1:[#6]:[#6]:[#6](-[#7]-[#6]-[#16](=[#8])=[#8]):[#6]:[#7]:1",
            "[#6]-[#6]-[#6]-[#6]-[#8]-[#6]1:[#6]:[#6]:[#6](-[#7]-[#6]-[#16](=[#8])=[#8]):[#6]:[#7]:1 |^1:11|",
            "c1ccncc1",
        ),
        (
            "Cn1cccc1CC1CCC1",
            "[#6]-[#7]1:[#6]:[#6]:[#6]:[#6]:1-[#6]-[#6]1-[#6]-[#6]-[#6]-1",
            "[#6]-[#7]1:[#6]:[#6]:[#6]:[#6]:1-[#6]-[#6]1-[#6]-[#6]-[#6]-1",
            "c1c[nH]c(CC2CCC2)c1",
        ),
        (
            "O=C1CC[C@H](C)N1CCc1ccccc1",
            "[#8]=[#6]1-[#6]-[#6]-[#6@H](-[#6])-[#7]-1-[#6]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1",
            "[#8]=[#6]1-[#6]-[#6]-[#6@H](-[#6])-[#7]-1-[#6]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1",
            "O=C1CCCN1CCc1ccccc1",
        ),
    ];
    for &(input, smarts, cx, scaffold) in cases {
        let m = p(input);
        assert_eq!(
            super::rdkit_smarts(&m, true, None).unwrap(),
            smarts,
            "{input}"
        );
        assert_eq!(super::rdkit_cx_smarts(&m).unwrap(), cx, "{input}");
        assert_eq!(
            super::rdkit_murcko_scaffold(&m).unwrap(),
            scaffold,
            "{input}"
        );
    }
    assert_eq!(
        super::rdkit_chiral_centers(&p("O=C1CC[C@H](C)N1CCc1ccccc1"), true).unwrap(),
        [(4, "S".to_string())]
    );
    assert_eq!(
        super::rdkit_stereoisomer_count(&p("CC(F)C(Cl)Br")).unwrap(),
        4
    );
    assert_eq!(
        super::rdkit_pdb_block(&p("CC(=O)[O-]"), None).unwrap(),
        "HETATM    1  C1  UNL     1       0.000   0.000   0.000  1.00  0.00           C  \n\
         HETATM    2  C2  UNL     1       0.000   0.000   0.000  1.00  0.00           C  \n\
         HETATM    3  O1  UNL     1       0.000   0.000   0.000  1.00  0.00           O  \n\
         HETATM    4  O2  UNL     1       0.000   0.000   0.000  1.00  0.00           O1-\n\
         CONECT    1    2\nCONECT    2    3    3    4\nEND\n"
    );
}

/// Keep the Rust surface for query and reaction SMARTS aligned with the
/// Python binding fixtures. These paths are substantial parsers/writers in
/// their own right and must not rely on binding-only coverage.
#[test]
fn rdkit_query_and_reaction_smarts_match_rdkit() {
    for (input, want) in [
        ("[N;H2,H1;!$(NC=O)]", "[N;H2,H1;!$(NC=O)]"),
        ("[#6;$([C;H3,H2])]-[O,N]", "[#6&$([C;H3,H2])]-[O,N]"),
    ] {
        assert_eq!(super::rdkit_smarts_to_smarts(input).as_deref(), Ok(want));
    }

    for (input, want) in [
        (
            "[C@@H:1]([NH2:2])([#6:3])[C:4]=[O:5]>>[C@@H:1]([NH:2]C(C)=O)([#6:3])[C:4]=[O:5]",
            "[N&H2:2][C@&H1:1]([#6:3])[C:4]=[O:5]>>[N&H1:2]([C@&H1:1]([#6:3])[C:4]=[O:5])C(C)=O",
        ),
        (
            "([C:1](=[O:2])[OH].[NH2:3])>>[C:1](=[O:2])[N:3]",
            "([C:1](=[O:2])[O&H1].[N&H2:3])>>[C:1](=[O:2])[N:3]",
        ),
    ] {
        assert_eq!(
            super::rdkit_reaction_to_smarts(input).as_deref(),
            Ok(want),
            "{input}"
        );
    }

    for invalid in ["([C:1].[O:2]>>[C:1][O:2]", "[C:1]>[O:2]"] {
        assert!(
            super::rdkit_reaction_to_smarts(invalid).is_err(),
            "{invalid}"
        );
    }
}

/// Exercise the pure-Rust half of `MolFromInchi` without requiring the
/// optional native IUPAC library. The adapter consumes exactly the plain data
/// returned by `GetStructFromINCHI`.
#[test]
fn rdkit_inchi_output_adapter_builds_and_validates_molecules() {
    use super::{InchiOutputAtom, rdkit_molecule_from_inchi_output};

    let atoms = vec![
        InchiOutputAtom {
            element: "C".into(),
            bonds: vec![(1, 1, 0)],
            num_iso_h: [3, 0, 0, 0],
            ..Default::default()
        },
        InchiOutputAtom {
            element: "C".into(),
            bonds: vec![(0, 1, 0), (2, 1, 0)],
            num_iso_h: [2, 0, 0, 0],
            ..Default::default()
        },
        InchiOutputAtom {
            element: "O".into(),
            bonds: vec![(1, 1, 0)],
            num_iso_h: [1, 0, 0, 0],
            ..Default::default()
        },
    ];
    let mol = rdkit_molecule_from_inchi_output(&atoms, &[]).expect("ethanol output");
    assert_eq!(super::rdkit_canonical_smiles(&mol).as_deref(), Ok("CCO"));

    let unknown = [InchiOutputAtom {
        element: "Xx".into(),
        ..Default::default()
    }];
    assert!(rdkit_molecule_from_inchi_output(&unknown, &[]).is_err());

    let bad_bond = [
        InchiOutputAtom {
            element: "C".into(),
            bonds: vec![(1, 9, 0)],
            ..Default::default()
        },
        InchiOutputAtom {
            element: "C".into(),
            ..Default::default()
        },
    ];
    assert!(rdkit_molecule_from_inchi_output(&bad_bond, &[]).is_err());
}

/// `GetStereoisomerCount` cases that need `FindPotentialStereo`'s
/// dependent ("possible") stereo (RDKit 2026.03.1).
#[test]
fn stereoisomer_count_uses_find_potential_stereo() {
    for (smiles, want) in [
        ("CC(C(=O)O)=C1CCC(C)CC1", 4),
        ("O[As]=O", 2),
        ("ON=C1C=CC(C=C1)=NO", 4),
        ("C1C[S+]2CC[S+]1CC2", 4),
        ("C12C3=C4C5=C1[Fe]23456789C%10C6=C7C8=C9%10", 1),
        ("BrCC(Br)COP(=O)(OCC(Br)CBr)OCC(Br)CBr", 16),
    ] {
        let mol = crate::parse(smiles).expect("parses");
        assert_eq!(
            super::rdkit_stereoisomer_count(&mol).unwrap(),
            want,
            "{smiles}"
        );
    }
}

/// RDKit 2026.03.1 `MolToMolBlock` switches to V3000 for a dative bond.
#[test]
fn mol_block_2d_writes_v3000_for_dative_bonds() {
    let mol = crate::parse("C[13CH2]O->[Fe]").expect("parses");
    let block = super::rdkit_mol_block_2d(&mol).unwrap();
    let want = "\n     RDKit          2D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n\
                M  V30 BEGIN CTAB\nM  V30 COUNTS 4 3 0 0 0\nM  V30 BEGIN ATOM\n\
                M  V30 1 C -1.979613 -0.136500 0.000000 0\n\
                M  V30 2 C -0.599379 0.450827 0.000000 0 MASS=13\n\
                M  V30 3 O 0.599379 -0.450827 0.000000 0\n\
                M  V30 4 Fe 1.979613 0.136500 0.000000 0 VAL=1\n\
                M  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 2 1 2 3\n\
                M  V30 3 9 3 4\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n";
    assert_eq!(block, want);
}

#[test]
fn allene_chirality_is_read_and_dropped() {
    // Chem.MolToSmiles(Chem.MolFromSmiles(s)), RDKit 2026.03.1.
    for (s, want) in [
        ("OC=[C@AL1]=CC", "CC=C=CO"),
        ("OC=[C@AL2]=CC", "CC=C=CO"),
        ("OC=[C@AL]=CC", "CC=C=CO"),
        ("C[C@AL1](F)Cl", "C[C](F)Cl"),
        ("[C@AL1H2]", "[CH2]"),
    ] {
        assert_eq!(rd(s).unwrap(), want, "input {s}");
    }
    for s in ["OC=[C@AL0]=CC", "[C@AL3]", "C[C@@AL1H](F)Cl", "[CH2@AL1]"] {
        assert!(crate::parse(s).is_err(), "input {s}");
    }
}

#[test]
fn mol_hash_matches_rdkit() {
    use super::RdkitHashFunction;
    let cases: &[(&str, &str, &str, Option<&str>)] = &[
        // rdMolHash.MolHash(Chem.MolFromSmiles(s), f[, True]), RDKit 2026.03.1.
        (
            r"C[C@H](N)C(=O)O",
            "AnonymousGraph",
            r"**(*)*(*)*",
            Some(r"**(*)*(*)*"),
        ),
        (
            r"C[C@H](N)C(=O)O",
            "ElementGraph",
            r"C[C@H](N)C(O)O",
            Some(r"C[C@H](N)C(O)O"),
        ),
        (
            r"C[C@H](N)C(=O)O",
            "CanonicalSmiles",
            r"C[C@H](N)C(=O)O",
            Some(r"C[C@H](N)C(=O)O"),
        ),
        (r"C[C@H](N)C(=O)O", "MurckoScaffold", r"", None),
        (r"C[C@H](N)C(=O)O", "ExtendedMurcko", r"", None),
        (
            r"C[C@H](N)C(=O)O",
            "MolFormula",
            r"C3H7NO2",
            Some(r"C3H7NO2"),
        ),
        (r"C[C@H](N)C(=O)O", "AtomBondCounts", r"6,5", Some(r"6,5")),
        (
            r"C[C@H](N)C(=O)O",
            "DegreeVector",
            r"0,2,0,4",
            Some(r"0,2,0,4"),
        ),
        (
            r"C[C@H](N)C(=O)O",
            "Mesomer",
            r"C[C@H](N)[C]([O])O_0",
            Some(r"C[C@H](N)[C]([O])O_0"),
        ),
        (
            r"C[C@H](N)C(=O)O",
            "HetAtomTautomer",
            r"C[C@H]([N])[C]([O])[O]_3_0",
            Some(r"C[C@H]([N])[C]([O])[O]_3_0"),
        ),
        (
            r"C[C@H](N)C(=O)O",
            "HetAtomProtomer",
            r"C[C@H]([N])[C]([O])[O]_3",
            Some(r"C[C@H]([N])[C]([O])[O]_3"),
        ),
        (
            r"C[C@H](N)C(=O)O",
            "RedoxPair",
            r"C[C@H](N)[C]([O])O",
            Some(r"C[C@H](N)[C]([O])O"),
        ),
        (
            r"C[C@H](N)C(=O)O",
            "Regioisomer",
            r"*N.CCC(=O)O",
            Some(r"*N.CCC(=O)O"),
        ),
        (r"C[C@H](N)C(=O)O", "NetCharge", r"0", Some(r"0")),
        (
            r"C[C@H](N)C(=O)O",
            "SmallWorldIndexBR",
            r"B5R0",
            Some(r"B5R0"),
        ),
        (
            r"C[C@H](N)C(=O)O",
            "SmallWorldIndexBRL",
            r"B5R0L0",
            Some(r"B5R0L0"),
        ),
        (
            r"C[C@H](N)C(=O)O",
            "ArthorSubstructureOrder",
            r"000600050100030003000029000000",
            Some(r"000600050100030003000029000000"),
        ),
        (
            r"C[C@H](N)C(=O)O",
            "HetAtomTautomerv2",
            r"[CH3]-[C@H](-[NH2])-[C](:[O]):[O]_1_0",
            Some(r"[CH3]-[C@H](-[NH2])-[C](:[O]):[O]_1_0"),
        ),
        (
            r"C[C@H](N)C(=O)O",
            "HetAtomProtomerv2",
            r"[CH3]-[C@H](-[NH2])-[C](:[O]):[O]_1",
            Some(r"[CH3]-[C@H](-[NH2])-[C](:[O]):[O]_1"),
        ),
        (
            r"Oc1ccccn1",
            "AnonymousGraph",
            r"**1*****1",
            Some(r"**1*****1"),
        ),
        (
            r"Oc1ccccn1",
            "ElementGraph",
            r"OC1CCCCN1",
            Some(r"OC1CCCCN1"),
        ),
        (
            r"Oc1ccccn1",
            "CanonicalSmiles",
            r"Oc1ccccn1",
            Some(r"Oc1ccccn1"),
        ),
        (
            r"Oc1ccccn1",
            "MurckoScaffold",
            r"c1ccncc1",
            Some(r"c1ccncc1"),
        ),
        (
            r"Oc1ccccn1",
            "ExtendedMurcko",
            r"*c1ccccn1",
            Some(r"*c1ccccn1"),
        ),
        (r"Oc1ccccn1", "MolFormula", r"C5H5NO", Some(r"C5H5NO")),
        (r"Oc1ccccn1", "AtomBondCounts", r"7,7", Some(r"7,7")),
        (r"Oc1ccccn1", "DegreeVector", r"0,1,5,1", Some(r"0,1,5,1")),
        (
            r"Oc1ccccn1",
            "Mesomer",
            r"O[C]1[CH][CH][CH][CH][N]1_0",
            Some(r"O[C]1[CH][CH][CH][CH][N]1_0"),
        ),
        (
            r"Oc1ccccn1",
            "HetAtomTautomer",
            r"[O][C]1[CH][CH][CH][CH][N]1_1_0",
            Some(r"[O][C]1[CH][CH][CH][CH][N]1_1_0"),
        ),
        (
            r"Oc1ccccn1",
            "HetAtomProtomer",
            r"[O][C]1[CH][CH][CH][CH][N]1_1",
            Some(r"[O][C]1[CH][CH][CH][CH][N]1_1"),
        ),
        (
            r"Oc1ccccn1",
            "RedoxPair",
            r"O[C]1[CH][CH][CH][CH][N]1",
            Some(r"O[C]1[CH][CH][CH][CH][N]1"),
        ),
        (
            r"Oc1ccccn1",
            "Regioisomer",
            r"*O.c1ccncc1",
            Some(r"*O.c1ccncc1"),
        ),
        (r"Oc1ccccn1", "NetCharge", r"0", Some(r"0")),
        (r"Oc1ccccn1", "SmallWorldIndexBR", r"B7R1", Some(r"B7R1")),
        (
            r"Oc1ccccn1",
            "SmallWorldIndexBRL",
            r"B7R1L5",
            Some(r"B7R1L5"),
        ),
        (
            r"Oc1ccccn1",
            "ArthorSubstructureOrder",
            r"00070007010005000200002d000000",
            Some(r"00070007010005000200002d000000"),
        ),
        (
            r"Oc1ccccn1",
            "HetAtomTautomerv2",
            r"[O]:[C]1:[C]:[C]:[C]:[C]:[N]:1_5_0",
            Some(r"[O]:[C]1:[C]:[C]:[C]:[C]:[N]:1_5_0"),
        ),
        (
            r"Oc1ccccn1",
            "HetAtomProtomerv2",
            r"[O]:[C]1:[C]:[C]:[C]:[C]:[N]:1_5",
            Some(r"[O]:[C]1:[C]:[C]:[C]:[C]:[N]:1_5"),
        ),
        (
            r"CC(=O)/C=C(\O)C",
            "AnonymousGraph",
            r"**(*)**(*)*",
            Some(r"**(*)**(*)*"),
        ),
        (
            r"CC(=O)/C=C(\O)C",
            "ElementGraph",
            r"CC(O)CC(C)O",
            Some(r"CC(O)CC(C)O"),
        ),
        (
            r"CC(=O)/C=C(\O)C",
            "CanonicalSmiles",
            r"CC(=O)/C=C(/C)O",
            Some(r"CC(=O)/C=C(/C)O"),
        ),
        (r"CC(=O)/C=C(\O)C", "MurckoScaffold", r"", None),
        (r"CC(=O)/C=C(\O)C", "ExtendedMurcko", r"", None),
        (r"CC(=O)/C=C(\O)C", "MolFormula", r"C5H8O2", Some(r"C5H8O2")),
        (r"CC(=O)/C=C(\O)C", "AtomBondCounts", r"7,6", Some(r"7,6")),
        (
            r"CC(=O)/C=C(\O)C",
            "DegreeVector",
            r"0,2,1,4",
            Some(r"0,2,1,4"),
        ),
        (
            r"CC(=O)/C=C(\O)C",
            "Mesomer",
            r"C[C]([O])[CH][C](C)O_0",
            Some(r"C[C]([O])[CH][C](C)O_0"),
        ),
        (
            r"CC(=O)/C=C(\O)C",
            "HetAtomTautomer",
            r"C[C]([O])[CH][C](C)[O]_1_0",
            Some(r"C[C]([O])[CH][C](C)[O]_1_0"),
        ),
        (
            r"CC(=O)/C=C(\O)C",
            "HetAtomProtomer",
            r"C[C]([O])[CH][C](C)[O]_1",
            Some(r"C[C]([O])[CH][C](C)[O]_1"),
        ),
        (
            r"CC(=O)/C=C(\O)C",
            "RedoxPair",
            r"C[C]([O])[CH][C](C)O",
            Some(r"C[C]([O])[CH][C](C)O"),
        ),
        (
            r"CC(=O)/C=C(\O)C",
            "Regioisomer",
            r"CC(=O)/C=C(/C)O",
            Some(r"CC(=O)/C=C(/C)O"),
        ),
        (r"CC(=O)/C=C(\O)C", "NetCharge", r"0", Some(r"0")),
        (
            r"CC(=O)/C=C(\O)C",
            "SmallWorldIndexBR",
            r"B6R0",
            Some(r"B6R0"),
        ),
        (
            r"CC(=O)/C=C(\O)C",
            "SmallWorldIndexBRL",
            r"B6R0L1",
            Some(r"B6R0L1"),
        ),
        (
            r"CC(=O)/C=C(\O)C",
            "ArthorSubstructureOrder",
            r"00070006010005000200002e000000",
            Some(r"00070006010005000200002e000000"),
        ),
        (
            r"CC(=O)/C=C(\O)C",
            "HetAtomTautomerv2",
            r"[C]:[C](:[O]):[C]:[C](:[C]):[O]_8_0",
            Some(r"[C]:[C](:[O]):[C]:[C](:[C]):[O]_8_0"),
        ),
        (
            r"CC(=O)/C=C(\O)C",
            "HetAtomProtomerv2",
            r"[C]:[C](:[O]):[C]:[C](:[C]):[O]_8",
            Some(r"[C]:[C](:[O]):[C]:[C](:[C]):[O]_8"),
        ),
        (
            r"[NH3+]CC([O-])=O.[Na+]",
            "AnonymousGraph",
            r"*.***(*)*",
            Some(r"*.***(*)*"),
        ),
        (
            r"[NH3+]CC([O-])=O.[Na+]",
            "ElementGraph",
            r"NCC(O)O.[Na]",
            Some(r"NCC(O)O.[Na]"),
        ),
        (
            r"[NH3+]CC([O-])=O.[Na+]",
            "CanonicalSmiles",
            r"[NH3+]CC(=O)[O-].[Na+]",
            Some(r"[NH3+]CC(=O)[O-].[Na+]"),
        ),
        (r"[NH3+]CC([O-])=O.[Na+]", "MurckoScaffold", r"", None),
        (r"[NH3+]CC([O-])=O.[Na+]", "ExtendedMurcko", r"", None),
        (
            r"[NH3+]CC([O-])=O.[Na+]",
            "MolFormula",
            r"C2H5NNaO2+",
            Some(r"C2H5NNaO2+"),
        ),
        (
            r"[NH3+]CC([O-])=O.[Na+]",
            "AtomBondCounts",
            r"6,4",
            Some(r"6,4"),
        ),
        (
            r"[NH3+]CC([O-])=O.[Na+]",
            "DegreeVector",
            r"0,1,1,3",
            Some(r"0,1,1,3"),
        ),
        (
            r"[NH3+]CC([O-])=O.[Na+]",
            "Mesomer",
            r"[NH3]C[C]([O])[O].[Na]_1",
            Some(r"[NH3]C[C]([O])[O].[Na]_1"),
        ),
        (
            r"[NH3+]CC([O-])=O.[Na+]",
            "HetAtomTautomer",
            r"[N]C[C]([O])[O].[Na]_3_1",
            Some(r"[N]C[C]([O])[O].[Na]_3_1"),
        ),
        (
            r"[NH3+]CC([O-])=O.[Na+]",
            "HetAtomProtomer",
            r"[N]C[C]([O])[O].[Na]_2",
            Some(r"[N]C[C]([O])[O].[Na]_2"),
        ),
        (
            r"[NH3+]CC([O-])=O.[Na+]",
            "RedoxPair",
            r"[NH3]C[C]([O])[O].[Na]",
            Some(r"[NH3]C[C]([O])[O].[Na]"),
        ),
        (
            r"[NH3+]CC([O-])=O.[Na+]",
            "Regioisomer",
            r"*[NH3+].CC(=O)[O-].[Na+]",
            Some(r"*[NH3+].CC(=O)[O-].[Na+]"),
        ),
        (r"[NH3+]CC([O-])=O.[Na+]", "NetCharge", r"1", Some(r"1")),
        (
            r"[NH3+]CC([O-])=O.[Na+]",
            "SmallWorldIndexBR",
            r"B4R4294967295",
            Some(r"B4R4294967295"),
        ),
        (
            r"[NH3+]CC([O-])=O.[Na+]",
            "SmallWorldIndexBRL",
            r"B4R4294967295L1",
            Some(r"B4R4294967295L1"),
        ),
        (
            r"[NH3+]CC([O-])=O.[Na+]",
            "ArthorSubstructureOrder",
            r"00060004020002000300002e010300",
            Some(r"00060004020002000300002e010300"),
        ),
        (
            r"[NH3+]CC([O-])=O.[Na+]",
            "HetAtomTautomerv2",
            r"[NH3+]-[CH2]-[C](=[O])-[O-].[Na+]_0_0",
            Some(r"[NH3+]-[CH2]-[C](=[O])-[O-].[Na+]_0_0"),
        ),
        (
            r"[NH3+]CC([O-])=O.[Na+]",
            "HetAtomProtomerv2",
            r"[NH3+]-[CH2]-[C](=[O])-[O-].[Na+]_0",
            Some(r"[NH3+]-[CH2]-[C](=[O])-[O-].[Na+]_0"),
        ),
        (
            r"C1CC[C@H]2CCCC[C@@H]2C1C/C=C/c1ccc(F)cc1",
            "AnonymousGraph",
            r"**1***(****2****3*****32)**1",
            Some(r"**1***(****2****3*****32)**1"),
        ),
        (
            r"C1CC[C@H]2CCCC[C@@H]2C1C/C=C/c1ccc(F)cc1",
            "ElementGraph",
            r"FC1CCC(CCCC2CCC[C@H]3CCCC[C@H]23)CC1",
            Some(r"FC1CCC(CCCC2CCC[C@H]3CCCC[C@H]23)CC1"),
        ),
        (
            r"C1CC[C@H]2CCCC[C@@H]2C1C/C=C/c1ccc(F)cc1",
            "CanonicalSmiles",
            r"Fc1ccc(/C=C/CC2CCC[C@H]3CCCC[C@H]23)cc1",
            Some(r"Fc1ccc(/C=C/CC2CCC[C@H]3CCCC[C@H]23)cc1"),
        ),
        (
            r"C1CC[C@H]2CCCC[C@@H]2C1C/C=C/c1ccc(F)cc1",
            "MurckoScaffold",
            r"C(=C\c1ccccc1)/CC1CCC[C@H]2CCCC[C@H]12",
            Some(r"C(=C\c1ccccc1)/CC1CCC[C@H]2CCCC[C@H]12"),
        ),
        (
            r"C1CC[C@H]2CCCC[C@@H]2C1C/C=C/c1ccc(F)cc1",
            "ExtendedMurcko",
            r"*c1ccc(/C=C/CC2CCC[C@H]3CCCC[C@H]23)cc1",
            Some(r"*c1ccc(/C=C/CC2CCC[C@H]3CCCC[C@H]23)cc1"),
        ),
        (
            r"C1CC[C@H]2CCCC[C@@H]2C1C/C=C/c1ccc(F)cc1",
            "MolFormula",
            r"C19H25F",
            Some(r"C19H25F"),
        ),
        (
            r"C1CC[C@H]2CCCC[C@@H]2C1C/C=C/c1ccc(F)cc1",
            "AtomBondCounts",
            r"20,22",
            Some(r"20,22"),
        ),
        (
            r"C1CC[C@H]2CCCC[C@@H]2C1C/C=C/c1ccc(F)cc1",
            "DegreeVector",
            r"0,5,14,1",
            Some(r"0,5,14,1"),
        ),
        (
            r"C1CC[C@H]2CCCC[C@@H]2C1C/C=C/c1ccc(F)cc1",
            "Mesomer",
            r"F[C]1[CH][CH][C]([CH][CH]CC2CCC[C@H]3CCCC[C@H]23)[CH][CH]1_0",
            Some(r"F[C]1[CH][CH][C]([CH][CH]CC2CCC[C@H]3CCCC[C@H]23)[CH][CH]1_0"),
        ),
        (
            r"C1CC[C@H]2CCCC[C@@H]2C1C/C=C/c1ccc(F)cc1",
            "HetAtomTautomer",
            r"F[C]1[CH][CH][C]([CH][CH]CC2CCC[C@H]3CCCC[C@H]23)[CH][CH]1_0_0",
            Some(r"F[C]1[CH][CH][C]([CH][CH]CC2CCC[C@H]3CCCC[C@H]23)[CH][CH]1_0_0"),
        ),
        (
            r"C1CC[C@H]2CCCC[C@@H]2C1C/C=C/c1ccc(F)cc1",
            "HetAtomProtomer",
            r"F[C]1[CH][CH][C]([CH][CH]CC2CCC[C@H]3CCCC[C@H]23)[CH][CH]1_0",
            Some(r"F[C]1[CH][CH][C]([CH][CH]CC2CCC[C@H]3CCCC[C@H]23)[CH][CH]1_0"),
        ),
        (
            r"C1CC[C@H]2CCCC[C@@H]2C1C/C=C/c1ccc(F)cc1",
            "RedoxPair",
            r"F[C]1[CH][CH][C]([CH][CH]CC2CCC[C@H]3CCCC[C@H]23)[CH][CH]1",
            Some(r"F[C]1[CH][CH][C]([CH][CH]CC2CCC[C@H]3CCCC[C@H]23)[CH][CH]1"),
        ),
        (
            r"C1CC[C@H]2CCCC[C@@H]2C1C/C=C/c1ccc(F)cc1",
            "Regioisomer",
            r"*C=CC*.*F.C1CC[C@H]2CCCC[C@@H]2C1.c1ccccc1",
            Some(r"*C=CC*.*F.C1CC[C@H]2CCCC[C@@H]2C1.c1ccccc1"),
        ),
        (
            r"C1CC[C@H]2CCCC[C@@H]2C1C/C=C/c1ccc(F)cc1",
            "NetCharge",
            r"0",
            Some(r"0"),
        ),
        (
            r"C1CC[C@H]2CCCC[C@@H]2C1C/C=C/c1ccc(F)cc1",
            "SmallWorldIndexBR",
            r"B22R3",
            Some(r"B22R3"),
        ),
        (
            r"C1CC[C@H]2CCCC[C@@H]2C1C/C=C/c1ccc(F)cc1",
            "SmallWorldIndexBRL",
            r"B22R3L14",
            Some(r"B22R3L14"),
        ),
        (
            r"C1CC[C@H]2CCCC[C@@H]2C1C/C=C/c1ccc(F)cc1",
            "ArthorSubstructureOrder",
            r"00140016010013000100007b000000",
            Some(r"00140016010013000100007b000000"),
        ),
        (
            r"C1CC[C@H]2CCCC[C@@H]2C1C/C=C/c1ccc(F)cc1",
            "HetAtomTautomerv2",
            r"[F]-[c]1:[cH]:[cH]:[c](/[CH]=[CH]/[CH2]-[CH]2-[CH2]-[CH2]-[CH2]-[C@H]3-[CH2]-[CH2]-[CH2]-[CH2]-[C@H]-2-3):[cH]:[cH]:1_0_0",
            Some(
                r"[F]-[c]1:[cH]:[cH]:[c](/[CH]=[CH]/[CH2]-[CH]2-[CH2]-[CH2]-[CH2]-[C@H]3-[CH2]-[CH2]-[CH2]-[CH2]-[C@H]-2-3):[cH]:[cH]:1_0_0",
            ),
        ),
        (
            r"C1CC[C@H]2CCCC[C@@H]2C1C/C=C/c1ccc(F)cc1",
            "HetAtomProtomerv2",
            r"[F]-[c]1:[cH]:[cH]:[c](/[CH]=[CH]/[CH2]-[CH]2-[CH2]-[CH2]-[CH2]-[C@H]3-[CH2]-[CH2]-[CH2]-[CH2]-[C@H]-2-3):[cH]:[cH]:1_0",
            Some(
                r"[F]-[c]1:[cH]:[cH]:[c](/[CH]=[CH]/[CH2]-[CH]2-[CH2]-[CH2]-[CH2]-[C@H]3-[CH2]-[CH2]-[CH2]-[CH2]-[C@H]-2-3):[cH]:[cH]:1_0",
            ),
        ),
        (r"[CH2]CC", "AnonymousGraph", r"[*]**", Some(r"[*]**")),
        (r"[CH2]CC", "ElementGraph", r"CCC", Some(r"CCC")),
        (
            r"[CH2]CC",
            "CanonicalSmiles",
            r"[CH2]CC",
            Some(r"[CH2]CC |^1:0|"),
        ),
        (r"[CH2]CC", "MurckoScaffold", r"", None),
        (r"[CH2]CC", "ExtendedMurcko", r"", None),
        (r"[CH2]CC", "MolFormula", r"C3H7", Some(r"C3H7")),
        (r"[CH2]CC", "AtomBondCounts", r"3,2", Some(r"3,2")),
        (r"[CH2]CC", "DegreeVector", r"0,0,1,2", Some(r"0,0,1,2")),
        (r"[CH2]CC", "Mesomer", r"[CH2]CC_0", Some(r"[CH2]CC_0")),
        (
            r"[CH2]CC",
            "HetAtomTautomer",
            r"[CH2]CC_0_0",
            Some(r"[CH2]CC_0_0"),
        ),
        (
            r"[CH2]CC",
            "HetAtomProtomer",
            r"[CH2]CC_0",
            Some(r"[CH2]CC_0"),
        ),
        (r"[CH2]CC", "RedoxPair", r"[CH2]CC", Some(r"[CH2]CC")),
        (
            r"[CH2]CC",
            "Regioisomer",
            r"[CH2]CC",
            Some(r"[CH2]CC |^1:0|"),
        ),
        (r"[CH2]CC", "NetCharge", r"0", Some(r"0")),
        (r"[CH2]CC", "SmallWorldIndexBR", r"B2R0", Some(r"B2R0")),
        (r"[CH2]CC", "SmallWorldIndexBRL", r"B2R0L1", Some(r"B2R0L1")),
        (
            r"[CH2]CC",
            "ArthorSubstructureOrder",
            r"000300020100030000000012010000",
            Some(r"000300020100030000000012010000"),
        ),
        (
            r"[CH2]CC",
            "HetAtomTautomerv2",
            r"[CH2]-[CH2]-[CH3]_0_0",
            Some(r"[CH2]-[CH2]-[CH3]_0_0"),
        ),
        (
            r"[CH2]CC",
            "HetAtomProtomerv2",
            r"[CH2]-[CH2]-[CH3]_0",
            Some(r"[CH2]-[CH2]-[CH3]_0"),
        ),
        (
            r"CCOC(=O)c1ccc2[nH]ccc2c1",
            "AnonymousGraph",
            r"****(*)*1***2****2*1",
            Some(r"****(*)*1***2****2*1"),
        ),
        (
            r"CCOC(=O)c1ccc2[nH]ccc2c1",
            "ElementGraph",
            r"CCOC(O)C1CCC2NCCC2C1",
            Some(r"CCOC(O)C1CCC2NCCC2C1"),
        ),
        (
            r"CCOC(=O)c1ccc2[nH]ccc2c1",
            "CanonicalSmiles",
            r"CCOC(=O)c1ccc2[nH]ccc2c1",
            Some(r"CCOC(=O)c1ccc2[nH]ccc2c1"),
        ),
        (
            r"CCOC(=O)c1ccc2[nH]ccc2c1",
            "MurckoScaffold",
            r"c1ccc2[nH]ccc2c1",
            Some(r"c1ccc2[nH]ccc2c1"),
        ),
        (
            r"CCOC(=O)c1ccc2[nH]ccc2c1",
            "ExtendedMurcko",
            r"*c1ccc2[nH]ccc2c1",
            Some(r"*c1ccc2[nH]ccc2c1"),
        ),
        (
            r"CCOC(=O)c1ccc2[nH]ccc2c1",
            "MolFormula",
            r"C11H11NO2",
            Some(r"C11H11NO2"),
        ),
        (
            r"CCOC(=O)c1ccc2[nH]ccc2c1",
            "AtomBondCounts",
            r"14,15",
            Some(r"14,15"),
        ),
        (
            r"CCOC(=O)c1ccc2[nH]ccc2c1",
            "DegreeVector",
            r"0,4,8,2",
            Some(r"0,4,8,2"),
        ),
        (
            r"CCOC(=O)c1ccc2[nH]ccc2c1",
            "Mesomer",
            r"CCO[C]([O])[C]1[CH][CH][C]2N[CH][CH][C]2[CH]1_0",
            Some(r"CCO[C]([O])[C]1[CH][CH][C]2N[CH][CH][C]2[CH]1_0"),
        ),
        (
            r"CCOC(=O)c1ccc2[nH]ccc2c1",
            "HetAtomTautomer",
            r"CCO[C]([O])[C]1[CH][CH][C]2[N][CH][CH][C]2[CH]1_1_0",
            Some(r"CCO[C]([O])[C]1[CH][CH][C]2[N][CH][CH][C]2[CH]1_1_0"),
        ),
        (
            r"CCOC(=O)c1ccc2[nH]ccc2c1",
            "HetAtomProtomer",
            r"CCO[C]([O])[C]1[CH][CH][C]2[N][CH][CH][C]2[CH]1_1",
            Some(r"CCO[C]([O])[C]1[CH][CH][C]2[N][CH][CH][C]2[CH]1_1"),
        ),
        (
            r"CCOC(=O)c1ccc2[nH]ccc2c1",
            "RedoxPair",
            r"CCO[C]([O])[C]1[CH][CH][C]2N[CH][CH][C]2[CH]1",
            Some(r"CCO[C]([O])[C]1[CH][CH][C]2N[CH][CH][C]2[CH]1"),
        ),
        (
            r"CCOC(=O)c1ccc2[nH]ccc2c1",
            "Regioisomer",
            r"*OC(*)=O.CC.c1ccc2[nH]ccc2c1",
            Some(r"*OC(*)=O.CC.c1ccc2[nH]ccc2c1"),
        ),
        (r"CCOC(=O)c1ccc2[nH]ccc2c1", "NetCharge", r"0", Some(r"0")),
        (
            r"CCOC(=O)c1ccc2[nH]ccc2c1",
            "SmallWorldIndexBR",
            r"B15R2",
            Some(r"B15R2"),
        ),
        (
            r"CCOC(=O)c1ccc2[nH]ccc2c1",
            "SmallWorldIndexBRL",
            r"B15R2L8",
            Some(r"B15R2L8"),
        ),
        (
            r"CCOC(=O)c1ccc2[nH]ccc2c1",
            "ArthorSubstructureOrder",
            r"000e000f01000b0003000059000000",
            Some(r"000e000f01000b0003000059000000"),
        ),
        (
            r"CCOC(=O)c1ccc2[nH]ccc2c1",
            "HetAtomTautomerv2",
            r"[CH3]-[CH2]-[O]-[C](:[O]):[C]1:[C]:[C]:[C]2:[N]:[C]:[C]:[C]:2:[C]:1_6_0",
            Some(r"[CH3]-[CH2]-[O]-[C](:[O]):[C]1:[C]:[C]:[C]2:[N]:[C]:[C]:[C]:2:[C]:1_6_0"),
        ),
        (
            r"CCOC(=O)c1ccc2[nH]ccc2c1",
            "HetAtomProtomerv2",
            r"[CH3]-[CH2]-[O]-[C](:[O]):[C]1:[C]:[C]:[C]2:[N]:[C]:[C]:[C]:2:[C]:1_6",
            Some(r"[CH3]-[CH2]-[O]-[C](:[O]):[C]1:[C]:[C]:[C]2:[N]:[C]:[C]:[C]:2:[C]:1_6"),
        ),
        (
            r"[Pt@SP1](Cl)(Cl)([NH3])[NH3]",
            "AnonymousGraph",
            r"*[*@SP1](*)(*)*",
            Some(r"*[*@SP1](*)(*)*"),
        ),
        (
            r"[Pt@SP1](Cl)(Cl)([NH3])[NH3]",
            "ElementGraph",
            r"[NH2][Pt@SP1]([NH2])([Cl])[Cl]",
            Some(r"[NH2][Pt@SP1]([NH2])([Cl])[Cl]"),
        ),
        (
            r"[Pt@SP1](Cl)(Cl)([NH3])[NH3]",
            "CanonicalSmiles",
            r"[NH3]->[Pt@SP1](<-[NH3])([Cl])[Cl]",
            Some(r"[NH3][Pt@SP1]([NH3])([Cl])[Cl] |C:0.0,2.1|"),
        ),
        (r"[Pt@SP1](Cl)(Cl)([NH3])[NH3]", "MurckoScaffold", r"", None),
        (r"[Pt@SP1](Cl)(Cl)([NH3])[NH3]", "ExtendedMurcko", r"", None),
        (
            r"[Pt@SP1](Cl)(Cl)([NH3])[NH3]",
            "MolFormula",
            r"Cl2H6N2Pt",
            Some(r"Cl2H6N2Pt"),
        ),
        (
            r"[Pt@SP1](Cl)(Cl)([NH3])[NH3]",
            "AtomBondCounts",
            r"5,4",
            Some(r"5,4"),
        ),
        (
            r"[Pt@SP1](Cl)(Cl)([NH3])[NH3]",
            "DegreeVector",
            r"1,0,0,4",
            Some(r"1,0,0,4"),
        ),
        (
            r"[Pt@SP1](Cl)(Cl)([NH3])[NH3]",
            "Mesomer",
            r"[NH3][Pt@SP1]([NH3])([Cl])[Cl]_0",
            Some(r"[NH3][Pt@SP1]([NH3])([Cl])[Cl]_0"),
        ),
        (
            r"[Pt@SP1](Cl)(Cl)([NH3])[NH3]",
            "HetAtomTautomer",
            r"[N][Pt@SP1]([N])([Cl])[Cl]_6_0",
            Some(r"[N][Pt@SP1]([N])([Cl])[Cl]_6_0"),
        ),
        (
            r"[Pt@SP1](Cl)(Cl)([NH3])[NH3]",
            "HetAtomProtomer",
            r"[N][Pt@SP1]([N])([Cl])[Cl]_6",
            Some(r"[N][Pt@SP1]([N])([Cl])[Cl]_6"),
        ),
        (
            r"[Pt@SP1](Cl)(Cl)([NH3])[NH3]",
            "RedoxPair",
            r"[NH3][Pt@SP1]([NH3])([Cl])[Cl]",
            Some(r"[NH3][Pt@SP1]([NH3])([Cl])[Cl]"),
        ),
        (
            r"[Pt@SP1](Cl)(Cl)([NH3])[NH3]",
            "Regioisomer",
            r"[NH3]->[Pt@SP1](<-[NH3])([Cl])[Cl]",
            Some(r"[NH3][Pt@SP1]([NH3])([Cl])[Cl] |C:0.0,2.1|"),
        ),
        (
            r"[Pt@SP1](Cl)(Cl)([NH3])[NH3]",
            "NetCharge",
            r"0",
            Some(r"0"),
        ),
        (
            r"[Pt@SP1](Cl)(Cl)([NH3])[NH3]",
            "SmallWorldIndexBR",
            r"B4R0",
            Some(r"B4R0"),
        ),
        (
            r"[Pt@SP1](Cl)(Cl)([NH3])[NH3]",
            "SmallWorldIndexBRL",
            r"B4R0L0",
            Some(r"B4R0L0"),
        ),
        (
            r"[Pt@SP1](Cl)(Cl)([NH3])[NH3]",
            "ArthorSubstructureOrder",
            r"00050004010000000400007e000000",
            Some(r"00050004010000000400007e000000"),
        ),
        (
            r"[Pt@SP1](Cl)(Cl)([NH3])[NH3]",
            "HetAtomTautomerv2",
            r"[NH3]->[Pt@SP1](<-[NH3])(-[Cl])-[Cl]_0_0",
            Some(r"[NH3]-[Pt@SP1](-[NH3])(-[Cl])-[Cl]_0_0 |C:0.0,2.1|"),
        ),
        (
            r"[Pt@SP1](Cl)(Cl)([NH3])[NH3]",
            "HetAtomProtomerv2",
            r"[NH3]->[Pt@SP1](<-[NH3])(-[Cl])-[Cl]_0",
            Some(r"[NH3]-[Pt@SP1](-[NH3])(-[Cl])-[Cl]_0 |C:0.0,2.1|"),
        ),
        (r"CC", "AnonymousGraph", r"**", Some(r"**")),
        (r"CC", "ElementGraph", r"CC", Some(r"CC")),
        (r"CC", "CanonicalSmiles", r"CC", Some(r"CC")),
        (r"CC", "MurckoScaffold", r"", None),
        (r"CC", "ExtendedMurcko", r"", None),
        (r"CC", "MolFormula", r"C2H6", Some(r"C2H6")),
        (r"CC", "AtomBondCounts", r"2,1", Some(r"2,1")),
        (r"CC", "DegreeVector", r"0,0,0,2", Some(r"0,0,0,2")),
        (r"CC", "Mesomer", r"CC_0", Some(r"CC_0")),
        (r"CC", "HetAtomTautomer", r"CC_0_0", Some(r"CC_0_0")),
        (r"CC", "HetAtomProtomer", r"CC_0", Some(r"CC_0")),
        (r"CC", "RedoxPair", r"CC", Some(r"CC")),
        (r"CC", "Regioisomer", r"CC", Some(r"CC")),
        (r"CC", "NetCharge", r"0", Some(r"0")),
        (r"CC", "SmallWorldIndexBR", r"B1R0", Some(r"B1R0")),
        (r"CC", "SmallWorldIndexBRL", r"B1R0L0", Some(r"B1R0L0")),
        (
            r"CC",
            "ArthorSubstructureOrder",
            r"00020001010002000000000c000000",
            Some(r"00020001010002000000000c000000"),
        ),
        (
            r"CC",
            "HetAtomTautomerv2",
            r"[CH3]-[CH3]_0_0",
            Some(r"[CH3]-[CH3]_0_0"),
        ),
        (
            r"CC",
            "HetAtomProtomerv2",
            r"[CH3]-[CH3]_0",
            Some(r"[CH3]-[CH3]_0"),
        ),
    ];
    for &(s, f, want, want_cx) in cases {
        let mol = crate::parse(s).expect("parses");
        let func = RdkitHashFunction::from_name(f).expect("function");
        assert_eq!(
            super::rdkit_mol_hash(&mol, func, false).unwrap(),
            want,
            "{f} {s}"
        );
        // `None`: RDKit raises (an empty scaffold has no output order).
        assert_eq!(
            super::rdkit_mol_hash(&mol, func, true).ok().as_deref(),
            want_cx,
            "{f} {s} (CX)"
        );
    }
}
