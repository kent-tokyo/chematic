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
