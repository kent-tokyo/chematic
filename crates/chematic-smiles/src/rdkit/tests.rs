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
fn unsupported_inputs_are_errors() {
    assert!(matches!(
        rd("F[Pt@SP1](Cl)(Br)I"),
        Err(RdkitSmilesError::Unsupported(_))
    ));
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
