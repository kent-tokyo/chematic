use chematic_smiles::rdkit_mol_from_mol2_block;
fn make(types: &[&str]) -> String {
    let n = types.len();
    let mut s = format!(
        "@<TRIPOS>MOLECULE\nfixture\n{n} {} 0 0 0\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n",
        n - 1
    );
    for (i, t) in types.iter().enumerate() {
        s.push_str(&format!(
            "{} A{} {:.3} {:.3} {:.3} {t} 1 RES 0\n",
            i + 1,
            i + 1,
            i as f64 * 1.2,
            (i % 2) as f64,
            (i % 3) as f64 * 0.25
        ));
    }
    s.push_str("@<TRIPOS>BOND\n");
    for i in 1..n {
        s.push_str(&format!("{i} 1 {} 1\n", i + 1));
    }
    s
}
#[test]
fn tripos_substructure_cleanup_matches_rdkit_reference() {
    // Independently generated with RDKit 2025.09.3 on these minimal graphs.
    let cases: &[(&[&str], &str)] = &[
        (
            &["P.3", "O.co2", "O.co2", "O.3", "O.3"],
            "O=P([O-])([O-])[O-]",
        ),
        (&["S.o2", "O.co2", "O.co2", "C.3"], "[C][S@@](=O)[O-]"),
        (&["C.2", "O.co2", "O.co2", "C.3"], "[C]C(=O)[O-]"),
        (&["C.cat", "N.pl3", "N.pl3", "N.pl3"], "[N+]=C([N-2])[N-2]"),
        (&["C.cat", "N.pl3", "N.pl3"], "[N+]=[C][N-2]"),
    ];
    for &(types, expected) in cases {
        let result = rdkit_mol_from_mol2_block(&make(types), true, true, true).unwrap();
        assert_eq!(result.smiles, expected, "{types:?}");
    }
}
