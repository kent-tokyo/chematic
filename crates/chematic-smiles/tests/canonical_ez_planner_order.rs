//! The complete aromatic E/Z planner checks its candidate spellings in
//! ascending order and skips plans that leave a stereo end without a
//! directional token; the output is the minimum of the exhaustive check
//! (strings recorded before that change).

const CASES: &[(&str, &str)] = &[
    (
        "CC(C)(C)/N=c1\\c(O)c(O)\\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
        "CC(C)(/N=c/1\\c(c(\\c1=N/[C@@H](Cc3ccc(cc3)NC(c2c(Cl)cncc2Cl)=O)C(=O)O)O)O)C",
    ),
    (
        "c1(O[H])c(=N\\C(C([H])([H])[H])(C([H])([H])[H])C([H])([H])[H])/c(=N/[C@@]([H])(C(c2c(c([H])c(c([H])c2[H])N([H])C(=O)c2c(c([H])nc([H])c2Cl)Cl)[H])([H])[H])C(O[H])=O)c1O[H]",
        "[H]C(C(C([H])([H])[H])(/N=c/3\\c(O[H])c(\\c3=N/[C@](C(O[H])=O)(C(c2c([H])c(c(c(c2[H])[H])N(C(=O)c1c(c([H])nc(c1Cl)[H])Cl)[H])[H])([H])[H])[H])O[H])C([H])([H])[H])([H])[H]",
    ),
    (
        "C/N=c1\\c(O)c(O)\\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
        "[C@H](Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)(C(O)=O)/N=c/1\\c(O)c(O)\\c1=N/C",
    ),
    (
        "CC(C)(C)/N=C(\\N)c1ccc(N2CCCN(c3ccc(/C(N)=N/C(C)(C)C)cc3)CC2)cc1",
        "C2CN(CCCN2c3ccc(/C(=N/C(C)(C)C)N)cc3)c1ccc(cc1)/C(=N/C(C)(C)C)N",
    ),
];

#[test]
fn planner_output_is_unchanged() {
    for &(input, want) in CASES {
        let mol = chematic_smiles::parse(input).unwrap();
        assert_eq!(chematic_smiles::canonical_smiles(&mol), want, "{input}");
    }
}
