use chematic_smiles::{parse, rdkit_smiles};

#[test]
fn branched_and_conflicting_alkene_markers_follow_rdkit() {
    // Canonical references generated independently with RDKit 2025.09.3.
    let cases = [
        (r"F/C(Cl)=C(/Br)\I", r"F/C(Cl)=C(/Br)I"),
        (r"F/C(/Cl)=C(/Br)\I", r"F/C(Cl)=C(/Br)I"),
        (r"F/C(\Cl)=C(/Br)\I", "FC(Cl)=C(Br)I"),
        (r"F/C=C/C=C/F", r"F/C=C/C=C/F"),
        (r"F/C(Cl)=C(/Br)/C=C(/Cl)\F", r"FC(Cl)=C(Br)/C=C(\F)Cl"),
        (r"C1/C=C/C=C/1", "C1=CCC=C1"),
        (
            r"F/C(Cl)=C(/Br)\C(=C(/F)Cl)/I",
            r"FC(Cl)=C(I)/C(Br)=C(\F)Cl",
        ),
    ];
    for (input, expected) in cases {
        let molecule = parse(input).unwrap();
        assert_eq!(
            rdkit_smiles(&molecule, &Default::default()).unwrap(),
            expected,
            "{input}"
        );
    }
}
