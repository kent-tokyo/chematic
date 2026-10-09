//! Expected values are RDKit 2026.03.1's
//! `rdMolStandardize.TautomerEnumerator()` results for
//! `m = Chem.MolFromSmiles(input)`: `Chem.MolToSmiles(te.Canonicalize(m))`,
//! `sorted(Chem.MolToSmiles(t) for t in te.Enumerate(m))`, the result
//! status and `TautomerEnumerator.ScoreTautomer(m)`.

use super::{
    RdkitTautomerStatus, rdkit_canonical_tautomer, rdkit_enumerate_tautomers, rdkit_tautomer_score,
};

fn check(input: &str, canon: &str, score: i32, tautomers: &[&str]) {
    let mol = crate::parse(input).expect("parses");
    assert_eq!(
        rdkit_canonical_tautomer(&mol).as_deref(),
        Ok(canon),
        "canonical tautomer of {input}"
    );
    assert_eq!(rdkit_tautomer_score(&mol), Ok(score), "score of {input}");
    let res = rdkit_enumerate_tautomers(&mol).expect("enumerates");
    assert_eq!(
        res.status,
        RdkitTautomerStatus::Completed,
        "status of {input}"
    );
    let mut got = res.smiles.clone();
    got.sort();
    assert_eq!(got, tautomers, "tautomers of {input}");
}

#[test]
fn enone_phenol_with_bond_stereo() {
    check(
        "CC(=O)/C=C/c1ccc(O)cc1",
        "CC(=O)C=Cc1ccc(O)cc1",
        254,
        &[
            "C=C(O)C=CC1=CCC(=O)C=C1",
            "C=C(O)C=CC1C=CC(=O)C=C1",
            "C=C(O)C=Cc1ccc(O)cc1",
            "C=C(O)CC=C1C=CC(=O)C=C1",
            "CC(=O)C=CC1=CCC(=O)C=C1",
            "CC(=O)C=CC1C=CC(=O)C=C1",
            "CC(=O)C=Cc1ccc(O)cc1",
            "CC(=O)CC=C1C=CC(=O)C=C1",
            "CC(O)=CC=C1C=CC(=O)C=C1",
        ],
    );
}

#[test]
fn ring_enone_keeps_unmodified_bond_stereo() {
    check(
        r"O=C1CCCCC/C1=C\c1ccccc1",
        r"O=C1CCCCC/C1=C\c1ccccc1",
        253,
        &[r"O=C1CCCCC/C1=C\c1ccccc1", r"OC1=CCCCC/C1=C\c1ccccc1"],
    );
}

#[test]
fn guanidine_bond_stereo_removed() {
    check(
        r"C/N=C(\NC)NCc1ccc(I)cc1",
        "CN=C(NC)NCc1ccc(I)cc1",
        257,
        &["CN=C(NC)NCc1ccc(I)cc1", "CNC(=NCc1ccc(I)cc1)NC"],
    );
}

#[test]
fn amino_acid_loses_sp3_stereo() {
    check(
        "N[C@H](CC(=O)O)C(=O)O",
        "NC(CC(=O)O)C(=O)O",
        10,
        &[
            "N=C(C=C(O)O)C(O)O",
            "N=C(CC(=O)O)C(O)O",
            "NC(=CC(=O)O)C(O)O",
            "NC(C=C(O)O)=C(O)O",
            "NC(C=C(O)O)C(=O)O",
            "NC(CC(=O)O)=C(O)O",
            "NC(CC(=O)O)C(=O)O",
        ],
    );
}

#[test]
fn phosphonate() {
    check(
        "N[C@@H](CNCP(=O)(O)O)C(=O)O",
        "NC(CNCP(=O)(O)O)C(=O)O",
        7,
        &[
            "N=C(CNCP(=O)(O)O)C(O)O",
            "NC(=CNCP(=O)(O)O)C(O)O",
            "NC(CNCP(=O)(O)O)=C(O)O",
            "NC(CNCP(=O)(O)O)C(=O)O",
        ],
    );
}

#[test]
fn adenine_nitrate() {
    check(
        "Nc1ncnc2c1ncn2CCO[N+](=O)[O-]",
        "Nc1ncnc2c1ncn2CCO[N+](=O)[O-]",
        202,
        &[
            "N=c1[nH]cnc2c1ncn2CCO[N+](=O)[O-]",
            "N=c1nc[nH]c2c1ncn2CCO[N+](=O)[O-]",
            "Nc1ncnc2c1ncn2CCO[N+](=O)[O-]",
        ],
    );
}

#[test]
fn phosphite_to_h_phosphonate() {
    check(
        "NCCCC(N)OP(O)O",
        "NCCCC(N)O[PH](=O)O",
        0,
        &["NCCCC(N)OP(O)O", "NCCCC(N)O[PH](=O)O"],
    );
}

#[test]
fn pyridopyrazinediol() {
    check(
        "Oc1nc2ccncc2nc1O",
        "O=c1[nH]c2ccncc2[nH]c1=O",
        200,
        &[
            "O=c1[nH]c2ccncc2[nH]c1=O",
            "O=c1[nH]c2ccncc2nc1O",
            "O=c1[nH]c2cnccc2nc1O",
            "O=c1nc2cc[nH]cc-2[nH]c1=O",
            "O=c1nc2cc[nH]cc-2nc1O",
            "Oc1nc2ccncc2nc1O",
        ],
    );
}

#[test]
fn biguanide() {
    check(
        "CCCCNC(=N)NC(=N)N",
        "CCCCN=C(N)N=C(N)N",
        11,
        &[
            "CCCCN=C(N)N=C(N)N",
            "CCCCN=C(N)NC(=N)N",
            "CCCCNC(=N)N=C(N)N",
            "CCCCNC(=N)NC(=N)N",
            "CCCCNC(N)=NC(=N)N",
        ],
    );
}

#[test]
fn aminoisoxazolol() {
    check(
        "NC(C(=O)O)c1cc(O)no1",
        "NC(C(=O)O)c1cc(=O)[nH]o1",
        109,
        &[
            "N=C(C(=O)O)C1C=C(O)NO1",
            "N=C(C(=O)O)C1CC(=O)NO1",
            "N=C(C(=O)O)C1CC(O)=NO1",
            "N=C(c1cc(=O)[nH]o1)C(O)O",
            "N=C(c1cc(O)no1)C(O)O",
            "NC(=C(O)O)c1cc(=O)[nH]o1",
            "NC(=C(O)O)c1cc(O)no1",
            "NC(C(=O)O)=C1C=C(O)NO1",
            "NC(C(=O)O)=C1CC(=O)NO1",
            "NC(C(=O)O)=C1CC(O)=NO1",
            "NC(C(=O)O)c1cc(=O)[nH]o1",
            "NC(C(=O)O)c1cc(O)no1",
        ],
    );
}

#[test]
fn cresol() {
    check(
        "Cc1ccc(C(C)C)cc1O",
        "Cc1ccc(C(C)C)cc1O",
        253,
        &[
            "CC(C)=C1C=CC(C)=C(O)C1",
            "CC(C)=C1C=CC(C)C(=O)C1",
            "CC(C)=C1C=CC(C)C(O)=C1",
            "CC(C)C1=CC(=O)C(C)C=C1",
            "CC1=CC=C(C(C)C)CC1=O",
            "CC1=CCC(=C(C)C)C=C1O",
            "CC1=CCC(=C(C)C)CC1=O",
            "CC1=CCC(C(C)C)=CC1=O",
            "Cc1ccc(C(C)C)cc1O",
        ],
    );
}

#[test]
fn cyclic_urea_keeps_distant_stereo() {
    check(
        "O=C1N[C@H](CO)C[C@@H](C(=O)O)N1",
        "O=C1NC(C(=O)O)C[C@@H](CO)N1",
        12,
        &[
            "O=C(O)C1C[C@@H](CO)N=C(O)N1",
            "O=C(O)C1C[C@@H](CO)NC(O)=N1",
            "O=C1NC(=C(O)O)C[C@@H](CO)N1",
            "O=C1NC(C(=O)O)C[C@@H](CO)N1",
            "OC[C@@H]1CC(=C(O)O)N=C(O)N1",
            "OC[C@@H]1CC(=C(O)O)NC(O)=N1",
        ],
    );
}

#[test]
fn max_transforms_reached() {
    let mol = crate::parse("CC1CC(=O)C2=C(O)c3c(O)ccc(O)c3CC2C1").expect("parses");
    let res = rdkit_enumerate_tautomers(&mol).expect("enumerates");
    assert_eq!(res.status, RdkitTautomerStatus::MaxTransformsReached);
    assert_eq!(res.smiles.len(), 237);
    assert_eq!(
        rdkit_canonical_tautomer(&mol).as_deref(),
        Ok("CC1CC(=O)C2C(=O)c3c(O)ccc(O)c3CC2C1")
    );
}
