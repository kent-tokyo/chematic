//! Regression tests for RDKit Normalize follow-up patterns.

use chematic_chem::normalize_groups;
use chematic_smiles::parse;

#[test]
fn normalize_groups_preserves_labelled_nitro_oxygen() {
    for smiles in [
        "c1ccccc1[N+](=O)[O-:28]",
        "c1ccccc1[N+]([O-])=[O:28]",
        "c1ccccc1[N+](=O)[18O-]",
        "C[N+]([O-])C",
    ] {
        let mol = parse(smiles).unwrap();
        assert_eq!(
            chematic_smiles::rdkit_canonical_smiles(&normalize_groups(&mol)).unwrap(),
            chematic_smiles::rdkit_canonical_smiles(&mol).unwrap(),
            "{smiles}"
        );
    }
}
#[test]
fn normalize_groups_preserves_oxygen_substituted_sulfur() {
    for smiles in ["O=S1OCCO1", "COS(=O)OC", "CS(=O)O", "CS(=O)[O-]"] {
        let mol = parse(smiles).unwrap();
        assert_eq!(
            chematic_smiles::rdkit_canonical_smiles(&normalize_groups(&mol)).unwrap(),
            chematic_smiles::rdkit_canonical_smiles(&mol).unwrap(),
            "{smiles}"
        );
    }
}

#[test]
fn normalize_groups_recombines_aromatic_sulfur_charge() {
    for (input, expected) in [
        ("[S-]c1[n+](C)cccc1", "Cn1ccccc1=S"),
        ("[S-]c1nc2ccccc2c2cccc[n+]12", "S=c1nc2ccccc2c2ccccn12"),
    ] {
        let mol = parse(input).unwrap();
        let normalized = normalize_groups(&mol);
        assert_eq!(
            chematic_smiles::rdkit_canonical_smiles(&normalized).unwrap(),
            expected,
            "{input}"
        );
        assert_eq!(
            chematic_smiles::rdkit_canonical_smiles(&normalize_groups(&normalized)).unwrap(),
            expected,
            "idempotence: {input}"
        );
    }
    let mapped = parse("[S-:1]c1[n+](C)cccc1").unwrap();
    assert_eq!(
        chematic_smiles::rdkit_canonical_smiles(&normalize_groups(&mapped)).unwrap(),
        chematic_smiles::rdkit_canonical_smiles(&mapped).unwrap()
    );
}

#[test]
fn normalize_groups_shifts_conjugated_iminium_charge() {
    for (input, expected) in [
        ("CNC=CC=[N+](C)C", "C[NH+]=CC=CN(C)C"),
        (
            "COC1=CC(=C(NC=CC2=[N+](C)C3=C(C=CC=C3)C2(C)C)C=C1)OC",
            "COc1ccc([NH+]=CC=C2N(C)c3ccccc3C2(C)C)c(OC)c1",
        ),
    ] {
        let mol = parse(input).unwrap();
        let normalized = normalize_groups(&mol);
        assert_eq!(
            chematic_smiles::rdkit_canonical_smiles(&normalized).unwrap(),
            expected,
            "{input}"
        );
        assert_eq!(
            chematic_smiles::rdkit_canonical_smiles(&normalize_groups(&normalized)).unwrap(),
            expected,
            "idempotence: {input}"
        );
    }
    for input in [
        "CN(C)C=CC=[N+](C)C",
        "CNC=CC=[N+](C)[O-]",
        "CN/C=C/C=[N+](C)C",
    ] {
        let mol = parse(input).unwrap();
        assert_eq!(
            chematic_smiles::rdkit_canonical_smiles(&normalize_groups(&mol)).unwrap(),
            chematic_smiles::rdkit_canonical_smiles(&mol).unwrap(),
            "{input}"
        );
    }
}
