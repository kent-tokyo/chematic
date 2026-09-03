//! Regression fixtures for cosmetic E/Z carrier placement.

use chematic_smiles::{canonical_smiles, canonical_smiles_stable_key, parse};

#[test]
fn trisubstituted_ez_carrier_spelling_converges() {
    let output_a = "O=C(c2cc(OC)c(OC)c(c2)OC)/C(=C/c1ccc(c(O)c1)OC)C";
    let output_b = "O=C(c2cc(OC)c(OC)c(c2)OC)C(=C/c1ccc(c(O)c1)OC)/C";
    let a = canonical_smiles(&parse(output_a).expect("fixture A parses"));
    let b = canonical_smiles(&parse(output_b).expect("fixture B parses"));
    assert_eq!(a, b, "E/Z carrier choice must not depend on input flank");
}

#[test]
fn stable_key_accepts_converged_ez_spelling() {
    let mol = parse("C/C=C/C").expect("fixture parses");
    assert_eq!(
        canonical_smiles_stable_key(&mol),
        Some(r"C(=C/C)\C".to_string())
    );
}
