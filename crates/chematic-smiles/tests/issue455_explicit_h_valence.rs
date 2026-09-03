//! Regression coverage for explicit bracket H counts at the validation boundary.

use chematic_core::validate_valence;
use chematic_smiles::parse;

#[test]
fn issue455_rejects_stale_h_on_neutral_ring_nitrogen() {
    let mol = parse("C(C[NH]1CCCC1)OC").expect("parser accepts the issue reproducer");
    assert!(
        !validate_valence(&mol).is_empty(),
        "explicit [NH] must not be accepted when the graph leaves no N-H valence"
    );
}

#[test]
fn explicit_h_regressions_preserve_valid_charged_and_aromatic_atoms() {
    for smiles in ["[NH4+]", "c1cc[nH]c1"] {
        let mol = parse(smiles).expect("valid regression fixture must parse");
        assert!(
            validate_valence(&mol).is_empty(),
            "{smiles} must remain valid"
        );
    }
}

#[test]
fn issue455_rejects_stale_h_radical_spelling() {
    let mol = parse("[CH]C").expect("parser accepts the stale-H fixture");
    assert!(
        !validate_valence(&mol).is_empty(),
        "[CH]C must be rejected when explicit H disagrees with graph valence"
    );
}
