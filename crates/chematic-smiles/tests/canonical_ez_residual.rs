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

#[test]
fn stable_key_fails_closed_for_coupled_ez_systems() {
    let mol = parse(r"CC1CNC(/C=C\C=C/C=C\C2(C)C(=C(O)C(C2)=O)C(/C=C\C=C/C=C\C=C/C=C1)=O)=O")
        .expect("coupled E/Z fixture parses");
    assert_eq!(canonical_smiles_stable_key(&mol), None);
}

#[test]
fn stable_key_fails_closed_for_aromatic_stash_residuals() {
    for input in [
        r"CC/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
        r"CCCC(C)/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
        r"COCC/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
    ] {
        let mol = parse(input).expect("fixture parses");
        assert_eq!(canonical_smiles_stable_key(&mol), None, "{input}");
    }
}

#[test]
fn aromatic_stash_residual_outputs_are_fail_closed_variants() {
    let fixtures = [
        (
            r"CC/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
            [
                r"c/3(c(/c(c3=N\CC)=N\[C@@H](Cc1ccc(NC(=O)c2c(cncc2Cl)Cl)cc1)C(O)=O)O)O",
                r"c3(c(c(/c3=N/CC)=N\[C@@H](Cc1ccc(NC(=O)c2c(cncc2Cl)Cl)cc1)C(O)=O)O)O",
            ],
        ),
        (
            r"CCCC(C)/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
            [
                r"c/1(c(/c(c1=N\[C@H](C(O)=O)Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)=N\C(C)CCC)O)O",
                r"c1(c(c(/c1=N/[C@H](C(O)=O)Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)=N\C(C)CCC)O)O",
            ],
        ),
        (
            r"COCC/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
            [
                r"c/1(O)c(O)/c(=N\[C@@H](Cc3ccc(cc3)NC(=O)c2c(Cl)cncc2Cl)C(=O)O)c1=N\CCOC",
                r"c1(O)c(O)c(=N/[C@@H](Cc3ccc(cc3)NC(=O)c2c(Cl)cncc2Cl)C(=O)O)\c1=N\CCOC",
            ],
        ),
    ];

    for (input, expected_variants) in fixtures {
        let mol = parse(input).expect("fixture parses");
        let output = canonical_smiles(&mol);
        assert!(
            expected_variants.contains(&output.as_str()),
            "unexpected residual output: {output}"
        );
        assert_eq!(
            canonical_smiles_stable_key(&mol),
            None,
            "residual must stay fail-closed: {input}"
        );
    }
}

#[test]
fn canonical_writer_keeps_ez_when_carrier_is_a_ring_edge() {
    // The carrier selected for this fused ring alkene can be a DFS ring
    // closure.  The close side of a SMILES ring edge is intentionally emitted
    // without a direction marker, so selecting that occurrence used to erase
    // the E/Z specification on re-parse.
    let input = "C1=C(/C=C/c2ccccc2)N2CCCN=C2c2ccccc21";
    let mol = parse(input).expect("fixture parses");
    let output = canonical_smiles(&mol);
    let reparsed = parse(&output)
        .unwrap_or_else(|e| panic!("canonical output '{output}' does not parse: {e}"));

    let directional_bonds = |m: &chematic_core::Molecule| {
        m.bonds()
            .filter(|(_, bond)| {
                matches!(
                    bond.order,
                    chematic_core::BondOrder::Up | chematic_core::BondOrder::Down
                )
            })
            .count()
    };
    assert_eq!(
        directional_bonds(&reparsed),
        directional_bonds(&mol),
        "canonicalization must not erase ring-adjacent E/Z markers: {output}"
    );
}
