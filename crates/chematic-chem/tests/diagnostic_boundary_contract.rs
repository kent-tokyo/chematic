//! Public errors must retain the input details needed by binding users.
use chematic_chem::{LocalResonanceLimits, enrichment_factors, enumerate_local_cation_resonance};
use chematic_core::{BondIdx, BondOrder, MoleculeBuilder};
use chematic_smiles::parse;

#[test]
fn resonance_diagnostics_identify_invalid_limits_and_unsupported_input() {
    let mut query_bond = parse("CC").unwrap();
    query_bond.set_bond_order(BondIdx(0), BondOrder::QueryAny);
    let cases = [
        (
            parse("CCO").unwrap(),
            LocalResonanceLimits {
                max_forms: 0,
                ..Default::default()
            },
            vec!["max_forms", "1..=1024", "got 0"],
        ),
        (
            parse("CCO").unwrap(),
            LocalResonanceLimits {
                max_atoms: 2,
                ..Default::default()
            },
            vec!["atom count 3", "limit 2"],
        ),
        (
            MoleculeBuilder::new().build(),
            LocalResonanceLimits::default(),
            vec!["non-empty molecule"],
        ),
        (
            parse("C[C@H](F)Cl").unwrap(),
            LocalResonanceLimits::default(),
            vec!["stereochemical annotations"],
        ),
        (
            parse("*C").unwrap(),
            LocalResonanceLimits::default(),
            vec!["query atom 0"],
        ),
        (
            query_bond,
            LocalResonanceLimits::default(),
            vec!["query bond 0"],
        ),
    ];
    for (mol, limits, details) in cases {
        let error = match enumerate_local_cation_resonance(&mol, &limits) {
            Err(error) => error,
            Ok(_) => panic!("unsupported input must fail before enumeration"),
        };
        let message = error.to_string();
        for detail in details {
            assert!(
                message.contains(detail),
                "{detail:?} missing from {message:?}"
            );
        }
    }
}

#[test]
fn screening_diagnostics_distinguish_missing_data_from_undefined_enrichment() {
    for (labels, fractions, reason) in [
        (&[][..], &[1.0][..], "ranking is empty"),
        (&[true][..], &[][..], "fraction list is empty"),
        (&[false][..], &[1.0][..], "no actives"),
    ] {
        assert!(
            enrichment_factors(labels, fractions)
                .unwrap_err()
                .to_string()
                .contains(reason)
        );
    }
}
