//! Tests extracted from `transform.rs`; kept as a child module so private helpers remain testable.

use super::*;
use crate::reaction::ReactionParseLimits;

#[test]
fn reactant_template_limits_accept_the_boundary_and_report_excess_counts() {
    for rdkit in [false, true] {
        let limits = ReactionParseLimits {
            max_components_per_side: 2,
            max_atoms_per_molecule: 3,
            max_bonds_per_molecule: 2,
            ..Default::default()
        };
        let (mols, queries) =
            parse_reactant_templates_with_limits("CCC.C", rdkit, &limits).unwrap();
        assert_eq!(mols.len(), 2);
        assert_eq!(queries.len(), 2);
        for (text, limits, resource, actual, limit) in [
            ("C.C.C", limits, "reactants", 3, 2),
            ("CCCC", limits, "atoms per molecule", 4, 3),
            ("C1CC1", limits, "bonds per molecule", 3, 2),
        ] {
            let error = match parse_reactant_templates_with_limits(text, rdkit, &limits) {
                Ok(_) => panic!("limit was not enforced: {text}"),
                Err(e) => e,
            };
            let message = error.to_string();
            assert!(
                message.contains(resource)
                    && message.contains(&format!("limit {limit}"))
                    && message.contains(&format!("got {actual}")),
                "{message}"
            );
            assert!(
                matches!(error,TransformError::SmirksParse(RxnError::ResourceLimit{resource:r,actual:a,limit:l}) if r==resource && a==actual && l==limit)
            );
        }
    }
}

#[test]
fn match_combination_limit_is_applied_before_product_generation() {
    let mol = chematic_smiles::parse("CC").unwrap();
    let smirks = "[C:1].[C:2]>>[C:1][C:2]";
    let error = match run_reactants_strict_with_limits(
        smirks,
        &[&mol, &mol],
        &ReactionTransformLimits { max_matches: 2 },
    ) {
        Ok(_) => panic!("four combinations must exceed limit two"),
        Err(e) => e,
    };
    assert!(
        error
            .to_string()
            .contains("reaction match combinations exceeds limit 2 (got 4)")
    );
    assert!(matches!(
        error,
        TransformError::ResourceLimit {
            resource: "reaction match combinations",
            actual: 4,
            limit: 2
        }
    ));
    let accepted = run_reactants_strict_with_limits(
        smirks,
        &[&mol, &mol],
        &ReactionTransformLimits { max_matches: 4 },
    )
    .unwrap();
    assert_eq!(accepted.len(), 4);
    for set in accepted {
        assert_eq!(set.len(), 1);
        assert_eq!(chematic_smiles::canonical_smiles(&set[0]), "CC");
    }
}

#[test]
fn component_grouping_keeps_nested_branches_and_recursive_queries_together() {
    for (text, expected) in [
        ("([C:1].[O:2]).[N:3]", vec!["[C:1].[O:2]", "[N:3]"]),
        ("(C(C)C.O).N", vec!["C(C)C.O", "N"]),
        (
            "([$([C](=O)O):1].[N:2]).O",
            vec!["[$([C](=O)O):1].[N:2]", "O"],
        ),
        ("(C.O)", vec!["C.O"]),
    ] {
        assert_eq!(split_components(text), expected);
    }
    let mol = chematic_smiles::parse("CO").unwrap();
    let out = run_reactants_strict("([C:1].[O:2])>>[C:1][O:2]", &[&mol]).unwrap();
    assert_eq!(out.len(), 1);
    assert_eq!(
        chematic_smiles::canonical_smiles(&out[0][0]),
        chematic_smiles::canonical_smiles(&mol)
    );
}

#[test]
fn explicit_h_product_spelling_is_normalized_and_other_query_primitives_are_rejected() {
    for (written, expected) in [
        ("[C;H2:1]O", "[CH2:1]O"),
        ("[N;H]", "[NH]"),
        ("[#6:1]", "[#6:1]"),
    ] {
        assert_eq!(normalize_product_templates(written).unwrap(), expected);
    }
    for text in [
        "[C;X4:1]",
        "[C;H123:1]",
        "[C;H1:x]",
        "[C,N:1]",
        "[!C:1]",
        "[$(C):1]",
    ] {
        let error = normalize_product_templates(text).unwrap_err();
        assert!(error.to_string().contains(text), "{text}");
        assert!(matches!(
            error,
            TransformError::SmirksParse(RxnError::UnsupportedProductPrimitive { .. })
        ));
    }
}
