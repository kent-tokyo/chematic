use chematic_rxn::{
    ComponentRole, ReactionDocumentError, ReactionLoss, ReactionQueryError, StoichiometryComponent,
    StoichiometryError, StoichiometryStep, analyze_components,
};

#[test]
fn public_reaction_errors_preserve_paths_limits_and_loss_counts() {
    let query_errors = [
        ReactionQueryError::ResourceLimit {
            resource: "rows",
            actual: 3,
            limit: 2,
        },
        ReactionQueryError::SmartsParseError {
            smarts: "[".into(),
            source: "unclosed atom".into(),
        },
        ReactionQueryError::MissingArrowDelimiter,
        ReactionQueryError::InvalidAgentsSection,
        ReactionQueryError::MapNumberMismatch {
            map_num: 7,
            message: "missing reactant".into(),
        },
    ];
    for error in query_errors {
        assert!(!error.to_string().is_empty());
    }

    let errors = [
        ReactionDocumentError::parse_message("broken RXN"),
        ReactionDocumentError::InvalidDocument("duplicate id".into()),
        ReactionDocumentError::Losses(vec![ReactionLoss {
            field: "conditions".into(),
            detail: "legacy format".into(),
        }]),
    ];
    assert!(errors[0].to_string().contains("broken RXN"));
    assert!(errors[1].to_string().contains("duplicate id"));
    assert!(errors[2].to_string().contains("1 field"));
}

#[test]
fn stoichiometry_rejects_invalid_coefficients_empty_components_and_bad_smiles() {
    let component = |id: &str, smiles: &str, coefficient: f64| StoichiometryComponent {
        id: id.into(),
        role: ComponentRole::Reactant,
        smiles: smiles.into(),
        coefficient,
    };
    for (component, expected) in [
        (component("negative", "C", -1.0), "invalid coefficient"),
        (component("nan", "C", f64::NAN), "invalid coefficient"),
        (component("empty", "", 1.0), "empty component"),
        (component("bad", "C1", 1.0), "failed to parse"),
    ] {
        let step = StoichiometryStep {
            id: "step".into(),
            components: vec![component],
        };
        let error = analyze_components(&step).unwrap_err();
        assert!(error.to_string().contains(expected));
        assert!(matches!(
            error,
            StoichiometryError::InvalidCoefficient { .. }
                | StoichiometryError::EmptyComponent { .. }
                | StoichiometryError::Parse { .. }
        ));
    }
}
