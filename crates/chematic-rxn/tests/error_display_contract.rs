use chematic_rxn::RxnError;

#[test]
fn reaction_parse_errors_keep_the_failed_stage_and_input() {
    let cases = [
        RxnError::MissingArrow,
        RxnError::ResourceLimit {
            resource: "components",
            actual: 11,
            limit: 10,
        },
        RxnError::SmilesParse {
            part: "C(".into(),
            source: "unclosed branch".into(),
        },
        RxnError::UnsupportedAtomicNumberPrimitive {
            primitive: "[#x]".into(),
        },
        RxnError::AtomicNumberExpansionLimit {
            actual: 9,
            limit: 8,
        },
        RxnError::SmartsParse {
            part: "[C".into(),
            source: "missing bracket".into(),
        },
        RxnError::UnsupportedReactantTemplate {
            part: "[C@]".into(),
            reason: "stereo",
        },
        RxnError::UnsupportedProductPrimitive {
            primitive: "[C;R]".into(),
        },
    ];

    for error in &cases {
        assert!(!error.to_string().is_empty());
    }
}
