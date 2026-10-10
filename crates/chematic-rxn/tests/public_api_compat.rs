use chematic_rxn::{ReactionTransformDiagnostics, TracedReactionTransformReport};

#[test]
fn traced_report_keeps_the_v1_struct_literal_shape() {
    let _report = TracedReactionTransformReport {
        products: Vec::new(),
        diagnostics: ReactionTransformDiagnostics {
            accepted_matches: 0,
            applied_products: 0,
            valence_rejected_matches: 0,
            truncated_matches: false,
        },
    };
}
