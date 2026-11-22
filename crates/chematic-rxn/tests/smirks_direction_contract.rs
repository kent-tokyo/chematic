use chematic_rxn::{PreparedReaction, RxnError, TransformError};

#[test]
fn smarts_only_directional_reactant_is_rejected_in_default_profile() {
    let smirks = "[#6;X3:1]=[#6:2](/[#6:3])>>[#6:1]1[#6:2]([#6:3])O1";
    let error = match PreparedReaction::new(smirks) {
        Ok(_) => panic!("unsupported stereo must not be silently discarded"),
        Err(error) => error,
    };
    assert!(
        matches!(error, TransformError::SmirksParse(RxnError::UnsupportedReactantTemplate { reason, .. }) if reason.contains("double-bond stereo"))
    );
}

#[test]
fn oversized_smirks_returns_a_resource_error() {
    let limit = chematic_rxn::ReactionParseLimits::default().max_input_bytes;
    let smirks = "C".repeat(limit + 1);
    let error = match PreparedReaction::new(&smirks) {
        Ok(_) => panic!("oversized template must be rejected"),
        Err(error) => error,
    };
    assert!(
        matches!(error, TransformError::SmirksParse(RxnError::ResourceLimit { resource: "input bytes", actual, limit: actual_limit }) if actual == limit + 1 && actual_limit == limit)
    );
}
