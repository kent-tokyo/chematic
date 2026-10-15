use chematic_smarts::QueryHydrogenError;

#[test]
fn query_hydrogen_errors_have_stable_messages() {
    let cases = [
        QueryHydrogenError::InvalidGraph,
        QueryHydrogenError::ResourceLimit,
        QueryHydrogenError::RecursiveQuery,
        QueryHydrogenError::UnsupportedStereo { atom: 7 },
    ];
    let messages: Vec<_> = cases.iter().map(ToString::to_string).collect();
    assert_eq!(
        messages,
        [
            "query has inconsistent bonds or adjacency",
            "query-H folding exceeds its graph or expression limits",
            "query-H folding does not support recursive SMARTS",
            "query atom 7 needs one removable H and a complete tetrahedral frame",
        ]
    );
}
