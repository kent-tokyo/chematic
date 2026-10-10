use crate::transform::*;
use chematic_core::Molecule;
use chematic_perception::find_sssr;
use chematic_smiles::{canonical_smiles, parse};
fn products(rows: Vec<Vec<Molecule>>) -> Vec<Vec<String>> {
    rows.into_iter()
        .map(|set| set.iter().map(canonical_smiles).collect())
        .collect()
}
#[test]
fn strict_variant_application_and_precomputed_rings_preserve_product_order() {
    let mol = parse("CCNC=O").unwrap();
    let rings = find_sssr(&mol);
    let limits = ReactionTransformLimits::default();
    for template in [
        "[#7:1][C:2](=[O:3])>>[#7:1][C:2](=[O:3])",
        "[N:1][C:2](=[O:3])>>[N:1][C:2](=[O:3])",
    ] {
        let prepared = PreparedReaction::new(template).unwrap();
        assert!(!prepared.has_tetrahedral_reactant_stereo());
        let strict = products(prepared.run_reactants_strict(&[&mol]).unwrap());
        assert_eq!(strict.len(), 1);
        assert_eq!(strict[0], vec![canonical_smiles(&parse("NC=O").unwrap())]);
        assert_eq!(
            products(
                prepared
                    .run_reactants_strict_with_limits(&[&mol], &limits)
                    .unwrap()
            ),
            strict
        );
        assert_eq!(
            products(
                prepared
                    .run_reactants_strict_with_rings(&[&mol], &[&rings])
                    .unwrap()
            ),
            strict
        );
        let ordinary = products(prepared.run_reactants(&[&mol]).unwrap());
        assert_eq!(ordinary, vec![vec![canonical_smiles(&mol)]]);
        let report = prepared
            .run_reactants_with_rings_and_limits_with_diagnostics(&[&mol], &[&rings], &limits)
            .unwrap();
        assert_eq!(report.diagnostics.accepted_matches, 1);
        assert_eq!(report.diagnostics.applied_products, 1);
        assert_eq!(products(report.products), ordinary);
        assert_eq!(
            products(
                prepared
                    .run_reactants_with_rings_and_limits(&[&mol], &[&rings], &limits)
                    .unwrap()
            ),
            ordinary
        );
        let expected = prepared
            .run_reactants_with_variant_diagnostics(&[&mol], &limits)
            .unwrap();
        let cached = prepared
            .run_reactants_with_rings_and_limits_with_variant_diagnostics(
                &[&mol],
                &[&rings],
                &limits,
            )
            .unwrap();
        assert_eq!(expected, cached);
        let traced = prepared
            .run_reactants_traced_with_diagnostics(&[&mol], &limits)
            .unwrap();
        assert_eq!(traced.diagnostics.accepted_matches, 1);
        assert!(traced.rejected_products.is_empty());
    }
}
#[test]
fn reaction_context_mismatches_return_typed_count_errors() {
    let mol = parse("CO").unwrap();
    let prepared = PreparedReaction::new("[C:1]>>[C:1]").unwrap();
    let matches = prepared.find_matches(&[&mol]).unwrap();
    let first = &matches[0];
    for error in [
        first.atom_map_positions("no reaction").unwrap_err(),
        first
            .atom_map_positions("[C:1].[O:2]>>[C:1][O:2]")
            .unwrap_err(),
        prepared.apply_match_traced(&[], first, true).err().unwrap(),
        prepared
            .apply_match_traced(
                &[&mol],
                &ReactionMatch {
                    per_reactant: vec![],
                },
                true,
            )
            .err()
            .unwrap(),
        prepared
            .find_matches_with_rings(&[&mol], &[])
            .err()
            .unwrap(),
    ] {
        let message = error.to_string();
        assert!(
            message.contains("SMIRKS") || message.contains("count mismatch"),
            "{message}"
        );
    }
    let error = prepared
        .run_reactants_with_limits(&[&mol], &ReactionTransformLimits { max_matches: 0 })
        .err()
        .unwrap();
    assert!(matches!(error, TransformError::ResourceLimit { .. }));
    assert!(error.to_string().contains("limit 0"));
}
#[test]
fn malformed_product_templates_refuse_unsupported_query_primitives() {
    for text in [
        "[C:1]>>[C;Hbogus:1]",
        "[C:1]>>[C;H3:x]",
        "[C:1]>>[C;H3:]",
        "[C:1]>>[C;H999:1]",
        "[C:1]>>(C)C",
    ] {
        let error = PreparedReaction::new(text).err().expect(text);
        assert!(error.to_string().contains("SMIRKS"), "{text}: {error}");
    }
}
#[test]
fn compatibility_refusal_codes_are_stable_for_bindings() {
    for (reason, code) in [
        (
            ReactionCompatibilityUnsupported::ChiralReactantTemplateSemantics,
            "chiral_reactant_template_semantics",
        ),
        (
            ReactionCompatibilityUnsupported::AmbiguousStereoBondOrder,
            "ambiguous_stereo_bond_order",
        ),
        (
            ReactionCompatibilityUnsupported::EzReactantTemplateSemantics,
            "ez_reactant_template_semantics",
        ),
        (
            ReactionCompatibilityUnsupported::RingModelUnresolved,
            "ring_model_unresolved",
        ),
    ] {
        assert_eq!(reason.reason_code(), code);
    }
}

#[test]
fn reaction_tetrahedral_permutations_and_carried_centers_match_pinned_rdkit() {
    let fixture: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.03.6-reaction-boundary.json"
    ))
    .unwrap();
    assert_eq!(fixture["rdkit_version"], "2026.03.6");
    let mut failures = Vec::new();
    for row in fixture["rows"].as_array().unwrap() {
        let source = row["reactant"].as_str().unwrap();
        let template = row["template"].as_str().unwrap();
        let mol = parse(source).unwrap();
        let outcome =
            run_reactants_traced_rdkit_2026_03_6(template, &[&mol], &Default::default()).unwrap();
        let RdkitProfileOutcome::Report(report) = outcome else {
            panic!("unexpected refusal: {template}");
        };
        let actual = report
            .products
            .iter()
            .map(|set| {
                set.iter()
                    .map(|p| chematic_smiles::rdkit_canonical_smiles(&p.molecule).unwrap())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let expected: Vec<Vec<String>> = serde_json::from_value(row["products"].clone()).unwrap();
        if actual != expected {
            failures.push(format!(
                "{source} with {template}: {actual:?} != {expected:?}"
            ));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn document_boundary_errors_preserve_the_authored_document() {
    use crate::document::*;
    use serde_json::json;
    let doc = ReactionDocument::from_reaction_smiles("[CH3:1][OH:2]>>[CH2:1]=[O:2]").unwrap();
    let original = serde_json::to_value(&doc).unwrap();
    let mutations = [
        ("/steps/0/id", json!("")),
        ("/steps/0/conditions", json!([{"key":"","value":"x"}])),
        (
            "/steps/0/conditions",
            json!([{"key":"pH","value":"1"},{"key":"pH","value":"2"}]),
        ),
        ("/steps/0/provenance", json!([{"source":"","kind":"test"}])),
        ("/provenance", json!([{"source":"test","kind":""}])),
        ("/steps/0/components/0/coefficient", json!(0)),
        ("/steps/0/components/0/smiles", json!("invalid")),
        (
            "/steps/0/components/0/atom_maps",
            json!([{"map_number":1,"atom_index":99}]),
        ),
        (
            "/steps/0/components/0/atom_maps",
            json!([{"map_number":1,"atom_index":0},{"map_number":1,"atom_index":0}]),
        ),
        (
            "/steps/0/components/0/atom_maps",
            json!([{"map_number":9,"atom_index":0}]),
        ),
        (
            "/steps/0/components/0/atom_maps",
            json!([{"map_number":1,"atom_index":0}]),
        ),
    ];
    for (path, value) in mutations {
        let mut changed = original.clone();
        *changed.pointer_mut(path).unwrap() = value;
        let error = ReactionDocument::from_json_str(&changed.to_string()).unwrap_err();
        assert!(!error.to_string().is_empty(), "{path}");
        assert_eq!(serde_json::to_value(&doc).unwrap(), original);
    }
    for edit in [
        json!({"op":"set_step_condition","step_id":"missing","key":"pH","value":"7"}),
        json!({"op":"set_component_coefficient","component_id":"missing","coefficient":2}),
        json!({"op":"set_document_id","id":""}),
    ] {
        assert!(doc.apply_json_edit(&edit.to_string()).is_err());
        assert_eq!(serde_json::to_value(&doc).unwrap(), original);
    }
    assert!(doc.apply_json_edit("{").is_err());
    assert!(ReactionDocument::from_json_str("{").is_err());
    assert!(ReactionDocument::from_reaction_smiles("invalid").is_err());
}
