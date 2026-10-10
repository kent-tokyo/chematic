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
