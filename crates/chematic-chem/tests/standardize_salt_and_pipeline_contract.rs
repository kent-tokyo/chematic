//! Salt heuristics and pipeline stage names not reached by other tests.
use chematic_chem::standardize::{SaltCatalog, remove_salts_with_catalog};
use chematic_chem::{
    StandardizationPipeline, StandardizationStep, StandardizeOptions, standardize,
};
use chematic_smiles::{canonical_smiles, parse};

fn stripped(smi: &str) -> String {
    canonical_smiles(&remove_salts_with_catalog(
        &parse(smi).unwrap(),
        &SaltCatalog::new(),
    ))
}

#[test]
fn small_carbon_free_fragment_is_treated_as_salt() {
    // Neutral, bonded `[Na]Cl` matches no catalog pattern (they are charged
    // ions); the carbon-free small-fragment heuristic still strips it.
    assert_eq!(
        stripped("[Na]Cl.CCCCO"),
        canonical_smiles(&parse("CCCCO").unwrap())
    );
    assert_eq!(
        stripped("ON=O.CCCCO"),
        canonical_smiles(&parse("CCCCO").unwrap())
    );
}

#[test]
fn carbon_containing_small_fragment_is_kept_when_largest() {
    assert_eq!(
        stripped("CO.[Na]Cl"),
        canonical_smiles(&parse("CO").unwrap())
    );
}

#[test]
fn step_names_are_unique_snake_case() {
    let steps = [
        StandardizationStep::DisconnectMetals,
        StandardizationStep::LargestFragment,
        StandardizationStep::NeutralizeCharges,
        StandardizationStep::NormalizeGroups,
        StandardizationStep::ZwitterionNormalization,
        StandardizationStep::RemoveExplicitHydrogens,
        StandardizationStep::CanonicalTautomer,
        StandardizationStep::FragmentParent,
        StandardizationStep::ChargeParent,
        StandardizationStep::IsotopeParent,
        StandardizationStep::StereoParent,
    ];
    let mut names: Vec<&str> = steps.iter().map(|s| s.as_str()).collect();
    assert!(
        names
            .iter()
            .all(|n| n.chars().all(|c| c.is_ascii_lowercase() || c == '_'))
    );
    names.sort_unstable();
    names.dedup();
    assert_eq!(names.len(), steps.len());
}

#[test]
fn disconnect_metals_keeps_directional_bonds() {
    let mol = parse("F/C=C/O[Na]").unwrap();
    let out = standardize(&mol, &StandardizeOptions::default());
    let smi = canonical_smiles(&out);
    assert!(smi.contains('/') || smi.contains('\\'), "{smi}");
    assert!(!smi.contains("Na") || smi.contains('.'), "{smi}");
    let (_, report) = StandardizationPipeline::new(StandardizeOptions::default()).run(&mol);
    assert!(
        report
            .steps
            .iter()
            .any(|s| s.step == StandardizationStep::DisconnectMetals)
    );
}
