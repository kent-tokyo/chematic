use chematic_rxn::{RdkitProfileOutcome, run_reactants_traced_rdkit_2026_03_6};

#[test]
fn disconnected_template_keeps_carried_chain_stereochemistry() {
    // RDKit 2025.09.3 RunReactants returns the same absolute configuration.
    let smiles = "FCCC[C@](Cl)(Br)I";
    let smirks = "([F:1].[C:2]([Cl:3])[Br:4])>>([F:1].[C:2]([Cl:3])[Br:4])";
    let molecule = chematic_smiles::parse(smiles).unwrap();
    let result =
        run_reactants_traced_rdkit_2026_03_6(smirks, &[&molecule], &Default::default()).unwrap();
    let RdkitProfileOutcome::Report(report) = result else {
        panic!("fixture must be supported")
    };
    assert_eq!(report.products.len(), 1);
    assert_eq!(report.products[0].len(), 1);
    let expected = chematic_smiles::canonical_smiles(&molecule);
    assert_eq!(
        chematic_smiles::canonical_smiles(&report.products[0][0].molecule),
        expected
    );
}
