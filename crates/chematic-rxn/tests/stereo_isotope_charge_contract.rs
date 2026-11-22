use chematic_rxn::run_reactants;
use chematic_smiles::parse;

#[test]
fn chiral_reaction_template_retains_isotope_and_charge_constraints() {
    // Match counts cross-checked with RDKit 2025.09.3.
    let rule = "[13C@@H:1]([F:2])([Cl:3])[O-:4]>>[13C@@H:1]([F:2])([Cl:3])[O-:4]";
    for (smiles, expected) in [
        ("[13C@@H](F)(Cl)[O-]", 1),
        ("[12C@@H](F)(Cl)[O-]", 0),
        ("[13C@@H](F)(Cl)O", 0),
    ] {
        let molecule = parse(smiles).unwrap();
        assert_eq!(
            run_reactants(rule, &[&molecule]).unwrap().len(),
            expected,
            "{smiles}"
        );
    }
}
