//! #769: the SMILES parser rejects a neutral oxygen with more than two
//! bonds, but a SMIRKS template is a pattern, not a molecule. BioTransformer
//! rule BTMR0839 writes its product oxygen as `[#8;A;H1X2:1]` between two
//! ring atoms; RDKit reads the template and finds no match on most
//! reactants. The template must still parse, and a product that would carry
//! the impossible oxygen is refused when the template is applied.

use chematic_rxn::PreparedReaction;

const BTMR0839: &str = "[#8;A;H1X2:1][#6:2]!@-[#6:3]@-[#6:4]!@-[#6:5][#8;A;H1X2]>>[#6:3]-1@-[#6:4]-[#6:5][#8;A;H1X2:1][#6:2]-1";

#[test]
fn an_impossible_product_atom_does_not_stop_the_template_parsing() {
    let prepared = PreparedReaction::shared(BTMR0839).expect("the template parses");
    let ethanol = chematic_smiles::parse("CCO").unwrap();
    assert!(prepared.run_reactants(&[&ethanol]).unwrap().is_empty());

    // A diol the reactant side matches: no product with a three-coordinate
    // neutral oxygen is returned.
    let diol = chematic_smiles::parse("OCC1CCCC1CO").unwrap();
    if let Ok(products) = prepared.run_reactants(&[&diol]) {
        for set in products {
            for product in set {
                let smiles = chematic_smiles::write(&product);
                assert!(chematic_smiles::parse(&smiles).is_ok(), "{smiles}");
            }
        }
    }
}

#[test]
fn reaction_smiles_of_molecules_keeps_the_valence_check() {
    assert!(chematic_rxn::parse_reaction("CO(C)C>>CO").is_err());
    assert!(chematic_rxn::parse_reaction("COC>>CO").is_ok());
}
