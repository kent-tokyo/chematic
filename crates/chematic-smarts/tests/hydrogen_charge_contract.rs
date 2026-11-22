use chematic_core::{Atom, Element, MoleculeBuilder};
use chematic_smarts::{find_matches, parse_smarts};

#[test]
fn hydrogen_queries_accept_numeric_and_repeated_charge_spellings() {
    // RDKit 2025.09.3 parses each pair as the same charged hydrogen query.
    for (numeric, repeated, charge) in [
        ("[H+2]", "[H++]", 2),
        ("[H-2]", "[H--]", -2),
        ("[2H+2]", "[2H++]", 2),
    ] {
        let isotope = numeric.starts_with("[2");
        for actual in [-2, 0, 2] {
            let mut atom = Atom::new(Element::H);
            atom.charge = actual;
            atom.isotope = isotope.then_some(2);
            let mut builder = MoleculeBuilder::new();
            builder.add_atom(atom);
            let molecule = builder.build();
            for pattern in [numeric, repeated] {
                let query = parse_smarts(pattern).unwrap();
                assert_eq!(
                    !find_matches(&query, &molecule).is_empty(),
                    actual == charge,
                    "{pattern}, charge {actual}"
                );
            }
        }
    }
}
