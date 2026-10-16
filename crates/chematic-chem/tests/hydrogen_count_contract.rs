//! Open Babel #3008/#2624 inspired total-H accounting, independent of notation.
use chematic_chem::{add_hydrogens, calc_mol_formula, num_hydrogens, remove_hydrogens};
use chematic_core::{Atom, AtomIdx, Element, MoleculeBuilder};

const CASES: &str = include_str!("../../../validation/openbabel_issues/hydrogen-counts/cases.tsv");

#[test]
fn explicit_bracket_inferred_isotopic_and_ionic_hydrogens_are_counted_once() {
    for row in CASES.lines().filter(|row| !row.starts_with('#')) {
        let (smiles, count) = row.split_once('\t').unwrap();
        let mol = chematic_smiles::parse(smiles).unwrap();
        let expected: usize = count.parse().unwrap();
        let before = chematic_smiles::canonical_smiles(&mol);
        assert_eq!(num_hydrogens(&mol), expected, "{smiles}");
        assert_eq!(num_hydrogens(&mol), expected, "repeat: {smiles}");
        assert_eq!(chematic_smiles::canonical_smiles(&mol), before);
        let again = chematic_smiles::parse(&before).unwrap();
        assert_eq!(
            num_hydrogens(&again),
            expected,
            "SMILES round trip: {smiles}"
        );
    }
}

#[test]
fn ordinary_hydrogen_expansion_and_folding_preserve_total_and_formula() {
    for smiles in [
        "C",
        "[CH4]",
        "[NH4+]",
        "[OH2].[OH2]",
        "CCO",
        "c1ccccc1",
        "[nH]1cccc1",
        "N[C@@H](C)O",
        "[2H]O[3H]",
        "[H-][Al+3]([H-])[H-]",
    ] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let expected = (num_hydrogens(&mol), calc_mol_formula(&mol));
        let expanded = add_hydrogens(&mol);
        let folded = remove_hydrogens(&expanded);
        for other in [&expanded, &folded] {
            assert_eq!(
                (num_hydrogens(other), calc_mol_formula(other)),
                expected,
                "{smiles}"
            );
        }
    }
}

#[test]
fn component_order_and_atom_permutations_preserve_total() {
    for smiles in ["[H].[OH2].[2H].[NH4+]", "[NH4+].[2H].[OH2].[H]"] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        assert_eq!(num_hydrogens(&mol), 8);
        let order: Vec<_> = (0..mol.atom_count())
            .rev()
            .map(|i| AtomIdx(i as u32))
            .collect();
        let reversed = mol.renumber_atoms(&order).unwrap();
        assert_eq!(num_hydrogens(&reversed), 8);
    }
}

#[test]
fn wildcard_backing_element_and_empty_graph_are_not_real_hydrogen_nodes() {
    let mut builder = MoleculeBuilder::new();
    assert_eq!(num_hydrogens(&MoleculeBuilder::new().build()), 0);
    let mut dummy = Atom::new(Element::H);
    dummy.wildcard = true;
    dummy.hydrogen_count = Some(255);
    builder.add_atom(dummy);
    assert_eq!(num_hydrogens(&builder.build()), 255);
}
