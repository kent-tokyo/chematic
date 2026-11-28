use chematic_core::{AtomIdx, BondOrder};
use chematic_smiles::parse;

#[test]
fn wildcard_brackets_keep_isotope_charge_hydrogen_and_map() {
    let mol = parse("[13*H2-:7]").unwrap();
    let atom = mol.atom(AtomIdx(0));
    assert_eq!(atom.isotope, Some(13));
    assert_eq!(atom.charge, -1);
    assert_eq!(atom.hydrogen_count, Some(2));
    assert_eq!(atom.atom_map, Some(7));
}

#[test]
fn dangling_bonds_fail_at_chain_boundaries() {
    for smiles in ["C-", "C-.", "C=.C"] {
        assert!(parse(smiles).is_err(), "{smiles}");
    }
}

#[test]
fn explicit_bond_before_branch_is_preserved() {
    let mol = parse("C=(C)C").unwrap();
    assert_eq!(mol.bond_count(), 2);
    assert_eq!(
        mol.bonds()
            .filter(|(_, b)| b.order == BondOrder::Double)
            .count(),
        1
    );
}

#[test]
fn malformed_atoms_bonds_and_charges_are_typed_errors() {
    for smiles in [
        "C(C=)", // explicit bond closed by a branch end
        "C()",   // empty branch
        "C(=)C", // explicit bond with nothing in the branch
        "[*:1",  // wildcard bracket without ']'
        "[Xx]",  // unknown bracket element
        "C1C1",  // ring closure duplicates the chain bond
    ] {
        assert!(parse(smiles).is_err(), "{smiles}");
    }
}

#[test]
fn double_minus_charge_parses() {
    assert_eq!(parse("[O--]").unwrap().atom(AtomIdx(0)).charge, -2);
}
