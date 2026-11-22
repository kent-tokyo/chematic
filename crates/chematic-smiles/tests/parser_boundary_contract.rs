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
