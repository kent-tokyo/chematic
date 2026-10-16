use chematic_core::{AtomIdx, implicit_hcount};
use chematic_mol::{parse_mol2_record, parse_mol2_source_record};

const PARTIAL_ETHANE: &str = "\
@<TRIPOS>MOLECULE
partial-ethane
7 6 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
1 C1 0 0 0 C.3 1 MOL 0
2 C2 1 0 0 C.3 1 MOL 0
3 H1 0 1 0 H 1 MOL 0
4 H2 0 -1 0 H 1 MOL 0
5 H3 1 1 0 H 1 MOL 0
6 H4 1 -1 0 H 1 MOL 0
7 H5 2 0 0 H 1 MOL 0
@<TRIPOS>BOND
1 1 2 1
2 1 3 1
3 1 4 1
4 2 5 1
5 2 6 1
6 2 7 1
";

const IMPLICIT_ETHANE: &str = "\
@<TRIPOS>MOLECULE
implicit-ethane
2 1 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
1 C1 0 0 0 C.3 1 MOL 0
2 C2 1 0 0 C.3 1 MOL 0
@<TRIPOS>BOND
1 1 2 1
";

#[test]
fn explicit_h_atoms_define_the_hydrogen_inventory_for_both_mol2_readers() {
    let ordinary = parse_mol2_record(PARTIAL_ETHANE).unwrap();
    let retained = parse_mol2_source_record(PARTIAL_ETHANE).unwrap();
    for mol in [&ordinary.molecule, retained.molecule()] {
        assert_eq!(mol.atom(AtomIdx(0)).hydrogen_count, Some(0));
        assert_eq!(mol.atom(AtomIdx(1)).hydrogen_count, Some(0));
        assert_eq!(implicit_hcount(mol, AtomIdx(0)), 0);
        assert_eq!(implicit_hcount(mol, AtomIdx(1)), 0);
    }
}

#[test]
fn no_explicit_h_atoms_retains_normal_valence_inference() {
    let ordinary = parse_mol2_record(IMPLICIT_ETHANE).unwrap();
    let retained = parse_mol2_source_record(IMPLICIT_ETHANE).unwrap();
    for mol in [&ordinary.molecule, retained.molecule()] {
        assert_eq!(mol.atom(AtomIdx(0)).hydrogen_count, None);
        assert_eq!(mol.atom(AtomIdx(1)).hydrogen_count, None);
        assert_eq!(implicit_hcount(mol, AtomIdx(0)), 3);
        assert_eq!(implicit_hcount(mol, AtomIdx(1)), 3);
    }
}
