use chematic_fp::rdkit_avalon::{
    AvalonAtom, AvalonBond, AvalonMolecule, RDKIT_AVALON_DEFAULT_BIT_FLAGS, avalon_fp_bytes,
};

#[test]
fn cyclic_fingerprint_is_stable_under_reordering_bond_records() {
    let atom = AvalonAtom {
        symbol: "C".into(),
        charge: 0,
        radical: 0,
        atext: String::new(),
    };
    let mut molecule = AvalonMolecule {
        atoms: vec![atom; 6],
        bonds: [[1, 2], [3, 4], [2, 3], [4, 5], [5, 6], [6, 1], [2, 5]]
            .into_iter()
            .map(|atoms| AvalonBond {
                atoms,
                bond_type: 1,
            })
            .collect(),
    };
    let first = avalon_fp_bytes(&molecule, 512, RDKIT_AVALON_DEFAULT_BIT_FLAGS);
    molecule.bonds.reverse();
    let reordered = avalon_fp_bytes(&molecule, 512, RDKIT_AVALON_DEFAULT_BIT_FLAGS);
    assert!(first.iter().any(|&byte| byte != 0));
    assert_eq!(first, reordered);
}
