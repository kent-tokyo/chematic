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

#[test]
fn shortcut_labels_affect_fingerprint_only_when_enabled() {
    let r = |label: &str| AvalonAtom {
        symbol: "R".into(),
        charge: 0,
        radical: 0,
        atext: label.into(),
    };
    let carbon = AvalonAtom {
        symbol: "C".into(),
        charge: 0,
        radical: 0,
        atext: String::new(),
    };
    let molecule = |label: &str| AvalonMolecule {
        atoms: vec![r(label), carbon.clone()],
        bonds: vec![AvalonBond {
            atoms: [1, 2],
            bond_type: 1,
        }],
    };
    // Avalon USE_SHORTCUT_LABELS bit, part of the bitFlags input.
    let without_labels = RDKIT_AVALON_DEFAULT_BIT_FLAGS & !0x800000;
    let labeled_a = avalon_fp_bytes(&molecule("Et"), 2048, RDKIT_AVALON_DEFAULT_BIT_FLAGS);
    let labeled_b = avalon_fp_bytes(&molecule("Me"), 2048, RDKIT_AVALON_DEFAULT_BIT_FLAGS);
    assert_ne!(labeled_a, labeled_b);
    assert_eq!(
        avalon_fp_bytes(&molecule("Et"), 2048, without_labels),
        avalon_fp_bytes(&molecule("Me"), 2048, without_labels)
    );
}
