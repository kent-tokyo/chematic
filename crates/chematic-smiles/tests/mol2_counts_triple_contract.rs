use chematic_core::{BondIdx, BondOrder};
use chematic_smiles::rdkit_mol_from_mol2_block;

#[test]
fn mol2_triple_bond_retains_order() {
    // RDKit 2025.09.3 returns two atoms connected by a triple bond.
    let block = "@<TRIPOS>MOLECULE\nfixture\n2 1\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.1 1 RES 0\n2 N1 1 0 0 N.1 1 RES 0\n@<TRIPOS>BOND\n1 1 2 3\n";
    let parsed = rdkit_mol_from_mol2_block(block, false, false, false).unwrap();
    assert_eq!(parsed.molecule.atom_count(), 2);
    assert_eq!(parsed.molecule.bond_count(), 1);
    assert_eq!(parsed.molecule.bond(BondIdx(0)).order, BondOrder::Triple);
}

#[test]
fn mol2_atom_count_without_bond_count_defaults_to_no_bonds() {
    // RDKit 2025.09.3 accepts this counts line as a one-atom molecule.
    let block =
        "@<TRIPOS>MOLECULE\nfixture\n1\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3 1 RES 0\n";
    let parsed = rdkit_mol_from_mol2_block(block, false, false, false).unwrap();
    assert_eq!(parsed.molecule.atom_count(), 1);
    assert_eq!(parsed.molecule.bond_count(), 0);
}

#[test]
fn mol2_rejects_empty_counts_and_truncated_atom_records() {
    for (block, expected) in [
        (
            "@<TRIPOS>MOLECULE\nfixture\n\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3\n",
            "Empty counts line",
        ),
        (
            "@<TRIPOS>MOLECULE\nfixture\n2 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3\n",
            "premature EOF",
        ),
    ] {
        let error = match rdkit_mol_from_mol2_block(block, false, false, false) {
            Ok(_) => panic!("malformed MOL2 must fail"),
            Err(error) => error,
        };
        assert!(error.to_string().contains(expected), "{error}");
    }
}
