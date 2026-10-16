//! Radical hydrogen counts in Open Babel 3.2.1 V3000 output.

use chematic_mol::read_mol_v3000_with_diagnostics;
use chematic_smiles::rdkit_canonical_smiles;

#[test]
fn v3000_radical_attributes_prevent_spurious_implicit_hydrogens() {
    for (name, block, expected) in [
        (
            "aromatic carbon",
            include_str!("fixtures/openbabel-v3000-radical-aromatic.mol"),
            r"O=C1C[18CH2]CCC/C1=C\c1[c]cc(F)cc1",
        ),
        (
            "oxygen",
            include_str!("fixtures/openbabel-v3000-radical-oxygen.mol"),
            "[O]C(=O)[C@H]1C[C@@H](/C=C/CP(=O)(O)O)CCN1",
        ),
        (
            "carbon",
            include_str!("fixtures/openbabel-v3000-radical-carbon.mol"),
            "[CH2]N(C)CCOc1ccc(/C(=C(/CC)c2ccccc2)c2ccc(O)cc2)cc1",
        ),
    ] {
        let report = read_mol_v3000_with_diagnostics(block).unwrap();
        assert!(
            report
                .metadata
                .v3000_atom_properties
                .iter()
                .any(|(_, attrs)| attrs.contains("RAD=2")),
            "{name}: RAD metadata"
        );
        assert_eq!(
            rdkit_canonical_smiles(&report.mol).unwrap(),
            expected,
            "{name}"
        );
    }
}
