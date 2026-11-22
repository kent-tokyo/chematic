use chematic_smiles::{RdkitSmilesParams, parse, rdkit_cx_smiles};

#[test]
fn cx_smiles_kekule_option_writes_aromatic_molecules_as_alternating_bonds() {
    // Fixed RDKit 2025.09.3 MolToCXSmiles references with doKekule enabled.
    let params = RdkitSmilesParams {
        kekule: true,
        ..Default::default()
    };
    for (input, expected) in [
        ("c1ccccc1", "C1=CC=CC=C1"),
        ("c1ccncc1", "C1=CC=NC=C1"),
        ("c1ccccc1O", "OC1=CC=CC=C1"),
    ] {
        assert_eq!(
            rdkit_cx_smiles(&parse(input).unwrap(), &params).unwrap(),
            expected
        );
    }
}
