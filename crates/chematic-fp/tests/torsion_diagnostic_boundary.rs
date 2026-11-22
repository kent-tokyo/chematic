use chematic_core::AtomIdx;
use chematic_fp::{
    RdkitTorsionDiagnosticError as Error, explain_rdkit_topological_torsion_code,
    rdkit_topological_torsion_path_code,
};

#[test]
fn heterogeneous_torsion_path_is_direction_independent() {
    let molecule = chematic_smiles::parse("CCNO").unwrap();
    let forward = [AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3)];
    let reverse = [AtomIdx(3), AtomIdx(2), AtomIdx(1), AtomIdx(0)];
    let code = rdkit_topological_torsion_path_code(&molecule, &forward).unwrap();
    assert_eq!(
        code,
        rdkit_topological_torsion_path_code(&molecule, &reverse).unwrap()
    );
    let explanation = explain_rdkit_topological_torsion_code(code, 4).unwrap();
    let symbols: Vec<_> = explanation.iter().map(|atom| atom.symbol).collect();
    assert!(symbols == ["C", "C", "N", "O"] || symbols == ["O", "N", "C", "C"]);
}

#[test]
fn diagnostic_path_lengths_are_bounded_before_packing_or_decoding() {
    let molecule = chematic_smiles::parse("CC").unwrap();
    assert!(matches!(
        rdkit_topological_torsion_path_code(&molecule, &[AtomIdx(0); 8]),
        Err(Error::PathTooLong {
            size: 8,
            maximum: 7
        })
    ));
    assert_eq!(
        explain_rdkit_topological_torsion_code(0, 0),
        Err(Error::EmptyPath)
    );
}
