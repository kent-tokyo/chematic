use chematic_core::BondOrder;
use chematic_smiles::{parse, rdkit_sanitized_model};

#[test]
fn sanitized_model_keeps_quadruple_metal_bonds() {
    // RDKit 2025.09.3 reads both as one QUADRUPLE bond.
    for smiles in ["[Re]$[Re]", "[Mo]$[Mo]"] {
        let molecule = parse(smiles).unwrap();
        let model = rdkit_sanitized_model(&molecule).unwrap();
        assert_eq!(model.bonds, [Some(BondOrder::Quadruple)], "{smiles}");
        assert_eq!(model.atoms.len(), 2);
    }
}
