use chematic_fp::rdkit_rdk_fp;
use chematic_smiles::parse;

#[test]
fn rdk_fingerprint_handles_quadruple_metal_bonds() {
    // RDKit 2025.09.3 RDKFingerprint defaults, 2048 bits.
    for (smiles, expected) in [
        ("[Re]$[Re]", vec![1211, 2012]),
        ("[Mo]$[Mo]", vec![510, 537]),
    ] {
        let bits = rdkit_rdk_fp(&parse(smiles).unwrap());
        let actual: Vec<_> = (0..2048).filter(|&bit| bits.get(bit)).collect();
        assert_eq!(actual, expected, "{smiles}");
    }
}
