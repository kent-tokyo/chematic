use chematic_core::AtomIdx;
use chematic_ff::{assign_mmff94_numeric_types, mmff94_charges_numeric};

// RDKit 2026.09.2: the first atom ring for each N5M is its five-membered
// pyrrole ring, even though a 16-membered macrocycle contains both nitrogens.
const PORPHYRIN: &str =
    "C=CC1=C(C)c2cc3[n-]c(cc4nc(cc5[n-]c(cc1n2)c(C)c5CCC(=O)[O-])C(CCC(=O)[O-])=C4C)c(C)c3C=C";

#[test]
fn fused_anionic_nitrogens_use_the_rdkit_first_ring() {
    let mol = chematic_chem::add_hydrogens(&chematic_smiles::parse(PORPHYRIN).unwrap());
    let rings = chematic_perception::rdkit_sssr_ring_order(&mol).unwrap_or_else(|| {
        chematic_perception::find_symmetrized_sssr(&mol)
            .rings()
            .to_vec()
    });
    for index in [8, 16] {
        let sizes: Vec<_> = rings
            .iter()
            .filter(|ring| ring.contains(&AtomIdx(index)))
            .map(Vec::len)
            .collect();
        assert_eq!(
            sizes.first(),
            Some(&5),
            "atom {index}, ring sizes {sizes:?}"
        );
    }
    let types = assign_mmff94_numeric_types(&mol).unwrap();
    let charges = mmff94_charges_numeric(&mol).unwrap();
    for index in [8, 16] {
        assert_eq!(types[index], 76);
        assert!(
            (charges[index] + 1.3).abs() < 1e-12,
            "atom {index}: {}",
            charges[index]
        );
    }
    assert!((charges.iter().sum::<f64>() + 6.0).abs() < 1e-12);
}
