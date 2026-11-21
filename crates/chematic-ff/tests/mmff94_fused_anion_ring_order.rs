use chematic_core::AtomIdx;
use chematic_ff::{assign_mmff94_numeric_types, mmff94_charges_numeric};

// RDKit 2026.09.2: the first atom ring for each N5M is its five-membered
// pyrrole ring, even though a 16-membered macrocycle contains both nitrogens.
const PORPHYRIN: &str =
    "C=CC1=C(C)c2cc3[n-]c(cc4nc(cc5[n-]c(cc1n2)c(C)c5CCC(=O)[O-])C(CCC(=O)[O-])=C4C)c(C)c3C=C";

#[test]
fn fused_anionic_nitrogens_use_the_rdkit_first_ring() {
    let mol = chematic_chem::add_hydrogens(&chematic_smiles::parse(PORPHYRIN).unwrap());
    let rings = chematic_perception::relevant_cycles(&mol, 10_000).unwrap();
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

#[test]
fn azolate_charge_sharing_matches_rdkit() {
    for (smiles, expected) in [
        ("c1cc[n-]c1", vec![-1.3]),
        ("c1c[n-]nc1", vec![-0.775, -0.775]),
        ("c1nn[n-]c1", vec![-0.65, -1.0 / 3.0, -0.65]),
    ] {
        let mol = chematic_chem::add_hydrogens(&chematic_smiles::parse(smiles).unwrap());
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        let charges = mmff94_charges_numeric(&mol).unwrap();
        let actual: Vec<_> = types
            .iter()
            .zip(charges)
            .filter_map(|(ty, q)| (*ty == 76).then_some(q))
            .collect();
        assert_eq!(actual.len(), expected.len(), "{smiles}");
        for (a, b) in actual.iter().zip(expected) {
            assert!((a - b).abs() < 1e-12, "{smiles}: {a} != {b}");
        }
    }
}

#[test]
fn fused_anion_energy_and_gradient_match_rdkit_on_common_coordinates() {
    use chematic_ff::rdkit_mmff::{MmffVariant, RdkitMmffField};
    let fixture: serde_json::Value =
        serde_json::from_str(include_str!("mmff94_fused_anion_fixture.json")).unwrap();
    let coords: Vec<[f64; 3]> = serde_json::from_value(fixture["coords"].clone()).unwrap();
    let mol = chematic_chem::add_hydrogens(&chematic_smiles::parse(PORPHYRIN).unwrap());
    for (name, variant) in [
        ("MMFF94", MmffVariant::Mmff94),
        ("MMFF94s", MmffVariant::Mmff94s),
    ] {
        let field =
            RdkitMmffField::with_terms(&mol, &coords, variant, 100.0, true, [true; 7]).unwrap();
        let expected: Vec<f64> =
            serde_json::from_value(fixture["variants"][name]["gradient"].clone()).unwrap();
        let mut gradient = vec![0.0; coords.len() * 3];
        field.gradient(coords.as_flattened(), &mut gradient);
        let energy = field.energy(coords.as_flattened());
        let reference = fixture["variants"][name]["energy"].as_f64().unwrap();
        assert!(
            (energy - reference).abs() < 1e-7,
            "{name}: {energy} != {reference}"
        );
        for (axis, (actual, expected)) in gradient.iter().zip(expected).enumerate() {
            assert!(
                (actual - expected).abs() < 1e-7,
                "{name} axis {axis}: {actual} != {expected}"
            );
        }
    }
}
