//! Exposed development macrocycles; oracle is RDKit 2026.09.2 on fixed coordinates.
use chematic_ff::rdkit_uff::RdkitUffField;
use std::collections::BTreeSet;

#[test]
fn aromatic_furans_in_macrocycles_match_rdkit_uff() {
    let fixture: serde_json::Value =
        serde_json::from_str(include_str!("uff_macrocycle_aromatic_fixture.json")).unwrap();
    for case in fixture["results"].as_array().unwrap() {
        let name = case["name"].as_str().unwrap();
        let smiles = case["smiles"].as_str().unwrap();
        let mol = chematic_chem::add_hydrogens(&chematic_smiles::parse(smiles).unwrap());
        let view = chematic_perception::apply_aromaticity_rdkit_parity_experimental(&mol).unwrap();
        let pairs: BTreeSet<_> = view
            .bonds()
            .filter_map(|(_, b)| {
                (b.order == chematic_core::BondOrder::Aromatic)
                    .then_some((b.atom1.0.min(b.atom2.0), b.atom1.0.max(b.atom2.0)))
            })
            .collect();
        let expected: BTreeSet<_> = case["aromatic_bonds"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| {
                let a = v[0].as_u64().unwrap() as u32;
                let b = v[1].as_u64().unwrap() as u32;
                (a.min(b), a.max(b))
            })
            .collect();
        assert_eq!(pairs, expected, "{name} aromatic bonds");
        let coords: Vec<[f64; 3]> = serde_json::from_value(case["coords"].clone()).unwrap();
        let pos: Vec<f64> = coords.iter().flatten().copied().collect();
        let expected_gradient: Vec<f64> = serde_json::from_value(case["gradient"].clone()).unwrap();
        let field = RdkitUffField::new(&mol, &coords, 10.0, true);
        let energy = field.energy(&pos);
        let want = case["energy"].as_f64().unwrap();
        assert!(
            (energy - want).abs() < 1e-6,
            "{name} energy {energy} vs {want}"
        );
        let mut grad = vec![0.0; pos.len()];
        field.gradient(&pos, &mut grad);
        for (i, (a, b)) in grad.iter().zip(&expected_gradient).enumerate() {
            assert!((a - b).abs() < 1e-5, "{name} gradient {i}: {a} vs {b}");
        }
    }
}
