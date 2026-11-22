use chematic_3d::etkdg_knowledge::get_torsion_preference;
use chematic_core::AtomIdx;

#[test]
fn legacy_torsion_preferences_preserve_nitrile_and_phosphate_rules() {
    // These are compatibility checks for the legacy heuristic, not measured energies.
    let cases = [
        ("CC[N+]#[C-]", [0, 1, 2, 3], 180., 0.15),
        ("c1ccccc1[N+]#[C-]", [4, 5, 6, 7], 0., 0.10),
        ("CC(=O)[N+]#[C-]", [0, 1, 3, 4], 180., 0.18),
        ("CP(=O)(O)OC", [0, 1, 4, 5], 180., 0.06),
        ("CP(=O)(O)OC", [5, 4, 1, 0], 180., 0.06),
        ("CC(=O)C(C)=O", [0, 1, 3, 4], 0., 0.10),
    ];
    for (smiles, indices, angle, penalty) in cases {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let [a, b, c, d] = indices.map(AtomIdx);
        let preference = get_torsion_preference(&mol, a, b, c, d).unwrap();
        assert_eq!(preference.angle_deg, angle, "{smiles}");
        assert_eq!(preference.penalty_per_degree, penalty, "{smiles}");
    }
}
