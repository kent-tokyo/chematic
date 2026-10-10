use crate::etkdg_knowledge::TorsionKnowledgeConfig;
use crate::etkdg_knowledge::build_torsion_knowledge;
use crate::etkdg_knowledge::{get_torsion_preference, score_torsion};
#[test]
fn torsion_rule_families_preserve_valid_paths_and_periodic_energy() {
    let fixture: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.03.1-descriptor-boundary.json"
    ))
    .unwrap();
    for text in fixture["rows"]
        .as_array()
        .unwrap()
        .iter()
        .step_by(2)
        .map(|r| r["smiles"].as_str().unwrap())
        .filter(|s| !s.is_empty())
        .chain([
            "CCOC(=S)NCC",
            "CCSSCC",
            "CCN=NC",
            "COOCC",
            "CNC(=S)NC",
            "CP(=O)(OC)OC",
            "C1CC2CCC1C2",
            "CC[Si](C)(C)OCC",
            "C=CC#CC",
            "CCCNC=O",
            "CCSC(=O)C",
            "CCN=S(=O)C",
            "[H]C([H])([H])C([H])([H])[H]",
            "[H]C(=O)NC",
            "[H]N(C)C(=O)C",
            "FC(=O)NC",
            "CC(=O)N(C)C",
            "CC(=O)NC(=O)C",
            "CC(=O)NOC",
            "CC(=O)NSC",
            "CC(=O)NNC",
            "NC(=O)NC",
            "NC(=S)NC",
            "CC(=O)Nc1ccccc1",
            "CC(=O)OC(=O)C",
            "COC(=O)OC",
            "COC(=S)SC",
            "CSC(=O)NC",
            "CS(=O)(=O)NC",
            "CS(=O)(=O)OC",
            "CS(=O)(=O)SC",
            "CSSSC",
            "NCCN",
            "OCCO",
            "NCCO",
            "CCN(O)O",
            "CCN=NC",
            "CN=NCC",
            "CCN=C(N)N",
            "CCN=C(C)C",
            "CC=NNC",
            "CC=NO",
            "CC=NSC",
            "CON=NC",
            "c1ccccc1NC",
            "c1ccccc1OC",
            "c1ccccc1SC",
            "c1ccccc1N(C)C",
            "c1ccccc1C(=O)C",
            "c1ccccc1C(=S)C",
            "c1ccccc1C=C(C)C",
            "c1ccccc1C#CC",
            "CC(C)=CC(=O)C",
            "CC(C)=CC(=S)C",
            "CC=CC=C",
            "FC=CF",
            "FC=CC",
            "FCCCCl",
            "ClCCBr",
            "COOC",
            "CNOC",
            "CNOO",
            "CNNC",
            "C1CCN(C)CC1",
            "C1CCOC1",
            "C1CCSC1",
            "CCP(=O)(O)OC",
            "CCP(=S)(O)OC",
            "COP(=O)(OC)OC",
        ])
    {
        let mol = chematic_smiles::parse(text).unwrap();
        for (b, _) in mol.atoms() {
            for (c, _) in mol.neighbors(b) {
                for (a, _) in mol.neighbors(b).filter(|(a, _)| *a != c) {
                    for (d, _) in mol.neighbors(c).filter(|(d, _)| *d != b && *d != a) {
                        if let Some(preference) = get_torsion_preference(&mol, a, b, c, d) {
                            assert!(preference.angle_deg.is_finite());
                            assert!(preference.penalty_per_degree > 0.0);
                            assert!(score_torsion(preference.angle_deg, &preference).abs() < 1e-12);
                            assert!(
                                (score_torsion(27.0, &preference)
                                    - score_torsion(387.0, &preference))
                                .abs()
                                    < 1e-9
                            );
                        }
                    }
                }
            }
        }
        for mask in 0..16 {
            let report = build_torsion_knowledge(
                &mol,
                &TorsionKnowledgeConfig {
                    use_exp_torsions: mask & 1 != 0,
                    use_small_ring_torsions: mask & 2 != 0,
                    use_macrocycle_torsions: mask & 4 != 0,
                    include_legacy_heuristic: mask & 8 != 0,
                    ..Default::default()
                },
            );
            let mut bonds = std::collections::HashSet::new();
            for p in &report.potentials {
                let (b, c) = p.central_bond;
                assert!(bonds.insert((b.min(c), b.max(c))));
                assert!(mol.bond_between(b, c).is_some());
                for pair in p.atoms.windows(2) {
                    assert!(mol.bond_between(pair[0], pair[1]).is_some());
                }
                assert!((p.energy(27.0) - p.energy(387.0)).abs() < 1e-8);
                let step = 1e-4;
                let numerical = (p.energy(27.0 + step) - p.energy(27.0 - step)) / (2.0 * step);
                assert!(
                    (numerical - p.d_energy_d_phi_deg(27.0)).abs() < 1e-7,
                    "{text}: {}",
                    p.rule_id
                );
            }
        }
    }
}
