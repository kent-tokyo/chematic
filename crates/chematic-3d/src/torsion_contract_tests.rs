//! These are contracts of the public legacy heuristic, not evidence that its
//! preferred angles reproduce measured conformer populations.
use crate::etkdg_knowledge::{build_smarts_torsion_map, get_torsion_preference, score_torsion};
use chematic_smiles::parse;
use std::collections::HashSet;

#[test]
fn torsion_scores_are_periodic_nonnegative_and_minimal_at_the_preference() {
    for smiles in [
        "CCCC",
        "CCOC",
        "CCNC",
        "CCN(C)C",
        "CC(=O)NC",
        "CC(=O)OC",
        "CSC",
        "CSSC",
        "CS(=O)C",
        "CS(=O)(=O)C",
        "CP(C)C",
        "CCC#N",
        "CC=CNC",
        "CC=CC(=O)C",
        "CC=CCO",
        "CC=NO",
        "CNNC",
        "CN=NC",
        "CONC",
        "CC(=O)NNC",
        "CNC(=O)NC",
        "CNC(=S)NC",
        "CC(=S)SC",
        "CN(C)N(C)C",
        "CNC(=O)OC",
        "CCSCC",
        "CCOC(=O)C",
        "CCNC(=O)C",
        "CCOC(=O)NC",
        "CCS(=O)(=O)NC",
        "CCOP(=O)(O)O",
        "CCN=CC",
        "CCN=C=O",
        "CC=CC=CC",
        "CC=CC#N",
        "CC=COC",
        "CC=CSC",
        "CC(=O)C(=O)C",
        "CC(=O)SC",
        "CC(=O)N(O)C",
        "c1ccccc1-c1ccccc1",
        "c1ccccc1CC",
        "c1ccccc1C(=O)C",
        "c1ccccc1OC",
        "c1ccccc1NC",
        "c1ccccc1SC",
        "c1ccccc1N(=O)=O",
        "c1ccccc1C=C",
        "c1ccccc1-c1ncccc1",
        "Cn1cccc1",
        "CCc1ncccc1",
        "CCc1ccoc1",
        "CCc1ccsc1",
        "CCc1cc[nH]c1",
        "c1ccoc1-c1ccccc1",
        "c1ccsc1-c1ccccc1",
        "CC(=O)n1cccc1",
        "COC",
        "C[Se]CC",
        "CC[Si](C)(C)C",
        "CC(F)(F)F",
    ] {
        let mol = parse(smiles).unwrap();
        for (_, bond) in mol.bonds() {
            for (b, c) in [(bond.atom1, bond.atom2), (bond.atom2, bond.atom1)] {
                for (a, _) in mol.neighbors(b).filter(|(a, _)| *a != c) {
                    for (d, _) in mol.neighbors(c).filter(|(d, _)| *d != b && *d != a) {
                        if let Some(pref) = get_torsion_preference(&mol, a, b, c, d) {
                            assert!(pref.angle_deg.is_finite());
                            assert!(
                                pref.penalty_per_degree.is_finite() && pref.penalty_per_degree > 0.
                            );
                            assert_eq!(score_torsion(pref.angle_deg, &pref), 0.);
                            for offset in [-180., -90., -1., 0., 1., 90., 180.] {
                                let angle = pref.angle_deg + offset;
                                let score = score_torsion(angle, &pref);
                                assert!(score >= 0.);
                                assert!((score - score_torsion(angle + 360., &pref)).abs() < 1e-10);
                                assert!((score - score_torsion(angle - 360., &pref)).abs() < 1e-10);
                            }
                        }
                    }
                }
            }
        }
        let map = build_smarts_torsion_map(&mol, &HashSet::new());
        for (&(a, b), pref) in &map {
            let reverse = map.get(&(b, a)).unwrap();
            assert_eq!(pref.angle_deg, reverse.angle_deg);
            assert_eq!(pref.penalty_per_degree, reverse.penalty_per_degree);
        }
        let excluded: HashSet<_> = map.keys().copied().collect();
        assert!(build_smarts_torsion_map(&mol, &excluded).is_empty());
    }
}
