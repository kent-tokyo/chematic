use chematic_core::{AtomIdx, Element};
use chematic_ff::mmff94_energy::mmff94_vdw_combined;
use chematic_ff::{assign_mmff94_numeric_types, mmff94_charges_numeric};
use chematic_smiles::parse;

// Oracle: RDKit 2026.03.6, MMFFGetMoleculeProperties on explicit-H molecules.
// Source: Release_2026_03_6/Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp,
// setMMFFHydrogenType. In particular, 29 is HOCC/HOCN and 33 is HOS.
#[test]
fn oxygen_bound_hydrogen_types_and_charges_match_their_environment() {
    for (smiles, expected_type, expected_charge) in [
        ("CCO[H]", 21, 0.4),
        ("[H]ON(C)C", 21, 0.4),
        ("c1ccccc1O[H]", 29, 0.45),
        ("[H]Oc1ccccc1", 29, 0.45),
        ("[H]OC=C", 29, 0.45),
        ("[H]OC=N", 29, 0.45),
        ("CC(=O)O[H]", 24, 0.5),
        ("[H]OP(=O)(O)O", 24, 0.5),
        ("[H]OS(=O)(=O)O", 33, 0.5),
    ] {
        let mol = parse(smiles).unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        let charges = mmff94_charges_numeric(&mol).unwrap();
        let mut checked = 0;
        for i in 0..mol.atom_count() {
            let idx = AtomIdx(i as u32);
            if mol.atom(idx).element == Element::H
                && mol
                    .neighbors(idx)
                    .any(|(nb, _)| mol.atom(nb).element == Element::O)
            {
                assert_eq!(types[i], expected_type, "{smiles}: H atom {i}");
                assert!(
                    (charges[i] - expected_charge).abs() < 1e-12,
                    "{smiles}: H atom {i} charge {} != {expected_charge}",
                    charges[i]
                );
                checked += 1;
            }
        }
        assert!(checked > 0, "fixture must contain explicit O-H: {smiles}");
    }
}

// Oracle: RDKit 2026.03.6 GetMMFFVdWParams (scaled radius/well depth).
// Source: Release_2026_03_6/Code/ForceField/MMFF/Nonbonded.cpp,
// calcUnscaledVdWMinimum: ANY donor suppresses the asymmetry correction,
// including donor/nonacceptor pairs. Donor/acceptor scaling follows afterward.
#[test]
fn donor_vdw_pairs_use_the_rdkit_combining_rule_in_both_orders() {
    for (a, b, expected_radius, expected_epsilon) in [
        (3, 23, 3.2988003775750503, 0.03275316289992711),
        (1, 23, 3.275768748270286, 0.03306703754626389),
        (7, 23, 2.44265298442308, 0.019039914981539462),
        (8, 23, 2.6568477686322076, 0.016843816258122306),
        (5, 5, 2.9698484809835, 0.021574173607307283),
        (1, 5, 3.598787981153042, 0.02807765636226562),
    ] {
        for (i, j) in [(a, b), (b, a)] {
            let (radius, epsilon) = mmff94_vdw_combined(i, j).unwrap();
            assert!(
                (radius - expected_radius).abs() < 1e-12,
                "types {i}/{j}: radius {radius} != {expected_radius}"
            );
            assert!(
                (epsilon - expected_epsilon).abs() < 1e-12,
                "types {i}/{j}: epsilon {epsilon} != {expected_epsilon}"
            );
        }
    }
}
