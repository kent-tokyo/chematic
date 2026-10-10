//! Fixed-version references for less common elements, rings, and explicit H.
use crate::{descriptors::*, topo_descriptors::*};
use chematic_core::MoleculeBuilder;
use serde_json::Value;
fn near(actual: f64, expected: f64, label: &str) {
    assert!(
        (actual - expected).abs() < 1e-8,
        "{label}: actual={actual}, expected={expected}"
    );
}
#[test]
fn descriptor_boundary_references_match_rdkit_2026_03_1() {
    let fixture: Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.03.1-descriptor-boundary.json"
    ))
    .unwrap();
    assert_eq!(fixture["rdkit_version"], "2026.03.1");
    for row in fixture["rows"].as_array().unwrap() {
        let text = row["smiles"].as_str().unwrap();
        let mol = if text.is_empty() {
            MoleculeBuilder::new().build()
        } else {
            chematic_smiles::parse(text).unwrap()
        };
        let mq: Vec<u32> = row["mqn"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_u64().unwrap() as u32)
            .collect();
        assert_eq!(mqn(&mol), mq, "MQN {text}");
        near(
            hall_kier_alpha(&mol),
            row["alpha"].as_f64().unwrap(),
            &format!("alpha {text}"),
        );
        for (i, f) in [chi0v, chi1v, chi2v, chi3v, chi4v].iter().enumerate() {
            near(
                f(&mol),
                row["chi_v"][i].as_f64().unwrap(),
                &format!("chi{i}v {text}"),
            );
        }
        assert_eq!(
            rdkit_num_rings(&mol),
            row["rings"].as_u64().unwrap() as usize,
            "rings {text}"
        );
        let c = carbon_types(&mol);
        let actual = [
            c.c1sp1, c.c2sp1, c.c1sp2, c.c2sp2, c.c3sp2, c.c1sp3, c.c2sp3, c.c3sp3,
        ];
        let expected: Vec<u32> = row["carbon_types"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_u64().unwrap() as u32)
            .collect();
        assert_eq!(actual.as_slice(), expected, "carbon types {text}");
        let counts = [
            num_carbons(&mol),
            num_nitrogens(&mol),
            num_oxygens(&mol),
            num_fluorines(&mol),
            num_chlorines(&mol),
            num_bromines(&mol),
            num_iodines(&mol),
            num_sulfurs(&mol),
            num_phosphorus(&mol),
        ];
        let expected: Vec<usize> = row["counts"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_u64().unwrap() as usize)
            .collect();
        assert_eq!(counts.as_slice(), expected, "element counts {text}");
    }
}
#[test]
fn empty_and_isolated_graph_descriptors_remain_defined() {
    for mol in [
        MoleculeBuilder::new().build(),
        chematic_smiles::parse("[H]").unwrap(),
        chematic_smiles::parse("[He]").unwrap(),
    ] {
        assert_eq!(fraction_rotatable_bonds(&mol), 0.0);
        assert_eq!(mde_carbon(&mol), [0.0; 10]);
        let bundle = distance_descriptor_bundle(&mol);
        assert_eq!(bundle.autocorr_2d, vec![0.0; 7]);
        assert_eq!(bundle.moran, vec![0.0; 7]);
        assert_eq!(bundle.geary, vec![1.0; 7]);
        assert_eq!(padmakar_ivan_index(&mol), 0);
        assert!(labute_asa_per_atom(&mol).iter().all(|x| x.is_finite()));
    }
}
#[test]
fn cns_mpo_piecewise_scores_respect_the_documented_boundaries() {
    // Methane has no ionizable center, so logD=logP at pH 7.4.
    let mol = chematic_smiles::parse("C").unwrap();
    near(
        cns_mpo_from_parts(&mol, 2.0, 40.0, 360.0, 0, 8.0),
        6.0,
        "maximum",
    );
    near(
        cns_mpo_from_parts(&mol, 5.0, 120.0, 500.0, 2, 10.0),
        0.0,
        "minimum",
    );
    near(
        cns_mpo_from_parts(&mol, 3.0, 20.0, 430.0, 1, 9.0),
        3.5,
        "interior",
    );
    for psa in [-1.0, 121.0] {
        near(
            cns_mpo_from_parts(&mol, 2.0, psa, 360.0, 0, 8.0),
            5.0,
            "PSA outside",
        );
    }
    near(
        cns_mpo_from_parts(&mol, 2.0, 105.0, 360.0, 0, 8.0),
        5.5,
        "PSA descending",
    );
}
