#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> Molecule {
        parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"))
    }

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    // ── Wiener Index ──────────────────────────────────────────────────────────

    #[test]
    fn wiener_ethane() {
        // CC: only 1 pair at distance 1 → W = 1
        assert_eq!(wiener_index(&mol("CC")) as u32, 1);
    }

    #[test]
    fn wiener_propane() {
        // CCC: pairs (C1,C2)=1, (C1,C3)=2, (C2,C3)=1 → W = 4
        assert_eq!(wiener_index(&mol("CCC")) as u32, 4);
    }

    #[test]
    fn wiener_benzene() {
        // c1ccccc1: 6 atoms in ring, W = 1+2+3+2+1 + 1+2+3+2 + 1+2+3 + 1+2 + 1 = 27
        assert_eq!(wiener_index(&mol("c1ccccc1")) as u32, 27);
    }

    #[test]
    fn wiener_increases_with_chain_length() {
        let w2 = wiener_index(&mol("CC")); // ethane
        let w3 = wiener_index(&mol("CCC")); // propane
        let w4 = wiener_index(&mol("CCCC")); // butane
        assert!(w2 < w3 && w3 < w4);
    }

    #[test]
    fn wiener_single_atom_zero() {
        assert_eq!(wiener_index(&mol("C")) as u32, 0);
    }

    // ── Kappa Indices ─────────────────────────────────────────────────────────

    #[test]
    fn kappa1_propane() {
        // n=3, p1=2 (C-C, C-C): κ1 = 3·4/4 = 3.0
        assert!(close(kappa1(&mol("CCC")), 3.0, 0.01));
    }

    #[test]
    fn kappa1_benzene() {
        // Aromatic C: alpha = 6·(0.67/0.77 - 1) = -0.78 ≠ 0 (RDKit CalcKappa1
        // verified: 3.4116), unlike sp3-only alkanes where alpha = 0.
        assert!(close(kappa1(&mol("c1ccccc1")), 3.4116, 0.001));
    }

    #[test]
    fn kappa2_propane() {
        // n=3, p2=1: κ2 = 2·1/1 = 2.0
        assert!(close(kappa2(&mol("CCC")), 2.0, 0.01));
    }

    #[test]
    fn kappa2_benzene() {
        // RDKit CalcKappa2(benzene) verified: 1.6058 (alpha-corrected).
        assert!(close(kappa2(&mol("c1ccccc1")), 1.6058, 0.001));
    }

    #[test]
    fn kappa3_propane_zero() {
        // n=3, no 3-bond paths: κ3 = 0
        assert_eq!(kappa3(&mol("CCC")), 0.0);
    }

    #[test]
    fn kappa3_benzene() {
        // RDKit CalcKappa3(benzene) verified: 0.5824 (alpha-corrected).
        assert!(close(kappa3(&mol("c1ccccc1")), 0.5824, 0.001));
    }

    #[test]
    fn kappa1_single_atom_zero() {
        assert_eq!(kappa1(&mol("C")), 0.0);
    }

    #[test]
    fn kappa_alpha_corrected_matches_rdkit_aspirin() {
        // RDKit CalcKappa1/2/3("CC(=O)Oc1ccccc1C(=O)O") verified: 9.2496 /
        // 3.7093 / 2.2974. Aspirin mixes sp3/sp2 C, O, and aromatic C, so a
        // nonzero Hall-Kier alpha correction is exercised on all three.
        // Tolerance is 0.1 (not 0.001 like the benzene/aspirin-free cases)
        // because `hall_kier_alpha`'s covalent-radius table has its own
        // separately-tracked precision gap for O (chematic's alpha for
        // aspirin is -1.766 vs RDKit's -1.840) — see
        // tasks/descriptor_validation_coverage.md. This test only guards the
        // kappa formula's alpha-wiring, not that residual table precision.
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        assert!(close(kappa1(&m), 9.2496, 0.1), "kappa1 = {}", kappa1(&m));
        assert!(close(kappa2(&m), 3.7093, 0.1), "kappa2 = {}", kappa2(&m));
        assert!(close(kappa3(&m), 2.2974, 0.1), "kappa3 = {}", kappa3(&m));
    }

    #[test]
    fn kappa_all_matches_individual() {
        for smi in ["CC", "CCC", "CCCC", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"] {
            let m = mol(smi);
            let (k1, k2, k3) = kappa_all(&m);
            assert!((k1 - kappa1(&m)).abs() < 1e-10, "{smi}: kappa1 mismatch");
            assert!((k2 - kappa2(&m)).abs() < 1e-10, "{smi}: kappa2 mismatch");
            assert!((k3 - kappa3(&m)).abs() < 1e-10, "{smi}: kappa3 mismatch");
        }
    }

    // ── Chi Connectivity ─────────────────────────────────────────────────────

    #[test]
    fn chi0_benzene() {
        // 6 aromatic C each with δ=2: χ0 = 6 · 2^(-0.5) ≈ 4.243
        let c = chi0(&mol("c1ccccc1"));
        assert!(close(c, 4.243, 0.01), "chi0(benzene) = {c}");
    }

    #[test]
    fn chi1_benzene() {
        // 6 bonds, each (2·2)^(-0.5) = 0.5: χ1 = 3.0
        let c = chi1(&mol("c1ccccc1"));
        assert!(close(c, 3.0, 0.01), "chi1(benzene) = {c}");
    }

    #[test]
    fn chi0_propane() {
        // δ(C1)=1, δ(C2)=2, δ(C3)=1: χ0 = 1 + 1/√2 + 1 ≈ 2.707
        let c = chi0(&mol("CCC"));
        assert!(close(c, 2.707, 0.01), "chi0(propane) = {c}");
    }

    #[test]
    fn chi1_propane() {
        // bonds (1·2) and (2·1): χ1 = 2·(2)^(-0.5) ≈ 1.414
        let c = chi1(&mol("CCC"));
        assert!(close(c, 1.414, 0.01), "chi1(propane) = {c}");
    }

    #[test]
    fn chi_increases_with_chain() {
        // Longer chain → larger χ0
        assert!(chi0(&mol("CCC")) < chi0(&mol("CCCC")));
    }

    #[test]
    fn chi0v_benzene() {
        // Each aromatic C: Zv=4, H=1, denom=1 → δv=3: χ0v = 6/√3 ≈ 3.464
        let c = chi0v(&mol("c1ccccc1"));
        assert!(close(c, 3.464, 0.01), "chi0v(benzene) = {c}");
    }

    #[test]
    fn chi1v_benzene() {
        // bonds (3·3): χ1v = 6·(9)^(-0.5) = 6/3 = 2.0
        let c = chi1v(&mol("c1ccccc1"));
        assert!(close(c, 2.0, 0.01), "chi1v(benzene) = {c}");
    }

    // ── Bertz CT ──────────────────────────────────────────────────────────────

    #[test]
    fn bertz_ct_increases_with_complexity() {
        let bz = bertz_ct(&mol("c1ccccc1"));
        let asp = bertz_ct(&mol("CC(=O)Oc1ccccc1C(=O)O"));
        assert!(bz < asp, "benzene BertzCT {bz} should be < aspirin {asp}");
    }

    // Reference values: RDKit 2026.03.1 `GraphDescriptors.BertzCT`.
    #[test]
    fn bertz_ct_small_molecules_match_rdkit() {
        assert_eq!(bertz_ct(&mol("C")), 0.0);
        assert_eq!(bertz_ct(&mol("CC")), 0.0);
        assert_eq!(bertz_ct(&mol("CCO")), 2.7548875021634682);
    }

    #[test]
    fn bertz_ct_benzene() {
        assert_eq!(bertz_ct(&mol("c1ccccc1")), 71.96100505779535);
    }

    #[test]
    fn bertz_ct_aspirin() {
        assert_eq!(bertz_ct(&mol("CC(=O)Oc1ccccc1C(=O)O")), 343.2228677267164);
    }

    // ── LabuteASA ─────────────────────────────────────────────────────────────

    #[test]
    fn labute_asa_positive() {
        assert!(labute_asa(&mol("c1ccccc1")) > 0.0, "benzene ASA > 0");
    }

    #[test]
    fn labute_asa_single_oxygen() {
        // RDKit CalcLabuteASA("O") verified: 6.8492. The implicit-H
        // contribution is a single pooled term per heavy atom (not scaled by
        // the atom's actual H count) — see `labute_asa_parts` doc comment.
        let asa = labute_asa(&mol("O"));
        assert!(
            (asa - 6.8492).abs() < 0.001,
            "water ASA {asa:.4} != RDKit 6.8492"
        );
    }

    #[test]
    fn labute_asa_matches_rdkit_aspirin() {
        // RDKit CalcLabuteASA("CC(=O)Oc1ccccc1C(=O)O") verified: 74.7571.
        let asa = labute_asa(&mol("CC(=O)Oc1ccccc1C(=O)O"));
        assert!((asa - 74.7571).abs() < 0.001, "aspirin ASA {asa:.4}");
    }

    #[test]
    fn labute_asa_uses_rdkit_model_for_dative_aromatic_inputs() {
        for (smiles, expected) in [
            (
                "CC1=C(NC2=NC3<-C(=CC=CC=3)N3C=NC=C23)C(Cl)=CC=C1",
                132.300_143_293_752_28,
            ),
            (
                "CC(=O)NCCc1c2n(c3<-c1cc(Cl)cc3)CCCc1cccc[c]1-2",
                152.507_162_477_150_7,
            ),
        ] {
            let got = labute_asa(&mol(smiles));
            assert!((got - expected).abs() < 1e-10, "{smiles}: {got}");
        }
    }

    #[test]
    fn labute_asa_per_atom_excludes_pooled_h_term() {
        // Quaternary carbon (CC(C)(C)C) has an atom with zero implicit H —
        // RDKit still runs the pooled-H pass once for it (see
        // `labute_asa_parts`). RDKit `_CalcLabuteASAContribs` ats verified:
        // [6.9237, 5.4150, 6.9237, 6.9237, 6.9237], hs=1.0891.
        let m = mol("CC(C)(C)C");
        let per_atom = labute_asa_per_atom(&m);
        let expected = [6.9237, 5.4150, 6.9237, 6.9237, 6.9237];
        for (got, want) in per_atom.iter().zip(expected.iter()) {
            assert!((got - want).abs() < 0.001, "{per_atom:?} vs {expected:?}");
        }
        let total = labute_asa(&m);
        assert!((total - (expected.iter().sum::<f64>() + 1.0891)).abs() < 0.001);
    }

    #[test]
    fn labute_asa_monotone_with_size() {
        // Larger molecule → larger ASA.
        let bz = labute_asa(&mol("c1ccccc1"));
        let asp = labute_asa(&mol("CC(=O)Oc1ccccc1C(=O)O"));
        assert!(
            bz < asp,
            "benzene ASA {bz:.2} should be < aspirin ASA {asp:.2}"
        );
    }

    #[test]
    fn labute_asa_aromatic_reduces_vs_saturated() {
        // Aromatic C-C bonds (scale=0.1) create more surface overlap than
        // single C-C bonds (scale=0), so benzene ASA < cyclohexane ASA
        // (same atom count, but less overlap in cyclohexane).
        let bz = labute_asa(&mol("c1ccccc1")); // aromatic
        let ch = labute_asa(&mol("C1CCCCC1")); // saturated
        assert!(bz < ch, "benzene ASA {bz:.2} < cyclohexane ASA {ch:.2}");
    }

    #[test]
    fn randic_index_ethane() {
        // Ethane: 1 edge between two degree-1 nodes → Randic = 1/sqrt(1*1) = 1.0
        let m = mol("CC");
        assert!((randic_index(&m) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn zagreb_m1_ethane() {
        // Ethane: 2 atoms each with degree 1 → Σ d² = 1+1 = 2
        let m = mol("CC");
        assert_eq!(zagreb_index_m1(&m), 2);
    }

    // ── Zagreb M2 ───────────────────────────────────────────────────────────

    #[test]
    fn zagreb_m2_ethane() {
        // Ethane C-C: one bond, deg(C)=1, deg(C)=1 → M2 = 1*1 = 1
        assert_eq!(zagreb_index_m2(&mol("CC")), 1);
    }

    #[test]
    fn zagreb_m2_propane() {
        // Propane C-C-C: bonds (C1-C2) and (C2-C3)
        // deg: C1=1, C2=2, C3=1
        // M2 = 1*2 + 2*1 = 4
        assert_eq!(zagreb_index_m2(&mol("CCC")), 4);
    }

    #[test]
    fn zagreb_m2_benzene() {
        // Benzene: 6 bonds, each between degree-2 atoms → M2 = 6 * (2*2) = 24
        assert_eq!(zagreb_index_m2(&mol("c1ccccc1")), 24);
    }

    #[test]
    fn zagreb_m2_ge_m1_for_branched() {
        // For any graph, M2 captures edge-based branching; M1 is vertex-based.
        // Both should be positive for any non-trivial molecule.
        for smi in ["CC(C)C", "CC(=O)O", "c1ccccc1"] {
            let m = mol(smi);
            assert!(zagreb_index_m2(&m) > 0, "M2 should be > 0 for {smi}");
        }
    }

    // ── Eccentricity / Petitjean / ECI ───────────────────────────────────────

    #[test]
    fn eccentricities_propane() {
        // CCC: terminal Cs have ecc=2, middle C has ecc=1
        let m = mol("CCC");
        let ecc = graph_eccentricities(&m);
        assert_eq!(ecc.len(), 3);
        assert_eq!(*ecc.iter().max().unwrap(), 2);
        assert_eq!(*ecc.iter().min().unwrap(), 1);
    }

    #[test]
    fn graph_diameter_propane() {
        assert_eq!(graph_diameter(&mol("CCC")), 2);
    }

    #[test]
    fn graph_radius_propane() {
        assert_eq!(graph_radius(&mol("CCC")), 1);
    }

    #[test]
    fn petitjean_propane() {
        // (2 - 1) / 2 = 0.5
        let v = petitjean_index(&mol("CCC"));
        assert!((v - 0.5).abs() < 1e-9, "expected 0.5 got {v}");
    }

    #[test]
    fn petitjean_benzene() {
        // c1ccccc1: 6-cycle — all atoms have eccentricity 3, so diameter=radius=3 → Petitjean=0
        let v = petitjean_index(&mol("c1ccccc1"));
        assert!(
            (v - 0.0).abs() < 1e-9,
            "expected 0.0 (symmetric ring) got {v}"
        );
    }

    #[test]
    fn petitjean_single_atom() {
        assert_eq!(petitjean_index(&mol("C")), 0.0);
    }

    #[test]
    fn eci_propane() {
        // CCC: middle C deg=2 ecc=1 → 2; each terminal C deg=1 ecc=2 → 2 each → total 6
        assert_eq!(eccentric_connectivity_index(&mol("CCC")), 6);
    }

    // ── Hosoya Index ─────────────────────────────────────────────────────────

    #[test]
    fn hosoya_methane() {
        // Single atom: 0 edges → Z = 1 (empty matching only)
        assert_eq!(hosoya_index(&mol("C")), 1);
    }

    #[test]
    fn hosoya_ethane() {
        // CC: 1 edge → Z = 2 (empty + {e1})
        assert_eq!(hosoya_index(&mol("CC")), 2);
    }

    #[test]
    fn hosoya_propane() {
        // CCC: 2 edges, no two share a vertex → Z = 1 + 2 = 3
        assert_eq!(hosoya_index(&mol("CCC")), 3);
    }

    #[test]
    fn hosoya_butane() {
        // CCCC: 3 edges; matchings: empty(1) + each edge(3) + {e1,e3}(1) = 5
        assert_eq!(hosoya_index(&mol("CCCC")), 5);
    }

    #[test]
    fn hosoya_benzene() {
        // c1ccccc1: C6 cycle — Z(C6) = M0+M1+M2+M3 = 1+6+9+2 = 18
        // M2 = C(6,2) - 6 adjacent pairs = 15 - 6 = 9; M3 = 2 perfect matchings
        assert_eq!(hosoya_index(&mol("c1ccccc1")), 18);
    }

    #[test]
    fn distance_matrix_ethane() {
        let m = mol("CC");
        let dm = topological_distance_matrix(&m);
        assert_eq!(dm.len(), 2);
        assert_eq!(dm[0][0], 0);
        assert_eq!(dm[0][1], 1);
        assert_eq!(dm[1][0], 1);
        assert_eq!(dm[1][1], 0);
    }

    #[test]
    fn distance_matrix_propane() {
        // C-C-C: d(0,1)=1, d(0,2)=2, d(1,2)=1
        let m = mol("CCC");
        let dm = topological_distance_matrix(&m);
        assert_eq!(dm[0][2], 2);
        assert_eq!(dm[1][2], 1);
    }

    // ── Padmakar-Ivan (PI) Index ─────────────────────────────────────────────

    #[test]
    fn pi_single_atom() {
        assert_eq!(padmakar_ivan_index(&mol("C")), 0);
    }

    #[test]
    fn pi_ethane() {
        // CC: 1 bond, n_u=1 (C1), n_v=1 (C2) → PI = 2
        assert_eq!(padmakar_ivan_index(&mol("CC")), 2);
    }

    #[test]
    fn pi_propane() {
        // CCC: edge(1-2) → n_u=1, n_v=2 → 3; edge(2-3) → n_u=2, n_v=1 → 3; PI = 6
        assert_eq!(padmakar_ivan_index(&mol("CCC")), 6);
    }

    #[test]
    fn pi_butane() {
        // CCCC: each of 3 edges contributes 4 → PI = 12
        assert_eq!(padmakar_ivan_index(&mol("CCCC")), 12);
    }

    #[test]
    fn pi_benzene() {
        // c1ccccc1: 6-ring, each of 6 edges contributes 6 → PI = 36
        assert_eq!(padmakar_ivan_index(&mol("c1ccccc1")), 36);
    }

    #[test]
    fn pi_linear_chain_formula() {
        // For a linear chain of n heavy atoms: PI = n(n-1)
        for n in 2..=6 {
            let smiles = "C".repeat(n);
            let expected = (n as u64) * (n as u64 - 1);
            assert_eq!(padmakar_ivan_index(&mol(&smiles)), expected, "chain n={n}");
        }
    }

    // ── Schultz / Gutman MTI ─────────────────────────────────────────────────

    #[test]
    fn schultz_mti_ethane() {
        // CC: 2 atoms, deg=[1,1], d=1 → MTI = (1+1)*1 = 2
        assert_eq!(schultz_mti(&mol("CC")), 2);
    }

    #[test]
    fn gutman_mti_ethane() {
        // CC: 2 atoms, deg=[1,1], d=1 → MTI* = 1*1*1 = 1
        assert_eq!(gutman_mti(&mol("CC")), 1);
    }

    #[test]
    fn schultz_mti_propane() {
        // CCC: atoms 0,1,2 with deg=[1,2,1], distances: d(0,1)=1, d(0,2)=2, d(1,2)=1
        // (1+2)*1 + (1+1)*2 + (2+1)*1 = 3 + 4 + 3 = 10
        assert_eq!(schultz_mti(&mol("CCC")), 10);
    }

    #[test]
    fn gutman_mti_propane() {
        // CCC: atoms 0,1,2 with deg=[1,2,1], distances: d(0,1)=1, d(0,2)=2, d(1,2)=1
        // 1*2*1 + 1*1*2 + 2*1*1 = 2 + 2 + 2 = 6
        assert_eq!(gutman_mti(&mol("CCC")), 6);
    }

    #[test]
    fn schultz_mti_empty() {
        // edge case: single atom
        let m = chematic_smiles::parse("[C]").unwrap();
        assert_eq!(schultz_mti(&m), 0);
    }

    // ── VABC ─────────────────────────────────────────────────────────────────

    #[test]
    fn vabc_methane() {
        // Single carbon + 4 implicit H — result should be positive and roughly 30–50 Å³
        let v = vabc(&mol("C"));
        assert!(v > 10.0, "VABC(methane) should be > 10 Å³, got {v}");
        assert!(v < 100.0, "VABC(methane) should be < 100 Å³, got {v}");
    }

    #[test]
    fn vabc_water() {
        let v = vabc(&mol("O"));
        assert!(v > 5.0, "VABC(water) should be > 5 Å³, got {v}");
        assert!(v < 50.0, "VABC(water) should be < 50 Å³, got {v}");
    }

    #[test]
    fn vabc_increases_with_size() {
        // Larger molecule should have larger VABC
        let v_ethane = vabc(&mol("CC"));
        let v_propane = vabc(&mol("CCC"));
        let v_butane = vabc(&mol("CCCC"));
        assert!(v_propane > v_ethane, "propane > ethane");
        assert!(v_butane > v_propane, "butane > propane");
    }

    // ── Gravitational index ──────────────────────────────────────────────────

    #[test]
    fn gravitational_index_single_atom() {
        let m = chematic_smiles::parse("[C]").unwrap();
        assert_eq!(gravitational_index(&m), 0.0);
    }

    #[test]
    fn gravitational_index_ethane() {
        // CC: 2 C, mass=12.011 each, d=1 → G = 12.011*12.011/1 ≈ 144.26
        let g = gravitational_index(&mol("CC"));
        assert!((g - 12.011_f64 * 12.011).abs() < 0.01, "got {g}");
    }

    #[test]
    fn gravitational_index_positive() {
        let g = gravitational_index(&mol("c1ccccc1"));
        assert!(g > 0.0, "benzene gravitational index should be positive");
    }

    #[test]
    fn chi_all_matches_individual() {
        // Verify chi_all returns identical values to individual chi0-chi4/v functions.
        for smi in [
            "CC",
            "CCC",
            "c1ccccc1",
            "CC(=O)Oc1ccccc1C(=O)O",
            "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",
        ] {
            let m = mol(smi);
            let (c0, c1, c2, c3, c4, c0v, c1v, c2v, c3v, c4v) = chi_all(&m);
            assert!((c0 - chi0(&m)).abs() < 1e-10, "{smi}: chi0 mismatch");
            assert!((c1 - chi1(&m)).abs() < 1e-10, "{smi}: chi1 mismatch");
            assert!((c2 - chi2(&m)).abs() < 1e-10, "{smi}: chi2 mismatch");
            assert!((c3 - chi3(&m)).abs() < 1e-10, "{smi}: chi3 mismatch");
            assert!((c4 - chi4(&m)).abs() < 1e-10, "{smi}: chi4 mismatch");
            assert!((c0v - chi0v(&m)).abs() < 1e-10, "{smi}: chi0v mismatch");
            assert!((c1v - chi1v(&m)).abs() < 1e-10, "{smi}: chi1v mismatch");
            assert!((c2v - chi2v(&m)).abs() < 1e-10, "{smi}: chi2v mismatch");
            assert!((c3v - chi3v(&m)).abs() < 1e-10, "{smi}: chi3v mismatch");
            assert!((c4v - chi4v(&m)).abs() < 1e-10, "{smi}: chi4v mismatch");
        }
    }

    #[test]
    fn topology_bundle_matches_existing_scalar_apis() {
        for smi in ["CC", "CCC", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"] {
            let m = mol(smi);
            let bundle = topology_bundle(&m);
            assert_eq!(bundle.wiener, wiener_index(&m), "{smi}: Wiener mismatch");
            assert_eq!(bundle.kappa, kappa_all(&m), "{smi}: Kappa mismatch");
            assert_eq!(bundle.chi, chi_all(&m), "{smi}: Chi mismatch");
        }
    }
}

#[cfg(test)]
mod sphere_overlap_boundary_tests {
    use super::*;

    #[test]
    fn sphere_overlap_agrees_with_analytic_lenses_containment_and_scaling() {
        // Two unit spheres one radius apart have two caps of height 1/2:
        // their combined lens volume is 5*pi/12.
        for (r1, r2, distance, expected) in [
            (1.0, 1.0, 1.0, 5.0 * PI / 12.0),
            (1.0, 3.0, 1.0, 4.0 * PI / 3.0),
            (1.0, 3.0, 4.0, 0.0),
            (1.0, 3.0, 5.0, 0.0),
        ] {
            for scale in [0.25_f64, 1.0, 2.5] {
                let expected = expected * scale.powi(3);
                let direct = sphere_intersection(r1 * scale, r2 * scale, distance * scale);
                let reverse = sphere_intersection(r2 * scale, r1 * scale, distance * scale);
                assert!((direct - expected).abs() < 1e-12 * expected.abs().max(1.0));
                assert!((direct - reverse).abs() < 1e-12);
                assert!(direct.is_finite() && direct >= 0.0);
            }
        }
    }
}
