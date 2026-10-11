use super::*;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mmff94_energy::{mmff94_bond_energy, mmff94_bond_energy_resolved};

    #[test]
    fn eq_level_resolves_type_via_ladder_and_falls_back_to_identity() {
        // Type 6's own eqLevel row: [6, 6, 6, 0] -- level4 (stage 2) is a
        // real substitution to 6 (no-op here since it's already 6), level5
        // (stage 3) is the universal wildcard 0.
        assert_eq!(eq_level(6, 0), 6);
        assert_eq!(eq_level(6, 3), 0);
        // A type with no row in MMFF94_EQ_LEVEL (e.g. an exotic/rare type
        // outside this 55-entry table) must fall back to identity at every
        // stage, never panic.
        assert_eq!(eq_level(99, 0), 99);
        assert_eq!(eq_level(99, 3), 99);
    }

    #[test]
    fn angle_energy_resolves_via_eq_level_ladder_when_exact_and_type0_both_miss() {
        // Issue #227 Stage B: (angle_type=0, ti=6, tj=40, tk=37) has no
        // exact-type row, but type 6's eqLevel table ([6, 6, 6, 0]) and
        // type 37's ([37, 2, 1, 0]) meet at stage 1 (level3: 6 and 2) --
        // live-RDKit-oracle-confirmed via
        // `scripts/mmff94_angle_bond_gap_classify.py` against
        // `chembl_tier_b_0164`'s real O-N(aromatic)-C(aromatic) triple.
        let p = mmff94_angle_energy(0, 6, 40, 37).expect("eqLevel ladder must resolve this");
        assert!((p.ka - 1.316).abs() < 1e-6, "ka={} should be 1.316", p.ka);
        assert!(
            (p.theta0 - 115.626).abs() < 1e-6,
            "theta0={} should be 115.626",
            p.theta0
        );
    }

    // ─── Issue #227 Stage C: empirical rule (eq. 20) ──────────────────────

    #[test]
    fn angle_empirical_reuses_wildcard_theta0_for_halogen_and_amine_c_substituents() {
        // angle_type=0, center=type1 (sp3 C, crd=4). All 7 confirmed against
        // the live RDKit oracle by `scripts/mmff94_angle_bond_gap_classify.py`
        // (issue #227 Stage C's own diagnosed table_gap residual). Each hits
        // the eqLevel ladder's level5 (0,0) substitution and lands on the
        // restored wildcard row `(0, 0, 1, 0)` -> theta0=108.9, ka=0.0
        // (`TableSearch::ZeroKa`) -- so theta0 must come back REUSED
        // verbatim, not recomputed by the from-scratch sub-rule (which would
        // also give 109.45 here anyway, since crd=4 -- the two are
        // numerically close for this central type, which is exactly why this
        // reuse-vs-recompute distinction was easy to miss).
        let cases: &[(u8, u8, f64)] = &[
            (11, 12, 1.2566039721725888), // C(F)(Cl)
            (11, 13, 1.199132258341949),  // C(F)(Br)
            (12, 13, 1.1553717810569428), // C(Cl)(Br)
            (11, 14, 1.146032669253874),  // C(F)(I)
            (12, 14, 1.1284626065931707), // C(Cl)(I)
            (13, 14, 1.1317750264511846), // C(Br)(I)
            (8, 11, 1.348550207611195),   // C(N)(F)
        ];
        for &(ti, tk, expected_ka) in cases {
            let r0_i = mmff94_bond_energy(0, 1, ti)
                .unwrap_or_else(|| panic!("C-{ti} bond must be tabulated"))
                .r0;
            let r0_k = mmff94_bond_energy(0, 1, tk)
                .unwrap_or_else(|| panic!("C-{tk} bond must be tabulated"))
                .r0;
            let (p, kind) = mmff94_angle_energy_resolved(0, ti, 1, tk, r0_i, r0_k, 0)
                .unwrap_or_else(|| panic!("({ti},1,{tk}) must resolve via empirical reuse"));
            assert_eq!(kind, Mmff94Resolution::EmpiricalAngle);
            assert!(
                (p.ka - expected_ka).abs() < 1e-3,
                "ti={ti} tk={tk} ka={} expected {expected_ka}",
                p.ka
            );
            assert!(
                (p.theta0 - 108.9).abs() < 1e-9,
                "theta0 must be the wildcard row's own value ({}), reused verbatim -- got {}",
                108.9,
                p.theta0
            );
            // The ka==0.0 filter's whole point: the plain function must
            // never surface this placeholder row as a real parameter.
            assert!(
                mmff94_angle_energy(0, ti, 1, tk).is_none(),
                "plain mmff94_angle_energy must stay None for a ka==0.0 row"
            );
        }
    }

    #[test]
    fn angle_type_3_2_64_resolves_via_new_c5b_eq_level_row() {
        // RDKit/COSMolKit/OpenEye advantage directive, Phase 4 (3D coverage
        // gap): 11 of the 265-molecule strict corpus's coverage failures
        // were exactly this triple (`chembl_tier_b_0180` and 10 others,
        // "Angle C(3)-C(2)-C(64)" -- an indole-ring C5B carbon beta to the
        // ring N, exocyclic to a maleimide-type ring). Confirmed
        // independently via a live RDKit oracle query on
        // `chembl_tier_b_0180`'s own SMILES
        // (`Cn1cc(C2=C(c3cn(CCCSC(=N)N)c4ccccc34)C(=O)NC2=O)c2ccccc21`):
        // `MMFFGetMoleculeProperties(...).GetMMFFAngleBendParams` returns
        // `(angleType=2, ka=0.893, theta0=118.456)` for this exact atom-type
        // triple -- both values reproduced exactly here, resolved via the
        // eqLevel ladder now that type 64 (C5B) has a row.
        let r0_23 = mmff94_bond_energy_resolved(0, 2, 3).unwrap().0.r0;
        let r0_264 = mmff94_bond_energy_resolved(0, 2, 64).unwrap().0.r0;
        let (params, kind) = mmff94_angle_energy_resolved(2, 3, 2, 64, r0_23, r0_264, 0)
            .expect("(angle_type=2, 3, 2, 64) must resolve now that type 64 has an eq_level row");
        assert_eq!(kind, Mmff94Resolution::EquivalentType { level: 3 });
        assert!(
            (params.ka - 0.893).abs() < 1e-9,
            "ka={} should be 0.893 (RDKit oracle value)",
            params.ka
        );
        assert!(
            (params.theta0 - 118.456).abs() < 1e-9,
            "theta0={} should be 118.456 (RDKit oracle value)",
            params.theta0
        );
        // Reversed i/k order must agree (angle bend is symmetric in i/k).
        let (params_rev, _) = mmff94_angle_energy_resolved(2, 64, 2, 3, r0_264, r0_23, 0)
            .expect("reversed (64, 2, 3) must also resolve");
        assert_eq!(params_rev.ka, params.ka);
        assert_eq!(params_rev.theta0, params.theta0);
    }

    #[test]
    fn angle_empirical_reuses_sulfur_wildcard_theta0_for_c5a() {
        // RDKit resolves this tuple through C5A's canonical ladder to the
        // wildcard row (0, 0, 18, 0), whose ka is zero and theta0 is 104.6.
        // The empirical force constant is then computed while reusing that
        // table theta0 (MMFF.V eq. 20), rather than applying the generic
        // tetrahedral 109.45-degree rule.
        let (params, resolution) = mmff94_angle_energy_resolved(0, 43, 18, 63, 1.71, 1.749, 0)
            .expect("C5A ladder should reach the sulfur wildcard fallback");
        assert_eq!(resolution, Mmff94Resolution::EmpiricalAngle);
        assert!((params.ka - 1.281584861919745).abs() < 1e-9);
        assert!((params.theta0 - 104.6).abs() < 1e-9);
    }
}
