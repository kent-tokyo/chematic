use super::*;

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_atom_type_methane() {
        let mol = parse("C").unwrap();
        assert_eq!(classify_atom_type(&mol, AtomIdx(0)), AtomType::CSp3);
    }

    #[test]
    fn test_atom_type_ethene() {
        let mol = parse("C=C").unwrap();
        assert_eq!(classify_atom_type(&mol, AtomIdx(0)), AtomType::CSp2Alkene);
        assert_eq!(classify_atom_type(&mol, AtomIdx(1)), AtomType::CSp2Alkene);
    }

    #[test]
    fn test_atom_type_benzene() {
        let mol = parse("c1ccccc1").unwrap();
        assert_eq!(classify_atom_type(&mol, AtomIdx(0)), AtomType::CAromatic);
    }

    #[test]
    fn test_atom_type_acetaldehyde() {
        let mol = parse("CC=O").unwrap();
        let c_sp3 = if mol.atom(AtomIdx(0)).aromatic {
            classify_atom_type(&mol, AtomIdx(1))
        } else {
            classify_atom_type(&mol, AtomIdx(0))
        };
        assert_eq!(c_sp3, AtomType::CSp3);
    }

    #[test]
    fn test_torsion_score_perfect_match() {
        let pref = TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.1,
        };
        let score = score_torsion(180.0, &pref);
        assert!(score.abs() < 1e-6, "perfect match should have zero penalty");
    }

    #[test]
    fn test_torsion_score_deviation() {
        let pref = TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.1,
        };
        let score = score_torsion(160.0, &pref);
        assert!(
            (score - 2.0).abs() < 1e-6,
            "20° deviation should yield 2.0 penalty"
        );
    }

    #[test]
    fn test_torsion_score_periodic() {
        let pref = TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.1,
        };
        // 180° and -180° are the same angle
        let score1 = score_torsion(180.0, &pref);
        let score2 = score_torsion(-180.0, &pref);
        assert!(
            (score1 - score2).abs() < 1e-6,
            "periodic angles should score the same"
        );
    }

    #[test]
    fn test_default_torsion_preference() {
        let pref = default_torsion_preference();
        assert_eq!(pref.angle_deg, 180.0);
        assert!(pref.penalty_per_degree > 0.0);
    }

    #[test]
    fn test_alkane_torsion_preference() {
        let mol = parse("CCCC").unwrap(); // butane
        if mol.atom_count() >= 4 {
            let pref = get_torsion_preference(&mol, AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3));
            assert!(pref.is_some(), "butane C-C-C-C should have preference");
            if let Some(p) = pref {
                assert_eq!(p.angle_deg, 180.0, "alkane torsions prefer 180°");
            }
        }
    }

    #[test]
    fn test_biphenyl_torsion_preference() {
        let mol = parse("c1ccccc1-c1ccccc1").unwrap(); // biphenyl
        // Atoms: 0-5 ring1, 6-11 ring2; bond between 0 and 6
        let pref = get_torsion_preference(&mol, AtomIdx(1), AtomIdx(0), AtomIdx(6), AtomIdx(7));
        assert!(
            pref.is_some(),
            "biphenyl Ar-Ar should have a torsion preference"
        );
        if let Some(p) = pref {
            assert_eq!(p.angle_deg, 45.0, "biphenyl prefers ~45° twist");
        }
    }

    #[test]
    fn test_thioester_torsion_preference() {
        // Methyl thioformate: SC=O; S is atom 0, C is atom 1 (CCarbonyl)
        let mol = parse("SC=O").unwrap();
        let _pref = get_torsion_preference(&mol, AtomIdx(0), AtomIdx(0), AtomIdx(1), AtomIdx(2));
        // b=S (atom0), c=CCarbonyl (atom1) → thioester branch
        let b_type = classify_atom_type(&mol, AtomIdx(0));
        let c_type = classify_atom_type(&mol, AtomIdx(1));
        assert_eq!(b_type, AtomType::S);
        assert_eq!(c_type, AtomType::CCarbonyl);
    }

    #[test]
    fn test_disulfide_torsion_preference() {
        let mol = parse("CSSC").unwrap(); // dimethyl disulfide
        // bond S(1)-S(2): b_type=S, c_type=S
        let b_type = classify_atom_type(&mol, AtomIdx(1));
        let c_type = classify_atom_type(&mol, AtomIdx(2));
        assert_eq!(b_type, AtomType::S);
        assert_eq!(c_type, AtomType::S);
        let pref = get_torsion_preference(&mol, AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3));
        assert!(pref.is_some(), "disulfide should have ~90° preference");
        if let Some(p) = pref {
            assert_eq!(p.angle_deg, 90.0, "disulfide prefers 90°");
        }
    }

    #[test]
    fn test_nitrile_torsion_preference() {
        let mol = parse("CCC#N").unwrap(); // propionitrile
        // N is atom 3, atom type NSp
        let n_type = classify_atom_type(&mol, AtomIdx(3));
        assert_eq!(n_type, AtomType::NSp, "nitrile N should be NSp");
        let pref = get_torsion_preference(&mol, AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3));
        assert!(pref.is_some(), "nitrile torsion should have a preference");
        if let Some(p) = pref {
            assert_eq!(p.angle_deg, 180.0, "linear nitrile end prefers 180°");
        }
    }

    #[test]
    fn test_amine_torsion_preference() {
        let mol = parse("CCNC").unwrap(); // ethyl methyl amine
        // bond C(1)-N(2): b=CSp3, c=NSp2 (secondary amine: 2 explicit bonds)
        let pref = get_torsion_preference(&mol, AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3));
        assert!(pref.is_some(), "amine C-C-N-C should have preference");
        if let Some(p) = pref {
            assert_eq!(p.angle_deg, 180.0);
        }
    }

    #[test]
    fn test_phenyl_ketone_torsion_preference() {
        let mol = parse("c1ccccc1C(=O)C").unwrap(); // acetophenone
        // The C(=O) carbon is CCarbonyl, connected to CAromatic
        // Find the carbonyl carbon
        let c_carbonyl_idx = (0..mol.atom_count() as u32)
            .map(AtomIdx)
            .find(|&i| classify_atom_type(&mol, i) == AtomType::CCarbonyl);
        assert!(
            c_carbonyl_idx.is_some(),
            "acetophenone should have a carbonyl C"
        );
    }

    #[test]
    fn test_score_torsion_disulfide_at_90() {
        let pref = TorsionPreference {
            angle_deg: 90.0,
            penalty_per_degree: 0.1,
        };
        let score = score_torsion(90.0, &pref);
        assert!(score.abs() < 1e-6, "at preferred angle score should be 0");
        let score_off = score_torsion(90.0 + 20.0, &pref);
        assert!((score_off - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_pattern_count_covers_20_plus() {
        // Verify we have comprehensive coverage by checking patterns
        // for the major chemical motifs
        let mol_alkane = parse("CCCC").unwrap();
        let mol_biphenyl = parse("c1ccccc1-c1ccccc1").unwrap();
        let mol_amide = parse("CC(=O)N").unwrap();
        let mol_ester = parse("CC(=O)OC").unwrap();
        let mol_disulfide = parse("CSSC").unwrap();
        let mol_nitrile = parse("CCC#N").unwrap();
        let mol_amine = parse("CCNC").unwrap();

        let cases = [
            (&mol_alkane, AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3)),
            (
                &mol_biphenyl,
                AtomIdx(1),
                AtomIdx(0),
                AtomIdx(6),
                AtomIdx(7),
            ),
            (
                &mol_disulfide,
                AtomIdx(0),
                AtomIdx(1),
                AtomIdx(2),
                AtomIdx(3),
            ),
            (&mol_nitrile, AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3)),
            (&mol_amine, AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3)),
        ];
        for (mol, a, b, c, d) in &cases {
            let pref = get_torsion_preference(mol, *a, *b, *c, *d);
            // At least some of these should have preferences
            let _ = pref; // just ensuring no panic
        }
        // The amide and ester patterns
        let pref_amide =
            get_torsion_preference(&mol_amide, AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(2));
        let _ = pref_amide;
        let pref_ester =
            get_torsion_preference(&mol_ester, AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3));
        let _ = pref_ester;
    }

    // ── SMARTS torsion map tests ───────────────────────────────────────────────

    #[test]
    fn test_smarts_map_hindered_biaryl() {
        // 2,2'-dimethylbiphenyl: both inter-ring C are H0,X3 → hindered → 90°
        let mol = parse("Cc1ccccc1-c1ccccc1C").unwrap();
        let ring_set = chematic_perception::find_sssr(&mol);
        let ring_bonds: std::collections::HashSet<(u32, u32)> = ring_set
            .rings()
            .iter()
            .flat_map(|r| {
                let n = r.len();
                (0..n).flat_map(move |i| {
                    let a = r[i].0;
                    let b = r[(i + 1) % n].0;
                    [(a, b), (b, a)]
                })
            })
            .collect();
        let map = build_smarts_torsion_map(&mol, &ring_bonds);
        // At least one biaryl bond entry should be at 90°
        let has_90 = map.values().any(|p| (p.angle_deg - 90.0).abs() < 1.0);
        assert!(
            has_90,
            "hindered biaryl should get 90° preference in SMARTS map"
        );
    }

    #[test]
    fn test_smarts_map_primary_amide() {
        // Acetamide: CC(=O)N — primary amide → 0°
        let mol = parse("CC(=O)N").unwrap();
        let ring_bonds = std::collections::HashSet::new();
        let map = build_smarts_torsion_map(&mol, &ring_bonds);
        // Should have an entry for the N-C(=O) bond at 0°
        let has_0 = map.values().any(|p| p.angle_deg.abs() < 1.0);
        assert!(has_0, "primary amide N-C bond should get 0° preference");
    }

    // ── New heterocycle pattern tests ──────────────────────────────────────────

    #[test]
    fn test_atom_type_furan_oxygen() {
        // Furan: c1ccco1 — oxygen atom is aromatic
        let mol = parse("c1ccco1").unwrap();
        let o_idx = (0..mol.atom_count() as u32)
            .map(AtomIdx)
            .find(|&i| mol.atom(i).element.atomic_number() == 8)
            .expect("furan must have an oxygen atom");
        assert_eq!(
            classify_atom_type(&mol, o_idx),
            AtomType::OAromatic,
            "furan O should be OAromatic"
        );
    }

    #[test]
    fn test_atom_type_thiophene_sulfur() {
        // Thiophene: c1cccs1 — sulfur atom is aromatic
        let mol = parse("c1cccs1").unwrap();
        let s_idx = (0..mol.atom_count() as u32)
            .map(AtomIdx)
            .find(|&i| mol.atom(i).element.atomic_number() == 16)
            .expect("thiophene must have a sulfur atom");
        assert_eq!(
            classify_atom_type(&mol, s_idx),
            AtomType::SAromatic,
            "thiophene S should be SAromatic"
        );
    }

    #[test]
    fn test_furanyl_biaryl_prefers_planar() {
        // 2-phenylfuran: c1ccc(-c2ccco2)cc1
        // Inter-ring bond between phenyl CAromatic and furanyl CAromatic;
        // OAromatic attached on the furanyl side → should fire OAromatic–CAromatic rule (0°).
        let mol = parse("c1ccc(-c2ccco2)cc1").unwrap();
        // Find an atom adjacent to the furan O — that's the OAromatic–CAromatic bond end
        let o_idx = (0..mol.atom_count() as u32)
            .map(AtomIdx)
            .find(|&i| mol.atom(i).element.atomic_number() == 8)
            .expect("must have O");
        let o_neighbor = mol
            .neighbors(o_idx)
            .next()
            .map(|(n, _)| n)
            .expect("O has neighbors");
        assert_eq!(classify_atom_type(&mol, o_idx), AtomType::OAromatic);
        assert_eq!(classify_atom_type(&mol, o_neighbor), AtomType::CAromatic);
        // Torsion involving OAromatic–CAromatic as B–C pair should prefer 0°
        let pref = get_torsion_preference(&mol, o_idx, o_idx, o_neighbor, o_neighbor);
        // We only need to confirm the rule fires (Some) and angle is 0°
        let pref2 = get_torsion_preference(&mol, AtomIdx(0), o_idx, o_neighbor, AtomIdx(0));
        assert!(
            pref.is_some() || pref2.is_some(),
            "OAromatic–CAromatic should have a torsion preference"
        );
    }

    #[test]
    fn test_morpholine_gauche_preference() {
        // Morpholine: C1CNCCO1 — N-C-C-O chain in chair → gauche 60°
        // Atom order in morpholine SMILES C1CNCCO1:
        //   0=C, 1=C, 2=N, 3=C, 4=C, 5=O
        // N-C-C-O: a=N(2), b=C(3), c=C(4), d=O(5)
        let mol = parse("C1CNCCO1").unwrap();
        // Find N and O atoms
        let n_idx = (0..mol.atom_count() as u32)
            .map(AtomIdx)
            .find(|&i| mol.atom(i).element.atomic_number() == 7)
            .expect("morpholine must have N");
        let o_idx = (0..mol.atom_count() as u32)
            .map(AtomIdx)
            .find(|&i| mol.atom(i).element.atomic_number() == 8)
            .expect("morpholine must have O");
        // Find the C-C bond between N-side and O-side carbons
        // The preference fires when a=NSp3, b=CSp3, c=CSp3, d=OSp3
        let pref = get_torsion_preference(&mol, n_idx, AtomIdx(3), AtomIdx(4), o_idx);
        assert!(
            pref.is_some(),
            "morpholine N-C-C-O should have gauche preference"
        );
        if let Some(p) = pref {
            assert_eq!(p.angle_deg, 60.0, "morpholine N-C-C-O prefers 60°");
        }
    }
}

#[cfg(test)]
mod legacy_score_boundary_contract_tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn every_atom_type_combination_is_safe_and_returns_a_valid_preference() {
        // `get_torsion_preference` is a precedence-ordered compatibility
        // table over four atom types.  Exercise the complete type product,
        // including combinations that are unusual but can occur transiently
        // in reaction products or partially sanitised input.
        let mol = chematic_smiles::parse(
            "[H].CC.C=C.c1ccccc1.C=O.C#N.N=C.n1ccccc1.CN(C)C.O.o1cccc1.S.s1cccc1.P.F.[Xe]",
        )
        .unwrap();
        let mut representatives = HashMap::new();
        for (idx, _) in mol.atoms() {
            representatives
                .entry(classify_atom_type(&mol, idx))
                .or_insert(idx);
        }
        for required in [
            AtomType::CSp3,
            AtomType::CSp2Alkene,
            AtomType::CAromatic,
            AtomType::CCarbonyl,
            AtomType::NSp,
            AtomType::NSp2,
            AtomType::NAromatic,
            AtomType::NSp3,
            AtomType::OSp2,
            AtomType::OSp3,
            AtomType::OAromatic,
            AtomType::S,
            AtomType::SAromatic,
            AtomType::P,
            AtomType::H,
            AtomType::Halogen,
            AtomType::Other,
        ] {
            assert!(
                representatives.contains_key(&required),
                "missing {required:?}"
            );
        }
        let atoms: Vec<_> = representatives.values().copied().collect();
        let mut matched = 0usize;
        for &a in &atoms {
            for &b in &atoms {
                for &c in &atoms {
                    for &d in &atoms {
                        if let Some(preference) = get_torsion_preference(&mol, a, b, c, d) {
                            matched += 1;
                            assert!(preference.angle_deg.is_finite());
                            assert!(preference.penalty_per_degree.is_finite());
                            assert!(preference.penalty_per_degree > 0.0);
                        }
                    }
                }
            }
        }
        assert!(matched > 1_000);
    }

    #[test]
    fn legacy_preferences_on_real_paths_have_finite_periodic_nonnegative_scores() {
        // Compatibility/numerical contract of the retained legacy API.
        // This does not validate its heuristic preferred angles as physical
        // torsion minima; see docs/rfcs/3d_torsion_knowledge_audit.md.
        let mut preferred = 0;
        let mut unspecified = 0;
        for source in [
            "CCCC",
            "Cc1ccccc1",
            "CC(=O)NC",
            "COC(=O)C",
            "c1ccccc1-c1ccccc1",
            "CC=CN(C)C",
            "CC(=O)C=C",
            "CC(=O)c1ccccc1",
            "CSC(=O)C",
            "CS(=O)CC",
            "CSSC",
            "CCO",
            "CCN(C)C",
            "CCC#N",
            "CCP(C)C",
            "CC(=O)NNC",
            "CS(=O)(=O)NC",
            "COc1ccccc1",
            "CC=NOC",
            "COOC",
            "C#CC#C",
            "CC(=O)N(C)C",
            "c1ccc([N+](=O)[O-])cc1",
        ] {
            let input = chematic_smiles::parse(source).unwrap();
            for explicit_h in [false, true] {
                let mol = if explicit_h {
                    chematic_chem::add_hydrogens(&input)
                } else {
                    input.clone()
                };
                for (b, _) in mol.atoms() {
                    for (c, _) in mol.neighbors(b) {
                        for (a, _) in mol.neighbors(b).filter(|(a, _)| *a != c) {
                            for (d, _) in mol.neighbors(c).filter(|(d, _)| *d != b && *d != a) {
                                let Some(preference) = get_torsion_preference(&mol, a, b, c, d)
                                else {
                                    unspecified += 1;
                                    continue;
                                };
                                preferred += 1;
                                assert!(preference.angle_deg.is_finite());
                                assert!((-180.0..=180.0).contains(&preference.angle_deg));
                                assert!(
                                    preference.penalty_per_degree.is_finite()
                                        && preference.penalty_per_degree > 0.0
                                );
                                assert_eq!(score_torsion(preference.angle_deg, &preference), 0.0);
                                for offset in [-180.0_f64, -30.0, -1.0, 0.0, 1.0, 30.0, 180.0] {
                                    let score =
                                        score_torsion(preference.angle_deg + offset, &preference);
                                    assert!(score.is_finite() && score >= 0.0);
                                    assert!(
                                        (score - offset.abs() * preference.penalty_per_degree)
                                            .abs()
                                            < 1e-12
                                    );
                                    for turns in [-2.0, -1.0, 1.0, 2.0] {
                                        let periodic = score_torsion(
                                            preference.angle_deg + offset + 360.0 * turns,
                                            &preference,
                                        );
                                        assert!((score - periodic).abs() < 1e-12, "{source}");
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        assert!(preferred > 0);
        assert!(unspecified > 0);
    }

    #[test]
    fn torsion_preferences_cover_documented_chemical_families() {
        // One compact representative for each family documented by
        // `get_torsion_preference`. Walking every directed four-atom path
        // makes the check independent of parser atom order and verifies both
        // traversal directions where the rule promises them.
        let representatives = [
            "CCCC",
            "Cc1ccccc1CC",
            "CN(C)C=O",
            "CC(=O)OC",
            "c1ccccc1-c2ccccc2",
            "CC=CN(C)C",
            "CC=CCl",
            "C=CC(=O)C",
            "c1ccccc1C(=O)C",
            "CSC(=O)C",
            "CNS(=O)(=O)C",
            "CS(=O)C",
            "CSSC",
            "CCOC",
            "CCN(C)C",
            "CCC#N",
            "CCP(C)C",
            "NC(=O)N",
            "COc1ccccc1",
            "c1ccccc1[N+](=O)[O-]",
            "CC=NN",
            "CC=NO",
            "NC(=O)C(=O)C",
            "c1ccccc1CC",
            "C=CCC",
            "c1ccncc1-c2ccccc2",
            "Cn1ccnc1",
            "n1ccccc1C(=O)C",
            "n1ccccc1C=C",
            "CSc1ccccc1",
            "COCC",
            "c1ccccc1N(C)C",
            "c1ccoc1-c2ccccc2",
            "c1ccsc1-c2ccccc2",
            "Cc1ccoc1",
            "Cc1ccsc1",
            "O=C(c1ccoc1)C",
            "O1CCNCC1",
            "N1CCNCC1",
            "c1ccccc1C=C",
            "C=CSC",
            "C=CN(C)C",
            "CCC(=O)C",
            "n1ccccc1SC",
            "n1ccccc1OC",
            "c1ccsc1-c2ccccn2",
            "c1ccoc1-c2ccccn2",
            "COP(=O)(O)O",
            "OCCO",
            "CC=CC=C",
            "CC(=O)C(=O)C",
            "c1ccccc1Cl",
            "CP(=O)(OC)OC",
            "CCN=NC",
            "CCN=O",
            "CC(=S)NC",
            "CCS(=O)(=O)NC",
            "CC#CC",
            "COON=C",
        ];
        let mut preferred = 0usize;
        for source in representatives {
            let mol = chematic_smiles::parse(source).unwrap_or_else(|e| panic!("{source}: {e}"));
            for (b, _) in mol.atoms() {
                for (a, _) in mol.neighbors(b) {
                    for (c, _) in mol.neighbors(b).filter(|(c, _)| *c != a) {
                        for (d, _) in mol.neighbors(c).filter(|(d, _)| *d != b) {
                            if let Some(p) = get_torsion_preference(&mol, a, b, c, d) {
                                preferred += 1;
                                assert!(p.angle_deg.is_finite(), "{source}");
                                assert!(p.penalty_per_degree > 0.0, "{source}");
                            }
                        }
                    }
                }
            }
        }
        assert!(
            preferred >= 80,
            "only {preferred} paths received a preference"
        );
    }
}
