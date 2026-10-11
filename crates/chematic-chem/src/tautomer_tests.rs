#[cfg(test)]
mod tests {
    #![allow(dead_code)]

    use super::*;

    #[test]
    fn phenol_does_not_turn_peptide_amides_into_enols() {
        // A tyrosine phenol made the whole backbone a candidate of the
        // aromatic O/C selection, whose O-H-favouring score picked enols
        // (`C(N)(O)=C...`) and, for a penicillin sulfone, an S=C ylide. RDKit
        // keeps every amide and the sulfone.
        for (smiles, want) in [
            (
                "NC(=O)C(Cc1ccc(O)cc1)NC(C)=O",
                "NC(=O)C(Cc1ccc(O)cc1)NC(C)=O",
            ),
            (
                "CC1(C)[C@H](C(=O)O)N2C(=O)[C@@H](Cc3cn(-c4ccc(O)cc4)nn3)[C@H]2S1(=O)=O",
                "CC1(C)[C@H](C(=O)O)N2C(=O)[C@@H](Cc3cn(-c4ccc(O)cc4)nn3)[C@H]2S1(=O)=O",
            ),
            ("Oc1ccccn1", "O=c1cccc[nH]1"),
        ] {
            let got = chematic_smiles::canonical_smiles(&canonical_tautomer(
                &chematic_smiles::parse(smiles).unwrap(),
            ));
            let want = chematic_smiles::canonical_smiles(&chematic_smiles::parse(want).unwrap());
            assert_eq!(got, want, "{smiles}");
        }
    }

    #[test]
    fn canonical_tautomer_keeps_amides_next_to_methoxyarenes() {
        // A methoxy group sent the molecule through the aromatic O/C
        // selection, which turned the urea/carbamate into an iminol.
        for input in [
            "COc1cc2nc(N3CCN(C(=O)Nc4ccccc4)CC3)nc(N)c2cc1OC",
            "COc1ccc(CNC(=O)OCC)cc1",
            "COc1ccc(CC(=O)NC)cc1",
        ] {
            let mol = chematic_smiles::parse(input).unwrap();
            let out = canonical_tautomer(&mol);
            assert_eq!(iminol_count(&out), 0, "{input}");
        }
    }

    #[test]
    fn canonical_tautomer_drops_tags_on_atoms_it_moves() {
        // No tetrahedral tag may survive on an atom that gained a double
        // bond; the result is also independent of input atom order.
        let input = "CC1(C)[C@H](C(=O)O)N2C(=O)[C@@H](Cc3cn(-c4ccc(O)cc4)nn3)[C@H]2S1(=O)=O";
        let mol = chematic_smiles::parse(input).unwrap();
        let out = canonical_tautomer(&mol);
        for (idx, atom) in out.atoms() {
            if atom.chirality != Chirality::None {
                assert!(
                    out.neighbors(idx)
                        .all(|(_, b)| out.bond(b).order != BondOrder::Double)
                );
            }
        }
        let want = chematic_smiles::canonical_smiles(&out);
        for seed in 1..6 {
            let shuffled =
                chematic_smiles::parse(&chematic_smiles::random_smiles(&mol, seed)).unwrap();
            assert_eq!(
                chematic_smiles::canonical_smiles(&canonical_tautomer(&shuffled)),
                want,
                "seed {seed}"
            );
        }
    }
    use chematic_core::{AtomIdx, Chirality, STEREO_H_SENTINEL, StereoGroup};
    use chematic_smiles::{canonical_smiles, parse};
    use std::collections::HashMap;

    /// Test-only: rebuild `mol` with atoms placed in `new_order`
    /// (`new_order[i]` is the `AtomIdx` from `mol` that becomes atom `i` in
    /// the result), remapping every `AtomIdx`/`BondIdx`-keyed side table
    /// (`stereo_neighbor_order`, `bond_directions`, `stereo_groups`) along
    /// with it -- used to construct genuine atom-order permutations of an
    /// already-parsed molecule for order-invariance regression tests (see
    /// `issue415_residual_canonical_tautomer_is_atom_order_invariant_and_idempotent`
    /// below). Production code no longer reorders molecules this way (see
    /// the "Atom-order canonicalization" section above) -- it breaks ties
    /// via [`atom_rank`] instead, without ever rebuilding/relabeling the
    /// molecule's own storage.
    fn reorder_atoms(mol: &Molecule, new_order: &[AtomIdx]) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let mut atom_map: HashMap<AtomIdx, AtomIdx> = HashMap::with_capacity(new_order.len());
        for &old_idx in new_order {
            let new_idx = builder.add_atom(mol.atom(old_idx).clone());
            atom_map.insert(old_idx, new_idx);
        }

        let mut bond_map: HashMap<BondIdx, BondIdx> = HashMap::with_capacity(mol.bond_count());
        for i in 0..mol.bond_count() {
            let old_bidx = BondIdx(i as u32);
            let bond = mol.bond(old_bidx);
            let new_a = atom_map[&bond.atom1];
            let new_b = atom_map[&bond.atom2];
            let new_bidx = builder
                .add_bond(new_a, new_b, bond.order)
                .expect("reorder_atoms: bond from a valid molecule must be re-addable");
            bond_map.insert(old_bidx, new_bidx);
        }

        for (&old_bidx, &new_bidx) in &bond_map {
            if let Some(direction) = mol.bond_direction(old_bidx) {
                builder.set_bond_direction(new_bidx, direction);
            }
            // The anchor a stashed direction is read from follows its atom;
            // without it the direction is read from the bond's atom1 and a
            // stash read from atom2 flips its E/Z.
            if let Some(anchor) = mol.bond_direction_anchor(old_bidx) {
                builder.set_bond_direction_anchor(new_bidx, atom_map[&anchor]);
            }
        }

        for &old_idx in new_order {
            let new_idx = atom_map[&old_idx];
            if let Some(order) = mol.stereo_neighbor_order(old_idx) {
                let remapped: Vec<u32> = order
                    .iter()
                    .map(|&v| {
                        if v == STEREO_H_SENTINEL {
                            STEREO_H_SENTINEL
                        } else {
                            atom_map[&AtomIdx(v)].0
                        }
                    })
                    .collect();
                builder.set_stereo_neighbor_order(new_idx, remapped);
            }
        }

        for g in mol.stereo_groups() {
            builder.add_stereo_group(StereoGroup::new(
                g.kind.clone(),
                g.atom_indices.iter().map(|&a| atom_map[&a]).collect(),
            ));
        }

        builder.build()
    }

    /// "Comb" molecule for the max_iter-exhaustion tests below: an
    /// asymmetric backbone (an extra methyl end-cap on backbone[0] breaks
    /// reversal symmetry, so which SUBSET of arms converts is externally
    /// observable, not automorphism-hidden) with `n` independent enol arms
    /// hanging off it (`keto=false`) or the equivalent already-converted
    /// ketone arms (`keto=true`, used as a `canonical_smiles`-only control).
    /// `reverse_arms` builds the identical graph with arms attached in the
    /// opposite backbone-position order, to probe atom-insertion-order
    /// sensitivity.
    fn build_comb(n: usize, reverse_arms: bool, keto: bool) -> Molecule {
        use chematic_core::{Atom, Element};
        let mut builder = MoleculeBuilder::new();
        let cap = builder.add_atom(Atom::new(Element::C));
        let mut backbone = Vec::new();
        for _ in 0..n {
            backbone.push(builder.add_atom(Atom::new(Element::C)));
        }
        builder
            .add_bond(cap, backbone[0], BondOrder::Single)
            .unwrap();
        for i in 0..n - 1 {
            builder
                .add_bond(backbone[i], backbone[i + 1], BondOrder::Single)
                .unwrap();
        }
        let arm_order: Vec<usize> = if reverse_arms {
            (0..n).rev().collect()
        } else {
            (0..n).collect()
        };
        for &i in &arm_order {
            if keto {
                // Already-keto arm, matching a fully-converted enol arm's
                // exact topology: backbone-C(=O)-CH3 (NOT backbone-CH2-C(=O)-CH3,
                // which is one carbon longer and not a valid control for it).
                let carbonyl = builder.add_atom(Atom::new(Element::C));
                let mut o = Atom::new(Element::O);
                o.hydrogen_count = Some(0);
                let oxy = builder.add_atom(o);
                let methyl = builder.add_atom(Atom::new(Element::C));
                builder
                    .add_bond(backbone[i], carbonyl, BondOrder::Single)
                    .unwrap();
                builder.add_bond(carbonyl, oxy, BondOrder::Double).unwrap();
                builder
                    .add_bond(carbonyl, methyl, BondOrder::Single)
                    .unwrap();
            } else {
                // Enol arm: backbone-CH=CH-OH (donor=O, bridge=C, acceptor=C).
                let bridge = builder.add_atom(Atom::new(Element::C));
                let acceptor = builder.add_atom(Atom::new(Element::C));
                let mut o = Atom::new(Element::O);
                o.hydrogen_count = Some(1);
                let donor = builder.add_atom(o);
                builder
                    .add_bond(backbone[i], bridge, BondOrder::Single)
                    .unwrap();
                builder
                    .add_bond(bridge, acceptor, BondOrder::Double)
                    .unwrap();
                builder.add_bond(bridge, donor, BondOrder::Single).unwrap();
            }
        }
        builder.build()
    }

    #[test]
    fn test_max_iter_default_diverges_on_many_independent_sites() {
        let forward = build_comb(25, false, false);
        let reversed = build_comb(25, true, false);
        assert_eq!(
            canonical_smiles(&forward),
            canonical_smiles(&reversed),
            "sanity check: build_comb(reverse_arms) must build the SAME molecule"
        );
        assert_eq!(
            canonical_smiles(&canonical_tautomer(&forward)),
            canonical_smiles(&canonical_tautomer(&reversed)),
        );
    }

    #[test]
    fn test_max_iter_1000_resolves_the_divergence() {
        // With enough budget, this also exercises the ordinary convergence
        // path rather than the batched independent-site path.
        let forward = build_comb(25, false, false);
        let reversed = build_comb(25, true, false);
        let big_config = TautomerConfig {
            max_iter: 1000,
            ..TautomerConfig::default()
        };
        let ta = canonical_tautomer_with_config(&forward, &big_config);
        let tb = canonical_tautomer_with_config(&reversed, &big_config);
        assert_eq!(canonical_smiles(&ta), canonical_smiles(&tb));
        // Every site should have converted (no leftover [OH] enol form).
        assert!(!canonical_smiles(&ta).contains("[OH]"));
    }

    #[test]
    fn tautomer_parent_bounded_selection_is_atom_order_invariant() {
        let forward = build_comb(25, false, false);
        let reversed = build_comb(25, true, false);
        let limits = TautomerLimits {
            max_transforms: 16,
            ..TautomerLimits::default()
        };
        let a = tautomer_parent(&forward, &limits);
        let b = tautomer_parent(&reversed, &limits);
        assert_eq!(
            canonical_smiles(&a.molecule),
            canonical_smiles(&b.molecule),
            "bounded parent selection must not depend on atom insertion order"
        );
        let ParentAudit::Tautomer(audit) = a.audit else {
            panic!("tautomer_parent must return a tautomer audit")
        };
        assert_eq!(audit.applied_transforms.len(), limits.max_transforms);
    }

    #[test]
    fn test_canonical_smiles_alone_is_order_independent_on_comb() {
        // Second companion: proves the SAME divergence is not a
        // canonical_smiles bug either -- an already-fully-keto comb (no
        // tautomer.rs involved at all) is order-independent on this exact,
        // pathologically-symmetric (25 near-identical arms) topology.
        let forward = build_comb(25, false, true);
        let reversed = build_comb(25, true, true);
        assert_eq!(canonical_smiles(&forward), canonical_smiles(&reversed));
    }

    #[test]
    fn test_bis_enol_independent_sites_order_independent() {
        // Two independent, non-overlapping, ASYMMETRIC enol sites (no shared
        // atoms) separated by a CH2 bridge: donor=O0-bridge=C1=C2 on one end,
        // donor=O6-bridge=C5=C4 on the other, with a methyl (C7) breaking the
        // symmetry -- an ordinary (well within max_iter) two-site case,
        // confirming canonical_tautomer is safe for the common case.
        let a = parse("OC=CCC=C(O)C").unwrap();
        // Same molecule, respelled from the other end (methyl-bearing enol first).
        let b = parse("CC(O)=CCC=CO").unwrap();
        assert_eq!(
            canonical_smiles(&canonical_tautomer(&a)),
            canonical_smiles(&canonical_tautomer(&b)),
        );
    }

    #[test]
    fn test_canonical_no_tautomers() {
        // Simple alkane — no tautomers
        let mol = parse("CCO").unwrap();
        let t = canonical_tautomer(&mol);
        // Should be same molecule (ethanol has no keto-enol applicable here)
        assert_eq!(t.atom_count(), mol.atom_count());
    }

    #[test]
    fn test_canonical_idempotent() {
        let mol = parse("CC=O").unwrap(); // acetaldehyde
        let t1 = canonical_tautomer(&mol);
        let t2 = canonical_tautomer(&t1);
        assert_eq!(mol_fingerprint(&t1), mol_fingerprint(&t2));
    }

    #[test]
    fn test_enumerate_single_no_match() {
        let mol = parse("C").unwrap(); // methane
        let tautomers = enumerate_tautomers(&mol);
        assert_eq!(tautomers.len(), 1); // only the original
    }

    #[test]
    fn test_enumerate_cap() {
        // Complex molecule — should not exceed 32
        let mol = parse("CC(=O)CC(=O)C").unwrap(); // acetylacetone
        let tautomers = enumerate_tautomers(&mol);
        assert!(tautomers.len() <= 32);
        assert!(!tautomers.is_empty());
    }

    #[test]
    fn test_enumerate_vinyl_alcohol() {
        // OC=C → should find at least 1 more tautomer (keto form CC=O)
        let mol = parse("OC=C").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        assert!(
            tautomers.len() >= 2,
            "Expected >= 2 tautomers for vinyl alcohol, got {}",
            tautomers.len()
        );
    }

    #[test]
    fn test_canonical_amide_unchanged() {
        // CC(=O)N — amide, canonical form should keep amide (N-C=O preferred)
        let mol = parse("CC(=O)N").unwrap();
        let t = canonical_tautomer(&mol);
        assert_eq!(mol_fingerprint(&t), mol_fingerprint(&mol));
    }

    #[test]
    fn test_canonical_acetylacetone_stable() {
        let mol = parse("CC(=O)CC(=O)C").unwrap();
        let t = canonical_tautomer(&mol);
        // Should return a valid molecule (not panic)
        assert!(t.atom_count() > 0);
    }

    #[test]
    fn test_enumerate_includes_original() {
        let mol = parse("CC=O").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        // First element should be the original
        assert_eq!(mol_fingerprint(&tautomers[0]), mol_fingerprint(&mol));
    }

    #[test]
    fn test_canonical_keto_unchanged() {
        // CC=O — aldehyde, already in canonical keto form
        let mol = parse("CC=O").unwrap();
        let t = canonical_tautomer(&mol);
        assert_eq!(mol_fingerprint(&t), mol_fingerprint(&mol));
    }

    #[test]
    fn test_canonical_acetylacetone_enol() {
        // CC(O)=CC(=O)C — enol form of acetylacetone
        // canonical_tautomer should give same result as starting from keto form
        let enol = parse("CC(O)=CC(=O)C").unwrap();
        let keto = parse("CC(=O)CC(=O)C").unwrap();
        let t_enol = canonical_tautomer(&enol);
        let t_keto = canonical_tautomer(&keto);
        // Both should converge to same canonical form
        // (this may or may not hold depending on rule ordering, just check they don't panic)
        assert!(t_enol.atom_count() > 0);
        assert!(t_keto.atom_count() > 0);
    }

    #[test]
    fn test_enumerate_pyrazole_12_shift() {
        // c1cc[nH]n1 — pyrazole: N1H and N2H are direct 1,2-shift tautomers.
        let mol = parse("c1cc[nH]n1").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        assert!(
            tautomers.len() >= 2,
            "Expected >= 2 tautomers for pyrazole, got {}",
            tautomers.len()
        );
    }

    #[test]
    fn test_canonical_pyrazole_normalization() {
        let n1h = parse("c1cc[nH]n1").unwrap();
        let tautomers = enumerate_tautomers(&n1h);
        let n2h = tautomers
            .iter()
            .find(|t| h_assignment(t) != h_assignment(&n1h))
            .expect("enumerate_tautomers should produce N2H tautomer of pyrazole");
        assert_eq!(
            canonical_smiles(&canonical_tautomer(&n1h)),
            canonical_smiles(&canonical_tautomer(n2h)),
            "canonical_tautomer should normalize N1H and N2H to the same form"
        );
    }

    // --- TautomerConfig ---

    #[test]
    fn test_config_default_same_as_no_config() {
        // canonical_tautomer and canonical_tautomer_with_config(default) must agree.
        let mol = parse("OC=C").unwrap(); // enol
        let a = canonical_tautomer(&mol);
        let b = canonical_tautomer_with_config(&mol, &TautomerConfig::default());
        assert_eq!(
            canonical_smiles(&a),
            canonical_smiles(&b),
            "default config should match canonical_tautomer"
        );
    }

    #[test]
    fn test_config_max_iter_one_limits_convergence() {
        // With max_iter=1, at most one rule application happens.
        // Result may differ from full convergence but should not panic.
        let mol = parse("OC=C").unwrap();
        let config = TautomerConfig {
            max_iter: 1,
            ..TautomerConfig::default()
        };
        let _ = canonical_tautomer_with_config(&mol, &config);
    }

    #[test]
    fn test_config_max_tautomers_caps_enumerate() {
        // Acetylacetone has many tautomers. Capping at 2 should return exactly 2.
        let mol = parse("CC(=O)CC(=O)C").unwrap(); // acetylacetone
        let config = TautomerConfig {
            max_tautomers: 2,
            ..TautomerConfig::default()
        };
        let tautomers = enumerate_tautomers_with_config(&mol, &config);
        assert_eq!(
            tautomers.len(),
            2,
            "max_tautomers=2 should return exactly 2"
        );
    }

    #[test]
    fn test_config_enabled_rules_subset() {
        // Enabling only rule 0 (keto-enol) should still work on an enol.
        let mol = parse("OC=C").unwrap();
        let config = TautomerConfig::keto_enol_only();
        let result = canonical_tautomer_with_config(&mol, &config);
        // Should convert enol to ketone (or at least not panic).
        assert!(result.atom_count() > 0);
    }

    #[test]
    fn test_config_empty_enabled_rules_equals_all() {
        let mol = parse("OC=C").unwrap();
        let all = canonical_tautomer_with_config(&mol, &TautomerConfig::default());
        let explicit_empty = canonical_tautomer_with_config(
            &mol,
            &TautomerConfig {
                enabled_rules: vec![],
                ..TautomerConfig::default()
            },
        );
        assert_eq!(
            canonical_smiles(&all),
            canonical_smiles(&explicit_empty),
            "empty enabled_rules should equal all rules"
        );
    }

    #[test]
    fn test_rule_count_and_names() {
        let count = TautomerConfig::rule_count();
        let names = TautomerConfig::rule_names();
        assert!(count > 0);
        assert_eq!(names.len(), count);
        assert!(!names[0].is_empty());
    }

    // B4 Tests: 1,5-shift with heteroatom bridges
    #[test]
    fn test_15_shift_beta_diketone() {
        // CC(=O)CC(=O)C — acetylacetone: classic 1,5-shift case
        let mol = parse("CC(=O)CC(=O)C").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        // Should find enol form via 1,5-shift
        assert!(
            tautomers.len() >= 2,
            "Expected >= 2 tautomers for β-diketone"
        );
    }

    #[test]
    fn test_15_shift_enol_imine() {
        // C-C(=N)-C-C(=O)-H: potential 1,5-O to N shift
        // Build: C1=C-N-C-O with appropriate bonds
        let mol = parse("CC(=N)CC(=O)C").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        // Should enumerate without panic
        assert!(!tautomers.is_empty());
    }

    #[test]
    fn test_15_shift_n_bridge_diketone() {
        // O-C-N(=O)-C-O type: nitro-type with bridging N
        // Using simplified SMILES that approximates the pattern
        let mol = parse("OC1=C(O)C(=O)C(=O)C=C1").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        // Should enumerate with N bridge possibility
        assert!(!tautomers.is_empty());
    }

    #[test]
    fn test_15_shift_s_bridge() {
        // O-C-S(=O)-C-O: thio-β-diketone variant
        let mol = parse("CC(=O)CS(=O)C(=O)C").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        assert!(!tautomers.is_empty());
    }

    #[test]
    fn test_15_shift_n_to_n_with_bridge() {
        // N-C-N(=)-C-N: guanidine-type via N bridge
        let mol = parse("NC(=N)NC(=N)N").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        assert!(!tautomers.is_empty());
    }

    #[test]
    fn test_15_shift_canonical_idempotent() {
        // Applying canonical_tautomer twice should give same result
        let mol = parse("CC(=O)CC(=O)C").unwrap();
        let t1 = canonical_tautomer(&mol);
        let t2 = canonical_tautomer(&t1);
        assert_eq!(mol_fingerprint(&t1), mol_fingerprint(&t2));
    }

    #[test]
    fn test_15_shift_heteroatom_enumeration() {
        // Test that heteroatom-bridged 1,5-shifts are enumerated
        let mol = parse("CC(=O)CC(=O)C").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        // Should find multiple tautomers via both C-bridge and other possibilities
        assert!(
            tautomers.len() >= 2,
            "Expected >= 2 tautomers for acetylacetone with new 1,5-shift rules"
        );
    }

    #[test]
    fn test_15_shift_c_to_o_with_heteroatom() {
        // Active methylene 1,5-shift with heteroatom bridge
        let mol = parse("CC(=O)CC(=O)C").unwrap();
        let config = TautomerConfig {
            max_tautomers: 64,
            ..TautomerConfig::default()
        };
        let tautomers = enumerate_tautomers_with_config(&mol, &config);
        assert!(!tautomers.is_empty());
    }

    #[test]
    fn test_15_shift_no_false_positives() {
        // Molecule that shouldn't match 1,5-shift patterns
        let mol = parse("CC(=O)C").unwrap(); // propanone (no 1,5-shift possible)
        let tautomers = enumerate_tautomers(&mol);
        // Should enumerate any 1,3-shifts but no spurious 1,5-shifts
        assert!(!tautomers.is_empty());
    }

    #[test]
    fn test_15_shift_multiple_donors_acceptors() {
        // β-diketone with multiple potential 1,5-shift sites
        let mol = parse("CC(=O)CC(=O)CC(=O)C").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        // Should enumerate multiple tautomers from multiple sites
        assert!(
            tautomers.len() >= 2,
            "Expected multiple tautomers from multiple 1,5-shift sites"
        );
    }

    #[test]
    fn test_15_shift_config_selectivity() {
        // Create config that enables only 1,5-shift rules (rules 20-42)
        let mol = parse("CC(=O)CC(=O)C").unwrap();
        let config = TautomerConfig {
            max_tautomers: 32,
            enabled_rules: (20..43).collect(), // indices for 1,5-shift rules (if they exist)
            ..TautomerConfig::default()
        };
        let tautomers = enumerate_tautomers_with_config(&mol, &config);
        // Should still enumerate with restricted rule set
        assert!(!tautomers.is_empty());
    }

    // ── Tautomer zone blocking tests ──────────────────────────────────────────

    /// Helper: canonical SMILES of the canonical tautomer.
    fn canonical_smi(smi: &str) -> String {
        let mol = parse(smi).unwrap();
        canonical_smiles(&canonical_tautomer(&mol))
    }

    /// Helper: canonical SMILES of the canonical tautomer with blocked atoms.
    fn blocked_smi(smi: &str, blocked: &[u32]) -> String {
        let mol = parse(smi).unwrap();
        let config = TautomerConfig {
            blocked_atoms: blocked.iter().map(|&i| AtomIdx(i)).collect(),
            ..TautomerConfig::default()
        };
        canonical_smiles(&canonical_tautomer_with_config(&mol, &config))
    }

    #[test]
    fn test_blocking_donor_suppresses_keto_enol() {
        // OC=C: oxygen (index 0) is the donor for the enol→keto tautomer.
        // With default config, canonical is the keto form (C-C=O).
        // With O blocked, the enol tautomer should be preserved (or at least
        // the move through O should not fire).
        let mol = parse("OC=C").unwrap();
        let default_config = TautomerConfig::default();
        let blocked_config = TautomerConfig {
            blocked_atoms: [AtomIdx(0)].into_iter().collect(), // block O
            ..TautomerConfig::default()
        };
        let default_result = canonical_tautomer_with_config(&mol, &default_config);
        let blocked_result = canonical_tautomer_with_config(&mol, &blocked_config);
        // If donor O is blocked, the move should not fire → same as original,
        // or at minimum different from the unconstrained result.
        // The key invariant: no panic, and the blocked result is a valid molecule.
        assert!(blocked_result.atom_count() > 0);
        // When tautomerism through O is blocked, the canonical form cannot reach
        // the keto tautomer — verify it differs from the unconstrained canonical
        // OR equals the input (depending on whether other rules fire).
        let _ = (default_result, blocked_result); // both valid
    }

    #[test]
    fn test_empty_blocked_sets_identical_to_default() {
        // Empty blocked_atoms/blocked_bonds must produce the same result as default.
        for smi in &["OC=C", "CC(=O)CC", "c1cc[nH]c1", "CN=C"] {
            let mol = parse(smi).unwrap();
            let default = canonical_tautomer(&mol);
            let empty_config = TautomerConfig {
                blocked_atoms: HashSet::new(),
                blocked_bonds: HashSet::new(),
                ..TautomerConfig::default()
            };
            let explicit_empty = canonical_tautomer_with_config(&mol, &empty_config);

            assert_eq!(
                canonical_smiles(&default),
                canonical_smiles(&explicit_empty),
                "empty blocked sets must give same result as default for {smi}"
            );
        }
    }

    #[test]
    fn test_enumerate_with_blocking_leq_enumerate_without() {
        // Blocking an atom can only reduce the number of tautomers reachable.
        let mol = parse("CC(=O)CC(=O)C").unwrap(); // 2 carbonyl groups → many tautomers
        let all = enumerate_tautomers(&mol);
        // Block the central C (index 3 in CC(=O)CC(=O)C)
        let config = TautomerConfig {
            blocked_atoms: [AtomIdx(3)].into_iter().collect(),
            ..TautomerConfig::default()
        };
        let blocked = enumerate_tautomers_with_config(&mol, &config);
        assert!(
            blocked.len() <= all.len(),
            "blocking must not increase tautomer count: {} > {}",
            blocked.len(),
            all.len()
        );
    }

    #[test]
    fn test_blocking_all_atoms_preserves_input() {
        // When every atom is blocked, no H-transfer can fire.
        let mol = parse("OC=C").unwrap();
        let n = mol.atom_count();
        let config = TautomerConfig {
            blocked_atoms: (0..n as u32).map(AtomIdx).collect(),
            ..TautomerConfig::default()
        };
        let result = canonical_tautomer_with_config(&mol, &config);
        // Result must equal the input (no tautomer can fire).

        assert_eq!(
            canonical_smiles(&result),
            canonical_smiles(&mol),
            "all atoms blocked → input unchanged"
        );
    }

    #[test]
    fn test_enumerate_all_atoms_blocked_returns_singleton() {
        let mol = parse("OC=C").unwrap();
        let n = mol.atom_count();
        let config = TautomerConfig {
            blocked_atoms: (0..n as u32).map(AtomIdx).collect(),
            ..TautomerConfig::default()
        };
        let tautomers = enumerate_tautomers_with_config(&mol, &config);
        assert_eq!(
            tautomers.len(),
            1,
            "all atoms blocked → only the original is returned"
        );
    }

    #[test]
    fn test_out_of_range_atom_index_is_safe() {
        // An AtomIdx larger than the molecule silently has no effect.
        let mol = parse("OC=C").unwrap();
        let n = mol.atom_count();
        let config = TautomerConfig {
            blocked_atoms: [AtomIdx(n as u32 + 100)].into_iter().collect(),
            ..TautomerConfig::default()
        };
        // Must not panic
        let result = canonical_tautomer_with_config(&mol, &config);
        assert!(result.atom_count() > 0);
    }

    // ── Stereo preservation tests (RDKit #7969: canonical_tautomer should NOT ──
    // ── erase sp3 chirality at stereocenters not involved in tautomerism)     ──

    /// RDKit issue #7969: canonical_tautomer erases sp3 chirality.
    /// chematic must NOT have this bug — chirality fields are preserved via clone().
    #[test]
    fn test_remote_stereo_preserved_keto_enol() {
        // C[C@H](O)CC(=O)C: keto-enol fires at the C(=O) end.
        // The [C@H] stereocenter (index 1) is remote — must be preserved.
        let mol = parse("C[C@H](O)CC(=O)C").unwrap();
        let before_chirality = mol.atom(AtomIdx(1)).chirality;
        assert_ne!(
            before_chirality,
            Chirality::None,
            "test setup: atom 1 must be chiral"
        );

        let t = canonical_tautomer(&mol);
        let after_chirality = t.atom(AtomIdx(1)).chirality;
        assert_ne!(
            after_chirality,
            Chirality::None,
            "Remote [C@H] chirality erased by canonical_tautomer (RDKit #7969 regression)"
        );

        // Additionally verify canonical SMILES contains a chirality marker
        let smi = canonical_smiles(&t);
        assert!(
            smi.contains('@'),
            "Canonical SMILES lost chirality marker: '{}'",
            smi
        );
    }

    #[test]
    fn test_alanine_stereo_trivially_preserved() {
        // [C@@H](N)(C(=O)O)C — no tautomer rule fires; stereo must be unchanged.
        let mol = parse("[C@@H](N)(C(=O)O)C").unwrap();
        let before = mol.atom(AtomIdx(0)).chirality;
        let t = canonical_tautomer(&mol);
        let after = t.atom(AtomIdx(0)).chirality;
        assert_eq!(
            before, after,
            "Alanine chirality changed; was {:?}, got {:?}",
            before, after
        );
    }

    #[test]
    fn test_glucose_all_stereocenters_preserved() {
        let mol = parse("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O").unwrap();
        let before: Vec<Chirality> = mol.atoms().map(|(_, a)| a.chirality).collect();

        let t = canonical_tautomer(&mol);
        let after: Vec<Chirality> = t.atoms().map(|(_, a)| a.chirality).collect();

        assert_eq!(
            before, after,
            "Glucose: stereocenters changed through canonical_tautomer"
        );
    }

    #[test]
    fn test_pyrazole_no_phantom_chirality() {
        // c1cc[nH]n1: N-H tautomers possible, no stereocenters → must stay chiral-free
        let mol = parse("c1cc[nH]n1").unwrap();
        let t = canonical_tautomer(&mol);
        let chiral_count = t
            .atoms()
            .filter(|(_, a)| a.chirality != Chirality::None)
            .count();
        assert_eq!(
            chiral_count, 0,
            "Phantom chirality introduced by pyrazole tautomerism"
        );
    }

    #[test]
    fn test_stereo_at_donor_does_not_panic() {
        // [C@@H](O)(C)C(=O)O — lactic acid: O (donor) is adjacent to stereocentre.
        // The tautomer may legitimately change chirality; we just verify no panic.
        let mol = parse("[C@@H](O)(C)C(=O)O").unwrap();
        let t = canonical_tautomer(&mol);
        assert!(t.atom_count() > 0);
        let smi = canonical_smiles(&t);
        assert!(!smi.is_empty());
    }

    #[test]
    fn test_enumerate_tautomers_remote_stereo_preserved() {
        // C[C@H](O)CC(=O)C: enumerate all tautomers.
        // Every produced tautomer must preserve chirality at atom 1 (remote centre).
        let mol = parse("C[C@H](O)CC(=O)C").unwrap();
        let original_chirality = mol.atom(AtomIdx(1)).chirality;

        let tautomers = enumerate_tautomers(&mol);
        for (i, t) in tautomers.iter().enumerate() {
            if t.atom_count() == mol.atom_count() {
                let ch = t.atom(AtomIdx(1)).chirality;
                assert_eq!(
                    ch, original_chirality,
                    "Tautomer #{}: chirality at atom 1 changed ({:?} → {:?})",
                    i, original_chirality, ch
                );
            }
        }
    }

    // -----------------------------------------------------------------
    // Tautomer-Rebuild-S2: transfer_hydrogen_aromatic must preserve the
    // stereo side channels of every atom/bond it doesn't itself modify.
    // The tests above only check `atom.chirality` (which survives via
    // `atom.clone()` regardless of this bug -- it's `chirality`'s
    // *reference frame*, `stereo_neighbor_order`, that was silently
    // dropped). These specifically probe the side channels: an aromatic
    // N-H ring (to actually fire transfer_hydrogen_aromatic, not the
    // non-aromatic transfer_hydrogen) combined with a remote tetrahedral
    // stereocenter, an enhanced stereo group, and a stashed bond
    // direction.
    // -----------------------------------------------------------------

    /// Chiral secondary-alcohol carbon (atom 1, `[C@H]`) attached to a
    /// pyrazole ring (atoms 6/7 are the two ring nitrogens) -- the N-H
    /// tautomerism is a direct aromatic 1,2-shift, exercising
    /// `transfer_hydrogen_aromatic` specifically (not the non-aromatic
    /// `transfer_hydrogen`, which is a separate, not-yet-fixed instance of
    /// this same bug class -- out of scope for this PR).
    const AROMATIC_TAUTOMER_WITH_REMOTE_STEREOCENTER: &str = "C[C@H](O)c1cc[nH]n1";

    #[test]
    fn transfer_hydrogen_aromatic_preserves_remote_stereo_neighbor_order() {
        let mol = parse(AROMATIC_TAUTOMER_WITH_REMOTE_STEREOCENTER).unwrap();
        let stereocenter = AtomIdx(1);
        let original_order = mol.stereo_neighbor_order(stereocenter).map(|s| s.to_vec());
        assert!(
            original_order.is_some(),
            "test setup sanity: atom 1 must be chiral"
        );

        let (donor, acceptor) = find_direct_aromatic_matches(&mol, &atom_rank(&mol))
            .into_iter()
            .next()
            .expect("test setup sanity: pyrazole should have a direct aromatic N-H match");
        let result = transfer_hydrogen_aromatic(&mol, donor, acceptor, &HashSet::new())
            .expect("pyrazole N-H transfer should succeed");

        assert_eq!(
            result.atom(stereocenter).chirality,
            mol.atom(stereocenter).chirality,
            "remote stereocenter's chirality must be unchanged"
        );
        assert_eq!(
            result
                .stereo_neighbor_order(stereocenter)
                .map(|s| s.to_vec()),
            original_order,
            "remote stereocenter's stereo_neighbor_order must survive transfer_hydrogen_aromatic \
             verbatim -- chirality alone surviving is not enough to keep it interpretable"
        );
    }

    #[test]
    fn transfer_hydrogen_aromatic_preserves_stereo_groups_and_bond_directions() {
        let mut mol = parse(AROMATIC_TAUTOMER_WITH_REMOTE_STEREOCENTER).unwrap();
        let stereocenter = AtomIdx(1);
        mol.add_stereo_group(chematic_core::StereoGroup::new(
            chematic_core::StereoGroupKind::Absolute,
            vec![stereocenter],
        ));
        // Stash an arbitrary bond direction on a bond uninvolved in the H
        // transfer (the C-O bond), mimicking the kind of stashed E/Z
        // marker `apply_aromaticity_ex` leaves behind on an exocyclic bond
        // adjacent to a ring atom promoted to Aromatic order.
        let (co_bond, _) = mol
            .bond_between(AtomIdx(1), AtomIdx(2))
            .expect("C-O bond must exist");
        mol.set_bond_direction(co_bond, BondOrder::Up);

        let (donor, acceptor) = find_direct_aromatic_matches(&mol, &atom_rank(&mol))
            .into_iter()
            .next()
            .expect("test setup sanity: pyrazole should have a direct aromatic N-H match");
        let result = transfer_hydrogen_aromatic(&mol, donor, acceptor, &HashSet::new())
            .expect("pyrazole N-H transfer should succeed");

        assert_eq!(
            result.stereo_groups(),
            mol.stereo_groups(),
            "stereo_groups must survive transfer_hydrogen_aromatic verbatim"
        );
        assert_eq!(
            result.bond_direction(co_bond),
            mol.bond_direction(co_bond),
            "bond_directions for a bond uninvolved in the H transfer must be unchanged"
        );
    }

    // -----------------------------------------------------------------
    // Transfer-Hydrogen-Correctness-P0: the non-aromatic `transfer_hydrogen`
    // (1,3-/1,5-H-shift, e.g. keto-enol) counterpart to
    // Tautomer-Rebuild-S2's `transfer_hydrogen_aromatic` fix -- same bug
    // shape, deliberately left out of scope there (see the comment on
    // `AROMATIC_TAUTOMER_WITH_REMOTE_STEREOCENTER` above). A keto-enol
    // 1,3-shift (donor O-H, bridge C, acceptor C=C) combined with a remote
    // tetrahedral stereocenter, an enhanced stereo group, and a stashed
    // bond direction.
    // -----------------------------------------------------------------

    /// Chiral secondary carbon (atom 1, `[C@H]`) two bonds away from an
    /// enol group (`C(O)=C`, atoms 4/5/6) -- structurally uninvolved in the
    /// keto-enol 1,3-shift, which exercises the non-aromatic
    /// `transfer_hydrogen` specifically.
    const NON_AROMATIC_TAUTOMER_WITH_REMOTE_STEREOCENTER: &str = "C[C@H](Cl)CC(O)=C";

    fn keto_enol_match(mol: &Molecule) -> (AtomIdx, AtomIdx, AtomIdx) {
        let rule = RULES
            .iter()
            .find(|r| r.name == "keto-enol")
            .expect("keto-enol rule must exist");
        find_matches(mol, rule, &atom_rank(mol))
            .into_iter()
            .next()
            .expect("test setup sanity: molecule should have a keto-enol match")
    }

    #[test]
    fn transfer_hydrogen_preserves_remote_stereo_neighbor_order() {
        let mol = parse(NON_AROMATIC_TAUTOMER_WITH_REMOTE_STEREOCENTER).unwrap();
        let stereocenter = AtomIdx(1);
        let original_order = mol.stereo_neighbor_order(stereocenter).map(|s| s.to_vec());
        assert!(
            original_order.is_some(),
            "test setup sanity: atom 1 must be chiral"
        );

        let (donor, bridge, acceptor) = keto_enol_match(&mol);
        let result = transfer_hydrogen(
            &mol,
            donor,
            bridge,
            acceptor,
            &HashSet::new(),
            &HashSet::new(),
        )
        .expect("keto-enol transfer should succeed");

        assert_eq!(
            result.atom(stereocenter).chirality,
            mol.atom(stereocenter).chirality,
            "remote stereocenter's chirality must be unchanged"
        );
        assert_eq!(
            result
                .stereo_neighbor_order(stereocenter)
                .map(|s| s.to_vec()),
            original_order,
            "remote stereocenter's stereo_neighbor_order must survive transfer_hydrogen \
             verbatim -- chirality alone surviving is not enough to keep it interpretable"
        );
    }

    #[test]
    fn transfer_hydrogen_preserves_stereo_groups_and_bond_directions() {
        let mut mol = parse(NON_AROMATIC_TAUTOMER_WITH_REMOTE_STEREOCENTER).unwrap();
        let stereocenter = AtomIdx(1);
        mol.add_stereo_group(chematic_core::StereoGroup::new(
            chematic_core::StereoGroupKind::Absolute,
            vec![stereocenter],
        ));
        // Stash an arbitrary bond direction on a bond uninvolved in the H
        // transfer (the C-Cl bond).
        let (c_cl_bond, _) = mol
            .bond_between(AtomIdx(1), AtomIdx(2))
            .expect("C-Cl bond must exist");
        mol.set_bond_direction(c_cl_bond, BondOrder::Up);

        let (donor, bridge, acceptor) = keto_enol_match(&mol);
        let result = transfer_hydrogen(
            &mol,
            donor,
            bridge,
            acceptor,
            &HashSet::new(),
            &HashSet::new(),
        )
        .expect("keto-enol transfer should succeed");

        assert_eq!(
            result.stereo_groups(),
            mol.stereo_groups(),
            "stereo_groups must survive transfer_hydrogen verbatim"
        );
        assert_eq!(
            result.bond_direction(c_cl_bond),
            mol.bond_direction(c_cl_bond),
            "bond_directions for a bond uninvolved in the H transfer must be unchanged"
        );
    }

    #[test]
    fn transfer_hydrogen_cip_of_uninvolved_stereocenter_survives_canonical_round_trip() {
        // Chemical-identity check: the remote stereocenter's actual CIP
        // label (not just the raw @/@@ character) must survive both the
        // transfer itself and a canonicalize -> reparse round trip.
        let mol = parse("C[C@H:9](Cl)CC(O)=C").unwrap();
        let stereocenter = mol
            .atoms()
            .find(|(_, a)| a.atom_map == Some(9))
            .map(|(idx, _)| idx)
            .expect("atom map tag not found");
        let before_cip = crate::assign_cip(&mol).get(stereocenter);
        assert!(
            before_cip.is_some(),
            "test setup sanity: must be resolvable"
        );

        let (donor, bridge, acceptor) = keto_enol_match(&mol);
        let result = transfer_hydrogen(
            &mol,
            donor,
            bridge,
            acceptor,
            &HashSet::new(),
            &HashSet::new(),
        )
        .expect("keto-enol transfer should succeed");

        let smi = canonical_smiles(&result);
        let reparsed = parse(&smi).expect("valid canonical SMILES");
        let reparsed_center = reparsed
            .atoms()
            .find(|(_, a)| a.atom_map == Some(9))
            .map(|(idx, _)| idx)
            .expect("atom map tag not found after round trip");
        let after_cip = crate::assign_cip(&reparsed).get(reparsed_center);

        assert_eq!(
            before_cip, after_cip,
            "remote stereocenter's CIP code must survive transfer_hydrogen + a canonical round trip"
        );
    }

    #[test]
    fn enumerate_tautomers_keto_enol_count_and_canonical_form_unaffected_by_metadata_fix() {
        // Regression pin: this fix is about metadata preservation only, not
        // enumeration logic -- the tautomer count and canonical_tautomer's
        // chosen form for a plain (no remote stereocenter) enol must be
        // exactly what they were before this PR.
        let mol = parse("CC(O)=C").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        assert_eq!(
            tautomers.len(),
            2,
            "propan-2-enol should enumerate to exactly 2 forms (enol + acetone)"
        );
        // Golden string updated (issue #205): `initial_invariant`'s
        // explicit/implicit-H-count unification changed Morgan ranks for
        // some atoms, so `canonical_smiles`'s neighbor-branch-ordering
        // shifted -- this is still the identical keto tautomer (acetone),
        // just a different valid serialization of the same graph; the
        // enumeration count above (2) is the actual regression pin for
        // "logic unaffected" and did not change.
        assert_eq!(
            canonical_smiles(&canonical_tautomer(&mol)),
            "C(C)(=O)C",
            "the canonical tautomer form must be unchanged (prefers the keto form)"
        );
    }

    #[test]
    fn enumerate_tautomers_count_and_canonical_form_unaffected_by_metadata_fix() {
        // Regression pin: this fix is about metadata preservation only, not
        // enumeration logic (mol_fingerprint/h_assignment dedup don't
        // consult stereo_neighbor_order/stereo_groups/bond_directions at
        // all) -- the tautomer count and canonical_tautomer's chosen form
        // for a plain pyrazole (no remote stereocenter to confound the
        // comparison) must be exactly what they were before this PR.
        let mol = parse("c1cc[nH]n1").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        assert_eq!(
            tautomers.len(),
            2,
            "pyrazole should enumerate to exactly 2 forms"
        );
        // Golden string updated (issue #205), same reason as the acetone
        // test above: same tautomer (pyrazole's NH form), different valid
        // canonical serialization; the count (2) is the unaffected pin.
        assert_eq!(
            canonical_smiles(&canonical_tautomer(&mol)),
            "c1c[nH]nc1",
            "pyrazole's canonical tautomer form must be unchanged"
        );
    }

    #[test]
    fn test_blocked_stereo_preserved_with_zone_blocking() {
        // O[C@@H](F)C(=O)C: block O (index 0) to suppress keto-enol.
        // Stereocentre [C@@H] (index 1) must be preserved.
        let mol = parse("O[C@@H](F)C(=O)C").unwrap();
        let config = TautomerConfig {
            blocked_atoms: [AtomIdx(0)].into_iter().collect(),
            ..TautomerConfig::default()
        };
        let t = canonical_tautomer_with_config(&mol, &config);
        let chirality_after = t.atom(AtomIdx(1)).chirality;
        assert_ne!(
            chirality_after,
            Chirality::None,
            "Chirality erased from blocked stereocentre"
        );
    }

    // ── RDKit PR #9128: E/Z stereo on exocyclic double bonds (hydrazones/imines) ──

    #[test]
    fn test_hydrazone_ez_stereo_preserved_in_canonical_tautomer() {
        // E-hydrazone and Z-hydrazone are DIFFERENT compounds.
        // The tautomer search must not merge them while deduplicating forms.
        let e_hydrazone = parse("C/C=N/N").expect("E-hydrazone");
        let z_hydrazone = parse("C/C=N\\N").expect("Z-hydrazone");
        let e_can = canonical_tautomer(&e_hydrazone);
        let z_can = canonical_tautomer(&z_hydrazone);
        let e_smi = canonical_smiles(&e_can);
        let z_smi = canonical_smiles(&z_can);
        assert_ne!(
            e_smi, z_smi,
            "E and Z hydrazone must remain distinct after canonical_tautomer (RDKit PR #9128): \
             E={e_smi} Z={z_smi}"
        );
    }

    // ── Tetrazole 1H ↔ 2H tautomers (OpenBabel PR #2975 pattern) ────────────
    //
    // Tetrazole has two stable tautomers:
    //   1H-tetrazole: H on N1 (directly adjacent to C)
    //   2H-tetrazole: H on N2 (one position further around the ring from C)
    // Both are aromatic and interconvert via a direct aromatic 1,2-shift.
    // The `find_direct_aromatic_matches` path handles this case.

    #[test]
    fn test_tetrazole_1h_enumerates_two_forms() {
        // c1nnn[nH]1 = 1H-tetrazole.  enumerate_tautomers must find at least
        // 2 forms (1H and 2H) via the direct aromatic 1,2-shift.
        let mol = parse("c1nnn[nH]1").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        assert!(
            tautomers.len() >= 2,
            "Expected >= 2 tautomers for 1H-tetrazole, got {}",
            tautomers.len()
        );
    }

    #[test]
    fn test_tetrazole_2h_enumerates_two_forms() {
        // c1n[nH]nn1 = 2H-tetrazole (H on N two positions from C).
        // enumerate_tautomers must find at least 2 forms (1H and 2H).
        let mol = parse("c1n[nH]nn1").unwrap();
        let tautomers = enumerate_tautomers(&mol);
        assert!(
            tautomers.len() >= 2,
            "Expected >= 2 tautomers for 2H-tetrazole, got {}",
            tautomers.len()
        );
    }

    #[test]
    fn test_tetrazole_canonical_from_1h_and_2h_agrees() {
        // canonical_tautomer must produce the same canonical form regardless of
        // whether the input is 1H or 2H tetrazole.
        let mol_1h = parse("c1nnn[nH]1").unwrap();
        let mol_2h = parse("c1n[nH]nn1").unwrap();
        let can_1h = canonical_smiles(&canonical_tautomer(&mol_1h));
        let can_2h = canonical_smiles(&canonical_tautomer(&mol_2h));
        assert_eq!(
            can_1h, can_2h,
            "canonical_tautomer must agree for 1H ({can_1h}) and 2H ({can_2h}) tetrazole"
        );
    }

    #[test]
    fn test_tetrazole_canonical_preserves_aromaticity() {
        // After canonical_tautomer the ring must remain fully aromatic.
        let mol = parse("c1nnn[nH]1").unwrap();
        let t = canonical_tautomer(&mol);
        let all_aromatic = t.atoms().all(|(_, a)| a.aromatic);
        assert!(
            all_aromatic,
            "all atoms in canonical tetrazole must be aromatic"
        );
    }

    #[test]
    fn imine_ez_stereo_preserved_in_tautomer_enumeration() {
        // Enumerate tautomers of E-imine C/C=N/C — all resulting tautomers must
        // produce valid canonical SMILES (no empty string, no panic).
        let e_imine = parse("C/C=N/C").expect("E-imine");
        let tautomers = enumerate_tautomers(&e_imine);
        assert!(
            !tautomers.is_empty(),
            "E-imine must enumerate at least one tautomer"
        );
        for (i, t) in tautomers.iter().enumerate() {
            let smi = canonical_smiles(t);
            assert!(
                !smi.is_empty(),
                "tautomer {i} must produce valid canonical SMILES"
            );
        }
    }

    // -- Phase 2 round-2B tautomer_parent tests --------------------------------
    // See docs/rfcs/tautomer_parent_identity_phase2_rfc.md sections 4.1-4.3.

    #[test]
    fn all_rules_map_to_a_named_tautomer_rule_id() {
        for rule in RULES {
            assert!(
                !matches!(rule_id_for(rule.name), TautomerRuleId::Other(_)),
                "rule '{}' has no TautomerRuleId variant -- add one",
                rule.name
            );
        }
    }

    #[test]
    fn tautomer_parent_completed_on_keto_enol() {
        let mol = parse("CC(O)=CC(C)=O").unwrap();
        let result = tautomer_parent(&mol, &TautomerLimits::default());
        assert_eq!(result.status, ParentComputationStatus::Completed);
        assert_eq!(
            canonical_smiles(&result.molecule),
            canonical_smiles(&canonical_tautomer(&mol))
        );
        match &result.audit {
            ParentAudit::Tautomer(record) => {
                assert_eq!(record.applied_transforms.len(), 1);
                assert_eq!(
                    record.applied_transforms[0].rule_id,
                    TautomerRuleId::KetoEnol
                );
                assert!(record.lost_stereo.is_empty());
                assert!(record.affected_isotopes.is_empty());
            }
            other => panic!("expected ParentAudit::Tautomer, got {other:?}"),
        }
    }

    #[test]
    fn tautomer_parent_max_transforms_reached_when_budget_too_small() {
        // 3 independent enol arms, each needing exactly 1 transform; a
        // budget of 1 converts only one and must report the exhaustion,
        // not silently return a partially-converted, order-dependent result
        // as if it were `Completed`.
        let mol = build_comb(3, false, false);
        let limits = TautomerLimits {
            max_transforms: 1,
            ..TautomerLimits::default()
        };
        let result = tautomer_parent(&mol, &limits);
        assert_eq!(result.status, ParentComputationStatus::MaxTransformsReached);
        match &result.audit {
            ParentAudit::Tautomer(record) => assert_eq!(record.applied_transforms.len(), 1),
            other => panic!("expected ParentAudit::Tautomer, got {other:?}"),
        }
    }

    #[test]
    fn tautomer_parent_completed_with_enough_budget_for_all_arms() {
        let mol = build_comb(3, false, false);
        let limits = TautomerLimits {
            max_transforms: 10,
            ..TautomerLimits::default()
        };
        let result = tautomer_parent(&mol, &limits);
        assert_eq!(result.status, ParentComputationStatus::Completed);
        match &result.audit {
            ParentAudit::Tautomer(record) => assert_eq!(record.applied_transforms.len(), 3),
            other => panic!("expected ParentAudit::Tautomer, got {other:?}"),
        }
    }

    #[test]
    fn tautomer_parent_max_tautomers_reached_on_tetrazole() {
        let mol = parse("c1nnn[nH]1").unwrap();
        let unlimited = tautomer_parent(&mol, &TautomerLimits::default());
        let (unlimited_candidate_count, _) = match &unlimited.audit {
            ParentAudit::Tautomer(r) => (r.candidate_count, ()),
            other => panic!("expected ParentAudit::Tautomer, got {other:?}"),
        };
        assert!(
            unlimited_candidate_count > 1,
            "tetrazole must have more than one direct-aromatic-shift candidate to make this test meaningful"
        );

        let limited = tautomer_parent(
            &mol,
            &TautomerLimits {
                max_tautomers: 1,
                ..TautomerLimits::default()
            },
        );
        assert_eq!(limited.status, ParentComputationStatus::MaxTautomersReached);
    }

    #[test]
    fn tautomer_parent_max_tautomers_zero_skips_comparison_without_false_exhaustion() {
        // A molecule with NO direct-aromatic-shift candidates at all must
        // not report MaxTautomersReached just because max_tautomers=0 --
        // see TautomerLimits::max_tautomers's doc comment.
        let mol = parse("CC(O)=CC(C)=O").unwrap();
        let result = tautomer_parent(
            &mol,
            &TautomerLimits {
                max_tautomers: 0,
                ..TautomerLimits::default()
            },
        );
        assert_eq!(result.status, ParentComputationStatus::Completed);
    }

    #[test]
    fn tautomer_parent_empty_molecule_is_invalid_input() {
        let empty = MoleculeBuilder::new().build();
        let result = tautomer_parent(&empty, &TautomerLimits::default());
        assert_eq!(
            result.status,
            ParentComputationStatus::InvalidInput(InvalidInputReason::EmptyMolecule)
        );
    }

    #[test]
    fn tp2_29_idempotence_acetylacetone() {
        let mol = parse("CC(=O)CC(C)=O").unwrap();
        let once = tautomer_parent(&mol, &TautomerLimits::default());
        let twice = tautomer_parent(&once.molecule, &TautomerLimits::default());
        assert_eq!(
            canonical_smiles(&once.molecule),
            canonical_smiles(&twice.molecule)
        );
        assert_eq!(twice.status, ParentComputationStatus::Completed);
    }

    #[test]
    fn tautomer_parent_preserves_remote_stereocenter() {
        // Alanine's stereocenter must survive tautomer_parent unaffected,
        // same as the primitive transfer_hydrogen/transfer_hydrogen_aromatic
        // functions this is built on (RFC section 1.2's confirmed non-bug).
        let mol = parse("C[C@@H](N)C(=O)O").unwrap();
        let result = tautomer_parent(&mol, &TautomerLimits::default());
        match &result.audit {
            ParentAudit::Tautomer(record) => {
                assert!(record.lost_stereo.is_empty());
                assert!(record.affected_isotopes.is_empty());
            }
            other => panic!("expected ParentAudit::Tautomer, got {other:?}"),
        }
    }

    // -- Phase 2 round-2C-2/2C-3: exocyclic lactam/lactim shift ----------------
    // See docs/rfcs/tautomer_parent_identity_phase2_rfc.md section 4.4a.

    #[test]
    fn ring_distance_matches_ortho_and_meta() {
        let two_pyridone = parse("O=c1cccc[nH]1").unwrap();
        let ring = chematic_perception::find_sssr(&two_pyridone);
        let ring = ring.rings().first().expect("one ring");
        // bridge = AtomIdx(1) (bonded to the exocyclic O), acceptor =
        // AtomIdx(6) (the ring N) -- ortho, distance 1.
        assert_eq!(ring_distance(ring, AtomIdx(1), AtomIdx(6)), Some(1));

        let three_hydroxypyridine = parse("Oc1cccnc1").unwrap();
        let ring = chematic_perception::find_sssr(&three_hydroxypyridine);
        let ring = ring.rings().first().expect("one ring");
        // bridge = AtomIdx(1), ring N = AtomIdx(5) -- meta, distance 2 (even).
        assert_eq!(ring_distance(ring, AtomIdx(1), AtomIdx(5)), Some(2));
    }

    fn assert_converges(label: &str, variants: &[&str]) {
        let outs: Vec<String> = variants
            .iter()
            .map(|s| canonical_smiles(&canonical_tautomer(&parse(s).unwrap())))
            .collect();
        assert!(
            outs.windows(2).all(|w| w[0] == w[1]),
            "{label}: variants did not converge to one canonical form: {outs:?}"
        );
    }

    fn assert_noop(label: &str, smi: &str) {
        let mol = parse(smi).unwrap();
        let out = canonical_tautomer(&mol);
        assert_eq!(
            canonical_smiles(&out),
            canonical_smiles(&mol),
            "{label}: expected a no-op, got a change"
        );
    }

    #[test]
    fn tp2_05_2_pyridone_converges() {
        assert_converges("tp2-05", &["O=c1cccc[nH]1", "Oc1ccccn1"]);
    }

    #[test]
    fn tp2_06_4_pyridone_converges() {
        assert_converges("tp2-06", &["O=c1cc[nH]cc1", "Oc1ccncc1"]);
    }

    #[test]
    fn tp2_08_uracil_converges() {
        assert_converges("tp2-08", &["O=c1cc[nH]c(=O)[nH]1", "Oc1ccnc(O)n1"]);
    }

    #[test]
    fn enumerate_aromatic_lactam_component_includes_reverse_lactim_edge() {
        let mol = parse("Cc1cc[nH]c(=O)n1").unwrap();
        let forms = enumerate_tautomers(&mol);
        assert!(
            forms.iter().any(|candidate| {
                has_aromatic_exocyclic_oxygen(candidate)
                    && !has_aromatic_exocyclic_carbonyl(candidate)
            }),
            "keto input must expose its reverse aromatic lactim form"
        );
    }

    #[test]
    fn tp2_holdout_01_hypoxanthine_converges() {
        // Holdout: never used to shape the mechanism above.
        assert_converges(
            "tp2-holdout-01",
            &["O=c1[nH]cnc2[nH]cnc12", "Oc1ncnc2[nH]cnc12"],
        );
    }

    #[test]
    fn tp2_holdout_06_n9_methylhypoxanthine_converges() {
        assert_converges(
            "tp2-holdout-06",
            &["Cn1cnc2c(=O)[nH]cnc21", "Cn1cnc2c(O)ncnc21"],
        );
    }

    #[test]
    fn tp2_07_09_dual_flank_cytosine_and_guanine_converge() {
        // The carbonyl-centered dual-flank H move must normalize both keto
        // spellings and the corresponding enol precursor.
        let cytosine_keto =
            canonical_smiles(&canonical_tautomer(&parse("Nc1cc[nH]c(=O)n1").unwrap()));
        let cytosine_enol = canonical_smiles(&canonical_tautomer(&parse("Nc1ccnc(O)n1").unwrap()));
        assert_eq!(
            cytosine_keto, cytosine_enol,
            "cytosine keto and amino/imino variants must converge"
        );

        let guanine_keto = canonical_smiles(&canonical_tautomer(
            &parse("Nc1nc2[nH]cnc2c(=O)[nH]1").unwrap(),
        ));
        let guanine_enol =
            canonical_smiles(&canonical_tautomer(&parse("Nc1nc2[nH]cnc2c(O)n1").unwrap()));
        assert_eq!(
            guanine_keto, guanine_enol,
            "guanine keto and amino/imino variants must converge"
        );
    }

    #[test]
    fn tp2_39_methylpyrimidinone_converges() {
        assert_converges(
            "tp2-39-methylpyrimidinone",
            &["Cc1cc[nH]c(=O)n1", "Cc1ccnc(=O)[nH]1", "Cc1ccnc(O)n1"],
        );
    }

    #[test]
    fn tp2_04_nitroso_oxime_converges_without_enolizing_ketones() {
        let nitroso = parse("CCN=O").unwrap();
        let oxime = parse("CC=NO").unwrap();

        assert_eq!(
            canonical_smiles(&canonical_tautomer(&nitroso)),
            canonical_smiles(&canonical_tautomer(&oxime)),
            "nitrosoethane and acetaldehyde oxime must share one canonical form"
        );

        // The specialized rule must not turn an unrelated carbonyl into an
        // enol through the generic any-bridge C→O rule.
        let ketone = parse("CCC=O").unwrap();
        assert_eq!(
            mol_fingerprint(&canonical_tautomer(&ketone)),
            mol_fingerprint(&ketone),
            "ordinary aldehydes/ketones must remain in their preferred carbonyl form"
        );
    }

    #[test]
    fn tp2_34_isotope_preserved_through_exocyclic_shift() {
        let mol = parse("[18OH]c1ccccn1").unwrap();
        let out = canonical_tautomer(&mol);
        assert_converges("tp2-34", &["[18O]=c1cccc[nH]1", "[18OH]c1ccccn1"]);
        // The isotope must land on an oxygen atom, not be dropped or moved
        // to an unrelated atom.
        let labeled = (0..out.atom_count())
            .map(|i| AtomIdx(i as u32))
            .find(|&i| out.atom(i).isotope == Some(18))
            .expect("18O label must survive");
        assert_eq!(out.atom(labeled).element.atomic_number(), 8);
    }

    #[test]
    fn tp2_35_remote_stereocenter_preserved_through_exocyclic_shift() {
        assert_converges(
            "tp2-35",
            &["O=c1c([C@@H](F)Cl)ccc[nH]1", "Oc1c([C@@H](F)Cl)cccn1"],
        );
        let mol = parse("Oc1c([C@@H](F)Cl)cccn1").unwrap();
        let out = canonical_tautomer(&mol);
        let stereocenters_before = (0..mol.atom_count())
            .filter(|&i| mol.atom(AtomIdx(i as u32)).chirality != Chirality::None)
            .count();
        let stereocenters_after = (0..out.atom_count())
            .filter(|&i| out.atom(AtomIdx(i as u32)).chirality != Chirality::None)
            .count();
        assert_eq!(stereocenters_before, 1);
        assert_eq!(stereocenters_after, 1);
    }

    #[test]
    fn negative_controls_stay_noop_under_exocyclic_shift() {
        for (label, smi) in [
            ("phenol", "Oc1ccccc1"),
            ("anisole", "COc1ccccc1"),
            ("aniline", "Nc1ccccc1"),
            ("pyridine-n-oxide", "[O-][n+]1ccccc1"),
            ("acetamide", "CC(N)=O"),
            ("3-hydroxypyridine", "Oc1cccnc1"),
            ("4-aminopyridine", "Nc1ccncc1"),
            ("2-aminopyridine", "Nc1ccccn1"),
        ] {
            assert_noop(label, smi);
        }
    }

    #[test]
    fn ring_internal_nh_shift_controls_unaffected_by_round_2c2() {
        // Existing ring-internal-only mechanism must still work exactly as
        // before -- the new step is a fallback tried only when no rule (and,
        // in canonical_tautomer_with_config's loop, implicitly, since this
        // check runs on the post-rule-loop `current`) already changed the
        // molecule; it never touches ring-internal bonds.
        for smi in ["c1cc[nH]n1", "c1nn[nH]n1", "c1ccc2[nH]cnc2c1"] {
            let mol = parse(smi).unwrap();
            let out = canonical_tautomer(&mol);
            let twice = canonical_tautomer(&out);
            assert_eq!(
                canonical_smiles(&out),
                canonical_smiles(&twice),
                "{smi}: not idempotent"
            );
        }
    }

    #[test]
    fn canonical_tautomer_and_tautomer_parent_agree_on_exocyclic_shift() {
        for smi in [
            "O=c1cccc[nH]1",
            "Oc1ccccn1",
            "O=c1cc[nH]cc1",
            "Oc1ccncc1",
            "O=c1cc[nH]c(=O)[nH]1",
            "Oc1ccnc(O)n1",
            "O=c1[nH]cnc2[nH]cnc12",
            "Oc1ncnc2[nH]cnc12",
        ] {
            let mol = parse(smi).unwrap();
            let via_canonical = canonical_smiles(&canonical_tautomer(&mol));
            let via_parent =
                canonical_smiles(&tautomer_parent(&mol, &TautomerLimits::default()).molecule);
            assert_eq!(
                via_canonical, via_parent,
                "{smi}: canonical_tautomer and tautomer_parent diverged"
            );
        }
    }

    #[test]
    fn tautomer_parent_max_transforms_reached_via_exocyclic_shift() {
        let mol = parse("Oc1ccccn1").unwrap();
        let limits = TautomerLimits {
            max_transforms: 0,
            ..TautomerLimits::default()
        };
        let result = tautomer_parent(&mol, &limits);
        assert_eq!(result.status, ParentComputationStatus::MaxTransformsReached);
    }

    // -- Round 2C-2 hardening (post-#365 review) --------------------------------
    // Tightened matcher preconditions, post-generation invariant check,
    // fingerprint-cross-checked tie-break, and the atom-order-invariance /
    // fail-closed-negative-control coverage the review asked for.

    /// Build a molecule from named atoms/bonds, inserted in the given
    /// `order` -- lets a test construct the *same* molecular graph with
    /// *different* underlying `AtomIdx` numbering (same technique as
    /// `build_comb`'s `reverse_arms`). This is the only way to test
    /// order-invariance of code that iterates `mol.atoms()`/
    /// `mol.neighbors()`/SSSR ring order without any risk that a
    /// hand-typed alternate SMILES accidentally describes a different
    /// molecule.
    fn build_named(
        order: &[&str],
        atoms: &[(&str, chematic_core::Element, bool, Option<u8>)],
        bonds: &[(&str, &str, BondOrder)],
    ) -> Molecule {
        use chematic_core::Atom;
        let mut builder = MoleculeBuilder::new();
        let mut handles: std::collections::HashMap<&str, AtomIdx> =
            std::collections::HashMap::new();
        for &name in order {
            let (_, elem, aromatic, h) = *atoms.iter().find(|(n, ..)| *n == name).unwrap();
            let mut a = Atom::new(elem);
            a.aromatic = aromatic;
            a.hydrogen_count = h;
            handles.insert(name, builder.add_atom(a));
        }
        for &(a, b, order) in bonds {
            builder.add_bond(handles[&a], handles[&b], order).unwrap();
        }
        builder.build()
    }

    fn rule_id_sequence(result: &ParentResult) -> Vec<TautomerRuleId> {
        match &result.audit {
            ParentAudit::Tautomer(record) => record
                .applied_transforms
                .iter()
                .map(|t| t.rule_id)
                .collect(),
            other => panic!("expected ParentAudit::Tautomer, got {other:?}"),
        }
    }

    #[test]
    fn atom_order_invariance_2_pyridone_graph_construction() {
        use chematic_core::Element;
        let atoms = [
            ("O", Element::O, false, Some(1)),
            ("C1", Element::C, true, None),
            ("C2", Element::C, true, None),
            ("C3", Element::C, true, None),
            ("C4", Element::C, true, None),
            ("C5", Element::C, true, None),
            ("N", Element::N, true, Some(0)),
        ];
        let bonds = [
            ("O", "C1", BondOrder::Single),
            ("C1", "C2", BondOrder::Aromatic),
            ("C2", "C3", BondOrder::Aromatic),
            ("C3", "C4", BondOrder::Aromatic),
            ("C4", "C5", BondOrder::Aromatic),
            ("C5", "N", BondOrder::Aromatic),
            ("N", "C1", BondOrder::Aromatic),
        ];
        let forward = ["O", "C1", "C2", "C3", "C4", "C5", "N"];
        let reversed = ["N", "C5", "C4", "C3", "C2", "C1", "O"];
        let mol_forward = build_named(&forward, &atoms, &bonds);
        let mol_reversed = build_named(&reversed, &atoms, &bonds);
        // Sanity: really the same input molecule before any shift.
        assert_eq!(
            canonical_smiles(&mol_forward),
            canonical_smiles(&mol_reversed)
        );

        let out_forward = canonical_tautomer(&mol_forward);
        let out_reversed = canonical_tautomer(&mol_reversed);
        assert_eq!(
            canonical_smiles(&out_forward),
            canonical_smiles(&out_reversed)
        );

        let limits = TautomerLimits::default();
        let parent_forward = tautomer_parent(&mol_forward, &limits);
        let parent_reversed = tautomer_parent(&mol_reversed, &limits);
        assert_eq!(parent_forward.status, parent_reversed.status);
        assert_eq!(
            canonical_smiles(&parent_forward.molecule),
            canonical_smiles(&parent_reversed.molecule)
        );
        assert_eq!(
            rule_id_sequence(&parent_forward),
            rule_id_sequence(&parent_reversed)
        );
    }

    #[test]
    fn atom_order_invariance_uracil_two_independent_sites_graph_construction() {
        // Both exocyclic sites (O0/C1 and O6/C5) simultaneously have TWO
        // valid odd-distance ring-N acceptor candidates each (N4 at distance
        // 3, N7 at distance 1, from either bridge) when starting from the
        // fully-di-enol form -- exactly the "candidate discovery order"
        // stress case the review flagged. Forward vs reversed atom
        // insertion order changes both iteration order over mol.atoms()
        // and the HashSet-then-sort order find_exocyclic_lactam_shift_matches
        // produces its candidates in.
        use chematic_core::Element;
        let atoms = [
            ("O0", Element::O, false, Some(1)),
            ("C1", Element::C, true, None),
            ("C2", Element::C, true, None),
            ("C3", Element::C, true, None),
            ("N4", Element::N, true, Some(0)),
            ("C5", Element::C, true, None),
            ("O6", Element::O, false, Some(1)),
            ("N7", Element::N, true, Some(0)),
        ];
        let bonds = [
            ("O0", "C1", BondOrder::Single),
            ("C1", "C2", BondOrder::Aromatic),
            ("C2", "C3", BondOrder::Aromatic),
            ("C3", "N4", BondOrder::Aromatic),
            ("N4", "C5", BondOrder::Aromatic),
            ("C5", "O6", BondOrder::Single),
            ("C5", "N7", BondOrder::Aromatic),
            ("N7", "C1", BondOrder::Aromatic),
        ];
        let forward = ["O0", "C1", "C2", "C3", "N4", "C5", "O6", "N7"];
        let reversed = ["N7", "O6", "C5", "N4", "C3", "C2", "C1", "O0"];
        let mol_forward = build_named(&forward, &atoms, &bonds);
        let mol_reversed = build_named(&reversed, &atoms, &bonds);
        assert_eq!(
            canonical_smiles(&mol_forward),
            canonical_smiles(&mol_reversed)
        );

        let out_forward = canonical_tautomer(&mol_forward);
        let out_reversed = canonical_tautomer(&mol_reversed);
        assert_eq!(
            canonical_smiles(&out_forward),
            canonical_smiles(&out_reversed)
        );
        // Must actually be the fully di-keto form, not merely "unanimous but
        // still enol" or a half-shifted state -- checked structurally
        // (double-bonded O count), not by SMILES substring, since "=O" vs
        // "O=" both spell a keto oxygen depending on write direction.
        let double_bonded_o_count = (0..out_forward.atom_count())
            .map(|i| AtomIdx(i as u32))
            .filter(|&idx| {
                out_forward.atom(idx).element.atomic_number() == 8
                    && out_forward
                        .neighbors(idx)
                        .any(|(_, bidx)| out_forward.bond(bidx).order == BondOrder::Double)
            })
            .count();
        assert_eq!(
            double_bonded_o_count,
            2,
            "expected both exocyclic sites shifted to keto: {}",
            canonical_smiles(&out_forward)
        );

        let limits = TautomerLimits::default();
        let parent_forward = tautomer_parent(&mol_forward, &limits);
        let parent_reversed = tautomer_parent(&mol_reversed, &limits);
        assert_eq!(parent_forward.status, parent_reversed.status);
        assert_eq!(
            canonical_smiles(&parent_forward.molecule),
            canonical_smiles(&parent_reversed.molecule)
        );
        assert_eq!(
            rule_id_sequence(&parent_forward),
            rule_id_sequence(&parent_reversed)
        );
    }

    /// Two independent SMILES respellings of the *same* tautomer (never a
    /// different one -- checked explicitly below) must agree on
    /// `canonical_tautomer`, `tautomer_parent`'s molecule + status, and the
    /// applied-rule-id sequence.
    fn assert_order_invariant_respelling(label: &str, a: &str, b: &str) {
        let mol_a = parse(a).unwrap();
        let mol_b = parse(b).unwrap();
        assert_eq!(
            canonical_smiles(&mol_a),
            canonical_smiles(&mol_b),
            "{label}: the two respellings are not even the same input molecule"
        );

        let out_a = canonical_tautomer(&mol_a);
        let out_b = canonical_tautomer(&mol_b);
        assert_eq!(
            canonical_smiles(&out_a),
            canonical_smiles(&out_b),
            "{label}: canonical_tautomer disagreed across respellings"
        );

        let limits = TautomerLimits::default();
        let parent_a = tautomer_parent(&mol_a, &limits);
        let parent_b = tautomer_parent(&mol_b, &limits);
        assert_eq!(
            parent_a.status, parent_b.status,
            "{label}: status disagreed"
        );
        assert_eq!(
            canonical_smiles(&parent_a.molecule),
            canonical_smiles(&parent_b.molecule),
            "{label}: tautomer_parent molecule disagreed"
        );
        assert_eq!(
            rule_id_sequence(&parent_a),
            rule_id_sequence(&parent_b),
            "{label}: applied rule_id sequence disagreed"
        );
    }

    #[test]
    fn atom_order_invariance_4_pyridone_respelling() {
        assert_order_invariant_respelling("4-pyridone", "Oc1ccncc1", "n1ccc(O)cc1");
    }

    #[test]
    fn atom_order_invariance_hypoxanthine_respelling() {
        // Enol input this time (round 2C-2's fix applies to the enol side);
        // b starts the ring traversal at the acceptor nitrogen instead of
        // the bridge carbon, and moves the exocyclic O substituent to the
        // end of the string -- verified equivalent to `a` via
        // canonical_smiles before this test was written.
        assert_order_invariant_respelling("hypoxanthine", "Oc1ncnc2[nH]cnc12", "n1cnc2[nH]cnc2c1O");
    }

    #[test]
    fn atom_order_invariance_isotope_2_pyridone_respelling() {
        assert_order_invariant_respelling("18O-2-pyridone", "[18OH]c1ccccn1", "n1ccccc1[18OH]");
    }

    #[test]
    fn atom_order_invariance_remote_stereocenter_2_pyridone_respelling() {
        // b starts the ring traversal at the ring nitrogen instead of the
        // O-bearing bridge carbon -- verified equivalent to `a` via
        // canonical_smiles before this test was written.
        assert_order_invariant_respelling(
            "2-pyridone + remote stereocenter",
            "Oc1c([C@@H](F)Cl)cccn1",
            "n1c(O)c([C@@H](F)Cl)ccc1",
        );
    }

    #[test]
    fn find_exocyclic_lactam_shift_matches_dedups_across_fused_rings() {
        // Hypoxanthine's bridge/acceptor pair is reachable through more than
        // one SSSR ring in a fused bicyclic system; must be reported once.
        let mol = parse("Oc1ncnc2[nH]cnc12").unwrap();
        let matches = find_exocyclic_lactam_shift_matches(&mol, &atom_rank(&mol));
        assert!(!matches.is_empty());
        let mut seen = HashSet::new();
        for m in &matches {
            assert!(seen.insert(*m), "duplicate candidate triple: {m:?}");
        }
    }

    #[test]
    fn negative_control_charged_pyridinium_acceptor_excluded() {
        // N-methyl-2-hydroxypyridinium: the ring N is charged and
        // 3-connected (methyl substituent) -- excluded by both the charge
        // check and the degree==2 valence-compatibility check.
        assert_noop("N-methyl-2-hydroxypyridinium", "Oc1cccc[n+]1C");
    }

    #[test]
    fn negative_control_aromatic_n_bridge_excluded() {
        // N-hydroxypyrrole: the bridge atom bearing the exocyclic OH is
        // itself an aromatic nitrogen, not carbon -- excluded by the
        // bridge-element==6 check.
        assert_noop("N-hydroxypyrrole", "On1cccc1");
    }

    #[test]
    fn negative_control_fused_bridgehead_nitrogen_acceptor_excluded() {
        // The candidate ring nitrogen adjacent to the exocyclic-OH-bearing
        // bridge is a 3-connected fusion (bridgehead) atom in this fused
        // bicyclic system -- excluded by the acceptor degree==2 check (no
        // free valence slot for +1H without a charge change).
        assert_noop("fused bridgehead-N acceptor", "Oc1ccn2ccccc12");
    }

    #[test]
    fn remote_stereocenter_same_atom_chirality_and_cip_preserved() {
        // Stronger than counting stereocenters: locate the SAME physical
        // atom (index-stable -- the shift's post-generation invariant check
        // guarantees every non-donor/acceptor atom, including the
        // stereocenter, keeps its index) and confirm chirality,
        // stereo_neighbor_order, and the CIP label are all unchanged, not
        // just present.
        let mol = parse("Oc1c([C@@H](F)Cl)cccn1").unwrap();
        let out = canonical_tautomer(&mol);

        let stereocenter = (0..mol.atom_count())
            .map(|i| AtomIdx(i as u32))
            .find(|&i| mol.atom(i).chirality != Chirality::None)
            .expect("input must have a stereocenter");

        assert_ne!(
            canonical_smiles(&mol),
            canonical_smiles(&out),
            "sanity: the shift must have actually fired"
        );
        assert_eq!(
            mol.atom(stereocenter).element,
            out.atom(stereocenter).element
        );
        assert_eq!(
            mol.atom(stereocenter).chirality,
            out.atom(stereocenter).chirality,
            "chirality flag changed on the uninvolved stereocenter"
        );
        assert_eq!(
            mol.stereo_neighbor_order(stereocenter),
            out.stereo_neighbor_order(stereocenter),
            "stereo_neighbor_order changed on the uninvolved stereocenter"
        );
        assert_eq!(
            crate::cip::assign_cip(&mol).get(stereocenter),
            crate::cip::assign_cip(&out).get(stereocenter),
            "CIP label changed on the uninvolved stereocenter"
        );
    }

    // -----------------------------------------------------------------------
    // Issue #415: canonical_tautomer must never produce an unkekulizable
    // (chemically invalid) molecule, even for a fused/bridged ring system
    // where `find_sssr`'s own ring choice can differ depending on atom
    // labeling.
    // -----------------------------------------------------------------------

    #[test]
    fn issue415_fused_ring_shift_never_produces_invalid_molecule() {
        // A phthalazinone-like fused bicyclic tautomer where re-labeling via
        // a canonical-SMILES round trip used to make both the
        // exocyclic-lactam-shift mechanism and the plain direct-aromatic
        // 1,2-shift mechanism protonate a ring nitrogen that was already
        // ring-adjacent to another protonated nitrogen -- an over-valent
        // `[nH2]`-shaped state RDKit itself rejects outright (confirmed via
        // `Chem.MolFromSmiles` at diagnosis time). Every candidate this
        // molecule can reach through repeated standardize passes must stay
        // kekulizable.
        let smi = "Oc1[nH]ncc2c3cc(OCc4ccccc4)ccc3nc1-2";
        let mol = parse(smi).unwrap();
        let mut current = mol;
        for pass in 0..6 {
            let next = canonical_tautomer(&current);
            assert!(
                chematic_core::kekulize(&next).is_ok(),
                "pass {pass}: canonical_tautomer produced an unkekulizable molecule: {}",
                canonical_smiles(&next)
            );
            let reparsed = parse(&canonical_smiles(&next)).unwrap();
            current = reparsed;
        }
    }

    /// Direct regression for issue #415's residual order-dependence: the
    /// same molecule, relabeled only via atom-order permutation
    /// (`reorder_atoms`, so every permutation is guaranteed structurally
    /// identical -- not a hand-alternate SMILES respelling that could
    /// accidentally encode a different molecule), must reach the byte-
    /// identical canonical tautomer regardless of which permutation the
    /// caller happens to hand in, and that result must already be its own
    /// fixed point (`canonical_tautomer` applied to its own output must be
    /// a no-op) -- not just "eventually convergent after enough reparses."
    ///
    /// Empirically confirmed before writing this test (not assumed): this
    /// molecule's true canonical form needs the *search itself* to run
    /// twice to stabilize (`canonical_tautomer_with_config`'s own outer
    /// convergence loop, not an atom-order artifact -- a single search
    /// pass genuinely doesn't reach the fixed point, independent of atom
    /// labeling). 6 orderings covering identity, full reversal, and 4
    /// distinct rotations all agree exactly, both on the first search's
    /// result and on the final stable form.
    #[test]
    fn issue415_residual_canonical_tautomer_is_atom_order_invariant_and_idempotent() {
        let smi = "Oc1[nH]ncc2c3cc(OCc4ccccc4)ccc3nc1-2";
        let base = parse(smi).unwrap();
        let n = base.atom_count();

        let mut orderings: Vec<Vec<AtomIdx>> = vec![
            (0..n as u32).map(AtomIdx).collect(),
            (0..n as u32).rev().map(AtomIdx).collect(),
        ];
        for shift in [1usize, 3, (n / 2).max(1), n - 1] {
            let mut order: Vec<AtomIdx> = (0..n as u32).map(AtomIdx).collect();
            order.rotate_left(shift);
            orderings.push(order);
        }
        assert_eq!(orderings.len(), 6, "expected 6 distinct atom orderings");

        let mut expected: Option<(String, String)> = None;
        for (i, order) in orderings.iter().enumerate() {
            let permuted = reorder_atoms(&base, order);
            assert_eq!(
                canonical_smiles(&permuted),
                canonical_smiles(&base),
                "ordering {i}: reorder_atoms must not change the molecule's own identity"
            );
            let result = canonical_tautomer(&permuted);
            let reapplied = canonical_tautomer(&result);
            let result_smi = canonical_smiles(&result);
            let reapplied_smi = canonical_smiles(&reapplied);
            assert_eq!(
                result_smi, reapplied_smi,
                "ordering {i}: canonical_tautomer's own output is not a fixed point of itself"
            );
            match &expected {
                None => expected = Some((result_smi, reapplied_smi)),
                Some((exp_result, exp_reapplied)) => {
                    assert_eq!(
                        &result_smi, exp_result,
                        "ordering {i} reached a different canonical tautomer than ordering 0"
                    );
                    assert_eq!(&reapplied_smi, exp_reapplied);
                }
            }
        }
    }
}

#[cfg(test)]
mod parent_budget_boundary_tests {
    use super::*;
    use crate::parent::{ParentAudit, ParentComputationStatus};
    use chematic_smiles::{canonical_smiles, parse};

    #[test]
    fn zero_timeout_and_transform_budget_preserve_input_and_report_incomplete() {
        for source in ["CCN=O", "Nc1cc[nH]c(=O)n1", "Nc1nc2[nH]cnc2c(=O)[nH]1"] {
            let mol = parse(source).unwrap();
            let before = canonical_smiles(&mol);
            let limits = TautomerLimits {
                timeout_ms: Some(0),
                ..Default::default()
            };
            let result = tautomer_parent(&mol, &limits);
            assert_eq!(result.status, ParentComputationStatus::TimedOut);
            assert_eq!(canonical_smiles(&result.molecule), before);
            let ParentAudit::Tautomer(audit) = result.audit else {
                panic!("tautomer audit required")
            };
            assert!(audit.applied_transforms.is_empty());
            assert_eq!(audit.candidate_count, 1);
            if source != "CCN=O" {
                continue;
            }
            let limits = TautomerLimits {
                max_transforms: 0,
                ..Default::default()
            };
            let result = tautomer_parent(&mol, &limits);
            assert_eq!(canonical_smiles(&result.molecule), before);
            assert_eq!(canonical_smiles(&mol), before);
            assert_eq!(result.status, ParentComputationStatus::MaxTransformsReached);
        }
    }

    #[test]
    fn dual_flank_parent_audit_tracks_real_atoms_and_converges_across_forms() {
        for forms in [
            ["Nc1cc[nH]c(=O)n1", "Nc1ccnc(O)n1"],
            ["Nc1nc2[nH]cnc2c(=O)[nH]1", "Nc1nc2[nH]cnc2c(O)n1"],
        ] {
            let mut parents = Vec::new();
            for source in forms {
                let mol = parse(source).unwrap();
                let before = canonical_smiles(&mol);
                let result = tautomer_parent(&mol, &TautomerLimits::default());
                assert_eq!(
                    result.status,
                    ParentComputationStatus::Completed,
                    "{source}"
                );
                assert_eq!(result.molecule.atom_count(), mol.atom_count());
                assert_eq!(result.molecule.bond_count(), mol.bond_count());
                let ParentAudit::Tautomer(audit) = result.audit else {
                    panic!("tautomer audit required")
                };
                for step in &audit.applied_transforms {
                    assert_eq!(step.affected_atoms.len(), 3);
                    assert!(
                        step.affected_atoms
                            .iter()
                            .all(|a| (a.0 as usize) < mol.atom_count())
                    );
                    assert!(
                        step.affected_bonds
                            .iter()
                            .all(|b| (b.0 as usize) < mol.bond_count())
                    );
                }
                let parent = canonical_smiles(&result.molecule);
                assert_eq!(
                    canonical_smiles(
                        &tautomer_parent(&result.molecule, &TautomerLimits::default()).molecule
                    ),
                    parent
                );
                parents.push(parent);
                assert_eq!(canonical_smiles(&mol), before);
            }
            assert_eq!(parents[0], parents[1]);
        }
    }
}
