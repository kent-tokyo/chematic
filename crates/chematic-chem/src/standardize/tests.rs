use super::*;
use crate::hash::mol_hash;
use chematic_smiles::parse;

#[test]
fn cleanup_keeps_neutral_acid_and_amine_unchanged() {
    for text in ["CC(=O)O", "CN"] {
        let mol = parse(text).unwrap();
        let cleaned = cleanup(&mol);
        assert_eq!(mol_hash(&cleaned), mol_hash(&mol), "{text}");
    }
}

#[test]
fn cleanup_preserves_aromatic_n_oxide_charge_and_valence() {
    for input in ["c1cc[n+]([O-])cc1", "[O-][n+]1ccccc1"] {
        let cleaned = cleanup(&parse(input).unwrap());
        assert_eq!(
            chematic_smiles::rdkit_canonical_smiles(&cleaned).as_deref(),
            Ok("[O-][n+]1ccccc1"),
            "{input}"
        );
    }
}

#[test]
fn cleanup_neutralizes_ortho_and_para_pyridinium_phenoxide() {
    assert_eq!(
        chematic_perception::assign_aromaticity(&parse("C[N+]1=CC=C(C=C1)[O-]").unwrap())
            .aromatic_atom_count(),
        6
    );
    for (input, expected) in [
        ("[O-]c1cc[n+](C)cc1", "Cn1ccc(=O)cc1"),
        ("C[N+]1=CC=C(C=C1)[O-]", "Cn1ccc(=O)cc1"),
        ("[O-]c1cccc[n+]1[O-]", "O=c1ccccn1[O-]"),
        ("c1c([O-])cc[n+](c1)C", "Cn1ccc(=O)cc1"),
        ("[O-]c1cc[nH+]cc1", "O=c1cc[nH]cc1"),
        ("[O-]c1cccc[nH+]1", "O=c1cccc[nH]1"),
        // The meta charge-separated form is outside this normalization.
        ("[O-]c1ccc[n+](C)c1", "C[n+]1cccc([O-])c1"),
        ("[O-]c1cc[n+](CC)cc1", "CCn1ccc(=O)cc1"),
        ("[O-]c1cc[n+](C)c(F)c1", "Cn1ccc(=O)cc1F"),
        ("[O-]c1cc[n+](C)c(Cl)c1", "Cn1ccc(=O)cc1Cl"),
        ("[O-]c1cc[n+](C)cc1C", "Cc1cn(C)ccc1=O"),
        // Fused systems are intentionally outside the simple-ring rule.
        ("[O-]c1ccc2cc[n+](C)cc2c1", "C[n+]1ccc2ccc([O-])cc2c1"),
        ("[O-]c1cc[n+](C)cc1[O-]", "Cn1ccc(=O)c([O-])c1"),
        ("[O-]c1cc[n+](C)cc1.[Na+]", "Cn1ccc(=O)cc1.[Na+]"),
        (
            "[O-]c1cc[n+](C)cc1.[O-]c1cc[n+](C)cc1",
            "Cn1ccc(=O)cc1.Cn1ccc(=O)cc1",
        ),
    ] {
        let cleaned = cleanup(&parse(input).unwrap());
        assert_eq!(
            chematic_smiles::rdkit_cleanup_canonical_smiles(&cleaned).as_deref(),
            Ok(expected),
            "{input}"
        );
        assert_eq!(
            chematic_smiles::rdkit_cleanup_canonical_smiles(&cleanup(&cleaned)).as_deref(),
            Ok(expected),
            "second Cleanup: {input}"
        );
    }
}

#[test]
fn cleanup_writer_keeps_isolated_metal_coordination_tag() {
    for (input, expected) in [
        ("F[Pt@SP1](Cl)(Br)I", "[Br-].[Cl-].[F-].[I-].[Pt@SP1+4]"),
        ("F[Pt@SP2](Cl)(Br)I", "[Br-].[Cl-].[F-].[I-].[Pt@SP2+4]"),
        (
            "C[Co@OH1](N)(N)(N)(N)Cl",
            "[CH3-].[Cl-].[Co@OH1+6].[NH2-].[NH2-].[NH2-].[NH2-]",
        ),
        (
            "C[Fe@TB2](F)(Cl)(Br)I",
            "[Br-].[CH3-].[Cl-].[F-].[Fe@TB2+5].[I-]",
        ),
    ] {
        let cleaned = cleanup(&parse(input).unwrap());
        assert_eq!(
            chematic_smiles::rdkit_cleanup_canonical_smiles(&cleaned).as_deref(),
            Ok(expected),
            "{input}"
        );
    }
    // The vacant coordination site changes the square-planar permutation
    // when the three ligands are removed.
    let partial = cleanup(&parse("F[Pt@SP1](Cl)Br").unwrap());
    assert_eq!(
        chematic_smiles::rdkit_cleanup_canonical_smiles(&partial).as_deref(),
        Ok("[Br-].[Cl-].[F-].[Pt@SP2+3]")
    );
}

#[test]
fn cleanup_disconnects_covalent_metal_ligands_with_charge_balance() {
    let mol = parse("Cl[Au]").unwrap();
    let cleaned = cleanup(&mol);
    assert_eq!(cleaned.bond_count(), 0);
    let charges: Vec<_> = cleaned.atoms().map(|(_, atom)| atom.charge).collect();
    assert_eq!(charges, [-1, 1]);
    assert_eq!(
        cleaned
            .atoms()
            .map(|(_, atom)| i32::from(atom.charge))
            .sum::<i32>(),
        0
    );
}

#[test]
fn cleanup_respects_default_metal_disconnector_element_scope() {
    for (input, expected) in [
        ("[Na]Cl", "[Na][Cl]"),
        ("[Li]C", "[Li][CH3]"),
        ("C[Mg]Br", "[CH3][Mg][Br]"),
        ("C[Hg]Cl", "[CH3][Hg][Cl]"),
        ("[Na]O", "[Na+].[OH-]"),
        ("[Fe]Cl", "[Cl-].[Fe+]"),
        ("[Sn]C", "[CH3][Sn]"),
    ] {
        let cleaned = cleanup(&parse(input).unwrap());
        assert_eq!(
            chematic_smiles::rdkit_canonical_smiles(&cleaned).as_deref(),
            Ok(expected),
            "{input}"
        );
    }
}

#[test]
fn cleanup_charges_bare_isolated_metals_like_rdkit() {
    // RDKit 2026.09.1 Cleanup changes only the bare isolated atoms in
    // this set. Explicit hydrides and bonded metals retain their form.
    for (input, expected) in [
        ("[Li]", "[Li+]"),
        ("[Na]", "[Na+]"),
        ("[K]", "[K+]"),
        ("[Mg]", "[Mg+2]"),
        ("[Ca]", "[Ca+2]"),
        ("[Rb]", "[Rb]"),
        ("[Cs]", "[Cs]"),
        ("[13Na:7]", "[13Na+:7]"),
        ("[LiH]", "[LiH]"),
        ("[MgH2]", "[MgH2]"),
        ("[H][Li]", "[LiH]"),
        ("[H][Mg][H]", "[MgH2]"),
        ("[Na]C", "[CH3][Na]"),
        ("[Na].c1ccccc1", "[Na+].c1ccccc1"),
        ("[Na].F[C@](Cl)(Br)I", "F[C@](Cl)(Br)I.[Na+]"),
    ] {
        let cleaned = cleanup(&parse(input).unwrap());
        assert_eq!(
            chematic_smiles::rdkit_cleanup_canonical_smiles(&cleaned).as_deref(),
            Ok(expected),
            "{input}"
        );
        assert_eq!(
            chematic_smiles::rdkit_cleanup_canonical_smiles(&cleanup(&cleaned)).as_deref(),
            Ok(expected),
            "second Cleanup: {input}"
        );
    }
}

#[test]
fn largest_fragment_two_fragments_picks_larger() {
    // "CC.CCC" — ethane (2 C) and propane (3 C)
    let mol = parse("CC.CCC").unwrap();
    let result = largest_fragment(&mol);
    assert_eq!(result.atom_count(), 3, "should keep propane (3 C)");
}

#[test]
fn largest_fragment_single_fragment_unchanged() {
    // "CC" — ethane, only one fragment
    let mol = parse("CC").unwrap();
    let result = largest_fragment(&mol);
    assert_eq!(result.atom_count(), 2);
}

#[test]
fn largest_fragment_keeps_benzene_over_ethane() {
    // "CC.c1ccccc1" — ethane (2 C) vs benzene (6 C)
    let mol = parse("CC.c1ccccc1").unwrap();
    let result = largest_fragment(&mol);
    assert_eq!(result.atom_count(), 6, "should keep benzene (6 atoms)");
}

#[test]
fn largest_fragment_ionic_pair_keeps_one_atom() {
    // "[Na+].[Cl-]" — both fragments are single atoms; either is fine
    let mol = parse("[Na+].[Cl-]").unwrap();
    let result = largest_fragment(&mol);
    assert_eq!(result.atom_count(), 1);
}

#[test]
fn neutralize_neutral_molecule_unchanged() {
    // "CC" is already neutral; no atom should gain/lose charge
    let mol = parse("CC").unwrap();
    let result = neutralize_charges(&mol);
    for i in 0..result.atom_count() {
        let atom = result.atom(AtomIdx(i as u32));
        assert_eq!(atom.charge, 0, "all atoms should remain neutral");
    }
}

#[test]
fn neutralize_acetate_oxygen() {
    // "CC(=O)[O-]" — acetate; the [O-] should become neutral with H added
    let mol = parse("CC(=O)[O-]").unwrap();
    let result = neutralize_charges(&mol);

    // Find the oxygen that was originally [O-]: it should now have charge 0
    // and hydrogen_count == Some(1).
    let neutralized_o = (0..result.atom_count())
        .map(|i| result.atom(AtomIdx(i as u32)))
        .find(|a| a.element == Element::O && a.hydrogen_count == Some(1));

    assert!(
        neutralized_o.is_some(),
        "neutralized [O-] should have hydrogen_count == Some(1)"
    );
    assert_eq!(
        neutralized_o.unwrap().charge,
        0,
        "neutralized [O-] should have charge == 0"
    );
}

#[test]
fn standardize_with_defaults() {
    // "CC(=O)[O-]" — acetate ion
    let mol = parse("CC(=O)[O-]").unwrap();
    let opts = StandardizeOptions::default();
    let result = standardize(&mol, &opts);

    // With default options, remove_explicit_h will be applied.
    // Check that [O-] was neutralized (charge should be 0).
    let has_neutral_o = (0..result.atom_count())
        .map(|i| result.atom(AtomIdx(i as u32)))
        .any(|a| a.element == Element::O && a.charge == 0);
    assert!(has_neutral_o, "acetate oxygen should be neutralized");

    // Should have at least 3 atoms (C, C, O) with no explicit H
    assert!(
        result.atom_count() >= 3,
        "should have at least 3 atoms after standardization"
    );
}

#[test]
fn standardize_skip_largest_fragment() {
    // "CC.CCC" — ethane and propane
    let mol = parse("CC.CCC").unwrap();
    let opts = StandardizeOptions {
        zwitterion_handling: ZwitterionHandling::Normalize,
        largest_fragment_only: false,
        ..Default::default()
    };
    let result = standardize(&mol, &opts);

    // Should keep both fragments
    assert_eq!(
        result.atom_count(),
        5,
        "should keep both fragments when largest_fragment_only=false"
    );
}

#[test]
fn pipeline_report_tracks_enabled_stage_changes() {
    let mol = parse("CC.CCC").unwrap();
    let pipeline = StandardizationPipeline::new(StandardizeOptions {
        largest_fragment_only: true,
        neutralize_charges: false,
        remove_explicit_h: false,
        canonical_tautomer: false,
        zwitterion_handling: ZwitterionHandling::Keep,
    });

    let (result, report) = pipeline.run(&mol);

    assert_eq!(result.atom_count(), 3);
    assert_eq!(report.status, PipelineStatus::Modified);
    assert!(report.changed());
    assert_eq!(report.steps.len(), 6);
    // Metal disconnection is the first (not enabled, so no change)
    assert_eq!(report.steps[0].step, StandardizationStep::DisconnectMetals);
    assert!(!report.steps[0].enabled);
    // NeutralizeCharges is applied second (not enabled, so no change)
    assert_eq!(report.steps[1].step, StandardizationStep::NeutralizeCharges);
    assert!(!report.steps[1].enabled);
    // LargestFragment is applied third and is enabled
    assert_eq!(report.steps[2].step, StandardizationStep::LargestFragment);
    assert!(report.steps[2].enabled);
    assert!(report.steps[2].changed);
}

#[test]
fn pipeline_report_marks_unchanged_clean_molecule() {
    let mol = parse("CC").unwrap();
    let pipeline = StandardizationPipeline::new(StandardizeOptions {
        canonical_tautomer: false,
        neutralize_charges: false,
        remove_explicit_h: false,
        largest_fragment_only: false,
        zwitterion_handling: ZwitterionHandling::Keep,
    });

    let (_result, report) = pipeline.run(&mol);

    assert_eq!(report.status, PipelineStatus::Unchanged);
    assert!(!report.changed());
    assert!(report.warnings.is_empty());
    assert!(report.steps.iter().all(|s| !s.enabled && !s.changed));
}

#[test]
fn pipeline_report_disconnects_metal_bonds() {
    let mol = parse("[Na]OC").unwrap();
    assert_eq!(mol.bond_count(), 2, "input has Na-O and O-C bonds");

    let pipeline = StandardizationPipeline::new(StandardizeOptions {
        canonical_tautomer: false,
        neutralize_charges: false,
        remove_explicit_h: false,
        largest_fragment_only: false,
        zwitterion_handling: ZwitterionHandling::Keep,
    });

    let (result, report) = pipeline.run(&mol);

    // Metal disconnection should run automatically, removing the Na-O bond
    assert_eq!(result.bond_count(), 1, "Na-O bond should be disconnected");
    // The remaining bond should be O-C
    assert!(
        result.bond(BondIdx(0)).atom1.0 < 3 && result.bond(BondIdx(0)).atom2.0 < 3,
        "remaining bond should connect organic atoms"
    );
    let metal_step = &report.steps[0];
    assert_eq!(metal_step.step, StandardizationStep::DisconnectMetals);
    assert!(metal_step.enabled);
    assert!(metal_step.changed);
    assert_eq!(metal_step.before.bonds, 2);
    assert_eq!(metal_step.after.bonds, 1);
}

#[test]
fn bug3_ionic_pair_neutralize_before_largest_fragment() {
    // BUG #3 Fix Verification: NeutralizeCharges must run BEFORE LargestFragment
    // Example: [NH4+].[OH-] (ammonium hydroxide)
    // NH4+ will be neutralized (reduce H by 1)
    // After neutralization: [NH3].[OH-] - still two fragments
    // LargestFragment then picks the larger one
    // Key: stage order affects which fragment is selected
    let mol = parse("[NH4+].[OH-]").unwrap();
    let pipeline = StandardizationPipeline::new(StandardizeOptions {
        largest_fragment_only: true,
        neutralize_charges: true,
        remove_explicit_h: false,
        canonical_tautomer: false,
        zwitterion_handling: ZwitterionHandling::Normalize,
    });

    let (_result, report) = pipeline.run(&mol);

    // Verify step order: NeutralizeCharges MUST come before LargestFragment
    assert_eq!(report.steps.len(), 6, "Should have 6 steps in pipeline");
    assert_eq!(
        report.steps[0].step,
        StandardizationStep::DisconnectMetals,
        "DisconnectMetals must be step 0"
    );
    assert_eq!(
        report.steps[1].step,
        StandardizationStep::NeutralizeCharges,
        "NeutralizeCharges must be step 1"
    );
    assert_eq!(
        report.steps[2].step,
        StandardizationStep::LargestFragment,
        "LargestFragment must be step 2"
    );

    // The test passes if step order is correct and the pipeline runs without error
    assert!(report.changed(), "Pipeline should report changes");
    assert_eq!(
        report.status,
        PipelineStatus::Modified,
        "Should be marked as Modified"
    );
    assert!(
        report.warnings.is_empty(),
        "valid ammonium hydroxide should not warn"
    );
}

// ── Parent structure extraction tests ────────────────────────────────────

#[test]
fn remove_isotopes_strips_isotope_labels() {
    // "[13C]CC" has one 13C isotope label
    let mol = parse("[13C]CC").unwrap();
    let result = remove_isotopes(&mol);
    for i in 0..result.atom_count() {
        assert_eq!(
            result.atom(chematic_core::AtomIdx(i as u32)).isotope,
            None,
            "atom {} should have no isotope",
            i
        );
    }
}

#[test]
fn remove_isotopes_preserves_structure() {
    // "[13C]CC" and "CC" should be structurally identical after isotope removal
    let mol = parse("[13C]CC").unwrap();
    let result = remove_isotopes(&mol);
    assert_eq!(result.atom_count(), 3, "atom count preserved");
    assert_eq!(result.bond_count(), 2, "bond count preserved");
}

#[test]
fn remove_stereo_strips_chirality() {
    // "N[C@@H](C)C(=O)O" — alanine with (S) stereochemistry
    let mol = parse("N[C@@H](C)C(=O)O").unwrap();
    let result = remove_stereo(&mol);
    for i in 0..result.atom_count() {
        use chematic_core::Chirality;
        assert_eq!(
            result.atom(chematic_core::AtomIdx(i as u32)).chirality,
            Chirality::None,
            "atom {} should have no chirality",
            i
        );
    }
}

#[test]
fn remove_stereo_converts_wedge_bonds_to_single() {
    // "C[C@H](O)C" — methylcarbinol with wedge stereochemistry
    // After remove_stereo, Up/Down bonds should become Single
    let mol = parse("C[C@H](O)C").unwrap();
    let result = remove_stereo(&mol);
    for i in 0..result.bond_count() {
        use chematic_core::BondOrder;
        let bond = result.bond(chematic_core::BondIdx(i as u32));
        assert_ne!(
            bond.order,
            BondOrder::Up,
            "bond {} should not be Up after stereo removal",
            i
        );
        assert_ne!(
            bond.order,
            BondOrder::Down,
            "bond {} should not be Down after stereo removal",
            i
        );
    }
}

#[test]
fn remove_stereo_preserves_structure() {
    // "N[C@@H](C)C(=O)O" (alanine) should keep 6 atoms, 5 bonds after stereo removal
    let mol = parse("N[C@@H](C)C(=O)O").unwrap();
    let result = remove_stereo(&mol);
    assert_eq!(result.atom_count(), 6, "atom count preserved");
    assert_eq!(result.bond_count(), 5, "bond count preserved");
}

#[test]
fn parent_variant_step_names_distinct() {
    // Verify that the 4 new parent variants have distinct step names
    let frag_parent = StandardizationStep::FragmentParent;
    let charge_parent = StandardizationStep::ChargeParent;
    let isotope_parent = StandardizationStep::IsotopeParent;
    let stereo_parent = StandardizationStep::StereoParent;

    assert_eq!(frag_parent.as_str(), "fragment_parent");
    assert_eq!(charge_parent.as_str(), "charge_parent");
    assert_eq!(isotope_parent.as_str(), "isotope_parent");
    assert_eq!(stereo_parent.as_str(), "stereo_parent");
}

// ── Cleanup transform tests (B3) ────────────────────────────────────────

#[test]
fn prefer_organic_removes_inorganic_salts() {
    // "CCO.[Na+].[Cl-]" — ethanol + sodium chloride
    let mol = parse("CCO.[Na+].[Cl-]").unwrap();
    assert_eq!(mol.atom_count(), 5, "input has CCO + Na + Cl");

    let result = prefer_organic(&mol);

    // Should keep only the organic ethanol fragment
    assert_eq!(result.atom_count(), 3, "should keep only ethanol (C, C, O)");
}

#[test]
fn prefer_organic_keeps_organic_if_no_inorganic() {
    // "CC" — ethane only
    let mol = parse("CC").unwrap();
    let result = prefer_organic(&mol);
    assert_eq!(result.atom_count(), 2, "ethane unchanged");
}

#[test]
fn prefer_organic_falls_back_to_largest() {
    // "C.C.C" — three separate carbons (all organic)
    let mol = parse("C.C.C").unwrap();
    let result = prefer_organic(&mol);
    // Should keep the largest organic fragment (any one of them, but all are size 1)
    assert_eq!(
        result.atom_count(),
        1,
        "falls back to largest fragment (one C)"
    );
}

#[test]
fn uncharge_neutralizes_all_charges() {
    // "[NH4+].[OH-]" — ammonium hydroxide
    let mol = parse("[NH4+].[OH-]").unwrap();
    assert!(
        mol.atoms().any(|(_, a)| a.charge != 0),
        "input has charged atoms"
    );

    let result = uncharge(&mol);

    // All atoms should be neutral
    for (_, atom) in result.atoms() {
        assert_eq!(atom.charge, 0, "all atoms should be neutral");
    }
}

#[test]
fn uncharge_follows_rdkit_uncharger() {
    let canon = |s: &str| chematic_smiles::canonical_smiles(&parse(s).unwrap());
    let run = |s: &str| chematic_smiles::canonical_smiles(&uncharge(&parse(s).unwrap()));
    for (input, expected) in [
        ("[NH4+].[OH-]", "N.O"),
        ("C[NH3+]", "CN"),
        ("CC(=O)[O-]", "CC(=O)O"),
        // Charge-separated groups stay (were radicals `[N](=O)[O]`).
        ("C[N+](=O)[O-]", "C[N+](=O)[O-]"),
        ("C[N+](C)(C)[O-]", "C[N+](C)(C)[O-]"),
        // A cation without H keeps its charge and one balancing anion.
        ("C[N+](C)(C)CC(=O)[O-]", "C[N+](C)(C)CC(=O)[O-]"),
        ("[Na+].CC(=O)[O-]", "[Na+].CC(=O)[O-]"),
        ("F[B-](F)(F)F", "F[B-](F)(F)F"),
    ] {
        assert_eq!(run(input), canon(expected), "{input}");
    }
}

#[test]
fn reionize_leaves_ethers_and_esters_alone() {
    // An ether/ester O was given a negative charge with two bonds,
    // depending on which C neighbour came first.
    for input in ["COc1ccccc1", "CC(=O)OC", "c1ccccc1OC(=O)C"] {
        let out = reionize(&parse(input).unwrap());
        assert!(out.atoms().all(|(_, a)| a.charge == 0), "{input}");
    }
    let canon = |s: &str| chematic_smiles::canonical_smiles(&parse(s).unwrap());
    let run = |s: &str| chematic_smiles::canonical_smiles(&reionize(&parse(s).unwrap()));
    assert_eq!(run("OC(=O)CCN"), canon("OC(=O)CCN"));
    // Neutral acid/base groups do not ionize independently.
    assert_eq!(run("Nc1ccccc1"), canon("Nc1ccccc1"));
    assert_eq!(run("CS(=O)(=O)NC"), canon("CS(=O)(=O)NC"));
}

#[test]
fn reionize_keeps_neutral_carboxylic_acids() {
    // "CC(=O)O" — acetic acid
    let mol = parse("CC(=O)O").unwrap();

    let result = reionize(&mol);

    // A neutral acid does not create an anion without an existing base.
    let has_negative_oxygen = result
        .atoms()
        .any(|(_, a)| a.element.atomic_number() == 8 && a.charge < 0);

    assert!(
        !has_negative_oxygen,
        "reionize must not independently deprotonate carboxylic acids"
    );
}

#[test]
fn reionize_keeps_neutral_amines() {
    // "CC(N)C" — secondary amine
    let mol = parse("CC(N)C").unwrap();

    let result = reionize(&mol);

    // Should have a positively charged nitrogen (ammonium)
    let has_positive_nitrogen = result
        .atoms()
        .any(|(_, a)| a.element.atomic_number() == 7 && a.charge > 0);

    assert!(
        !has_positive_nitrogen,
        "reionize must not independently protonate amines"
    );
}

#[test]
fn reionize_protects_amide_nitrogen() {
    // BUG FIX #2: Amide nitrogen should NOT be protonated
    // "CC(=O)N" — primary amide
    let mol = parse("CC(=O)N").unwrap();

    let result = reionize(&mol);

    // Should NOT have a positively charged nitrogen
    let has_positive_nitrogen = result
        .atoms()
        .any(|(_, a)| a.element.atomic_number() == 7 && a.charge > 0);

    assert!(
        !has_positive_nitrogen,
        "reionize should NOT protonate amide nitrogen"
    );
}

#[test]
fn reionize_protects_thioamide_nitrogen() {
    // BUG FIX #2: Thioamide nitrogen should also NOT be protonated
    // "CC(=S)N" — thioamide
    let mol = parse("CC(=S)N").unwrap();

    let result = reionize(&mol);

    // Should NOT have a positively charged nitrogen
    let has_positive_nitrogen = result
        .atoms()
        .any(|(_, a)| a.element.atomic_number() == 7 && a.charge > 0);

    assert!(
        !has_positive_nitrogen,
        "reionize should NOT protonate thioamide nitrogen (C=S conjugation)"
    );
}

// B2: normalize_groups expansion tests

#[test]
fn normalize_groups_nitro() {
    // Standard nitro group: [N+](=O)[O-]
    let mol = parse("C[N+](=O)[O-]").unwrap();
    let result = normalize_groups(&mol);

    // All atoms should be neutral after normalization
    let all_neutral = result.atoms().all(|(_, a)| a.charge == 0);
    assert!(all_neutral, "nitro group should be neutralized");

    // N-O bonds should be double
    let mut has_double_bond = false;
    for (_, bond) in result.bonds() {
        let a1 = result.atom(bond.atom1);
        let a2 = result.atom(bond.atom2);
        if ((a1.element.atomic_number() == 7 && a2.element.atomic_number() == 8)
            || (a1.element.atomic_number() == 8 && a2.element.atomic_number() == 7))
            && bond.order == chematic_core::BondOrder::Double
        {
            has_double_bond = true;
        }
    }
    assert!(has_double_bond, "nitro should have N=O double bond");
}

#[test]
fn normalize_groups_azide() {
    // Azide: [N-][N+]#N
    let mol = parse("[N-][N+]#N").unwrap();
    let result = normalize_groups(&mol);

    // All atoms should be neutral
    let all_neutral = result.atoms().all(|(_, a)| a.charge == 0);
    assert!(all_neutral, "azide should be neutralized");

    // Check for N=N bonds (converted from single)
    let mut has_double_bond_count = 0;
    for (_, bond) in result.bonds() {
        let a1 = result.atom(bond.atom1);
        let a2 = result.atom(bond.atom2);
        if a1.element.atomic_number() == 7
            && a2.element.atomic_number() == 7
            && bond.order == chematic_core::BondOrder::Double
        {
            has_double_bond_count += 1;
        }
    }
    assert!(
        has_double_bond_count > 0,
        "azide should have N=N double bonds after normalization"
    );
}

#[test]
fn normalize_groups_sulfoxide() {
    // Sulfoxide: S(=O)(C)(C)
    let mol = parse("C[S](=O)C").unwrap();
    let result = normalize_groups(&mol);

    let mut has_charge_separated_s_o = false;
    for (_, bond) in result.bonds() {
        let a1 = result.atom(bond.atom1);
        let a2 = result.atom(bond.atom2);
        has_charge_separated_s_o |= a1.element.atomic_number() == 16
            && a1.charge == 1
            && a2.element.atomic_number() == 8
            && a2.charge == -1
            && bond.order == chematic_core::BondOrder::Single
            || a2.element.atomic_number() == 16
                && a2.charge == 1
                && a1.element.atomic_number() == 8
                && a1.charge == -1
                && bond.order == chematic_core::BondOrder::Single;
    }
    assert!(
        has_charge_separated_s_o,
        "sulfoxide should normalize to [S+]-[O-]"
    );

    // Sulfones are four-coordinate and must not be rewritten by the
    // three-coordinate sulfoxide rule.
    let sulfone = parse("CS(=O)(=O)C").unwrap();
    assert_eq!(
        chematic_smiles::canonical_smiles(&normalize_groups(&sulfone)),
        chematic_smiles::canonical_smiles(&sulfone)
    );
}

// ── RDKit PR #9051: StereoGroup cleanup for non-chiral atoms ────────────

#[test]
fn remove_stereo_clears_stereo_groups() {
    // remove_stereo must produce a molecule with no stereo groups (RDKit PR #9051).
    // It rebuilds via MoleculeBuilder::new() so groups are already empty; this
    // test guards against regressions where from_molecule() is accidentally used.
    use chematic_core::{AtomIdx, StereoGroup, StereoGroupKind};
    let mut mol = parse("[C@@H](F)(Cl)Br").unwrap();
    mol.add_stereo_group(StereoGroup::new(
        StereoGroupKind::Absolute,
        vec![AtomIdx(0)],
    ));
    assert_eq!(
        mol.stereo_groups().len(),
        1,
        "precondition: group was added"
    );
    let stripped = remove_stereo(&mol);
    assert_eq!(
        stripped.stereo_groups().len(),
        0,
        "remove_stereo must clear stereo groups"
    );
}

#[test]
fn clean_stereo_groups_drops_non_chiral_atoms() {
    // clean_stereo_groups filters out atoms with Chirality::None (RDKit PR #9051).
    use chematic_core::{AtomIdx, StereoGroup, StereoGroupKind};
    // atom 0 = CH3 (not chiral), atom 1 = @@ chiral center
    let mut mol = parse("C[C@@H](F)Cl").unwrap();
    mol.add_stereo_group(StereoGroup::new(
        StereoGroupKind::Absolute,
        vec![AtomIdx(0), AtomIdx(1)], // atom 0 is NOT chiral
    ));
    let cleaned = clean_stereo_groups(&mol);
    assert_eq!(
        cleaned.stereo_groups().len(),
        1,
        "group must survive with 1 atom"
    );
    assert_eq!(
        cleaned.stereo_groups()[0].atom_indices,
        vec![AtomIdx(1)],
        "only chiral atom must remain in group"
    );
}

#[test]
fn clean_stereo_groups_drops_empty_groups() {
    // A group with no chiral atoms at all must be removed entirely.
    use chematic_core::{AtomIdx, StereoGroup, StereoGroupKind};
    let mut mol = parse("CC").unwrap(); // no chiral atoms
    mol.add_stereo_group(StereoGroup::new(
        StereoGroupKind::Absolute,
        vec![AtomIdx(0)],
    ));
    let cleaned = clean_stereo_groups(&mol);
    assert_eq!(
        cleaned.stereo_groups().len(),
        0,
        "empty group must be removed"
    );
}

#[test]
fn clean_stereo_groups_preserves_valid_groups() {
    // A group whose atoms are all chiral must be kept intact.
    use chematic_core::{AtomIdx, StereoGroup, StereoGroupKind};
    let mut mol = parse("[C@@H](F)(Cl)Br").unwrap();
    mol.add_stereo_group(StereoGroup::new(
        StereoGroupKind::Absolute,
        vec![AtomIdx(0)],
    ));
    let cleaned = clean_stereo_groups(&mol);
    assert_eq!(
        cleaned.stereo_groups().len(),
        1,
        "valid group must be preserved"
    );
    assert_eq!(cleaned.stereo_groups()[0].atom_indices, vec![AtomIdx(0)]);
}

#[test]
fn normalize_groups_mixed_nitro_and_azide() {
    // Molecule with both nitro and azide groups
    let mol = parse("C[N+](=O)[O-].N[N+](=O)[O-]").unwrap();
    let result = normalize_groups(&mol);

    // All atoms should be neutral
    let all_neutral = result.atoms().all(|(_, a)| a.charge == 0);
    assert!(all_neutral, "both nitro and azide should be neutralized");
}

// -- Phase 1 fragment-policy fixture tests --------------------------------
// Mirrors validation/standardization_phase1_fixtures.jsonl and
// validation/standardization_phase1_holdout.jsonl (ids referenced in each
// test name/comment). See docs/rfcs/explainable_standardization_phase1_rfc.md.

fn canon(s: &str) -> String {
    chematic_smiles::canonical_smiles(&parse(s).unwrap())
}

/// Assert `remove_salts(input)` keeps exactly the fragment `expected_kept`
/// describes (compared by canonical SMILES, not string equality).
fn assert_kept(id: &str, input: &str, expected_kept: &str) {
    let mol = parse(input).unwrap();
    let result = remove_salts(&mol);
    assert_eq!(
        chematic_smiles::canonical_smiles(&result),
        canon(expected_kept),
        "{id}: expected kept fragment {expected_kept:?} for input {input:?}"
    );
}

#[test]
fn phase1_std_p1_01_sodium_acetate() {
    assert_kept("std-p1-01", "CC(=O)[O-].[Na+]", "CC(=O)[O-]");
}

#[test]
fn phase1_std_p1_02_potassium_benzoate() {
    assert_kept("std-p1-02", "c1ccccc1C(=O)[O-].[K+]", "c1ccccc1C(=O)[O-]");
}

#[test]
fn phase1_std_p1_03_amine_hcl_neutral_form() {
    assert_kept(
        "std-p1-03",
        "CC(C)NCC(O)COc1ccccc1CCOC.Cl",
        "CC(C)NCC(O)COc1ccccc1CCOC",
    );
}

#[test]
fn phase1_std_p1_04_triethylamine_hbr() {
    assert_kept("std-p1-04", "CCN(CC)CC.Br", "CCN(CC)CC");
}

#[test]
fn phase1_std_p1_05_sodium_octanoate_long_chain() {
    assert_kept("std-p1-05", "CCCCCCCC(=O)[O-].[Na+]", "CCCCCCCC(=O)[O-]");
}

#[test]
fn phase1_std_p1_06_amine_hcl_hydrate_3_fragment() {
    assert_kept(
        "std-p1-06",
        "CC(C)NCC(O)c1ccc(O)c(CO)c1.Cl.O",
        "CC(C)NCC(O)c1ccc(O)c(CO)c1",
    );
}

#[test]
fn phase1_std_p1_07_10_zwitterions_are_noop() {
    for (id, smi) in [
        ("std-p1-07", "NCC(=O)[O-]"),
        ("std-p1-08", "C[C@@H](N)C(=O)[O-]"),
        ("std-p1-09", "NCCS(=O)(=O)[O-]"),
        ("std-p1-10", "NCCCC(=O)[O-]"),
    ] {
        assert_kept(id, smi, smi);
    }
}

#[test]
fn phase1_std_p1_11_13_hydrates() {
    assert_kept(
        "std-p1-11",
        "Cn1cnc2c1c(=O)n(C)c(=O)n2C.O",
        "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    );
    assert_kept(
        "std-p1-12",
        "OCC(O)C(O)C(O)C(O)CO.O",
        "OCC(O)C(O)C(O)C(O)CO",
    );
    assert_kept("std-p1-13", "CC(=O)Nc1ccc(O)cc1.O", "CC(=O)Nc1ccc(O)cc1");
}

#[test]
fn phase1_std_p1_14_anhydrous_negative_control() {
    // Same parent as std-p1-11 with no water fragment present: must be a
    // pure no-op, not strip anything from the one remaining fragment.
    assert_kept(
        "std-p1-14",
        "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
        "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    );
}

#[test]
fn phase1_std_p1_15_ferrocene_ionic_fragments() {
    // Fe2+ removed; either cyclopentadienide is an equally valid "kept"
    // choice since both are identical -- assert on the kept fragment's
    // OWN structure, not on which literal input position it came from.
    let mol = parse("[Fe+2].c1cc[cH-]c1.c1cc[cH-]c1").unwrap();
    let result = remove_salts(&mol);
    assert_eq!(
        chematic_smiles::canonical_smiles(&result),
        canon("c1cc[cH-]c1"),
        "std-p1-15: expected a cyclopentadienide fragment kept, Fe2+ removed"
    );
}

#[test]
fn phase1_std_p1_16_ethylmagnesium_bromide() {
    assert_kept("std-p1-16", "CC[Mg+].[Br-]", "CC[Mg+]");
}

#[test]
fn phase1_std_p1_17_cisplatin_single_fragment_noop() {
    // Bonded Cl/N ligands, not disconnected fragments -- must not be
    // misparsed as separate salt components.
    assert_kept("std-p1-17", "Cl[Pt](Cl)(N)N", "Cl[Pt](Cl)(N)N");
}

#[test]
fn phase1_std_p1_18_amine_tartrate_salt() {
    assert_kept(
        "std-p1-18",
        "CN1CCC(CC1)Nc1ccccc1.OC(=O)C(O)C(O)C(=O)O",
        "CN1CCC(CC1)Nc1ccccc1",
    );
}

#[test]
fn phase1_std_p1_19_two_api_cocrystal_ranked_by_size() {
    // Deliberate Phase 1 scope decision (see RFC section 6 / PR notes):
    // no ambiguity-margin detector is implemented -- a heavy-atom-count
    // margin alone cannot distinguish this genuinely-ambiguous cocrystal
    // from std-p1-holdout-10's unambiguous pentane/butane pair (both have
    // a 1-heavy-atom gap). This fixture's outcome is the same size-ranked
    // decision as any other multi-organic-fragment case: acetaminophen
    // (11 heavy atoms) outranks salicylic acid (10). The RFC's proposed
    // "flag close decisions" mechanism is a disclosed, not-yet-implemented
    // gap, not something this test asserts.
    assert_kept(
        "std-p1-19",
        "CC(=O)Nc1ccc(O)cc1.OC(=O)c1ccccc1O",
        "CC(=O)Nc1ccc(O)cc1",
    );
}

#[test]
fn phase1_std_p1_20_22_isotope_containing() {
    assert_kept("std-p1-20", "[13CH3]C(=O)[O-].[Na+]", "[13CH3]C(=O)[O-]");
    assert_kept(
        "std-p1-21",
        "OC[C@H]1O[C@@H]([13OH])[C@H](O)[C@@H](O)[C@@H]1O.O",
        "OC[C@H]1O[C@@H]([13OH])[C@H](O)[C@@H](O)[C@@H]1O",
    );
    assert_kept("std-p1-22", "[13CH3][NH3+].[Cl-]", "[13CH3][NH3+]");
}

#[test]
fn phase1_std_p1_23_tied_spelling_is_atom_order_invariant() {
    // REGRESSION for confirmed defect 1.1(a): pre-RFC largest_fragment()
    // picked a DIFFERENT fragment for these two spellings.
    let a = remove_salts(&parse("CCC.CCN").unwrap());
    let b = remove_salts(&parse("CCN.CCC").unwrap());
    assert_eq!(
        chematic_smiles::canonical_smiles(&a),
        chematic_smiles::canonical_smiles(&b),
        "std-p1-23: kept fragment must not depend on input spelling order"
    );
}

#[test]
fn phase1_std_p1_24_heavy_atom_count_not_explicit_h_count() {
    // REGRESSION for confirmed defect 1.1(b).
    assert_kept("std-p1-24", "CCC.[H]C([H])([H])[H]", "CCC");
}

#[test]
fn phase1_std_p1_25_benzene_hexane_tied_spelling_invariant() {
    let a = remove_salts(&parse("c1ccccc1.CCCCCC").unwrap());
    let b = remove_salts(&parse("CCCCCC.c1ccccc1").unwrap());
    assert_eq!(
        chematic_smiles::canonical_smiles(&a),
        chematic_smiles::canonical_smiles(&b),
        "std-p1-25: 6-heavy-atom tie must resolve the same way regardless of spelling"
    );
}

#[test]
fn phase1_std_p1_26_27_charged_fragments() {
    assert_kept("std-p1-26", "CC(C)(C)[NH3+].[Br-]", "CC(C)(C)[NH3+]");
    assert_kept("std-p1-27", "c1cc[nH+]cc1.[Cl-]", "c1cc[nH+]cc1");
}

#[test]
fn phase1_std_p1_28_choline_chloride_correct_for_the_right_reason() {
    // REGRESSION for confirmed defect 1.1(c): the legacy SaltCatalog's
    // "ammonium" SMARTS also classifies choline itself as salt (verified
    // empirically: SaltCatalog::default().is_salt() on choline alone
    // returns true), so the pre-RFC code only got this right via an
    // unrelated size-comparison fallback. The new structural policy
    // (no catalog consulted by default) must keep choline because it is
    // organic and heavy-atom-dominant, independent of any nitrogen-charge
    // pattern.
    assert_kept("std-p1-28", "C[N+](C)(C)CCO.[Cl-]", "C[N+](C)(C)CCO");
    assert!(
        SaltCatalog::default().is_salt(&parse("C[N+](C)(C)CCO").unwrap()),
        "sanity check: the legacy catalog's false positive on choline alone is still present \
             (expected -- SaltCatalog is intentionally untouched, kept as opt-in legacy behavior)"
    );
}

#[test]
fn phase1_std_p1_29_30_no_organic_fragment_abstains_with_warning() {
    for (id, smi) in [("std-p1-29", "[Na+].[Cl-]"), ("std-p1-30", "[Cl-]")] {
        let mol = parse(smi).unwrap();
        let (_output, record) = select_fragment(&mol, &FragmentPolicy::default());
        assert!(
            record.abstained.is_some(),
            "{id}: expected an abstain reason when no fragment classifies as Kept"
        );
        assert!(
            !record.warnings.is_empty(),
            "{id}: expected a warning surfaced alongside the abstain reason"
        );
    }
}

#[test]
fn phase1_std_p1_31_single_fragment_is_pure_noop() {
    let mol = parse("CCO").unwrap();
    let (output, record) = select_fragment(&mol, &FragmentPolicy::default());
    assert_eq!(chematic_smiles::canonical_smiles(&output), canon("CCO"));
    assert_eq!(record.fragments.len(), 1);
    assert!(record.abstained.is_none());
    assert!(matches!(
        record.fragments[0].decision,
        FragmentDecision::Kept { .. }
    ));
}

#[test]
fn phase1_std_p1_32_empty_molecule_completes_not_errors() {
    let empty = MoleculeBuilder::new().build();
    let (output, record) = select_fragment(&empty, &FragmentPolicy::default());
    assert_eq!(output.atom_count(), 0);
    assert!(record.fragments.is_empty());
    assert!(record.abstained.is_none());
}

#[test]
fn phase1_std_p1_33_duplicate_fragment_genuine_tie() {
    let mol = parse("CC(=O)[O-].[Na+].CC(=O)[O-].[Na+]").unwrap();
    let (output, record) = select_fragment(&mol, &FragmentPolicy::default());
    assert_eq!(
        chematic_smiles::canonical_smiles(&output),
        canon("CC(=O)[O-]")
    );
    let kept_count = record
        .fragments
        .iter()
        .filter(|f| matches!(f.decision, FragmentDecision::Kept { .. }))
        .count();
    assert_eq!(
        kept_count, 1,
        "exactly one fragment is reported Kept even though a duplicate exists"
    );
    let duplicate_marked = record.fragments.iter().any(|f| {
        matches!(
            &f.decision,
            FragmentDecision::Removed { rule_id, .. } if rule_id == "duplicate_of_kept_fragment"
        )
    });
    assert!(
        duplicate_marked,
        "the identical second acetate must be recorded as a duplicate, not a ranked-out loser"
    );
}

#[test]
fn phase1_std_p1_34_carbon_containing_solvate_known_gap() {
    // Disclosed Phase 1 limitation: DMSO contains carbon, so it does NOT
    // match the no-carbon structural salt heuristic and is kept (and
    // must be recorded as Kept-by-ranking, not Removed-as-solvent).
    let mol = parse("c1ccc2c(c1)cccc2.CS(=O)C").unwrap();
    let (output, record) = select_fragment(&mol, &FragmentPolicy::default());
    assert_eq!(
        chematic_smiles::canonical_smiles(&output),
        canon("c1ccc2c(c1)cccc2")
    );
    let dmso_decision = record
        .fragments
        .iter()
        .find(|f| f.snapshot.heavy_atom_count == 4)
        .map(|f| f.decision.clone())
        .expect("DMSO fragment present");
    assert!(
        matches!(
            dmso_decision,
            FragmentDecision::Removed { ref rule_id, .. } if rule_id == "ranked_out_by_size"
        ),
        "DMSO must be recorded as ranked-out-by-size (Kept-class, just smaller), not as a matched salt/solvent rule"
    );
}

#[test]
fn phase1_holdout_01_dopamine_wins_regardless_of_n_threshold() {
    // Dopamine (11 heavy atoms) vs H2SO4 (5 heavy atoms, S+4O -- corrected
    // from the RFC/holdout fixture's original "5 heavy atoms" note, which
    // is right) is not close enough for the open N threshold (RFC
    // section 6) to matter: dopamine wins by ranking regardless of
    // whether H2SO4 is classified StructuralSalt (excluded from
    // candidates) or Kept-but-outranked. This confirms the threshold
    // question only affects the audit log's rule_id here, not the
    // fragment-selection outcome.
    assert_kept(
        "std-p1-holdout-01",
        "CC(N)Cc1ccc(O)c(O)c1.OS(=O)(=O)O",
        "CC(N)Cc1ccc(O)c(O)c1",
    );
}

#[test]
fn phase1_holdout_02_phosphoric_acid_outranks_small_amine_disclosed_gap() {
    // Correction to the original holdout fixture's arithmetic: H3PO4 has
    // 4 oxygens (P + 4*O = 5 heavy atoms), not 3 -- it is NOT a boundary
    // case at N=4, it is squarely over it, same bucket as H2SO4. Unlike
    // holdout-01, ethanolamine (NCCO, 4 heavy atoms) is SMALLER than
    // phosphoric acid (5 heavy atoms), so under a pure heavy-atom-count
    // policy with no named-acid recognition, phosphoric acid outranks the
    // amine and is kept. This is a disclosed limitation of a purely
    // structural, no-named-list policy -- not a bug -- and is exactly the
    // kind of case a future, carefully-scoped small always-strip
    // extension (common mineral acids) would need to address, per RFC
    // section 6.
    assert_kept("std-p1-holdout-02", "NCCO.OP(=O)(O)O", "OP(=O)(O)O");
}

#[test]
fn phase1_holdout_03_calcium_diacetate_via_general_heuristic() {
    // Ca2+ is NOT on the tiny always-strip element list; it is still
    // removed, but via the general structural fallback rule, not
    // "monatomic_always_strip_ion" -- the audit log must say so.
    let mol = parse("CC(=O)[O-].CC(=O)[O-].[Ca+2]").unwrap();
    let (output, record) = select_fragment(&mol, &FragmentPolicy::default());
    assert_eq!(
        chematic_smiles::canonical_smiles(&output),
        canon("CC(=O)[O-]")
    );
    let ca_decision = record
        .fragments
        .iter()
        .find(|f| f.snapshot.heavy_atom_count == 1 && f.snapshot.formula.contains("Ca"))
        .map(|f| f.decision.clone())
        .expect("Ca2+ fragment present");
    assert!(
        matches!(
            ca_decision,
            FragmentDecision::Removed { ref rule_id, .. } if rule_id == "structural_no_carbon_small_fragment"
        ),
        "Ca2+ removal must be attributed to the general heuristic, not the always-strip-ion list"
    );
}

#[test]
fn phase1_holdout_04_05_always_strip_coverage() {
    assert_kept("std-p1-holdout-04", "CC(=O)[O-].[Li+]", "CC(=O)[O-]");
    assert_kept("std-p1-holdout-05", "c1ccc(cc1)CCN.F", "c1ccc(cc1)CCN");
}

#[test]
fn phase1_holdout_06_dihydrate_both_removed() {
    let mol = parse("Cn1cnc2c1c(=O)n(C)c(=O)n2C.O.O").unwrap();
    let (output, record) = select_fragment(&mol, &FragmentPolicy::default());
    assert_eq!(
        chematic_smiles::canonical_smiles(&output),
        canon("Cn1cnc2c1c(=O)n(C)c(=O)n2C")
    );
    let water_removed_count = record
            .fragments
            .iter()
            .filter(|f| matches!(&f.decision, FragmentDecision::Removed { rule_id, .. } if rule_id == "water_always_strip"))
            .count();
    assert_eq!(
        water_removed_count, 2,
        "both water fragments must be independently recorded"
    );
}

#[test]
fn phase1_holdout_07_three_fragment_mixed_rationale() {
    let mol = parse("CN1CCC(CC1)Nc1ccccc1.OC(=O)c1ccccc1O.Cl").unwrap();
    let (output, record) = select_fragment(&mol, &FragmentPolicy::default());
    assert_eq!(
        chematic_smiles::canonical_smiles(&output),
        canon("CN1CCC(CC1)Nc1ccccc1")
    );
    assert_eq!(record.fragments.len(), 3);
    let rule_ids: Vec<&str> = record
        .fragments
        .iter()
        .filter_map(|f| match &f.decision {
            FragmentDecision::Removed { rule_id, .. } => Some(rule_id.as_str()),
            FragmentDecision::Kept { .. } => None,
        })
        .collect();
    assert!(rule_ids.contains(&"ranked_out_by_size"));
    assert!(rule_ids.contains(&"structural_no_carbon_small_fragment"));
}

#[test]
fn phase1_holdout_08_isotope_sugar_hydrate_generalization() {
    assert_kept(
        "std-p1-holdout-08",
        "OC[C@H]1O[C@@H]([13OH])[C@H](O)[C@@H](O)[C@@H]1O.O",
        "OC[C@H]1O[C@@H]([13OH])[C@H](O)[C@@H](O)[C@@H]1O",
    );
}

#[test]
fn phase1_holdout_09_10_non_tie_sanity_checks() {
    assert_kept("std-p1-holdout-09", "CC.CCCCCCCCCC", "CCCCCCCCCC");
    assert_kept("std-p1-holdout-10", "CCCCC.CCCC", "CCCCC");
}

// -- Phase 2 round-2B Parent-function fixture tests -----------------------
// Mirrors validation/tautomer_parent_identity_phase2_fixtures.jsonl's
// tp2-17..22. See docs/rfcs/tautomer_parent_identity_phase2_rfc.md
// section 4.3.

#[test]
fn tp2_17_charge_parent_ammonium_acetate_single_fragment_result() {
    // charge_parent is NOT neutralize_charges: it selects the fragment
    // parent first (acetate: 4 heavy atoms, has carbon, over ammonium:
    // 1 heavy atom, no carbon), THEN neutralizes that one fragment.
    let mol = parse("CC(=O)[O-].[NH4+]").unwrap();
    let (result, record) = charge_parent(&mol);
    assert_eq!(chematic_smiles::canonical_smiles(&result), canon("CC(O)=O"));
    assert_eq!(record.rule_id, "charge_parent_v1");
    assert_eq!(
        record.fragments.len(),
        2,
        "inherited from fragment_parent's selection"
    );
}

#[test]
fn tp2_18_charge_parent_zwitterion_amino_acid() {
    let mol = parse("[NH3+]CC(=O)[O-]").unwrap();
    let (result, _) = charge_parent(&mol);
    assert_eq!(
        chematic_smiles::canonical_smiles(&result),
        canon("NCC(O)=O")
    );
}

#[test]
fn tp2_19_isotope_parent_deuterated_ethanol() {
    let mol = parse("[2H]C([2H])([2H])CO").unwrap();
    let (result, record) = isotope_parent(&mol);
    assert_eq!(
        chematic_smiles::canonical_smiles(&result),
        canon("[H]C(CO)([H])[H]")
    );
    assert!(record.fragments.is_empty());
}

#[test]
fn tp2_20_isotope_parent_preserves_stereo() {
    // Expected value corrected under issue #399: `remove_isotopes`
    // rebuilt the molecule via a bare `MoleculeBuilder` without carrying
    // `stereo_neighbor_order` forward, so this test's original expected
    // string had the anomeric-adjacent ring stereocenter silently
    // flipped -- baked into the "golden" value instead of being caught.
    // Independently verified via RDKit: parse the input, zero every
    // atom's isotope, `AssignStereochemistry(cleanIt=True, force=True)`,
    // and compare `MolToInchi` -- only the value below (not the old one)
    // matches the isotope-free reference's InChI.
    let mol = parse("OC[C@H]1O[C@@H]([13CH2]O)[C@H](O)[C@@H](O)[C@@H]1O").unwrap();
    let (result, _) = isotope_parent(&mol);
    assert_eq!(
        chematic_smiles::canonical_smiles(&result),
        canon("[C@@H]1(CO)O[C@@H]([C@@H](O)[C@@H]([C@H]1O)O)CO")
    );
}

#[test]
fn tp2_21_stereo_parent_alanine() {
    let mol = parse("C[C@@H](N)C(=O)O").unwrap();
    let (result, record) = stereo_parent(&mol);
    assert_eq!(
        chematic_smiles::canonical_smiles(&result),
        canon("C(C(O)=O)(C)N")
    );
    assert!(record.fragments.is_empty());
}

#[test]
fn tp2_22_fragment_parent_hcl_salt() {
    let mol = parse("CN1CCC(CC1)Nc1ccccc1.Cl").unwrap();
    let (result, _) = fragment_parent(&mol);
    assert_eq!(
        chematic_smiles::canonical_smiles(&result),
        canon("N(c2ccccc2)C1CCN(CC1)C")
    );
}

#[test]
fn holdout_04_parent_functions_noop_on_toluene() {
    let mol = parse("Cc1ccccc1").unwrap();
    let expected = chematic_smiles::canonical_smiles(&mol);
    assert_eq!(
        chematic_smiles::canonical_smiles(&fragment_parent(&mol).0),
        expected
    );
    assert_eq!(
        chematic_smiles::canonical_smiles(&charge_parent(&mol).0),
        expected
    );
    assert_eq!(
        chematic_smiles::canonical_smiles(&isotope_parent(&mol).0),
        expected
    );
    assert_eq!(
        chematic_smiles::canonical_smiles(&stereo_parent(&mol).0),
        expected
    );
}

// -----------------------------------------------------------------------
// Issue #407: normalize_zwitterion must never invent a proton. Every
// check here is a property (atom/element/H-count/charge conservation),
// not a spot check against one hardcoded expected string -- per the
// acceptance criteria, RDKit's MolStandardize output is deliberately
// never used as a blind oracle for these.
// -----------------------------------------------------------------------

/// Per-element heavy-atom counts plus total H (implicit `hydrogen_count`
/// summed across all atoms -- explicit H atoms, if any, are also heavy
/// atoms already counted in `counts` and contribute 0 to this sum, so
/// nothing is double-counted).
fn atom_and_h_census(mol: &Molecule) -> (std::collections::BTreeMap<&'static str, u32>, u32) {
    let mut counts: std::collections::BTreeMap<&'static str, u32> =
        std::collections::BTreeMap::new();
    let mut total_h: u32 = 0;
    for (_, atom) in mol.atoms() {
        *counts.entry(atom.element.symbol()).or_insert(0) += 1;
        total_h += atom.hydrogen_count.unwrap_or(0) as u32;
    }
    (counts, total_h)
}

fn net_charge(mol: &Molecule) -> i32 {
    mol.atoms().map(|(_, a)| a.charge as i32).sum()
}

/// Rebuild `mol` with atom indices reversed (last atom becomes index 0,
/// etc.), bonds remapped accordingly -- an isomorphic relabeling used to
/// check that `normalize_zwitterion`'s output doesn't depend on which
/// atom happens to be enumerated first.
fn reverse_atom_order(mol: &Molecule) -> Molecule {
    let n = mol.atom_count();
    // old atom i -> new position (n - 1 - i): the last original atom
    // becomes index 0, and so on.
    let remap: HashMap<AtomIdx, AtomIdx> = (0..n)
        .map(|i| (AtomIdx(i as u32), AtomIdx((n - 1 - i) as u32)))
        .collect();
    let mut slots: Vec<Option<chematic_core::Atom>> = vec![None; n];
    for i in 0..n {
        let old_idx = AtomIdx(i as u32);
        slots[remap[&old_idx].0 as usize] = Some(mol.atom(old_idx).clone());
    }
    let mut builder = MoleculeBuilder::new();
    for slot in slots {
        builder.add_atom(slot.expect("every slot filled by construction"));
    }
    copy_bonds(mol, &mut builder, &remap);
    builder.build()
}

fn assert_conserves_atoms_and_protons(label: &str, smi: &str) {
    let mol = parse(smi).unwrap();
    let result = normalize_zwitterion(&mol);
    let (before_counts, before_h) = atom_and_h_census(&mol);
    let (after_counts, after_h) = atom_and_h_census(&result);
    assert_eq!(
        before_counts, after_counts,
        "{label}: per-element heavy-atom counts must be conserved"
    );
    assert_eq!(
        before_h, after_h,
        "{label}: total H count (implicit) must be conserved -- a proton must never be invented or destroyed"
    );
    assert_eq!(
        result.atom_count(),
        mol.atom_count(),
        "{label}: atom count must be conserved"
    );
}

fn assert_idempotent(label: &str, smi: &str) {
    let mol = parse(smi).unwrap();
    let once = normalize_zwitterion(&mol);
    let once_canon = chematic_smiles::canonical_smiles(&once);
    let reparsed = parse(&once_canon).unwrap();
    let twice = normalize_zwitterion(&reparsed);
    let twice_canon = chematic_smiles::canonical_smiles(&twice);
    assert_eq!(
        once_canon, twice_canon,
        "{label}: normalize_zwitterion must be idempotent (re-running on its own output is a no-op)"
    );
}

#[test]
fn issue407_known_repro_molecules_conserve_atoms_and_are_idempotent() {
    // The exact 3 molecules from issue #407: diazo-N,N'-dioxide-like
    // fragments where neither charged nitrogen has an available proton.
    // Before the fix: the negative oxygen was unconditionally given a
    // proton it invented from nowhere, changing the formula and net
    // charge; the two nitrogens (which have no H to give) were untouched.
    for (label, smi) in [
        ("407-1", "CC12CCCCC1(Br)[N+]([O-])=[N+]2[O-]"),
        ("407-2", "CC12CCCCC1(Cl)[N+]([O-])=[N+]2[O-]"),
        ("407-3", "CC1(C)[N+]([O-])=[N+]([O-])C1(C)Br"),
    ] {
        assert_conserves_atoms_and_protons(label, smi);
        assert_idempotent(label, smi);
        let mol = parse(smi).unwrap();
        let result = normalize_zwitterion(&mol);
        assert_eq!(
            net_charge(&result),
            net_charge(&mol),
            "{label}: net charge must be unchanged when no atom has a transferable proton"
        );
        // Neither +N has an available H, so this pair must be left
        // completely untouched -- not just "charge-neutral overall".
        let canon_before = chematic_smiles::canonical_smiles(&mol);
        let canon_after = chematic_smiles::canonical_smiles(&result);
        assert_eq!(
            canon_before, canon_after,
            "{label}: with no transferable proton anywhere, normalize_zwitterion must be a full no-op"
        );
    }
}

#[test]
fn issue407_genuine_amino_acid_zwitterion_still_normalizes() {
    // Positive control: a real protonation-state zwitterion (donor has
    // an available H) must still be neutralized as before the fix.
    let mol = parse("[NH3+]CC(=O)[O-]").unwrap();
    let result = normalize_zwitterion(&mol);
    assert_eq!(
        net_charge(&result),
        0,
        "genuine zwitterion must reach net charge 0"
    );
    assert_eq!(
        chematic_smiles::canonical_smiles(&result),
        chematic_smiles::canonical_smiles(&parse("NCC(=O)O").unwrap()),
        "genuine zwitterion must still normalize to the neutral amino acid form"
    );
    assert_conserves_atoms_and_protons("genuine-zwitterion", "[NH3+]CC(=O)[O-]");
    assert_idempotent("genuine-zwitterion", "[NH3+]CC(=O)[O-]");
}

#[test]
fn issue407_negative_controls_no_regression() {
    // Nitro and pyridine-N-oxide: both have a +N/-O pair (has_zwitterion
    // fires) but the +N never has an available proton in either case --
    // must be a no-op, not a corruption, and must not regress into
    // "always neutralize nitro/N-oxide groups" (a different, wrong fix).
    for (label, smi) in [
        ("nitrobenzene", "c1ccccc1[N+](=O)[O-]"),
        ("pyridine-n-oxide", "[O-][n+]1ccccc1"),
    ] {
        assert_conserves_atoms_and_protons(label, smi);
        assert_idempotent(label, smi);
        let mol = parse(smi).unwrap();
        let result = normalize_zwitterion(&mol);
        assert_eq!(
            net_charge(&result),
            net_charge(&mol),
            "{label}: no transferable proton anywhere -- net charge must be unchanged"
        );
        assert_eq!(
            chematic_smiles::canonical_smiles(&mol),
            chematic_smiles::canonical_smiles(&result),
            "{label}: must be a full no-op, not a partial/silent mutation"
        );
    }
}

#[test]
fn issue407_atom_permutation_does_not_change_outcome() {
    // A symmetric case with two equidistant candidate positive atoms
    // (both with an available proton) for one negative atom: whichever
    // one the BFS-tie-break picks, the OUTCOME (net charge, formula,
    // canonical structure after normalization) must not depend on which
    // atom happened to be enumerated first.
    let mol = parse("[NH3+]CC([O-])C[NH3+]").unwrap();
    let reversed = reverse_atom_order(&mol);
    let result = normalize_zwitterion(&mol);
    let result_reversed = normalize_zwitterion(&reversed);
    assert_eq!(
        net_charge(&result),
        net_charge(&result_reversed),
        "atom-order permutation must not change the net charge outcome"
    );
    assert_eq!(
        chematic_smiles::canonical_smiles(&result),
        chematic_smiles::canonical_smiles(&result_reversed),
        "atom-order permutation must not change the normalized structure \
             (whichever of the two equidistant +N atoms the BFS tie-break \
             picks, the two are chemically interchangeable here, so the \
             canonical result must converge either way)"
    );
}
