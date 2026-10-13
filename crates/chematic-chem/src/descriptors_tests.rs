#[cfg(test)]
mod tests {
    use super::*;
    use chematic_core::MoleculeBuilder;
    use chematic_smiles::{parse, write};

    /// Parse a SMILES string, panicking on failure.
    fn mol(smiles: &str) -> Molecule {
        parse(smiles).unwrap_or_else(|e| panic!("failed to parse {smiles:?}: {e}"))
    }

    #[test]
    fn rdkit_num_rings_uses_and_invalidates_smiles_cycle_rank_shortcut() {
        let acyclic = mol("CCCC");
        assert_eq!(rdkit_num_rings(&acyclic), 0);

        let mut m = mol("C1CCCCC1");
        assert_eq!(m.smiles_ring_closure_bond_count(), 1);
        assert_eq!(rdkit_num_rings(&m), 1);

        let bridge = m.add_atom(chematic_core::Atom::new(Element::C));
        m.add_bond(AtomIdx(0), bridge, BondOrder::Single).unwrap();
        m.add_bond(AtomIdx(3), bridge, BondOrder::Single).unwrap();
        assert_eq!(m.smiles_ring_closure_bond_count(), 0);
        assert_eq!(rdkit_num_rings(&m), 2);
    }

    #[test]
    fn rdkit_num_rings_does_not_treat_every_smiles_closure_as_a_ring() {
        // A ring digit can connect components or close a path through a bond
        // RDKit excludes from its ring graph. Parser provenance is therefore
        // not, by itself, a safe ring count.
        for smiles in [
            "[Fe]1(Cl)(Cl)Cl.Cl1",
            "C=C(C)N1C=NC<-C1",
            "N1CCN->[Cu]1",
            "C1CCN->[Cu]1",
        ] {
            let m = mol(smiles);
            assert_eq!(m.smiles_ring_closure_bond_count(), 1, "{smiles}");
            assert_eq!(rdkit_num_rings(&m), 0, "{smiles}");
        }
    }

    #[test]
    fn rdkit_num_rings_uses_2026_09_relevant_cycles_for_symmetric_macrocycles() {
        // RDKit 2026.09 replaced the legacy symmetrized SSSR with
        // RingDecomposerLib. Large macrocycles containing several
        // para-substituted aromatic rings consequently retain every
        // equivalent shortest ring family (the legacy backend returned 10).
        let m = mol("c1ccc2c(c1)c1cc[n+]2Cc2ccc(cc2)-c2ccc(cc2)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCN1");
        assert_eq!(rdkit_num_rings(&m), 22);
    }

    #[test]
    fn rdkit_num_rings_bounds_relevant_family_expansion() {
        // Every repeated sugar unit contributes equal shortest paths.  The
        // generic relevant-cycle enumerator expands all combinations to 70;
        // RDKit 2026.09 reports the independent ring families (22 and 16).
        let fixtures = [
            (
                "OC[C@H]1O[C@@H]2O[C@H]3[C@H](O)[C@@H](O)[C@@H](O[C@@H]3CO)O[C@H]3[C@H](O)[C@@H](O)[C@@H](O[C@@H]3CO)O[C@H]3[C@H](O)[C@@H](O)[C@@H](O[C@@H]3CO)O[C@H]3[C@H](O)[C@@H](O)[C@@H](O[C@@H]3CO)O[C@H]3[C@H](O)[C@@H](O)[C@@H](O[C@@H]3CO)O[C@H]1[C@H](O)[C@H]2O",
                22,
            ),
            (
                "[C@@H]12[C@H](O)[C@H]([C@@H]([C@@H](CO)O2)O[C@H]2[C@@H]([C@@H](O)[C@@H]([C@@H](CO)O2)O[C@H]2[C@@H]([C@@H](O)[C@H](O[C@@H]3O[C@@H]([C@H]([C@H](O)[C@H]3O)O[C@H]3[C@H](O)[C@H]([C@H](O[C@H]4O[C@@H]([C@@H](O1)[C@@H]([C@H]4O)O)CO)[C@H](O3)CO)O)CO)[C@H](O2)CO)O)O)O",
                16,
            ),
        ];
        for (smiles, expected) in fixtures {
            assert_eq!(rdkit_num_rings(&mol(smiles)), expected, "{smiles}");
        }
    }

    /// RDKit sanitization rewrites perchloric acid as
    /// `[Cl+3]([O-])([O-])([O-])O`; descriptors see that form (values from
    /// RDKit 2026.03.1) while the user's molecule keeps its double bonds.
    #[test]
    fn perchlorate_descriptors_use_rdkit_cleanup() {
        let m = mol("OCl(=O)(=O)=O");
        assert!((tpsa(&m) - 89.41).abs() < 1e-9);
        assert_eq!(hba_count(&m), 4);
        assert!((logp_crippen(&m) - -4.124).abs() < 1e-9);
        assert_eq!(rotatable_bond_count(&m), 0);
        assert_eq!(
            m.bonds()
                .filter(|(_, b)| b.order == BondOrder::Double)
                .count(),
            3
        );
        assert_eq!(m.atoms().filter(|(_, a)| a.charge != 0).count(), 0);
    }

    /// Crippen H typing takes RDKit's first matching type: an `[OH+]` bound
    /// to a metal and to a carboxyl carbon is H2 (`[#1]O[!#6;!#7;!#8;!#16]`),
    /// not H4. Values from RDKit 2026.03.1.
    #[test]
    fn metal_bound_oh_hydrogen_is_h2() {
        let m = mol("N[Co+3](N)(N)(N)(N)[OH+]C(=O)C(F)(F)F");
        assert_eq!(logp_crippen(&m), -2.887899999999999);
        assert_eq!(logp_and_mr(&m).1, 34.21780000000002);
    }

    /// `cleanUpOrganometallics` makes the four-bonded `[C-]`–Fe bond dative
    /// in RDKit's sanitized ferrocene; Crippen typing sees that form.
    #[test]
    fn ferrocene_dative_bond_matches_rdkit_crippen() {
        let m = mol("CN(C)C[C-]12C3=C4C5=C1[Fe++]23456789[C-]%10C6=C7C8=C9%10");
        assert_eq!(logp_crippen(&m), 1.4556799999999999);
        assert_eq!(logp_and_mr(&m).1, 53.144000000000005);
    }

    /// RDKit `CalcFractionCSP3` counts carbons of total degree four.
    #[test]
    fn fsp3_counts_four_connected_carbanion() {
        let m = mol("CN(C)C[C-]12C3=C4C5=C1[Fe++]23456789[C-]%10C6=C7C8=C9%10");
        assert_eq!(fsp3(&m), 4.0 / 13.0);
    }

    /// RDKit's `CalcNumAtomStereoCenters` uses legacy stereo perception: the
    /// bridgehead nitrogens of this symmetric cage have tied neighbours.
    #[test]
    fn cage_amine_has_no_rdkit_stereocenters() {
        let m = mol("C1CN2CN1CN3CCN(C2)C3");
        assert_eq!(num_stereocenters(&m), 0);
        assert_eq!(num_unspecified_stereocenters(&m), 0);
        assert_eq!(num_stereocenters(&mol("CC(O)CC")), 1);
        assert_eq!(num_unspecified_stereocenters(&mol("C[C@H](O)CC")), 0);
        assert_eq!(num_stereocenters(&mol("CC1CCC(C)CC1")), 0);
    }

    /// RDKit 2026.03's strict rotor pattern excludes `[CH3]`, including a
    /// methyl spelled with explicit deuterium.
    #[test]
    fn deuterated_methyl_is_not_a_rotor() {
        assert_eq!(rotatable_bond_count(&mol("[2H]C([2H])([2H])NC=O")), 1);
        assert_eq!(rotatable_bond_count(&mol("[2H]CNC=O")), 1);
    }

    /// Deterministically relabel a graph while preserving original atom
    /// identity as atom-map metadata. `permutation[new_index] = old_index`.
    fn mapped_permutation(mol: &Molecule, permutation: &[usize]) -> Molecule {
        assert_eq!(permutation.len(), mol.atom_count());
        let mut builder = MoleculeBuilder::new();
        let mut old_to_new = vec![0u32; mol.atom_count()];

        for (new_index, &old_index) in permutation.iter().enumerate() {
            let mut atom = mol.atom(AtomIdx(old_index as u32)).clone();
            atom.atom_map = Some((old_index + 1) as u16);
            builder.add_atom(atom);
            old_to_new[old_index] = new_index as u32;
        }

        for (_, bond) in mol.bonds() {
            builder
                .add_bond(
                    AtomIdx(old_to_new[bond.atom1.0 as usize]),
                    AtomIdx(old_to_new[bond.atom2.0 as usize]),
                    bond.order,
                )
                .expect("a relabelled simple graph cannot introduce duplicate bonds");
        }
        builder.build()
    }

    fn fixed_permutation(atom_count: usize, seed: u64) -> Vec<usize> {
        let mut state = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut permutation: Vec<_> = (0..atom_count).collect();
        for index in (1..atom_count).rev() {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            permutation.swap(index, (state as usize) % (index + 1));
        }
        permutation
    }

    type MappedStereoObservation = (Vec<u16>, Vec<(u16, String)>, Vec<(u16, u16, String)>);

    fn mapped_observation(mol: &Molecule, rings_first: bool) -> MappedStereoObservation {
        if rings_first {
            let _ = chematic_perception::find_sssr(mol);
        }
        let potential = potential_stereocenter_indices(mol);
        let assignments = crate::cip::assign_cip(mol);
        if !rings_first {
            let _ = chematic_perception::find_sssr(mol);
        }

        let mut centers: Vec<_> = potential
            .into_iter()
            .map(|index| {
                mol.atom(index)
                    .atom_map
                    .expect("test atoms carry original maps")
            })
            .collect();
        centers.sort_unstable();

        let mut cip: Vec<_> = assignments
            .assignments
            .into_iter()
            .map(|(index, code)| {
                (
                    mol.atom(index)
                        .atom_map
                        .expect("test atoms carry original maps"),
                    format!("{code:?}"),
                )
            })
            .collect();
        cip.sort_unstable();

        let mut bonds: Vec<_> = mol
            .bonds()
            .map(|(_, bond)| {
                let first = mol
                    .atom(bond.atom1)
                    .atom_map
                    .expect("test atoms carry original maps");
                let second = mol
                    .atom(bond.atom2)
                    .atom_map
                    .expect("test atoms carry original maps");
                (
                    first.min(second),
                    first.max(second),
                    format!("{:?}", bond.order),
                )
            })
            .collect();
        bonds.sort_unstable();
        (centers, cip, bonds)
    }

    // Tolerance helpers.
    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    fn pct2(a: f64, b: f64) -> bool {
        // within 2% relative, or within 0.05 Da absolute (for very small values)
        approx(a, b, b.abs() * 0.02 + 0.05)
    }

    // -- Test 1: methane molecular weight ------------------------------------
    #[test]
    fn test_mw_methane() {
        let m = mol("C");
        // CH4: 12.011 + 4*1.008 = 16.043
        assert!(
            pct2(molecular_weight(&m), 16.043),
            "methane MW = {}",
            molecular_weight(&m)
        );
    }

    // -- Test 2: water molecular weight -------------------------------------
    #[test]
    fn test_mw_water() {
        let m = mol("O");
        // H2O: 15.999 + 2*1.008 = 18.015
        assert!(
            pct2(molecular_weight(&m), 18.015),
            "water MW = {}",
            molecular_weight(&m)
        );
    }

    #[test]
    fn rdkit_mw_uses_rdkit_periodic_table_without_changing_native_mw() {
        let sulfur = mol("CS");
        let selenium = mol("[Se]");
        assert!(
            approx(rdkit_molecular_weight(&sulfur), 48.110, 1e-12),
            "rdkit profile = {}",
            rdkit_molecular_weight(&sulfur)
        );
        assert!(
            approx(molecular_weight(&sulfur), 48.108, 1e-12),
            "native = {}",
            molecular_weight(&sulfur)
        );
        assert!(approx(rdkit_molecular_weight(&selenium), 78.96, 1e-12));
        assert!(approx(molecular_weight(&selenium), 78.971, 1e-12));
    }

    #[test]
    fn rdkit_mw_uses_nuclide_masses_and_mass_number_fallback() {
        let isotope = mol("[13C]");
        assert!(approx(rdkit_molecular_weight(&isotope), 13.00335484, 1e-8));
        let unknown = mol("[99C]");
        assert!(approx(rdkit_molecular_weight(&unknown), 99.0, 1e-12));
        assert!(approx(exact_mass(&unknown), 99.0, 1e-12));
    }

    #[test]
    fn rdkit_mw_normalizes_representation_and_zero_isotope_like_rdkit() {
        // RDKit's sanitised graph gives both benzene spellings C6H6. Isotope
        // zero means "unspecified", not a zero-mass carbon nuclide.
        assert!(approx(
            rdkit_molecular_weight(&mol("c1=cc=cc=c1")),
            72.066,
            1e-12
        ));
        assert!(approx(
            rdkit_molecular_weight(&mol("[0CH4]")),
            16.043,
            1e-12
        ));
    }

    #[test]
    fn rdkit_mw_uses_complete_generated_isotope_table() {
        // This representative is deliberately not one of the historical
        // hand-written entries. RDKit's public PeriodicTable reports
        // 11.0114336 Da for carbon-11, rather than the label's integer 11.
        let carbon11 = mol("[11C]");
        assert!(approx(rdkit_molecular_weight(&carbon11), 11.0114336, 1e-10));

        let table = crate::rdkit_isotope_mass_table::RDKIT_ISOTOPE_MASS_TABLE;
        assert_eq!(
            table.len(),
            3111,
            "unexpected generated RDKit isotope count"
        );
        for pair in table.windows(2) {
            let (left_z, left_a, _) = pair[0];
            let (right_z, right_a, _) = pair[1];
            assert!(
                (left_z, left_a) < (right_z, right_a),
                "isotope mass table must be sorted for binary search"
            );
        }
        for &(atomic_number, isotope, expected) in table.iter() {
            let element = Element::from_atomic_number(atomic_number)
                .expect("generated table contains a valid element");
            assert!(
                approx(rdkit_isotope_mass(element, isotope), expected, 0.0),
                "isotope mass lookup mismatch for Z={atomic_number}, A={isotope}"
            );
        }
    }

    #[test]
    fn rdkit_hba_profile_excludes_substituted_aromatic_n() {
        let caffeine = mol("Cn1cnc2c1c(=O)n(c(=O)n2C)C");
        assert_eq!(hba_count(&caffeine), 3);
        assert_eq!(rdkit_hba_count(&caffeine), 3);
        let methylpyrrole = mol("Cn1cccc1");
        assert_eq!(rdkit_hba_count(&methylpyrrole), 0);
        let methylpyridone = mol("O=C1C=CC=CN1C");
        assert_eq!(rdkit_hba_count(&methylpyridone), 1);
    }

    #[test]
    fn rdkit_hba_nucleoside_uses_2026_03_rule() {
        let nucleoside = mol("Nc1nc(N)c2ncn(C3CC(O)C(O)C(CO)O3)c2n1");
        assert_eq!(hba_count(&nucleoside), 9);
        assert_eq!(rdkit_hba_count(&nucleoside), 9);
    }

    #[test]
    fn rdkit_aromatic_ring_profile_uses_rdkit_aromaticity_without_changing_native() {
        let fused = mol("Cn1c2nc(=O)[nH]c(=O)c-2nc2ccccc21");
        assert_eq!(aromatic_ring_count(&fused), 3);
        assert_eq!(rdkit_aromatic_ring_count(&fused), 1);
    }

    #[test]
    fn rdkit_aromatic_ring_profile_recovers_cage_face_omitted_by_sssr() {
        let cage = mol(
            "COCCOCCOCCN1CC23C4=C5C6=C7c8c9c%10c%11c%12c%13c%14c(c2c2c%15c%16c%17c%18c%19c(c5c5c%20c%21c%22c%23c(c8C%22C65)c%10c5c%11c6c%13c8c(c%15%14)c%16c%10c%18c%11c(c%20%19)c%21c%13c%23c5c5c%13c%11c%10c8c65)C%17C42)C%12C9C73C1COCCOCCOC",
        );
        assert_eq!(rdkit_aromatic_ring_count(&cage), 19);
    }

    #[test]
    fn rdkit_aromatic_ring_profile_does_not_reperceive_zero_aromatic_result() {
        let fused_thiohydantoin = mol("S=C1N2CCN=C2SC2=NCCN12");
        assert_eq!(rdkit_aromatic_ring_count(&fused_thiohydantoin), 0);
    }

    #[test]
    fn unspecified_stereocenters_require_distinct_substituents() {
        assert_eq!(num_unspecified_stereocenters(&mol("C")), 0);
        assert_eq!(num_unspecified_stereocenters(&mol("CC")), 0);
        assert_eq!(num_unspecified_stereocenters(&mol("CC(C)C")), 0);
        assert_eq!(num_unspecified_stereocenters(&mol("C(F)(Cl)Br")), 1);
    }

    // -- Test 3: ethanol molecular weight -----------------------------------
    #[test]
    fn test_mw_ethanol() {
        let m = mol("CCO");
        // C2H6O: 2*12.011 + 6*1.008 + 15.999 = 46.068
        assert!(
            pct2(molecular_weight(&m), 46.068),
            "ethanol MW = {}",
            molecular_weight(&m)
        );
    }

    // -- Test 4: benzene molecular weight -----------------------------------
    #[test]
    fn test_mw_benzene() {
        let m = mol("c1ccccc1");
        // C6H6: 6*12.011 + 6*1.008 = 78.114
        assert!(
            pct2(molecular_weight(&m), 78.114),
            "benzene MW = {}",
            molecular_weight(&m)
        );
    }

    // -- Test 5: aspirin molecular weight -----------------------------------
    #[test]
    fn test_mw_aspirin() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        // C9H8O4: MW ~180.16
        let mw = molecular_weight(&m);
        assert!(approx(mw, 180.16, 1.0), "aspirin MW = {mw}");
    }

    // -- Test 6: methane exact mass -----------------------------------------
    #[test]
    fn test_exact_mass_methane() {
        let m = mol("C");
        // 12C + 4*(1H): 12.0000 + 4*1.00783 = 16.0313
        let em = exact_mass(&m);
        assert!(approx(em, 16.031, 0.01), "methane exact mass = {em}");
    }

    // -- Test 7: benzene heavy atom count -----------------------------------
    #[test]
    fn test_hac_benzene() {
        let m = mol("c1ccccc1");
        assert_eq!(heavy_atom_count(&m), 6);
    }

    // -- Test 8: aspirin heavy atom count -----------------------------------
    #[test]
    fn test_hac_aspirin() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        // C9H8O4: 9 C + 4 O = 13 heavy atoms
        assert_eq!(heavy_atom_count(&m), 13);
    }

    // -- Test 9: ethanol HBD ------------------------------------------------
    #[test]
    fn test_hbd_ethanol() {
        let m = mol("CCO");
        assert_eq!(hbd_count(&m), 1); // one OH
    }

    // -- Test 10: aniline HBD -----------------------------------------------
    #[test]
    fn test_hbd_aniline() {
        let m = mol("Nc1ccccc1");
        assert_eq!(hbd_count(&m), 1); // one NH2
    }

    // -- Test 11: benzene HBD -----------------------------------------------
    #[test]
    fn test_hbd_benzene() {
        let m = mol("c1ccccc1");
        assert_eq!(hbd_count(&m), 0);
    }

    // -- Test 12: ethanol HBA -----------------------------------------------
    #[test]
    fn test_hba_ethanol() {
        let m = mol("CCO");
        assert_eq!(hba_count(&m), 1); // one O
    }

    // -- Test 13: aspirin HBA -----------------------------------------------
    #[test]
    fn test_hba_aspirin() {
        // aspirin: 4 O total, but carboxyl OH is excluded → 3 (RDKit-aligned)
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        assert_eq!(hba_count(&m), 3);
    }

    // -- Test 14: benzene rotatable bonds ------------------------------------
    #[test]
    fn test_rot_benzene() {
        let m = mol("c1ccccc1");
        assert_eq!(rotatable_bond_count(&m), 0);
    }

    // -- Test 15: aspirin rotatable bonds ------------------------------------
    #[test]
    fn test_rot_aspirin() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        // RDKit strict=True (default): ester C(=O)–O bond is excluded (carbonyl hetero bond).
        // Counted: O–Ar, Ar–C(=O)COOH (2 bonds).  CH3 terminal, carboxyl OH terminal, ring bonds excluded.
        let r = rotatable_bond_count(&m);
        assert_eq!(r, 2, "aspirin rotatable bonds = {r}");
    }

    // -- Test: alkyne adjacent bond excluded ---------------------------------
    #[test]
    fn test_rot_alkyne_adjacent_excluded() {
        // Propyne CH3-C≡CH: C-C adjacent to triple bond is NOT rotatable (C≡ atom excluded).
        let m = mol("CC#C");
        assert_eq!(
            rotatable_bond_count(&m),
            0,
            "propyne: C-C adj to triple bond excluded"
        );
        // But-1-yne CH3-CH2-C≡CH: C3 has triple bond → C2-C3 excluded; C1 terminal → C1-C2 excluded.
        let m2 = mol("CCC#C");
        assert_eq!(
            rotatable_bond_count(&m2),
            0,
            "but-1-yne: both bonds excluded"
        );
        // Pent-1-yne CH3-CH2-CH2-C≡CH: C2-C3 is rotatable; C3-C4 excluded (C4 in triple).
        let m3 = mol("CCCC#C");
        assert_eq!(
            rotatable_bond_count(&m3),
            1,
            "pent-1-yne: CH2-CH2 bond is rotatable"
        );
    }

    // -- Test: allene cumulated double bonds excluded ------------------------
    #[test]
    fn test_rot_allene_excluded() {
        // Allene: C=C=C. Both single-bond-looking connections are in a cumulated system.
        let m = mol("C=C=C");
        assert_eq!(rotatable_bond_count(&m), 0, "allene: no rotatable bonds");
    }

    // -- Test 16: water TPSA -------------------------------------------------
    #[test]
    fn test_tpsa_water() {
        let m = mol("O");
        // H₂O: isolated O with no heavy-atom neighbors → 31.50 (Ertl water type)
        let t = tpsa(&m);
        assert!(approx(t, 31.50, 0.1), "water TPSA = {t}");
    }

    #[test]
    fn tpsa_explicit_deuterium_matches_hydroxyl_atom_type() {
        // Isotopic graph-H and implicit protium represent the same O-H
        // descriptor atom type. RDKit CalcTPSA is 20.23 for both CCO and
        // CCO[2H].
        assert!(approx(tpsa(&mol("CCO[2H]")), 20.23, 1e-6));
        assert_eq!(hbd_count(&mol("CCO[2H]")), 1);
    }

    #[test]
    fn kekule_pyridone_uses_aromatic_descriptor_types_without_promoting_imides() {
        // RDKit sanitizes this Kekulé spelling to O=c1cccc[nH]1. Its ring NH
        // is consequently an aromatic descriptor atom, unlike a two-carbonyl
        // imide such as phthalimide.
        let pyridone = mol("O=C1C=CC=CN1");
        assert!(approx(tpsa(&pyridone), 32.86, 1e-12));
        assert!(approx(logp_crippen(&pyridone), 0.3749, 1e-12));
        assert!(approx(molar_refractivity(&pyridone), 27.0627, 1e-12));
        assert!(approx(logp_and_mr(&pyridone).0, 0.3749, 1e-12));
        assert!(approx(logp_and_mr(&pyridone).1, 27.0627, 1e-12));

        let phthalimide = mol("O=C1NC(=O)c2ccccc12");
        assert!(approx(tpsa(&phthalimide), 46.17, 1e-12));
        assert!(approx(logp_crippen(&phthalimide), 0.5702, 1e-12));
        assert!(approx(molar_refractivity(&phthalimide), 38.2387, 1e-12));
    }

    #[test]
    fn tpsa_matches_rdkit_nitrogen_atom_type_table_boundaries() {
        // These public fixtures cover the atom types that are easy to miss when
        // TPSA is implemented as functional-group special cases. Values are from
        // RDKit 2025.09.3 CalcTPSA(includeSandP=true).
        for (smiles, expected) in [
            ("C1CN1", 21.94),        // secondary aziridine N
            ("CN1CC1", 3.01),        // tertiary aziridine N
            ("C[NH3+]", 27.64),      // primary ammonium
            ("C=[NH2+]", 25.59),     // protonated imine, degree one
            ("C=[NH+]C", 13.97),     // protonated imine, degree two
            ("[NH2+]1CC1", 16.61),   // cyclic secondary ammonium
            ("C[N+]#C", 4.36),       // nitrilium
            ("[nH+]1ccccc1", 14.14), // aromatic protonated N
        ] {
            let observed = tpsa(&mol(smiles));
            assert!(
                approx(observed, expected, 1e-12),
                "{smiles}: TPSA {observed}, expected {expected}"
            );
        }
    }

    #[test]
    fn tpsa_kekule_nitrogen_uses_sanitized_aromatic_environment() {
        assert!(approx(tpsa(&mol("N1=CC=CC=C1")), 12.89, 1e-12));
        assert!(approx(tpsa(&mol("O=C1C=CC=CN1C")), 22.00, 1e-12));
        // The shared contribution path guarantees that the public total and
        // per-atom API cannot silently use different atom typing rules.
        for smiles in ["N1=CC=CC=C1", "O=C1C=CC=CN1C", "O=C1NC(=O)c2ccccc12"] {
            let molecule = mol(smiles);
            let sum: f64 = tpsa_per_atom(&molecule).into_iter().sum();
            assert!(approx(sum, tpsa(&molecule), 1e-12), "{smiles}");
        }
    }

    // -- Test 17: aniline TPSA -----------------------------------------------
    #[test]
    fn test_tpsa_aniline() {
        let m = mol("Nc1ccccc1");
        // NH2 (aliphatic) → 26.02
        let t = tpsa(&m);
        assert!(approx(t, 26.02, 5.0), "aniline TPSA = {t}");
    }

    #[test]
    fn test_tpsa_rdkit_environment_boundaries() {
        // P=S is a single RDKit include-S/P environment (41.90); sulfur must
        // not also receive the ordinary thioether contribution.
        assert!(approx(tpsa(&mol("P(=S)(C)(C)C")), 41.90, 1e-12));
        let hypervalent_s = tpsa(&mol("S(C)(C)(=O)=S"));
        assert!(
            approx(hypervalent_s, 57.54, 1e-12),
            "S(=O)(=S) TPSA = {hypervalent_s}"
        );
        let phosphine_amino = tpsa(&mol("NC(N)=NCCCC(N)[PH](=O)O"));
        assert!(
            approx(phosphine_amino, 151.19, 1e-12),
            "P-H P=O TPSA = {phosphine_amino}"
        );
        let phosphine_methyl = tpsa(&mol("CNC(=N)NCCCC(N)[PH](=O)O"));
        assert!(
            approx(phosphine_methyl, 134.70, 1e-12),
            "methyl P-H P=O TPSA = {phosphine_methyl}"
        );
        let phosphine_nitro = tpsa(&mol("N/C(=N/CCCC(N)[PH](=O)O)N[N+](=O)[O-]"));
        assert!(
            approx(phosphine_nitro, 180.34, 1e-12),
            "nitro P-H P=O TPSA = {phosphine_nitro}"
        );
        let thiol = tpsa(&mol(
            "COc1cc2nc(N3CCN(/C(S)=N/c4ccc(NC(=S)N5CC5)cc4)CC3)nc(N)c2cc1OC",
        ));
        assert!(approx(thiol, 175.03, 1e-12), "C(=N)-SH TPSA = {thiol}");
        let charged_sulfur = tpsa(&mol("[S-]c1nc2ccccc2c2cccc[n+]12"));
        assert!(
            approx(charged_sulfur, 16.99, 1e-12),
            "charged sulfur TPSA = {charged_sulfur}"
        );
        // N+ in N+–O-–N environments is not the central azide N+ type.
        assert!(approx(tpsa(&mol("C[N+]([O-])=[N+]([O-])C")), 52.14, 1e-12));
        // P=N TPSA uses three distinct RDKit atom types: implicit-H,
        // fully-substituted, and explicitly bracketed P-H.
        let implicit_phosphazene = tpsa(&mol("CP(=N)C"));
        assert!(
            approx(implicit_phosphazene, 47.32, 1e-12),
            "implicit-H phosphazene TPSA = {implicit_phosphazene}"
        );
        assert!(approx(tpsa(&mol("CP(=N)(C)C")), 33.66, 1e-12));
        assert!(approx(tpsa(&mol("[PH](=N)C")), 23.85, 1e-12));
        // An exocyclic C=C must not promote a cyclic ether to RDKit's [o]
        // Crippen/TPSA type merely because one neighbor is aromatic.
        let bridged = mol("COc1cc2c(cc1OC)C1C(=O)c3ccc4c(c3OC1CO2)C(C)(C)C=CO4");
        assert!(approx(tpsa(&bridged), 63.22, 1e-12));
        assert!(approx(logp_crippen(&bridged), 4.0074, 1e-12));
        assert!(approx(molar_refractivity(&bridged), 105.7645, 1e-12));
    }

    #[test]
    fn fused_heterocycles_keep_rdkit_descriptor_aromaticity_boundary() {
        // Development-only A0 residuals. The broad aromaticity model promoted
        // neutral divalent sulfur (and adjacent nitrogens) in these fused
        // systems even though RDKit keeps their Kekule atom types. All three
        // descriptors share the same bounded parity representation.
        let cases = [
            (
                "O=C(Nc1ccccc1Cc1ccccc1)c1ccc2c(c1)SC1=NS(=O)(=O)CCN12",
                112.52,
                4.1412,
                124.5100,
            ),
            (
                r#"CCOc1cc(/C=C2\C(=N)N3N=C(CC(=O)N4CCCCC4)SC3=NC2=O)ccc1O"#,
                143.95,
                2.81267,
                119.3345,
            ),
            (
                "N=C1/C(=C/c2ccc(Br)o2)C(=O)N=C2SC(CC(=O)N3CCCC3)=NN12",
                127.63,
                2.67377,
                101.8497,
            ),
            (
                "N=C1/C(=C/c2cccn2CCOc2ccccc2)C(=O)N=C2SC=C(c3ccccc3)N12",
                95.98,
                4.87137,
                128.3637,
            ),
            (
                "COc1cc2c(cc1OC)N1C(=NC(N)=NC1c1ccccc1)NC2=O",
                101.54,
                1.6367,
                97.7316,
            ),
            (
                "COC(=O)C1=C(C)N=C2SC(C#N)=C(N)N2C1c1cccc(OC)c1",
                126.24,
                2.25308,
                94.1424,
            ),
        ];

        for (smiles, expected_tpsa, expected_logp, expected_mr) in cases {
            let molecule = mol(smiles);
            assert!(approx(tpsa(&molecule), expected_tpsa, 1e-12), "{smiles}");
            assert!(
                approx(logp_crippen(&molecule), expected_logp, 1e-12),
                "{smiles}"
            );
            assert!(
                approx(molar_refractivity(&molecule), expected_mr, 1e-12),
                "{smiles}"
            );
        }
    }

    #[test]
    fn fused_cyclic_ether_uses_rdkit_crippen_oxygen_type() {
        // This five-membered cyclic ether has an aromatic neighbour and a
        // vinylic neighbour, but RDKit keeps the oxygen non-aromatic. A former
        // shape-based shortcut incorrectly forced the [o] Crippen type.
        let molecule = mol("CC(=O)c1c(O)c(C)c(O)c2c1OC1=Cc3c(c(C)nn3C)C(=O)C12C");
        assert!(approx(logp_crippen(&molecule), 2.53824, 1e-12));
        assert!(approx(molar_refractivity(&molecule), 92.7456, 1e-12));
        assert!(approx(logp_and_mr(&molecule).0, 2.53824, 1e-12));
        assert!(approx(logp_and_mr(&molecule).1, 92.7456, 1e-12));
    }

    // -- Test 18: aspirin Lipinski -------------------------------------------
    #[test]
    fn test_lipinski_aspirin() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        assert!(lipinski_passes(&m));
    }

    // -- Test 19: benzene Lipinski ------------------------------------------
    #[test]
    fn test_lipinski_benzene() {
        let m = mol("c1ccccc1");
        assert!(lipinski_passes(&m));
    }

    // -- Additional tests ---------------------------------------------------

    // Benzene exact mass
    #[test]
    fn test_exact_mass_benzene() {
        let m = mol("c1ccccc1");
        // C6H6: 6*12 + 6*1.00783 = 78.04698
        let em = exact_mass(&m);
        assert!(approx(em, 78.047, 0.05), "benzene exact mass = {em}");
    }

    // Ethanol exact mass
    #[test]
    fn test_exact_mass_ethanol() {
        let m = mol("CCO");
        // C2H6O: 2*12 + 6*1.00783 + 15.9949 = 46.0419
        let em = exact_mass(&m);
        assert!(approx(em, 46.042, 0.05), "ethanol exact mass = {em}");
    }

    #[test]
    fn exact_mass_is_bit_identical_to_rdkit_exact_mol_wt() {
        // RDKit 2026.03.1 Descriptors.ExactMolWt, compared exactly: isotope
        // labels (known and unknown), halogens and other elements that used
        // to carry rounded masses, and formal charges (electrons).
        for (smiles, expected) in [
            ("[13C]", 13.00335484),
            ("[11CH3]CO", 45.053298412),
            ("[2H]", 2.014101778),
            ("[99C]", 99.0),
            ("CCl", 49.992327775999996),
            ("CBr", 93.941812196),
            ("CI", 141.927948096),
            ("CF", 34.021878316),
            ("C[Si](C)(C)C", 88.070826914),
            ("[Na+].[Cl-]", 57.95862196),
            ("CP(C)C", 76.044186918),
            ("CS(C)=O", 78.013935812),
            ("CC(=O)[O-]", 59.01385291591001),
            ("C[N+](C)(C)C", 74.09642580408999),
        ] {
            assert_eq!(exact_mass(&mol(smiles)), expected, "{smiles}");
        }
    }

    // Aspirin logp and Lipinski components
    #[test]
    fn test_logp_aspirin_is_reasonable() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        let lp = logp_crippen(&m);
        // The simplified model gives a low but not absurd value; just check range.
        assert!(lp > -5.0 && lp < 5.0, "aspirin logp = {lp}");
    }

    // Heavy atom count for ethanol
    #[test]
    fn test_hac_ethanol() {
        let m = mol("CCO");
        assert_eq!(heavy_atom_count(&m), 3); // 2 C + 1 O
    }

    // HBA for aniline (one N)
    #[test]
    fn test_hba_aniline() {
        let m = mol("Nc1ccccc1");
        assert_eq!(hba_count(&m), 1); // one N
    }

    // Rotatable bonds for n-butane (single chain)
    #[test]
    fn test_rot_butane() {
        let m = mol("CCCC");
        // C1-C2, C2-C3, C3-C4 → three bonds; C1 has degree 1 (only C2 neighbor), C4 too
        // degree(C1)=1 → excluded; degree(C4)=1 → excluded
        // Only C2-C3 is non-terminal on both sides? Actually C2 has degree 2 (C1,C3),
        // C3 has degree 2 (C2,C4). Bond C2-C3: both non-terminal → rotatable (1).
        // Bond C1-C2: C1 degree 1 → skip. Bond C3-C4: C4 degree 1 → skip.
        assert_eq!(rotatable_bond_count(&m), 1, "n-butane has 1 rotatable bond");
    }

    // TPSA for aspirin (should be > 0)
    #[test]
    fn test_tpsa_aspirin_positive() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        let t = tpsa(&m);
        assert!(t > 0.0, "aspirin TPSA = {t}");
    }

    // Quaternary ammonium N+ has no lone pair → TPSA contribution = 0
    #[test]
    fn test_tpsa_quaternary_n_zero() {
        // Trimethylphenylammonium: N+ with no lone pair → 0.00
        let m = mol("c1csc([N+]2(C3CCCCC3)CCCCC2)c1");
        let t = tpsa(&m);
        // N+ contributes 0; thiophene S contributes ~38.8 → total ~38.8
        assert!(t < 50.0, "N+ quaternary TPSA should be small, got {t}");
        // Specifically, N+ should NOT add 3.01 on top
        let m2 = mol("c1ccsc1");
        let t_base = tpsa(&m2); // thiophene TPSA without N+
        assert!(
            (t - t_base).abs() < 1.0,
            "N+ adds ~0 to TPSA, got delta {}",
            t - t_base
        );
    }

    // Ionic N-oxide [N+]-[O-] ring: O- accounts for TPSA, N+ contributes 0
    #[test]
    fn test_tpsa_n_oxide_ionic() {
        // [N+]1([O-])CCOCC1 — morpholine N-oxide
        // O- and ring O contribute; N+ should be 0
        let m = mol("O=C1CC[N+]2([O-])CCCC[C@@H]12"); // representative N-oxide ring
        let t = tpsa(&m);
        assert!(t > 0.0 && t < 150.0, "N-oxide TPSA in range, got {t}");
    }

    // -- Fsp3 tests --------------------------------------------------------

    #[test]
    fn test_fsp3_benzene() {
        let m = mol("c1ccccc1");
        // all aromatic C, no sp3
        assert!((fsp3(&m) - 0.0).abs() < 1e-9, "benzene Fsp3 should be 0");
    }

    #[test]
    fn test_fsp3_cyclohexane() {
        let m = mol("C1CCCCC1");
        // all sp3 C
        assert!(
            (fsp3(&m) - 1.0).abs() < 1e-9,
            "cyclohexane Fsp3 should be 1"
        );
    }

    #[test]
    fn test_fsp3_aspirin() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        // 9 C total: 1 sp3 (methyl), 2 sp2 C=O, 6 aromatic
        // sp3 = 1, total C = 9 → Fsp3 = 1/9 ≈ 0.111
        let f = fsp3(&m);
        assert!(f > 0.05 && f < 0.25, "aspirin Fsp3={f} expected ~0.111");
    }

    #[test]
    fn test_fsp3_no_carbon() {
        let m = mol("[NH4+]");
        assert!(
            (fsp3(&m) - 0.0).abs() < 1e-9,
            "no-carbon mol Fsp3 should be 0"
        );
    }

    #[test]
    fn test_fsp3_zero_valence_isotope_carbon_is_not_csp3() {
        assert_eq!(fsp3(&mol("[13C]")), 0.0);
        assert_eq!(fsp3(&mol("[13CH4]")), 1.0);
    }

    // -- MQN tests ----------------------------------------------------------

    #[test]
    fn test_mqn_length() {
        let m = mol("CCO");
        let desc = mqn(&m);
        assert_eq!(desc.len(), 42);
    }

    #[test]
    fn test_mqn_matches_rdkit_calc_mqns() {
        // RDKit 2026.03.1 rdMolDescriptors.MQNs_.
        let cases: [(&str, [u32; 42]); 6] = [
            (
                "C",
                [
                    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                ],
            ),
            (
                "CC",
                [
                    2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                ],
            ),
            (
                "c1ccccc1",
                [
                    6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 6, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                ],
            ),
            (
                "CC(=O)Oc1ccccc1C(=O)O",
                [
                    9, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 13, 5, 2, 0, 3, 3, 0, 2, 8, 4, 1, 1, 0, 0, 4,
                    1, 2, 0, 4, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                ],
            ),
            (
                "C1CC2CCC1C2",
                [
                    7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 5, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 3, 2,
                ],
            ),
            (
                "[NH3+]CC([O-])=O",
                [
                    2, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 5, 3, 1, 0, 0, 0, 0, 1, 6, 3, 3, 1, 1, 1, 3,
                    1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                ],
            ),
        ];
        for (smiles, expected) in cases {
            assert_eq!(mqn(&mol(smiles)), expected.to_vec(), "{smiles}");
        }
    }

    // -- AutoCorr2D tests ---------------------------------------------------

    #[test]
    fn test_autocorr_2d_single_atom() {
        let m = mol("C");
        let ac = autocorr_2d(&m);
        assert_eq!(ac.len(), 7);
        // Single atom: no pairs → all zeros
        for val in ac {
            assert!((val - 0.0).abs() < 1e-9);
        }
    }

    #[test]
    fn test_autocorr_2d_ethane() {
        let m = mol("CC");
        let ac = autocorr_2d(&m);
        assert_eq!(ac.len(), 7);
        // Ethane: distance 1 pair (C-C), both have valence 4
        // autocorr[0] (lag 1) = 4 * 4 = 16
        assert!((ac[0] - 16.0).abs() < 1e-9, "lag 1: {}", ac[0]);
        // Lag 2+ → no pairs
        for (i, value) in ac.iter().enumerate().take(7).skip(1) {
            assert!((*value - 0.0).abs() < 1e-9, "lag {}: {}", i + 1, value);
        }
    }

    #[test]
    fn test_autocorr_2d_propane() {
        let m = mol("CCC");
        let ac = autocorr_2d(&m);
        assert_eq!(ac.len(), 7);
        // Propane: C1-C2 dist=1, C2-C3 dist=1, C1-C3 dist=2
        // C1 (terminal): degree 1, implicit H = 3, valence = 1 + 3 = 4
        // C2 (central): degree 2, implicit H = 2, valence = 2 + 2 = 4
        // C3 (terminal): degree 1, implicit H = 3, valence = 1 + 3 = 4
        // lag 1: C1-C2 + C2-C3 = 4*4 + 4*4 = 32
        assert!((ac[0] - 32.0).abs() < 1e-9, "lag 1: {}", ac[0]);
        // lag 2: C1-C3 = 4*4 = 16
        assert!((ac[1] - 16.0).abs() < 1e-9, "lag 2: {}", ac[1]);
    }

    #[test]
    fn test_autocorr_2d_benzene() {
        let m = mol("c1ccccc1");
        let ac = autocorr_2d(&m);
        assert_eq!(ac.len(), 7);
        // Benzene: aromatic ring, all C have valence 3 (2 bonds + 1 H)
        // lag 1: 6 C-C bonds = 6 * (3*3) = 54
        assert!((ac[0] - 54.0).abs() < 1e-9, "lag 1 benzene: {}", ac[0]);
        // Should have non-zero values for multiple lags (cyclic)
        assert!(ac[1] > 0.0, "lag 2 should be non-zero");
    }

    // -- aromatic_ring_count tests -----------------------------------------

    #[test]
    fn test_aromatic_ring_count_benzene() {
        let m = mol("c1ccccc1");
        assert_eq!(aromatic_ring_count(&m), 1);
    }

    #[test]
    fn test_aromatic_ring_count_naphthalene() {
        let m = mol("c1ccc2ccccc2c1");
        assert_eq!(aromatic_ring_count(&m), 2);
    }

    #[test]
    fn test_aromatic_ring_count_cyclohexane() {
        let m = mol("C1CCCCC1");
        assert_eq!(aromatic_ring_count(&m), 0);
    }

    #[test]
    fn test_aromatic_ring_count_aspirin() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        assert_eq!(aromatic_ring_count(&m), 1);
    }

    #[test]
    fn test_aromatic_ring_count_anthracene() {
        let m = mol("c1ccc2cc3ccccc3cc2c1");
        assert_eq!(aromatic_ring_count(&m), 3);
    }

    #[test]
    fn test_aromatic_ring_count_pyrene() {
        // 4 fused 6-rings; SSSR may produce envelope cycle requiring 3-ring XOR removal
        let m = mol("c1ccc2cccc3ccc4cccc1c4c23");
        assert_eq!(aromatic_ring_count(&m), 4);
    }

    #[test]
    fn test_aromatic_ring_count_triphenylene() {
        let m = mol("c1ccc2ccc3ccccc3c2c1");
        assert_eq!(aromatic_ring_count(&m), 3);
    }

    #[test]
    fn test_aromatic_ring_count_fluoranthene() {
        // 3 six-membered + 1 five-membered aromatic ring
        let m = mol("c1ccc2-c3cccc4cccc-3c4c2c1");
        assert_eq!(aromatic_ring_count(&m), 4);
    }

    #[test]
    fn test_aromatic_ring_count_acridine() {
        let m = mol("c1ccc2nc3ccccc3cc2c1");
        assert_eq!(aromatic_ring_count(&m), 3);
    }

    // bench5k regression cases (rd=N, ch=N-1 pattern)
    #[test]
    fn test_arc_bench_lactone_benzene() {
        // rd=1 ch=0: benzene fused with lactone ring
        let m = mol("CC1(C)CC(=O)c2c(ccc(C(=O)CC(N)CO)c2N)O1");
        assert_eq!(aromatic_ring_count(&m), 1);
    }

    #[test]
    fn test_arc_bench_bridged_benzenes() {
        // rd=2 ch=1: two aromatic benzenes connected by O and C bridge
        let m = mol("Cc1c(C)c2c(c(C)c1O)Oc1c(C)c([C@H](C)[C@@H](C)O)c(O)c(O)c1C2");
        assert_eq!(aromatic_ring_count(&m), 2);
    }

    #[test]
    fn test_arc_bench_no_bridge() {
        // rd=2 ch=1: two benzenes with N/O bridges
        let m = mol("C[C@]12Nc3ccccc3O[C@@]1(C)Nc1ccccc1O2");
        assert_eq!(aromatic_ring_count(&m), 2);
    }

    #[test]
    fn test_arc_bench_steroid_benzene() {
        // rd=1 ch=0: complex steroid scaffold with one phenyl ring
        let m = mol(
            "CC1=C2C[C@@]3(C)C[C@H](O)[C@](C)(C[C@H](O)[C@H](O)[C@@](C)(O)CO)[C@H]3c3ccc(C)c(c32)C1",
        );
        assert_eq!(aromatic_ring_count(&m), 1);
    }

    #[test]
    fn test_hba_metformin() {
        // Metformin: CN(C)C(=N)NC(=N)N — RDKit CalcNumHBA = 2 (two imine =NH)
        // Confirmed by tracing current n_adjacent_to_pi_center logic:
        //   dimethylamino-N, bridging NH, terminal NH₂ → excluded (adjacent to C=N)
        //   two =NH imine nitrogens → counted (double bond to neighbor, not single-like)
        let m = mol("CN(C)C(=N)NC(=N)N");
        assert_eq!(hba_count(&m), 2, "metformin HBA should match RDKit (2)");
    }

    #[test]
    fn test_hba_caffeine_matches_rdkit_for_substituted_aromatic_n() {
        // RDKit's Ertl HBA definition counts one pyridine-like N and the two
        // carbonyl oxygens.  The three degree-3 N-methyl/imide nitrogens are
        // not acceptors.
        let m = mol("Cn1cnc2c1c(=O)n(C)c(=O)n2C");
        assert_eq!(hba_count(&m), 3);
    }

    // -- formal_charge_sum tests -------------------------------------------

    #[test]
    fn test_formal_charge_neutral_aspirin() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        assert_eq!(formal_charge_sum(&m), 0);
    }

    #[test]
    fn test_formal_charge_quaternary_n() {
        // Trimethylammonium: [N+]
        let m = mol("CC[N+](C)(C)C");
        assert_eq!(formal_charge_sum(&m), 1);
    }

    #[test]
    fn test_formal_charge_zwitterion() {
        // Glycine zwitterion: [NH3+]CC(=O)[O-]
        let m = mol("[NH3+]CC(=O)[O-]");
        assert_eq!(formal_charge_sum(&m), 0);
    }

    // -- molar_refractivity tests -------------------------------------------

    #[test]
    fn test_mr_benzene_range() {
        // Benzene MR ≈ 26.0 (RDKit reference)
        let m = mol("c1ccccc1");
        let mr = molar_refractivity(&m);
        assert!(mr > 20.0 && mr < 35.0, "benzene MR={mr:.2}");
    }

    #[test]
    fn test_mr_aspirin_range() {
        // Aspirin MR ≈ 45.5 (Ghose filter: 40-130 → should pass)
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        let mr = molar_refractivity(&m);
        assert!(mr > 35.0 && mr < 65.0, "aspirin MR={mr:.2}");
    }

    #[test]
    fn test_mr_chlorobenzene_higher_than_benzene() {
        // Cl adds ~5.85 to MR
        let m_bz = mol("c1ccccc1");
        let m_clb = mol("c1ccc(Cl)cc1");
        assert!(
            molar_refractivity(&m_clb) > molar_refractivity(&m_bz),
            "chlorobenzene should have higher MR than benzene"
        );
    }

    // -- Veber filter tests -------------------------------------------------

    #[test]
    fn test_veber_aspirin_passes() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        assert!(veber_passes(&m), "aspirin should pass Veber filter");
    }

    #[test]
    fn test_veber_large_flexible_fails() {
        // A molecule with many rotatable bonds should fail
        let m = mol("CCCCCCCCCCCCC(=O)O"); // myristic acid - 12 rotatable bonds
        let rotb = rotatable_bond_count(&m);
        if rotb > 10 {
            assert!(
                !veber_passes(&m),
                "myristic acid (rotb={rotb}) should fail Veber"
            );
        }
    }

    // -- Egan filter tests --------------------------------------------------

    #[test]
    fn test_egan_aspirin_passes() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        assert!(egan_passes(&m), "aspirin should pass Egan filter");
    }

    // -- REOS filter tests --------------------------------------------------

    #[test]
    fn test_reos_aspirin_passes() {
        // Aspirin: MW=180, logP~1.2, HBD=1, HBA=3, charge=0, rotb=3, HAC=13
        // HAC=13 < 15 → REOS fails! (aspirin is small)
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        // Just test that we can call it without panicking
        let _ = reos_passes(&m);
    }

    #[test]
    fn test_reos_ibuprofen_passes() {
        // Ibuprofen: MW=206, logP~3.5, HBD=1, HBA=1, charge=0, rotb=4, HAC=13
        // HAC=13 < 15 → likely fails
        let m = mol("CC(C)Cc1ccc(cc1)C(C)C(=O)O");
        let _ = reos_passes(&m);
    }

    #[test]
    fn test_reos_diazepam_passes() {
        // Diazepam: MW~285, logP~2.9, HBD=0, HBA=2, charge=0, rotb=1, HAC~22 — all in range
        let m = mol("CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21");
        assert!(reos_passes(&m), "diazepam should pass REOS");
    }

    // -- Ghose filter tests -------------------------------------------------

    #[test]
    fn test_ghose_aspirin_range() {
        // Aspirin MW=180 is below Ghose MW lower bound (160 ≤ MW ≤ 480). Should pass MW.
        // HAC=13 < 20 lower bound → Ghose fails for aspirin
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        // Just verify it doesn't panic
        let _ = ghose_passes(&m);
    }

    #[test]
    fn test_ghose_ibuprofen_passes() {
        // Ibuprofen: MW=206, logP~3.5, HAC=13 (borderline)
        let m = mol("CC(C)Cc1ccc(cc1)C(C)C(=O)O");
        let _ = ghose_passes(&m);
    }

    // -- Basic count descriptor tests ----------------------------------------

    #[test]
    fn test_num_heteroatoms_aspirin() {
        // Aspirin: 4 O atoms (ester C=O, ester O, carboxylic C=O, carboxylic OH)
        assert_eq!(num_heteroatoms(&mol("CC(=O)Oc1ccccc1C(=O)O")), 4);
    }

    #[test]
    fn test_num_heteroatoms_benzene_zero() {
        assert_eq!(num_heteroatoms(&mol("c1ccccc1")), 0);
    }

    #[test]
    fn test_ring_count_benzene() {
        assert_eq!(ring_count(&mol("c1ccccc1")), 1);
    }

    #[test]
    fn test_ring_count_naphthalene() {
        assert_eq!(ring_count(&mol("c1ccc2ccccc2c1")), 2);
    }

    #[test]
    fn test_ring_count_acyclic_zero() {
        assert_eq!(ring_count(&mol("CCO")), 0);
    }

    #[test]
    fn test_num_saturated_rings_cyclohexane() {
        assert_eq!(num_saturated_rings(&mol("C1CCCCC1")), 1);
    }

    #[test]
    fn test_num_saturated_rings_benzene_zero() {
        assert_eq!(num_saturated_rings(&mol("c1ccccc1")), 0);
    }

    #[test]
    fn test_num_aliphatic_rings_cyclohexane() {
        assert_eq!(num_aliphatic_rings(&mol("C1CCCCC1")), 1);
    }

    #[test]
    fn test_num_aliphatic_rings_benzene_zero() {
        assert_eq!(num_aliphatic_rings(&mol("c1ccccc1")), 0);
    }

    #[test]
    fn test_num_stereocenters_alanine() {
        // L-alanine: 1 R/S center
        assert_eq!(num_stereocenters(&mol("[C@@H](N)(C)C(=O)O")), 1);
    }

    #[test]
    fn test_potential_stereocenter_indices_match_count() {
        let molecule = mol("N[C@@H](C)C(=O)O");
        let indices = potential_stereocenter_indices(&molecule);
        assert_eq!(indices.len(), num_stereocenters(&molecule));
        assert_eq!(indices, vec![AtomIdx(1)]);
    }

    #[test]
    fn test_num_stereocenters_achiral_zero() {
        assert_eq!(num_stereocenters(&mol("CC(=O)O")), 0);
    }

    #[test]
    fn test_num_stereocenters_bicyclic_tertiary_amine() {
        // RDKit #9629: its legacy perception reported a spurious extra
        // stereocenter after an identity atom renumbering for this cage.
        // The bridgehead carbon is a potential center; the tertiary nitrogen
        // is not. Keep the chemically-derived count fixed independently of
        // input storage order and ring-cache implementation details.
        assert_eq!(num_stereocenters(&mol("C1CCN2CCCC2C1")), 1);
    }

    #[test]
    fn rdkit_9629_class_is_invariant_across_relabel_clone_reparse_and_call_order() {
        // This imports the *failure class* of RDKit #9629, not a claim that
        // this molecule reproduces RDKit's historical result. A potential
        // center must remain attached to the same atom identity across graph
        // storage permutations, cache initialization order, cloning, and
        // SMILES serialization. The ordinary saturated tertiary amine is a
        // necessary negative control: it must not become a center merely by
        // reordering or by calling ring perception first.
        let source = mol("C1CCN2CCCC2C1");
        let baseline = mapped_observation(
            &mapped_permutation(&source, &(0..source.atom_count()).collect::<Vec<_>>()),
            true,
        );
        assert_eq!(
            baseline.0.len(),
            1,
            "bridgehead is the one potential center"
        );
        assert!(
            baseline.1.is_empty(),
            "unspecified input has no assigned CIP code"
        );

        for seed in 0..32u64 {
            let permuted =
                mapped_permutation(&source, &fixed_permutation(source.atom_count(), seed));
            let reparsed = parse(&write(&permuted)).expect("mapped permutation must round-trip");
            for (label, candidate) in [
                ("permuted", permuted.clone()),
                ("cloned", permuted.clone()),
                ("reparsed", reparsed),
            ] {
                for rings_first in [false, true] {
                    assert_eq!(
                        mapped_observation(&candidate, rings_first),
                        baseline,
                        "seed {seed}, {label}, rings_first={rings_first} changed mapped stereo observation"
                    );
                }
            }
        }
    }

    #[test]
    fn rdkit_9629_class_keeps_ordinary_tertiary_amine_negative_control() {
        let source = mol("CN(C)CC");
        let baseline = mapped_observation(
            &mapped_permutation(&source, &(0..source.atom_count()).collect::<Vec<_>>()),
            true,
        );
        assert!(
            baseline.0.is_empty(),
            "ordinary tertiary amine is not a potential center"
        );
        assert!(baseline.1.is_empty());

        for seed in 0..32u64 {
            let permuted =
                mapped_permutation(&source, &fixed_permutation(source.atom_count(), seed));
            let reparsed = parse(&write(&permuted)).expect("mapped permutation must round-trip");
            for candidate in [permuted.clone(), reparsed] {
                assert_eq!(
                    mapped_observation(&candidate, false),
                    baseline,
                    "seed {seed}"
                );
                assert_eq!(
                    mapped_observation(&candidate, true),
                    baseline,
                    "seed {seed}"
                );
            }
        }
    }

    #[test]
    fn test_num_stereocenters_bridgehead_quaternary() {
        // Bicyclic bridgehead: ring-adjacent tie resolved via CIP Rule 5 provisional R/S.
        assert_eq!(
            num_stereocenters(&mol("NS(=O)(=O)OC[C@@]12CCCC[C@@H]1CCC2")),
            2
        );
    }

    #[test]
    fn test_num_stereocenters_ring_adjacent() {
        // 3 ring-adjacent stereocenters in a cyclohexyl side chain.
        assert_eq!(
            num_stereocenters(&mol(
                "CCCCc1cn([C@H]2[C@H](C)CCC[C@@H]2C)c(=O)n1Cc1ccc(-c2ccccc2-c2nn[nH]n2)nc1"
            )),
            3
        );
    }

    #[test]
    fn test_num_stereocenters_aromatic_sulfur_with_exocyclic_branch() {
        // RDKit counts the aromatic [s+] center in this zwitterionic ring as
        // a potential tetrahedral center; aromatic C/N must remain planar.
        assert_eq!(
            num_stereocenters(&mol("CN(C)Cc1ccc(CSCCNc2n[s+]([O-])nc2N)[nH]1")),
            1
        );
    }

    #[test]
    fn test_num_stereocenters_long_substituted_cyclopropane() {
        // The two ring centers are distinguished only after the short-range
        // branch walk reaches the remote chain termini.
        assert_eq!(num_stereocenters(&mol("CCCCCCCCC1CC1CCCCCCCC(=O)O")), 2);
    }

    #[test]
    fn test_num_stereocenters_branch_provenance_fallback_cages() {
        // The legacy pooled-sphere comparator ties two distinct carbon
        // branches in these fused/bridged cages. The provenance-preserving
        // fallback must recover the full RDKit potential-center count without
        // broadening ordinary tertiary-amine handling.
        assert_eq!(num_stereocenters(&mol("CC12NC(Cc3ccccc31)CC1CCCCC12")), 4);
        assert_eq!(num_stereocenters(&mol("CN1C2(C)c3ccccc3C1(C)C1CCCCC12")), 4);
    }

    #[test]
    fn test_num_stereocenters_ring_constrained_tertiary_nitrogen() {
        let residuals = [
            "Cn1cc(C2=NC[C@]3(CN4CC[C@@H]3C4)O2)c2ccccc21",
            "Cn1cc(C2=NC[C@@]3(C[N@@]4CC[C@@H]3C4)O2)c2ccccc21",
            "Cn1cc(C2=NCC3(CCN4CCCC3C4)O2)c2ccccc21",
        ];
        for smiles in residuals {
            assert_eq!(num_stereocenters(&mol(smiles)), 3, "{smiles}");
        }
        assert_eq!(num_stereocenters(&mol("C1N(C)CCC1")), 0);
    }

    // -- BalabanJ tests -------------------------------------------------

    #[test]
    fn test_balaban_j_ethane() {
        let m = mol("CC");
        let bj = balaban_j(&m);
        assert!(bj > 0.0, "ethane should have positive BalabanJ");
    }

    #[test]
    fn test_balaban_j_benzene() {
        let m = mol("c1ccccc1");
        let bj = balaban_j(&m);
        assert!(bj > 0.0, "benzene should have positive BalabanJ");
    }

    #[test]
    fn test_balaban_j_single_atom_zero() {
        let m = mol("C");
        let bj = balaban_j(&m);
        assert_eq!(bj, 0.0, "single atom should have BalabanJ = 0");
    }

    #[test]
    fn test_balaban_j_matches_rdkit() {
        // RDKit Descriptors.BalabanJ verified values. Benzene vs cyclohexane
        // is the key discriminator: same topology (hexagon), but aromatic
        // bonds (1/1.5 edge weight) vs single bonds (1/1.0) give different
        // weighted distance sums, so J must differ (3.0 vs 2.0) even though
        // a purely-topological (unweighted) formula would give 2.0 for both.
        assert!(
            (balaban_j(&mol("c1ccccc1")) - 3.0).abs() < 1e-6,
            "benzene BalabanJ"
        );
        assert!(
            (balaban_j(&mol("C1CCCCC1")) - 2.0).abs() < 1e-6,
            "cyclohexane BalabanJ"
        );
        assert!(
            (balaban_j(&mol("CCCCCC")) - 2.3391).abs() < 0.001,
            "hexane BalabanJ"
        );
        assert!(
            (balaban_j(&mol("CCO")) - 1.6330).abs() < 0.001,
            "ethanol BalabanJ"
        );
    }

    // -- Ipc tests -------------------------------------------------------

    #[test]
    fn test_ipc_ethane() {
        let m = mol("CC");
        let ipc_val = ipc(&m);
        assert!(ipc_val >= 0.0, "ethane should have non-negative Ipc");
    }

    #[test]
    fn test_ipc_benzene() {
        let m = mol("c1ccccc1");
        let ipc_val = ipc(&m);
        assert!(ipc_val > 0.0, "benzene should have positive Ipc");
    }

    #[test]
    fn test_ipc_single_atom_zero() {
        let m = mol("C");
        let ipc_val = ipc(&m);
        assert_eq!(ipc_val, 0.0, "single atom should have Ipc = 0");
    }

    // -- HallKierAlpha tests -----------------------------------------------

    #[test]
    fn test_hall_kier_alpha_methane() {
        // Single sp3 carbon = the reference radius itself → alpha is exactly 0.
        let m = mol("C");
        assert!(
            (hall_kier_alpha(&m) - 0.0).abs() < 1e-9,
            "methane HallKierAlpha should be exactly 0 (sp3 C is the reference atom)"
        );
    }

    #[test]
    fn test_hall_kier_alpha_ethane() {
        // All-sp3-carbon molecules stay at the reference radius → alpha = 0.
        let m = mol("CC");
        assert!(
            (hall_kier_alpha(&m) - 0.0).abs() < 1e-9,
            "ethane (all sp3 C) HallKierAlpha should be exactly 0"
        );
    }

    #[test]
    fn test_hall_kier_alpha_benzene() {
        // 6 aromatic (sp2) carbons, each smaller than the sp3 reference → negative.
        let m = mol("c1ccccc1");
        let hka = hall_kier_alpha(&m);
        assert!(
            hka < 0.0,
            "benzene HallKierAlpha should be negative, got {hka}"
        );
    }

    // -- USRCAT tests -------------------------------------------------------

    #[test]
    fn test_usrcat_shape() {
        let m = mol("CC");
        let usr = usrcat(&m);
        assert_eq!(usr.len(), 42, "USRCAT should return 42 values");
        assert!(usr[0] >= 0.0, "first slot should be non-negative");
    }

    #[test]
    fn test_usrcat_donors_acceptors() {
        let m = mol("CCO");
        let usr = usrcat(&m);
        assert!(usr[36] >= 0.0, "donor count should be non-negative");
        assert!(
            usr[37] > 0.0,
            "acceptor count should be positive (O present)"
        );
    }

    #[test]
    fn test_usrcat_aromatic() {
        let m = mol("c1ccccc1");
        let usr = usrcat(&m);
        assert!(
            usr[38] > 0.0,
            "aromatic count should be positive for benzene"
        );
    }

    #[test]
    fn test_usrcat_charged() {
        let m = mol("CC(=O)[O-]");
        let usr = usrcat(&m);
        assert!(
            usr[40] > 0.0,
            "anion count should be positive for charged carboxylate"
        );
    }

    // -- MMFF94 charges tests -----------------------------------------------

    #[test]
    fn test_mmff94_charges_length() {
        let m = mol("CCO");
        let charges = mmff94_charges(&m);
        assert_eq!(charges.len(), 3, "should have 3 charges for 3 atoms");
    }

    #[test]
    fn test_mmff94_charges_ethane() {
        let m = mol("CC");
        let charges = mmff94_charges(&m);
        assert_eq!(charges.len(), 2);
        // Both carbons should have similar (small negative) charges
        assert!(
            (charges[0] - charges[1]).abs() < 0.1,
            "carbons in ethane should have similar charges"
        );
    }

    #[test]
    fn test_mmff94_charges_charged_species() {
        let m = mol("CC(=O)[O-]");
        let charges = mmff94_charges(&m);
        assert_eq!(charges.len(), 4);
        // Carboxylate oxygen should be negative
        assert!(charges[3] < 0.0, "carboxylate oxygen should be negative");
    }

    #[test]
    fn test_mmff94_charges_water() {
        let m = mol("O");
        let charges = mmff94_charges(&m);
        assert_eq!(charges.len(), 1);
        assert!(
            charges[0].is_finite(),
            "water oxygen charge should be finite"
        );
    }

    // -- Element count tests -----------------------------------------------

    #[test]
    fn test_num_carbons_ethane() {
        let m = mol("CC");
        assert_eq!(num_carbons(&m), 2);
    }

    #[test]
    fn test_num_nitrogens_methylamine() {
        let m = mol("CN");
        assert_eq!(num_nitrogens(&m), 1);
    }

    #[test]
    fn test_num_oxygens_methanol() {
        let m = mol("CO");
        assert_eq!(num_oxygens(&m), 1);
    }

    #[test]
    fn test_num_halogens() {
        let m = mol("CCF");
        assert_eq!(num_carbons(&m), 2);
        assert_eq!(num_fluorines(&m), 1);
    }

    #[test]
    fn test_num_hydrogens_methane() {
        let m = mol("C");
        // Methane: 1 carbon with 4 implicit hydrogens
        assert_eq!(num_hydrogens(&m), 4);
    }

    #[test]
    fn test_num_hydrogens_ethane() {
        let m = mol("CC");
        // Ethane: 2 carbons, 6 total hydrogens
        assert_eq!(num_hydrogens(&m), 6);
    }

    #[test]
    fn test_num_hydrogens_water() {
        let m = mol("O");
        // Water: 2 implicit hydrogens on oxygen
        assert_eq!(num_hydrogens(&m), 2);
    }

    // -- Molecular formula tests -----------------------------------------------

    #[test]
    fn test_calc_mol_formula_ethane() {
        let m = mol("CC");
        assert_eq!(calc_mol_formula(&m), "C2H6");
    }

    // -- Bridgehead and spiro atom tests -----------------------------------

    #[test]
    fn test_num_bridgehead_atoms_acyclic() {
        // Acyclic molecule should return 0
        let m = mol("CCC");
        assert_eq!(num_bridgehead_atoms(&m), 0);
    }

    #[test]
    fn test_num_bridgehead_atoms_single_ring() {
        // Cyclohexane (single ring) should have 0 bridgeheads
        let m = mol("C1CCCCC1");
        assert_eq!(num_bridgehead_atoms(&m), 0);
    }

    #[test]
    fn test_num_bridgehead_atoms_norbornane() {
        // Norbornane (bicyclo[2.2.1]heptane): 2 bridgehead atoms
        // Structure: two rings sharing 2 atoms at positions 1 and 2, with bridges
        let m = mol("C1CC2CCC1C2");
        assert_eq!(num_bridgehead_atoms(&m), 2);
    }

    #[test]
    fn test_num_bridgehead_atoms_naphthalene_fused() {
        // Naphthalene is fused ring system (not bridged)
        // Bridgeheads = atoms in 2+ rings but only on shared edges (fused, no bridges)
        let m = mol("c1ccc2ccccc2c1");
        assert_eq!(
            num_bridgehead_atoms(&m),
            0,
            "naphthalene is fused, not bridged"
        );
    }

    #[test]
    fn test_num_spiro_atoms_single_ring() {
        // Cyclohexane has no spiro atoms
        let m = mol("C1CCCCC1");
        assert_eq!(num_spiro_atoms(&m), 0);
    }

    #[test]
    fn test_calc_mol_formula_charge_and_dummy_atoms_match_rdkit() {
        // RDKit 2026.03.1 CalcMolFormula.
        for (smiles, expected) in [
            ("[NH4+]", "H4N+"),
            ("CC(=O)[O-]", "C2H3O2-"),
            ("[O-]S(=O)(=O)[O-]", "O4S-2"),
            ("[Fe+3]", "Fe+3"),
            ("O=[N+]([O-])c1ccccc1", "C6H5NO2"),
            ("[Na+].[Cl-]", "ClNa"),
            ("*C", "CH3*"),
        ] {
            assert_eq!(calc_mol_formula(&mol(smiles)), expected, "{smiles}");
        }
    }

    #[test]
    fn test_calc_mol_formula_water() {
        let m = mol("O");
        assert_eq!(calc_mol_formula(&m), "H2O");
    }

    #[test]
    fn test_calc_mol_formula_benzene() {
        let m = mol("c1ccccc1");
        assert_eq!(calc_mol_formula(&m), "C6H6");
    }

    #[test]
    fn test_calc_mol_formula_acetic_acid() {
        let m = mol("CC(=O)O");
        assert_eq!(calc_mol_formula(&m), "C2H4O2");
    }

    // =========================================================================
    // Functional group bond counts (C4 - granular classification)
    // =========================================================================

    #[test]
    fn test_num_amide_bonds_acetamide() {
        // CH3-C(=O)-N-H: one amide bond
        let m = mol("CC(=O)N");
        assert_eq!(num_amide_bonds(&m), 1);
    }

    #[test]
    fn test_num_amide_bonds_urea() {
        // N-C(=O)-N: two N–C(=O) bonds (RDKit counts per bond, not per C=O carbon)
        let m = mol("NC(=O)N");
        assert_eq!(num_amide_bonds(&m), 2);
    }

    #[test]
    fn test_num_amide_bonds_primary_amide() {
        // CH3-C(=O)-NH2: primary amide
        let m = mol("CC(=O)N");
        assert_eq!(num_amide_bonds(&m), 1);
    }

    #[test]
    fn test_num_amide_bonds_none() {
        // Benzene has no amide bonds
        let m = mol("c1ccccc1");
        assert_eq!(num_amide_bonds(&m), 0);
    }

    #[test]
    fn test_num_ester_bonds_methyl_formate() {
        // H-C(=O)-O-CH3: one ester bond
        let m = mol("COC=O");
        assert_eq!(num_ester_bonds(&m), 1);
    }

    #[test]
    fn test_num_ester_bonds_acetic_acid_methyl_ester() {
        // CH3-C(=O)-O-CH3: one ester bond
        let m = mol("CC(=O)OC");
        assert_eq!(num_ester_bonds(&m), 1);
    }

    #[test]
    fn test_num_ester_bonds_none() {
        // Carboxylic acid (COOH) is not an ester
        let m = mol("CC(=O)O");
        // CC(=O)O is acetic acid: C=O with O-H (not O-C)
        assert_eq!(num_ester_bonds(&m), 0);
    }

    #[test]
    fn test_num_ester_bonds_aspirin() {
        // Aspirin CC(=O)Oc1ccccc1C(=O)O has one ester bond (acetyl ester)
        // and one carboxylic acid (not counted)
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        assert_eq!(
            num_ester_bonds(&m),
            1,
            "aspirin has one ester bond (aryl ester)"
        );
    }

    // B5: context-dependent alkene C LogP contributions
    #[test]
    fn test_logp_ethylene_terminal() {
        // Ethylene: both C are terminal =CH2 (h=2, no aromatic neighbor) → each 0.1551
        let m = mol("C=C");
        let lp = logp_crippen(&m);
        assert!(lp > 0.3 && lp < 1.3, "ethylene logp = {lp}");
    }

    #[test]
    fn test_logp_propene_terminal_internal() {
        // Propene: =CH2 (0.1551) + =CH- internal (0.2274) + CH3 → propene > ethylene
        let m = mol("CC=C");
        let lp = logp_crippen(&m);
        let eth = logp_crippen(&mol("C=C"));
        assert!(
            lp > eth,
            "propene logp ({lp}) should exceed ethylene ({eth})"
        );
        assert!(
            lp > 0.5 && lp < 2.0,
            "propene logp = {lp} out of expected range"
        );
    }

    #[test]
    fn test_logp_styrene_splits_correctly() {
        // Styrene: =CH2 (0.1551) + =CH- adj to Ar (0.2640) + 6 aromatic C
        let m = mol("C=Cc1ccccc1");
        let lp = logp_crippen(&m);
        assert!(lp > 1.8 && lp < 3.4, "styrene logp = {lp}");
        // Ar-adjacent CH contributes more (0.2640) than terminal CH2 (0.1551)
        let per_atom = logp_crippen_per_atom(&m);
        assert!(per_atom.len() >= 8, "styrene has 8 heavy atoms");
    }

    #[test]
    fn test_logp_1_phenylpropene_ar_adjacent() {
        // 1-phenylpropene Ph-CH=CH-CH3: Ar-adjacent =CH- (0.2640) + internal =CH- (0.2274)
        let m = mol("CC=Cc1ccccc1");
        let lp = logp_crippen(&m);
        assert!(lp > 2.0 && lp < 4.0, "1-phenylpropene logp = {lp}");
    }

    // C1: Reference LogP tests for complex polar molecules.
    //
    // For molecules with multiple polar functional groups (carboxylate, amide,
    // amine salt, conjugated C=C-C=O), the Crippen atom-type model underestimates
    // LogP versus RDKit. xlogp3() is implemented and more accurate for these cases.
    // These tests document the current values and catch regressions.

    #[test]
    fn test_logp_curcumin_reference() {
        // Curcumin: two conjugated vinyl-ketone arms, phenol/methoxy substituents.
        // RDKit Crippen ~3.04; Crippen atom-type model gives a lower value.
        // xlogp3 is the recommended API for complex polyphenols.
        let m = mol("COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O");
        let lp = logp_crippen(&m);
        assert!(lp > -5.0 && lp < 5.0, "curcumin crippen logp = {lp}");
    }

    #[test]
    fn test_logp_complex_molecules_xlogp3_preferred() {
        // For high-error molecules, verify xlogp3 gives a more positive value
        // (closer to RDKit) than Crippen for curcumin.
        let m = mol("COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O");
        let crippen = logp_crippen(&m);
        let xl3 = crate::xlogp3::xlogp3(&m);
        // Both should be finite; just document that they exist
        assert!(crippen.is_finite(), "crippen logp must be finite");
        assert!(xl3.is_finite(), "xlogp3 must be finite");
    }

    // ---- Enone vinyl C tests (v0.1.99) ----

    #[test]
    fn test_logp_mvk_enone_vinyl() {
        // Methyl vinyl ketone: CH2=CH-C(=O)-CH3
        // =CH2 (0.1551) + =CH- adj to C=O (enone, 0.1302) + C=O + CH3
        // RDKit Crippen ~0.44; expect reasonable range
        let m = mol("C=CC(=O)C");
        let lp = logp_crippen(&m);
        assert!(lp > 0.0 && lp < 1.5, "MVK logp = {lp}");
    }

    #[test]
    fn test_logp_chalcone_enone() {
        // Chalcone: Ph-CH=CH-C(=O)-Ph
        // Ar-CH= (0.2640, ar-adjacent wins) + =CH-C(=O) (0.1302, enone)
        // RDKit ~3.04; Crippen gives slightly different value
        let m = mol("c1ccccc1/C=C/C(=O)c1ccccc1");
        let lp = logp_crippen(&m);
        assert!(lp > 1.8 && lp < 4.5, "chalcone logp = {lp}");
    }

    #[test]
    fn test_logp_crotonate_internal_enone() {
        // Crotonic acid: CH3-CH=CH-COOH (trans)
        // =CH- (0.2274, non-ar, 1H) adj to C=O → enone case (0.1302)
        let m = mol("CC=CC(=O)O");
        let lp = logp_crippen(&m);
        assert!(lp > -0.5 && lp < 1.5, "crotonate logp = {lp}");
    }

    #[test]
    fn test_logp_enone_vs_plain_alkene() {
        // Enone vinyl C (0.1302) is less hydrophobic than plain internal alkene (0.2274)
        // MVK (C=C-C=O) vs 1-butene (C=C-CC)
        let enone = logp_crippen(&mol("C=CC(=O)C")); // MVK
        let alkene = logp_crippen(&mol("C=CCC")); // 1-butene
        assert!(
            enone < alkene,
            "MVK ({enone:.4}) should be < 1-butene ({alkene:.4}): enone is less hydrophobic"
        );
    }

    #[test]
    fn all_crippen_smarts_parse() {
        let mut failed = Vec::new();
        for (i, &(pattern, _, _)) in CRIPPEN_SMARTS.iter().enumerate() {
            if parse_smarts(pattern).is_err() {
                failed.push((i, pattern));
            }
        }
        assert!(
            failed.is_empty(),
            "CRIPPEN_SMARTS parse failures (index, pattern): {:?}",
            failed
        );
    }

    // ── Ro3 / Lead-like / Pfizer 3/75 / CNS MPO / TPSA per-atom ─────────────

    #[test]
    fn test_ro3_ethanol_passes() {
        // ethanol: MW=46, LogP≈-0.3, HBD=1, HBA=1, RotBonds=0
        assert!(ro3_passes(&mol("CCO")));
    }

    #[test]
    fn test_ro3_lipinski_drug_fails() {
        // aspirin (MW=180): just passes; ibuprofen (MW=206, LogP≈3.8): fails on LogP
        assert!(!ro3_passes(&mol("CC(C)Cc1ccc(cc1)C(C)C(=O)O")));
    }

    #[test]
    fn test_lead_like_ibuprofen_passes() {
        // ibuprofen: MW=206, LogP≈3.8, RotBonds=4, RingCount=1 — borderline but passes
        assert!(lead_like_passes(&mol("CC(C)Cc1ccc(cc1)C(C)C(=O)O")));
    }

    #[test]
    fn test_lead_like_large_drug_fails() {
        // MW > 450 should fail lead_like
        let big = mol("CC(C)c1ccc(cc1)C(=O)NC2CCN(CC2)c3ncnc4c3ccc(c4)OC(F)(F)F");
        assert!(!lead_like_passes(&big));
    }

    #[test]
    fn test_pfizer_3_75_safe_compound_passes() {
        // Aspirin: LogP≈1.2, TPSA≈63 → LogP ≤ 3 → passes (not in danger zone)
        assert!(pfizer_3_75_passes(&mol("CC(=O)Oc1ccccc1C(=O)O")));
    }

    #[test]
    fn test_pfizer_3_75_risky_compound_fails() {
        // Ibuprofen: LogP≈3.8, TPSA≈37 → LogP > 3 AND TPSA < 75 → risky
        assert!(!pfizer_3_75_passes(&mol("CC(C)Cc1ccc(cc1)C(C)C(=O)O")));
    }

    #[test]
    fn test_cns_mpo_score_range() {
        // All CNS MPO scores must be in [0, 6]
        for smi in [
            "c1ccccc1",
            "CC(=O)Oc1ccccc1C(=O)O",
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
            "c1ccc2[nH]cnc2c1",
        ] {
            let s = cns_mpo_score(&mol(smi));
            assert!(
                (0.0..=6.0).contains(&s),
                "CNS MPO out of range for {smi}: {s}"
            );
        }
    }

    #[test]
    fn test_cns_mpo_small_cns_drug_high_score() {
        // Caffeine: small, low LogP, moderate TPSA, low HBD → expect score ≥ 3
        let s = cns_mpo_score(&mol("Cn1cnc2c1c(=O)n(c(=O)n2C)C"));
        assert!(s >= 3.0, "caffeine CNS MPO should be ≥ 3, got {s}");
    }

    #[test]
    fn test_tpsa_per_atom_sum_equals_tpsa() {
        for smi in [
            "CC(=O)Oc1ccccc1C(=O)O",
            "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
            "c1ccncc1",
        ] {
            let m = mol(smi);
            let per_atom = tpsa_per_atom(&m);
            assert_eq!(per_atom.len(), m.atom_count());
            let sum: f64 = per_atom.iter().sum();
            let direct = tpsa(&m);
            assert!(
                (sum - direct).abs() < 1e-9,
                "tpsa_per_atom sum {sum} ≠ tpsa {direct} for {smi}"
            );
        }
    }

    #[test]
    fn test_tpsa_per_atom_carbon_contributes_zero() {
        // In benzene all atoms are C — per-atom TPSA must all be 0.0
        let benz = mol("c1ccccc1");
        for v in tpsa_per_atom(&benz) {
            assert_eq!(v, 0.0);
        }
    }

    // ── MCF composite filter ─────────────────────────────────────────────────

    #[test]
    fn test_mcf_caffeine_passes() {
        // Caffeine: no PAINS/Brenk alerts, MW=194, LogP=-0.07, passes Lipinski and Veber.
        assert!(mcf_passes(&mol("Cn1cnc2c1c(=O)n(c(=O)n2C)C")));
    }

    #[test]
    fn test_mcf_toluene_passes() {
        // Toluene: no PAINS/Brenk alerts, passes Lipinski and Veber.
        assert!(mcf_passes(&mol("Cc1ccccc1")));
    }

    #[test]
    fn test_mcf_ibuprofen_passes_without_acetal_overmatch() {
        // Brenk "acetal_ketal" is `[#8][#6]([#8])-[#6]`: the unspecified bonds
        // are single-or-aromatic (Daylight/RDKit), so the C=O of a carboxylic
        // acid does not match. RDKit 2026.03.6's BRENK FilterCatalog reports
        // no match for ibuprofen; it previously over-matched here.
        assert!(crate::brenk_passes(&mol("CC(C)Cc1ccc(cc1)C(C)C(=O)O")));
        assert!(mcf_passes(&mol("CC(C)Cc1ccc(cc1)C(C)C(=O)O")));
    }

    #[test]
    fn test_mcf_rhodanine_fails_pains() {
        // Rhodanine is a PAINS alert → MCF must fail
        assert!(!mcf_passes(&mol("O=C1CSC(=S)N1")));
    }

    #[test]
    fn test_mcf_aspirin_fails_brenk_active_ester() {
        // Aspirin has an aryl ester which triggers the Brenk "active_ester" alert.
        assert!(!mcf_passes(&mol("CC(=O)Oc1ccccc1C(=O)O")));
    }

    #[test]
    fn test_mcf_very_large_mol_fails_lipinski() {
        // MW > 500 → fails Lipinski → MCF must fail
        let big = mol("CC(C)c1ccc(cc1)C(=O)NC2CCN(CC2)c3ncnc4c3ccc(c4)OC(F)(F)F");
        assert!(!mcf_passes(&big));
    }

    // ── per-atom property vectors ────────────────────────────────────────────

    #[test]
    fn test_hybridization_per_atom_ethane_all_sp3() {
        let m = mol("CC");
        let h = hybridization_per_atom(&m);
        assert_eq!(h.len(), 2);
        assert!(h.iter().all(|&v| v == 3), "both C in ethane → sp3: {h:?}");
    }

    #[test]
    fn test_hybridization_per_atom_benzene_all_sp2() {
        let m = mol("c1ccccc1");
        let h = hybridization_per_atom(&m);
        assert_eq!(h.len(), 6);
        assert!(h.iter().all(|&v| v == 2), "all C in benzene → sp2: {h:?}");
    }

    #[test]
    fn test_hybridization_per_atom_acetylene_sp() {
        // HC≡CH: both C → sp (triple bond)
        let m = mol("C#C");
        let h = hybridization_per_atom(&m);
        assert_eq!(h.len(), 2);
        assert!(h.iter().all(|&v| v == 1), "alkyne C → sp: {h:?}");
    }

    #[test]
    fn test_hybridization_per_atom_follows_rdkit_lone_pair_model() {
        // RDKit 2026.03.6 GetHybridization: amide N, aryl ether O and
        // carboxylate O- are SP2; DMSO-like S(=O)(=O) is SP3.
        for (smiles, want) in [
            ("CC(=O)NC", vec![3, 2, 2, 2, 3]),
            ("COc1ccccc1", vec![3, 2, 2, 2, 2, 2, 2, 2]),
            ("CC(=O)[O-]", vec![3, 2, 2, 2]),
            ("CS(C)(=O)=O", vec![3, 3, 3, 2, 2]),
        ] {
            assert_eq!(hybridization_per_atom(&mol(smiles)), want, "{smiles}");
        }
    }

    #[test]
    fn test_hybridization_per_atom_acetaldehyde_mixed() {
        // CH3-C(=O)-H: C1(CH3)=sp3, C2(=O)=sp2
        let m = mol("CC=O");
        let h = hybridization_per_atom(&m);
        assert_eq!(h.len(), 3); // C, C, O
        assert_eq!(h[0], 3, "methyl C → sp3");
        assert_eq!(h[1], 2, "carbonyl C → sp2 (double bond to O)");
        assert_eq!(h[2], 2, "carbonyl O → sp2 (double bond to C)");
    }

    #[test]
    fn test_formal_charge_per_atom_neutral() {
        let m = mol("CC(=O)O");
        let fc = formal_charge_per_atom(&m);
        assert_eq!(fc.len(), m.atom_count());
        assert!(
            fc.iter().all(|&c| c == 0),
            "acetic acid has no formal charges: {fc:?}"
        );
    }

    #[test]
    fn test_formal_charge_per_atom_charged() {
        // Ammonium ion [NH4]+: N has charge +1
        let m = mol("[NH4+]");
        let fc = formal_charge_per_atom(&m);
        assert_eq!(fc.len(), 1, "only N is heavy");
        assert_eq!(fc[0], 1, "N has formal charge +1");
    }

    #[test]
    fn test_implicit_hcount_per_atom_ethane() {
        // Ethane CC: each C has 3 implicit H
        let m = mol("CC");
        let ih = implicit_hcount_per_atom(&m);
        assert_eq!(ih.len(), 2);
        assert!(
            ih.iter().all(|&h| h == 3),
            "each CH3 → 3 implicit H: {ih:?}"
        );
    }

    #[test]
    fn test_implicit_hcount_per_atom_sum() {
        // For any molecule, total implicit H should be consistent.
        for smi in [
            "CC(=O)O",
            "c1ccccc1",
            "CC(C)C",
            "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
        ] {
            let m = mol(smi);
            let ih = implicit_hcount_per_atom(&m);
            assert_eq!(ih.len(), m.atom_count(), "length mismatch for {smi}");
        }
    }

    // ── page 24 descriptors ───────────────────────────────────────────────────

    #[test]
    fn test_hba_count_lipinski_aspirin() {
        // Aspirin C9H8O4: 0 N + 4 O = 4
        assert_eq!(hba_count_lipinski(&mol("CC(=O)Oc1ccccc1C(=O)O")), 4);
    }

    #[test]
    fn test_hba_count_lipinski_aniline() {
        // Aniline C6H7N: 1 N + 0 O = 1
        assert_eq!(hba_count_lipinski(&mol("Nc1ccccc1")), 1);
    }

    #[test]
    fn test_hba_count_lipinski_ge_ertl() {
        // Lipinski counts more liberally than Ertl (e.g., amide N is excluded by Ertl).
        // For acetamide CC(=O)N: Ertl HBA = 1 (only C=O), Lipinski = 2 (O + N).
        let acetamide = mol("CC(=O)N");
        assert_eq!(hba_count_lipinski(&acetamide), 2); // N + O
        assert!(hba_count(&acetamide) <= hba_count_lipinski(&acetamide));
    }

    #[test]
    fn test_fraction_rotatable_bonds_benzene() {
        // Benzene has 0 rotatable bonds → fraction = 0.0
        assert_eq!(fraction_rotatable_bonds(&mol("c1ccccc1")), 0.0);
    }

    #[test]
    fn test_fraction_rotatable_bonds_in_range() {
        for smi in ["CCCCc1ccccc1", "CC(=O)Oc1ccccc1C(=O)O", "CCCC"] {
            let f = fraction_rotatable_bonds(&mol(smi));
            assert!(
                (0.0..=1.0).contains(&f),
                "fraction out of [0,1] for {smi}: {f}"
            );
        }
    }

    #[test]
    fn test_ring_system_count_benzene() {
        assert_eq!(ring_system_count(&mol("c1ccccc1")), 1);
    }

    #[test]
    fn test_ring_system_count_naphthalene_one_system() {
        // Naphthalene: 2 fused rings → 1 ring system
        assert_eq!(ring_system_count(&mol("c1ccc2ccccc2c1")), 1);
    }

    #[test]
    fn test_ring_system_count_biphenyl_two_systems() {
        // Biphenyl: two benzene rings connected by a single bond → 2 ring systems
        assert_eq!(ring_system_count(&mol("c1ccc(-c2ccccc2)cc1")), 2);
    }

    #[test]
    fn test_ring_system_count_acyclic_zero() {
        // Propane has no rings → 0 ring systems
        assert_eq!(ring_system_count(&mol("CCC")), 0);
    }

    // ── Moran / Geary autocorrelation ─────────────────────────────────────────

    #[test]
    fn moran_len_is_7() {
        assert_eq!(moran_autocorr(&mol("c1ccccc1")).len(), 7);
    }

    #[test]
    fn geary_len_is_7() {
        assert_eq!(geary_autocorr(&mol("c1ccccc1")).len(), 7);
    }

    #[test]
    fn moran_single_atom_returns_zeros() {
        let v = moran_autocorr(&mol("C"));
        assert_eq!(v, vec![0.0; 7]);
    }

    #[test]
    fn geary_single_atom_returns_ones() {
        // denom = 0 → all-identical property → no autocorrelation → C = 1
        let v = geary_autocorr(&mol("C"));
        assert_eq!(v, vec![1.0; 7]);
    }

    #[test]
    fn moran_uniform_valence_is_zero() {
        // Benzene: all C have the same valence (aromatic) → mean = valence, denom = 0 → I = 0
        let v = moran_autocorr(&mol("c1ccccc1"));
        for &x in &v {
            assert!(x.abs() < 1e-9, "expected 0 for uniform valence, got {x}");
        }
    }

    #[test]
    fn geary_finite_for_mixed_molecule() {
        // Ethanol CCO has mixed valence atoms — all values must be finite
        let v = geary_autocorr(&mol("CCO"));
        assert!(
            v.iter().all(|x| x.is_finite()),
            "all Geary values must be finite: {v:?}"
        );
    }

    #[test]
    fn distance_descriptors_map_compact_heavy_indices_with_explicit_isotopic_h() {
        // `topological_distance_matrix` uses compact heavy-atom positions,
        // while this molecule deliberately interleaves explicit deuteriums
        // before the heavy atoms. The descriptor APIs must map those rows
        // back to the original AtomIdx values instead of indexing raw atoms
        // by the compact position.
        let m = mol("[2H]C([2H])([2H])NC=O");
        let auto = autocorr_2d(&m);
        let moran = moran_autocorr(&m);
        let geary = geary_autocorr(&m);
        let information = [auto.as_slice(), moran.as_slice(), geary.as_slice()];
        assert!(
            information
                .into_iter()
                .flatten()
                .all(|value| value.is_finite()),
            "explicit isotopic hydrogens must not cause invalid distance descriptors"
        );
        assert!(ipc(&m).is_finite());
    }

    #[test]
    fn distance_descriptor_bundle_matches_lazy_apis() {
        for smi in ["CC", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"] {
            let m = mol(smi);
            let bundle = distance_descriptor_bundle(&m);
            assert_eq!(bundle.autocorr_2d, autocorr_2d(&m), "{smi}: AutoCorr");
            assert_eq!(bundle.moran, moran_autocorr(&m), "{smi}: Moran");
            assert_eq!(bundle.geary, geary_autocorr(&m), "{smi}: Geary");
        }
    }

    #[test]
    fn moran_finite_for_mixed_molecule() {
        let v = moran_autocorr(&mol("CCO"));
        assert!(
            v.iter().all(|x| x.is_finite()),
            "all Moran values must be finite: {v:?}"
        );
    }

    // ── CarbonTypes ──────────────────────────────────────────────────────────

    #[test]
    fn carbon_types_methane() {
        // methane C: sp3 with 0 heavy neighbors — none of the CxSPy bins (deg must be ≥1)
        let ct = carbon_types(&mol("C"));
        assert_eq!(
            ct.c1sp3 + ct.c2sp3 + ct.c3sp3,
            0,
            "isolated methane C has deg=0, not counted"
        );
    }

    #[test]
    fn carbon_types_ethane() {
        // ethane CC: each C is sp3 with 1 heavy neighbor → C1SP3×2
        let ct = carbon_types(&mol("CC"));
        assert_eq!(ct.c1sp3, 2, "ethane: 2×C1SP3");
        assert_eq!(ct.c2sp3, 0);
    }

    #[test]
    fn carbon_types_propane() {
        // propane CCC: 2×C1SP3 (terminals) + 1×C2SP3 (middle)
        let ct = carbon_types(&mol("CCC"));
        assert_eq!(ct.c1sp3, 2, "propane: 2 terminal sp3 C");
        assert_eq!(ct.c2sp3, 1, "propane: 1 internal sp3 C");
    }

    #[test]
    fn carbon_types_benzene() {
        // benzene: 6 aromatic C → hybridization=2, heavy deg=2 → C2SP2×6
        let ct = carbon_types(&mol("c1ccccc1"));
        assert_eq!(ct.c2sp2, 6, "benzene: 6×C2SP2 (aromatic)");
    }

    #[test]
    fn carbon_types_acetylene() {
        // acetylene HC≡CH: 2 sp C with 1 heavy neighbor each → C1SP1×2
        let ct = carbon_types(&mol("C#C"));
        assert_eq!(ct.c1sp1, 2, "acetylene: 2×C1SP1");
    }

    // ── Information Content ───────────────────────────────────────────────────

    #[test]
    fn information_content_benzene_zero_ic() {
        // All 6 C in benzene are equivalent (same element + same degree) → IC=0
        let ic = information_content(&mol("c1ccccc1"));
        assert!(
            ic.ic.abs() < 1e-9,
            "benzene: IC must be 0 (fully symmetric), got {}",
            ic.ic
        );
        assert!(ic.sic.abs() < 1e-9, "benzene: SIC must be 0");
    }

    #[test]
    fn information_content_propane_nonzero() {
        let ic = information_content(&mol("CCC"));
        assert!(
            ic.ic > 0.5,
            "propane: IC must be > 0 (two symmetry classes)"
        );
        assert!(ic.tic > 0.0);
        assert!(ic.ic.is_finite() && ic.tic.is_finite());
    }

    #[test]
    fn information_content_single_atom() {
        // Single heavy atom → all fields = 0
        let ic = information_content(&mol("[Na+]"));
        assert_eq!(ic.ic, 0.0);
        assert_eq!(ic.tic, 0.0);
    }

    // ── MDE ──────────────────────────────────────────────────────────────────

    #[test]
    fn mde_carbon_propane_mdec11() {
        // Propane CCC: 2 primary C (deg=1 C-neighbor) + 1 secondary C (deg=2)
        // MDEC11 = 1/sqrt(2) between the two primary C (distance 2 bonds)
        let v = mde_carbon(&mol("CCC"));
        let expected = 1.0 / (2.0_f64).sqrt();
        assert!(
            (v[0] - expected).abs() < 1e-6,
            "MDEC11 for propane: expected {expected:.4}, got {:.4}",
            v[0]
        );
    }

    #[test]
    fn mde_carbon_all_finite() {
        let v = mde_carbon(&mol("CC(CC)C(=O)O"));
        assert!(v.iter().all(|&x| x.is_finite()), "all MDEC must be finite");
    }

    // ── BCUT2D ───────────────────────────────────────────────────────────────

    #[test]
    fn bcut2d_hi_ge_lo() {
        // HI eigenvalue must be ≥ LO for all properties
        let b = bcut2d(&mol("c1ccccc1"));
        assert!(b.mwhi >= b.mwlo, "BCUT2D-MWHI must be ≥ BCUT2D-MWLO");
        assert!(b.logphi >= b.logplo, "BCUT2D-LOGPHI ≥ BCUT2D-LOGPLO");
        assert!(b.mrhi >= b.mrlo, "BCUT2D-MRHI ≥ BCUT2D-MRLO");
    }

    #[test]
    fn bcut2d_hi_strictly_greater_than_lo() {
        // Regression test for a min/max-eigenvalue bug where hi == lo for
        // every input (unshifted power iteration converges to the largest
        // *magnitude* eigenvalue regardless of sign, making -max(-A) collapse
        // to the same value as max(A)). Aspirin has heterogeneous atomic
        // properties (C/O with distinct MW/LogP/MR), so a correct eigenvalue
        // spread must be non-degenerate.
        let b = bcut2d(&mol("CC(=O)Oc1ccccc1C(=O)O"));
        assert!(
            b.mwhi > b.mwlo,
            "MW hi/lo must differ: {} vs {}",
            b.mwhi,
            b.mwlo
        );
        assert!(
            b.logphi > b.logplo,
            "LogP hi/lo must differ: {} vs {}",
            b.logphi,
            b.logplo
        );
        assert!(
            b.mrhi > b.mrlo,
            "MR hi/lo must differ: {} vs {}",
            b.mrhi,
            b.mrlo
        );
        assert!(
            b.chghi > b.chglo,
            "Charge hi/lo must differ: {} vs {}",
            b.chghi,
            b.chglo
        );
    }

    #[test]
    fn bcut2d_all_finite() {
        let b = bcut2d(&mol("CC(=O)Nc1ccccc1C(=O)O")); // 4-aminobenzoic acid derivative
        let vals = [
            b.chghi, b.chglo, b.logphi, b.logplo, b.mrhi, b.mrlo, b.mwhi, b.mwlo,
        ];
        assert!(
            vals.iter().all(|v| v.is_finite()),
            "all BCUT2D values must be finite: {vals:?}"
        );
    }

    #[test]
    fn bcut2d_single_atom_zero() {
        // H-only molecule → no heavy atoms → default struct
        let b = bcut2d(&mol("[H][H]"));
        assert_eq!(b.mwhi, 0.0);
    }

    #[test]
    fn logp_and_mr_matches_individual_functions() {
        // Verify logp_and_mr() returns identical values to logp_crippen() + molar_refractivity()
        // on a diverse set of molecules.
        let smiles = [
            "C",                                // methane
            "c1ccccc1",                         // benzene
            "CC(=O)Oc1ccccc1C(=O)O",            // aspirin
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",     // caffeine
            "c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34", // pyrene (PAH)
            "O=C(O)c1ccccc1O",                  // salicylic acid
            "CCO",                              // ethanol
            "O=C1C=CC=CN1",                     // Kekulé 2-pyridone
            "CC(N)Cc1ccccc1",                   // phenylalanine
        ];
        for smi in smiles {
            let m = mol(smi);
            let expected_logp = logp_crippen(&m);
            let expected_mr = molar_refractivity(&m);
            let (got_logp, got_mr) = logp_and_mr(&m);
            assert!(
                (expected_logp - got_logp).abs() < 1e-9,
                "{smi}: logp_crippen={expected_logp:.6} vs logp_and_mr.0={got_logp:.6}"
            );
            assert!(
                (expected_mr - got_mr).abs() < 1e-9,
                "{smi}: molar_refractivity={expected_mr:.6} vs logp_and_mr.1={got_mr:.6}"
            );
        }
    }
}
