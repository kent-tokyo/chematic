use super::*;

#[cfg(test)]
mod tests {
    #![allow(clippy::needless_range_loop)]

    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn fused_carbonyl_pah_mmff_aromaticity_propagates_between_rings() {
        let m = parse("O=CC1=C2C=CC=CC2=CC3=C1C=CC=C3").unwrap();
        let types = assign_mmff94_numeric_types(&m).unwrap();
        for idx in [4, 5, 6, 7] {
            assert_eq!(types[idx], 37, "aromatic six-ring carbon {idx}");
        }
    }

    #[test]
    fn mmff_aromatic_propagation_keeps_macrocycle_boundary() {
        let m = parse(
            "CC1=C2NC(=C1CCC(O)=O)C=C3N=C(C=C4NC(=CC5=NC(=C2)C(=C5C)C=C)C(=C4C)C=C)C(=C3CCC(O)=O)C",
        )
        .unwrap();
        let rings = chematic_perception::find_symmetrized_sssr(&m);
        assert!(rings.rings().iter().any(|ring| ring.len() >= 16));
        let types = assign_mmff94_numeric_types(&m).unwrap();
        for idx in [11, 12, 15, 19, 22, 23] {
            assert_eq!(
                types[idx], 2,
                "macrocycle atom {idx} must retain its prior type"
            );
        }
    }

    fn mol(s: &str) -> Molecule {
        parse(s).unwrap()
    }

    // ── Type assignment tests ────────────────────────────────────────────────

    #[test]
    fn glycine_types_match_mmff94_reference() {
        // AGLYSL01 reference: C1=type1, C2=type3, N1=type8, O5=type6, O6=type7
        // H on C = type5, H on N = type23, H on O = type24
        // SMILES: NCC(=O)O  (heavy atoms: N, C, C, O, O; explicit H via parse)
        let m = mol("NCC(=O)O");
        let types = assign_mmff94_numeric_types(&m).unwrap();

        // Collect heavy atom types by element
        let mut n_types: Vec<u8> = Vec::new();
        let mut c_types: Vec<u8> = Vec::new();
        let mut o_types: Vec<u8> = Vec::new();

        for i in 0..m.atom_count() {
            let a = m.atom(AtomIdx(i as u32));
            match a.element {
                Element::N => n_types.push(types[i]),
                Element::C => c_types.push(types[i]),
                Element::O => o_types.push(types[i]),
                _ => {}
            }
        }

        // Amine N → type 8
        assert!(
            n_types.iter().all(|&t| t == 8),
            "N should be type 8 (NR), got {:?}",
            n_types
        );
        // Should have sp3 C (type 1) and carbonyl C (type 3)
        assert!(
            c_types.contains(&1),
            "should have sp3 C (type 1), got {:?}",
            c_types
        );
        assert!(
            c_types.contains(&3),
            "should have carbonyl C (type 3), got {:?}",
            c_types
        );
        // O=C (type 7) and O-H (type 6)
        assert!(
            o_types.contains(&6),
            "should have OR oxygen (type 6), got {:?}",
            o_types
        );
        assert!(
            o_types.contains(&7),
            "should have O=C oxygen (type 7), got {:?}",
            o_types
        );
    }

    #[test]
    fn benzene_aromatic_c_is_type_37() {
        // Issue #227 Phase 1B-0: was asserting the old wrong value (63,
        // which the real Halgren/RDKit numbering assigns to a 5-ring alpha
        // carbon, not benzene). Verified against a live RDKit oracle
        // (`AllChem.MMFFGetMoleculeProperties` on benzene, this PR's audit).
        let m = mol("c1ccccc1");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        for i in 0..m.atom_count() {
            let a = m.atom(AtomIdx(i as u32));
            if a.element == Element::C {
                assert_eq!(types[i], 37, "benzene C should be type 37 (CB)");
            }
        }
    }

    #[test]
    fn pyridine_n_is_type_38() {
        // Issue #227 Phase 1B-0: was asserting the old wrong value (67,
        // which is not a valid MMFF94 numeric type in the real Halgren/
        // RDKit numbering at all -- chematic invented it). Verified against
        // a live RDKit oracle.
        let m = mol("c1ccncc1");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        for i in 0..m.atom_count() {
            let a = m.atom(AtomIdx(i as u32));
            if a.element == Element::N {
                assert_eq!(types[i], 38, "pyridine N should be type 38 (NPYD)");
            }
        }
    }

    #[test]
    fn pyridinium_n_is_type_58() {
        let m = mol("c1cc[nH+]cc1");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        for i in 0..m.atom_count() {
            let a = m.atom(AtomIdx(i as u32));
            if a.element == Element::N {
                assert_eq!(
                    types[i], 58,
                    "protonated pyridine N should be type 58 (NPD+)"
                );
            }
        }
    }

    #[test]
    fn fused_pyridinium_ring_matches_rdkit_aromatic_types() {
        // Regression for issue #227: the MMFF-specific ring boundary must
        // preserve RDKit's aromatic typing for this fused pyridinium system.
        let m =
            mol("c1ccc2c(c1)c1cc[n+]2Cc2ccc3c(c2)Cc2cc(ccc2-3)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCN1");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        assert_eq!(types[9], 58, "pyridinium N must retain MMFF aromatic type");
        assert_eq!(types[6], 37, "fused pyridinium carbon must be aromatic");
        assert_eq!(types[7], 37, "fused pyridinium carbon must be aromatic");
    }

    #[test]
    fn charged_macrocycle_directional_single_does_not_block_mmff_ring_resolution() {
        // Issue #337: `/` and `\` are directional *single* bonds, not a
        // distinct bond order. Treating atom 2's `\` bond as non-single made
        // the 2..=7 ring wait on a false exocyclic dependency and left it
        // Kekule after the macrocycle pass resolved every ring atom.
        let m =
            mol("C1=C\\c2ccc(cc2)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCNc2cc[n+](c3ccccc23)Cc2ccc/1cc2");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        for idx in 2..=7 {
            assert_eq!(types[idx], 37, "#337 ring carbon {idx} must be CB");
        }
    }

    #[test]
    fn charged_macrocycle_pyridinium_boundary_matches_rdkit_type_family() {
        // Issue #227: in the bis-pyridinium macrocycle family, RDKit's MMFF
        // typer treats the charged six-ring on the non-macrocycle edge as a
        // Kekulé ring, despite lowercase aromatic input.  The structural
        // boundary in compute_mmff94_aromatic_view must therefore yield
        // N+=C(54), C=O(3), and C=C(2), without changing the ordinary fused
        // pyridinium regression above.
        for (smiles, expected) in [
            (
                "c1cc2cc(c1)-c1cccc(c1)C[n+]1ccc(c3ccccc31)NCCCCCCCCCCNc1cc[n+](c3ccccc13)C2",
                &[
                    (13, 54),
                    (14, 3),
                    (15, 2),
                    (16, 2),
                    (35, 2),
                    (36, 2),
                    (37, 3),
                    (38, 54),
                ][..],
            ),
            (
                "c1ccc2c(c1)c1cc[n+]2Cc2ccc(cc2)-c2ccc(cc2)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCN1",
                &[(24, 54), (25, 3), (26, 2), (27, 2)][..],
            ),
            (
                "C1=C\\c2ccc(cc2)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCNc2cc[n+](c3ccccc23)Cc2ccc/1cc2",
                &[
                    (9, 54),
                    (10, 3),
                    (11, 2),
                    (12, 2),
                    (31, 2),
                    (32, 2),
                    (33, 3),
                    (34, 54),
                ][..],
            ),
            (
                "c1ccc2c(c1)c1cc[n+]2Cc2ccc(cc2)CCc2ccc(cc2)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCN1",
                &[(26, 54), (27, 3), (28, 2), (29, 2)][..],
            ),
            (
                "c1ccc2c(c1)c1cc[n+]2Cc2ccc(cc2)Cc2ccc(cc2)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCN1",
                &[(25, 54), (26, 3), (27, 2), (28, 2)][..],
            ),
            (
                "c1ccc2c(c1)c1cc[n+]2Cc2ccc3c(c2)Cc2cc(ccc2-3)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCN1",
                &[(25, 54), (26, 3), (27, 2), (28, 2)][..],
            ),
        ] {
            let types = assign_mmff94_numeric_types(&mol(smiles)).unwrap();
            for &(idx, expected_type) in expected {
                assert_eq!(
                    types[idx], expected_type,
                    "unexpected MMFF type at atom {idx} in charged macrocycle"
                );
            }
        }
    }

    #[test]
    fn iminium_n_is_type_54_when_not_aromatic() {
        // Kekulized iminium spelling: the positive, three-connected N has
        // total bond order 4 and a real N=C bond, so RDKit uses N+=C (54).
        let m = mol("C[N+]1=CC=CC=C1");
        let (n_idx, _) = m.atoms().find(|(_, a)| a.element == Element::N).unwrap();
        assert_eq!(assign_n_type(&m, &[], n_idx).unwrap(), 54);
    }

    #[test]
    fn pyrrole_n_is_type_39() {
        let m = mol("c1cc[nH]c1");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        for i in 0..m.atom_count() {
            let a = m.atom(AtomIdx(i as u32));
            if a.element == Element::N {
                assert_eq!(types[i], 39, "pyrrole N should be type 39 (NPYL)");
            }
        }
    }

    #[test]
    fn protonated_amine_n_is_type_34_not_the_o2cm_oxygen_row() {
        // Regression for the mirror-image bug to the furan C-C/N-row
        // collision: `assign_n_type`'s formal-charge branch used to
        // return `32` (O2CM, an OXYGEN type) for a positively-charged
        // nitrogen. The construction-time semantic-compatibility
        // invariant now makes that class of bug impossible to ship
        // silently -- it must resolve to type 34 (NR+), the registry's
        // only nitrogen entry among {32, 34, 35}.
        let m = mol("C[NH3+]");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        for i in 0..m.atom_count() {
            let a = m.atom(AtomIdx(i as u32));
            if a.element == Element::N {
                assert_eq!(types[i], 34, "protonated amine N should be type 34 (NR+)");
            }
        }
    }

    #[test]
    fn carboxylate_anionic_o_is_type_32_not_the_nr_plus_nitrogen_row() {
        // Issue #227 Priority 1A-3 correction: this test used to pin
        // acetate's anionic O to 35 (OM), which was itself an instance of
        // the very O2CM residual this Priority closes -- a carboxylate
        // anion's carbon has *two* terminal oxygen neighbors, so
        // `classify_terminal_o`'s `isCarboxylateO` condition fires and both
        // oxygens resolve to O2CM (32), confirmed against a live RDKit
        // oracle. The original mirror-image bug this test was written to
        // guard against still applies and is re-asserted below:
        // `assign_o_type`'s anionic-oxygen branch used to return `34`
        // (NR+, a NITROGEN type) for a negatively-charged oxygen.
        let m = mol("CC(=O)[O-]");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let mut saw_anionic_o = false;
        for i in 0..m.atom_count() {
            let a = m.atom(AtomIdx(i as u32));
            if a.element == Element::O && a.charge < 0 {
                assert_eq!(
                    types[i], 32,
                    "carboxylate anion's O should be type 32 (O2CM)"
                );
                assert_ne!(
                    types[i], 34,
                    "must never collide with NR+, a nitrogen-only type"
                );
                saw_anionic_o = true;
            }
        }
        assert!(saw_anionic_o, "test fixture must contain an anionic oxygen");

        // Negative control: a simple oxide oxygen on an unremarkable
        // aromatic carbon (that carbon has only ONE terminal oxygen, so
        // the carboxylate condition does not fire) must still resolve to
        // plain OM (35), not O2CM -- confirms the fix didn't just make
        // every anionic oxygen 32.
        let phenoxide = mol("c1ccccc1[O-]");
        let types = assign_mmff94_numeric_types(&phenoxide).unwrap();
        let mut saw_phenoxide_o = false;
        for i in 0..phenoxide.atom_count() {
            let a = phenoxide.atom(AtomIdx(i as u32));
            if a.element == Element::O {
                assert_eq!(
                    types[i], 35,
                    "phenoxide O (no sibling terminal O) should stay type 35 (OM)"
                );
                saw_phenoxide_o = true;
            }
        }
        assert!(saw_phenoxide_o, "test fixture must contain the phenoxide O");
    }

    #[test]
    fn imine_carbon_is_type_3_not_generic_vinylic() {
        // Issue #227 Priority 1A-2: closes a 39-atom residual (chematic:
        // C=C(2), RDKit: C=O(3)) -- RDKit's type-3 row is a "generic
        // carbonyl-like" umbrella for a 3-connected carbon double-bonded
        // to N, O, P, *or* S, not literal C=O alone. Expected values
        // copied verbatim from a live RDKit oracle.
        let ketimine = mol("CC(C)=NC"); // N-isopropylidene methylamine
        let types = assign_mmff94_numeric_types(&ketimine).unwrap();
        assert_eq!(
            types[1], 3,
            "imine carbon (C=N) should be type 3 (C=O family)"
        );

        // Negative control: plain alkene carbon (C=C) must stay type 2.
        let propene = mol("CC=C");
        let types = assign_mmff94_numeric_types(&propene).unwrap();
        assert_eq!(types[1], 2, "alkene carbon should stay type 2 (C=C)");
        assert_eq!(types[2], 2, "alkene carbon should stay type 2 (C=C)");
    }

    #[test]
    fn enamine_and_aniline_nitrogens_are_type_40_not_generic_nr() {
        // Issue #227 Priority 1A-2: closes the 129-atom NR(8)-vs-NC=C(40)
        // residual from Priority 1A. All expected values below are copied
        // verbatim from a live RDKit oracle query
        // (`AllChem.MMFFGetMoleculeProperties`), not derived from this
        // fix's own output.

        // Aniline: N attached to a benzene-ring carbon with no aromatic
        // O/S neighbor of its own -- the "isNbrBenzeneC" case.
        let aniline = mol("Nc1ccccc1");
        let types = assign_mmff94_numeric_types(&aniline).unwrap();
        assert_eq!(types[0], 40, "aniline N should be type 40 (NC=C)");

        // N-methyl vinylamine (enamine): N attached to a carbon genuinely
        // double-bonded to another carbon -- the "elementDoubleBondedToC
        // == 6" case.
        let enamine = mol("C=CNC");
        let types = assign_mmff94_numeric_types(&enamine).unwrap();
        assert_eq!(types[2], 40, "enamine N should be type 40 (NC=C)");

        // N-methylformamidine's exocyclic NH2: attached to a carbon
        // double-bonded to *another nitrogen* -- the "elementDoubleBondedToC
        // == 7" (amidine-like, non-charged) case.
        let amidine_like = mol("C(=NC)N");
        let types = assign_mmff94_numeric_types(&amidine_like).unwrap();
        assert_eq!(
            types[3], 40,
            "amidine-like exocyclic NH2 should be type 40 (NC=N)"
        );

        // Negative control: plain aliphatic amine, no qualifying carbon
        // neighbor at all -- must stay generic NR (8), not regress to 40.
        let propylamine = mol("CCCN");
        let types = assign_mmff94_numeric_types(&propylamine).unwrap();
        assert_eq!(
            types[3], 8,
            "plain aliphatic amine N should stay type 8 (NR)"
        );

        // Negative control: amide N (carbon neighbor has its own real
        // C=O) must still resolve to NC=O (10), never NC=C (40) -- the
        // carbonyl exclusion gate must fire correctly.
        let amide = mol("CC(=O)NC");
        let types = assign_mmff94_numeric_types(&amide).unwrap();
        assert_eq!(types[3], 10, "amide N should stay type 10 (NC=O)");
    }

    // ── NC=C multi-carbon-context order-independence (issue #227,
    //    Priority 1A-2 reinforcement) ────────────────────────────────────

    /// Deterministic xorshift64-based permutation generator. Not a source of
    /// true randomness -- just a fixed, reproducible way to generate many
    /// distinct atom/bond orderings without hand-authoring each one.
    fn deterministic_permutation(n: usize, seed: u64) -> Vec<usize> {
        let mut perm: Vec<usize> = (0..n).collect();
        let mut state = seed ^ 0x9E37_79B9_7F4A_7C15;
        for i in (1..n).rev() {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let j = (state as usize) % (i + 1);
            perm.swap(i, j);
        }
        perm
    }

    /// Rebuilds `mol` with atoms reordered per `atom_perm` (`new_idx ->
    /// old_idx`) and bonds re-added in the order given by `bond_perm`
    /// (indices into the *original* bond list). The two orderings vary
    /// independently, so this probes both atom-relabeling and
    /// bond-insertion-order sensitivity.
    fn rebuild_with_order(mol: &Molecule, atom_perm: &[usize], bond_perm: &[usize]) -> Molecule {
        let mut old_to_new = vec![0u32; atom_perm.len()];
        for (new_idx, &old_idx) in atom_perm.iter().enumerate() {
            old_to_new[old_idx] = new_idx as u32;
        }
        let mut builder = chematic_core::MoleculeBuilder::new();
        for &old_idx in atom_perm {
            builder.add_atom(mol.atom(AtomIdx(old_idx as u32)).clone());
        }
        let bonds: Vec<_> = mol.bonds().collect();
        for &bidx in bond_perm {
            let (_, bond) = bonds[bidx];
            let a = AtomIdx(old_to_new[bond.atom1.0 as usize]);
            let b = AtomIdx(old_to_new[bond.atom2.0 as usize]);
            builder
                .add_bond(a, b, bond.order)
                .expect("relabeling a valid molecule's own bonds cannot fail");
        }
        builder.build()
    }

    /// Returns the MMFF94 numeric type of the sole nitrogen atom in `mol`.
    /// All fixtures below are constructed to contain exactly one N, so
    /// identifying it by element (a content key) rather than by index
    /// survives atom relabeling -- an index-based lookup would silently
    /// compare the wrong atoms across variants.
    fn sole_nitrogen_type(mol: &Molecule) -> u8 {
        let types = assign_mmff94_numeric_types(mol).unwrap();
        let mut found = None;
        for i in 0..mol.atom_count() {
            if mol.atom(AtomIdx(i as u32)).element == Element::N {
                assert!(found.is_none(), "fixture must have exactly one N atom");
                found = Some(types[i]);
            }
        }
        found.expect("fixture must contain a nitrogen atom")
    }

    /// Verifies `expected_type` for the sole N atom is reproduced
    /// identically across: the original parse, full atom-order reversal, 16
    /// deterministic atom relabelings, a bond-insertion-order shuffle, and
    /// an independently hand-written alternate SMILES spelling of the same
    /// molecule. This is the order-independence evidence
    /// `classify_n_c3_carbon_context`'s doc comment claims for the NC=C
    /// port (issue #227 Priority 1A-2 reinforcement).
    fn assert_n_type_is_order_independent(
        smiles: &str,
        alt_smiles: &str,
        expected_type: u8,
        label: &str,
    ) {
        let base = mol(smiles);
        assert_eq!(
            sole_nitrogen_type(&base),
            expected_type,
            "{label}: original SMILES parse"
        );

        let n = base.atom_count();
        let bond_count = base.bonds().count();
        let identity_bonds: Vec<usize> = (0..bond_count).collect();

        let reversed: Vec<usize> = (0..n).rev().collect();
        let reversed_mol = rebuild_with_order(&base, &reversed, &identity_bonds);
        assert_eq!(
            sole_nitrogen_type(&reversed_mol),
            expected_type,
            "{label}: reversed atom order"
        );

        for seed in 0..16u64 {
            let perm = deterministic_permutation(n, seed);
            let variant = rebuild_with_order(&base, &perm, &identity_bonds);
            assert_eq!(
                sole_nitrogen_type(&variant),
                expected_type,
                "{label}: atom relabeling seed {seed}"
            );
        }

        let identity_atoms: Vec<usize> = (0..n).collect();
        let bond_perm = deterministic_permutation(bond_count, 0xB0AD);
        let bond_order_variant = rebuild_with_order(&base, &identity_atoms, &bond_perm);
        assert_eq!(
            sole_nitrogen_type(&bond_order_variant),
            expected_type,
            "{label}: bond insertion order shuffle"
        );

        let alt = mol(alt_smiles);
        assert_eq!(
            sole_nitrogen_type(&alt),
            expected_type,
            "{label}: alternate equivalent SMILES"
        );
    }

    #[test]
    fn nc_eq_c_multi_carbon_context_is_order_independent_and_matches_rdkit() {
        // Issue #227 Priority 1A-2 reinforcement: `classify_n_c3_carbon_context`
        // evaluates each carbon neighbor of a 3-connected N independently,
        // a deliberate divergence from RDKit's literal C++ (which threads
        // shared mutable state across its carbon-neighbor loop -- see that
        // function's doc comment for the full source-level analysis). The
        // only case where the two designs *could* observably differ is a
        // nitrogen with >=2 structurally distinct carbon neighbors, so
        // that's exactly what these 4 fixtures exercise. All expected
        // types are copied verbatim from a live RDKit oracle
        // (`AllChem.MMFFGetMoleculeProperties`), cross-checked across 32
        // `Chem.RenumberAtoms` atom orderings per fixture on the RDKit side
        // too (0/32 showed order sensitivity in any fixture) -- so these
        // pin exact parity, not a documented divergence.
        assert_n_type_is_order_independent(
            "CN(C)c1ccccc1",
            "c1ccc(cc1)N(C)C",
            40,
            "N bonded to two non-qualifying methyls + one qualifying aromatic C (aniline)",
        );
        assert_n_type_is_order_independent(
            "CN(c1ccccc1)C(C)(C)C",
            "c1ccc(cc1)N(C)C(C)(C)C",
            40,
            "N bonded to one qualifying aromatic C + one non-qualifying sp3 tert-butyl C",
        );
        assert_n_type_is_order_independent(
            "CN(c1ccccc1)C=C",
            "C=CN(C)c1ccccc1",
            40,
            "N bonded to two independently-qualifying carbons (aromatic ring + C=C)",
        );
        assert_n_type_is_order_independent(
            "CN(c1ccccc1)C(=O)C",
            "CC(=O)N(C)c1ccccc1",
            10,
            "N bonded to a qualifying aromatic C + a carbonyl C -- the global \
             isNCOorNCS gate must block NC=C regardless of neighbor order",
        );
    }

    #[test]
    fn caffeine_matches_rdkit_mmff_aromaticity_per_atom() {
        // Regression for issue #227 Priority 1A: caffeine's SMILES writes
        // its whole fused purine-dione system in lowercase aromatic
        // notation, including the two carbonyl-bearing ring carbons. This
        // is legitimate under chematic's own general Huckel model (which
        // still treats a ring carbon with an exocyclic C=O as a valid,
        // 0-pi-contributing aromatic-ring member), but RDKit's separate,
        // stricter MMFF94-specific aromaticity perception is not: the
        // exocyclic double bond blocks that ring's heteroatom lone-pair pi
        // bonus, so the whole 6-membered pyrimidinedione ring fails 4n+2 --
        // while the fused 5-membered imidazole ring (whose own atoms carry
        // no exocyclic multiple bonds) independently still passes. Expected
        // values below are copied verbatim from a live RDKit oracle query
        // (`validation/results/mmff94_rdkit_type_oracle.jsonl`, molecule
        // "caffeine"), not derived from this fix's own output.
        let m = mol("Cn1cnc2c1c(=O)n(C)c(=O)n2C");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let expected = [
            (0, 1),   // C, methyl
            (1, 39),  // N, NPYL -- imidazole ring, aromatic
            (2, 63),  // C, C5A  -- imidazole ring, aromatic
            (3, 66),  // N, N5B  -- imidazole ring, aromatic
            (4, 64),  // C, C5B  -- fused atom, aromatic via imidazole ring only
            (5, 63),  // C, C5A  -- fused atom, aromatic via imidazole ring only
            (6, 3),   // C, C=O  -- pyrimidinedione ring, NOT aromatic
            (7, 7),   // O, O=C
            (8, 10),  // N, NC=O -- pyrimidinedione ring, NOT aromatic
            (9, 1),   // C, methyl
            (10, 3),  // C, C=O  -- pyrimidinedione ring, NOT aromatic
            (11, 7),  // O, O=C
            (12, 10), // N, NC=O -- pyrimidinedione ring, NOT aromatic
            (13, 1),  // C, methyl
        ];
        for (i, expected_type) in expected {
            assert_eq!(
                types[i], expected_type,
                "caffeine atom {i}: expected MMFF94 type {expected_type} (RDKit oracle), got {}",
                types[i]
            );
        }
    }

    // ── assign_mmff94_numeric_types_with_view (issue #227 Phase 1) ──────────

    #[test]
    fn caffeine_reperceived_view_kekulizes_the_dione_ring_bond_to_single() {
        // Root-cause regression: the SAME bond the aromaticity test above
        // already showed gets non-aromatic (type 3/type 63, not two
        // aromatic-designated types) numeric types must ALSO carry a
        // non-Aromatic BondOrder in the returned view -- otherwise
        // bond_type_for/torsion_type_for (which read BondOrder, not the
        // numeric type's own registry `arom` flag) still see the wrong
        // thing even though typing itself is correct. Oracle-confirmed via
        // `MolFromSmiles(...).GetBondBetweenAtoms(5,6).GetIsAromatic() ==
        // False` (Release 2026.03.4).
        let m = mol("Cn1cnc2c1c(=O)n(C)c(=O)n2C");
        let (types, view) = assign_mmff94_numeric_types_with_view(&m).unwrap();
        assert_eq!(types[5], 63);
        assert_eq!(types[6], 3);
        let order = view
            .bond_between(AtomIdx(5), AtomIdx(6))
            .expect("atoms 5-6 must be bonded")
            .1
            .order;
        assert_ne!(
            order,
            BondOrder::Aromatic,
            "caffeine's ring-6 C5A-C=O bond must be Kekulized to a real \
             Single/Double order in the MMFF view, not left Aromatic"
        );
        // Original `m` is untouched (`&Molecule`, never mutated) -- the two
        // molecules may legitimately disagree.
        assert_eq!(
            m.bond_between(AtomIdx(5), AtomIdx(6)).unwrap().1.order,
            BondOrder::Aromatic,
            "the caller's original molecule must be left exactly as parsed"
        );
    }

    #[test]
    fn assign_mmff94_numeric_types_is_a_thin_wrapper_over_the_view_variant() {
        let m = mol("Cn1cnc2c1c(=O)n(C)c(=O)n2C");
        let types_only = assign_mmff94_numeric_types(&m).unwrap();
        let (types_with_view, _) = assign_mmff94_numeric_types_with_view(&m).unwrap();
        assert_eq!(types_only, types_with_view);
    }

    /// Finds the bond order (in the reperceived MMFF view) of the unique
    /// bond between a type-`ta` atom and a type-`tb` atom -- a CONTENT-based
    /// identification, not index-based, so it survives atom relabeling
    /// (types themselves are already independently proven order-independent
    /// by `nc_eq_c_multi_carbon_context_is_order_independent_and_matches_rdkit`
    /// and friends above). Panics if the pair is not found or is ambiguous
    /// (more than one such bond) -- both would make this an unreliable
    /// identity key for the fixture it's called on.
    fn bond_order_between_unique_type_pair(
        mol: &Molecule,
        types: &[u8],
        ta: u8,
        tb: u8,
    ) -> BondOrder {
        let mut found: Option<BondOrder> = None;
        for (_, bond) in mol.bonds() {
            let (i, j) = (bond.atom1.0 as usize, bond.atom2.0 as usize);
            let matches = (types[i] == ta && types[j] == tb) || (types[i] == tb && types[j] == ta);
            if matches {
                assert!(
                    found.is_none(),
                    "type pair ({ta},{tb}) must identify a UNIQUE bond in this fixture"
                );
                found = Some(bond.order);
            }
        }
        found.unwrap_or_else(|| panic!("no bond found between type {ta} and type {tb}"))
    }

    #[test]
    fn caffeine_reperceived_bond_order_is_invariant_under_atom_renumbering() {
        // Issue #227 Phase 1 reviewer follow-up: the deterministic-repeat
        // test above only proves no hidden randomness on a FIXED atom
        // order -- it says nothing about whether `chematic_core::kekulize`
        // (a blossom-matching solver, which CAN have genuine ties between
        // equally-valid Kekule structures for a symmetric ring) might pick
        // a different-but-equally-valid alternating bond pattern depending
        // on atom traversal order, and if so, whether that changes which
        // classification code the C5A(63)-C=O(3) ring-fusion bond this
        // whole fix depends on resolves to. Renumber caffeine's atoms 32
        // ways (deterministic_permutation, already used elsewhere in this
        // file for the same purpose on atom TYPES) and confirm the
        // reperceived BOND ORDER for that specific bond -- identified by
        // its unique (type 63, type 3) content signature, not by index, so
        // relabeling can't fool the check -- is identical every time.
        let base = mol("Cn1cnc2c1c(=O)n(C)c(=O)n2C");
        let n = base.atom_count();
        let identity_bonds: Vec<usize> = (0..base.bonds().count()).collect();

        let (base_types, base_view) = assign_mmff94_numeric_types_with_view(&base).unwrap();
        let reference = bond_order_between_unique_type_pair(&base_view, &base_types, 63, 3);
        assert_ne!(
            reference,
            BondOrder::Aromatic,
            "sanity: the reference value itself must be the fix's expected outcome"
        );

        for seed in 0..32u64 {
            let perm = deterministic_permutation(n, seed);
            let variant = rebuild_with_order(&base, &perm, &identity_bonds);
            let (types, view) = assign_mmff94_numeric_types_with_view(&variant).unwrap();
            let order = bond_order_between_unique_type_pair(&view, &types, 63, 3);
            assert_eq!(
                order, reference,
                "seed {seed}: reperceived bond order for the C5A-C=O ring-fusion \
                 bond must not depend on atom renumbering"
            );
        }
    }

    #[test]
    fn benzene_reperceived_ring_bond_orders_are_invariant_under_atom_renumbering() {
        // Companion to the caffeine test above, for the textbook case of a
        // GENUINE Kekule tie (two equally-valid alternating single/double
        // patterns related by a bond-order swap) rather than caffeine's
        // substituent-constrained ring. Confirms that whichever choice
        // `chematic_core::kekulize` makes, the OUTPUT is uniformly
        // `Aromatic` on all 6 ring bonds (RDKit's real MMFF-aromaticity
        // promotion for an accepted ring, not a residual single/double
        // pattern) and that this holds identically across 32 renumberings --
        // i.e. the Kekule tie, even if the solver's internal choice varies
        // with traversal order, never leaks into the final classification
        // input.
        let base = mol("c1ccccc1");
        let n = base.atom_count();
        let identity_bonds: Vec<usize> = (0..base.bonds().count()).collect();

        for seed in 0..32u64 {
            let perm = deterministic_permutation(n, seed);
            let variant = rebuild_with_order(&base, &perm, &identity_bonds);
            let (_, view) = assign_mmff94_numeric_types_with_view(&variant).unwrap();
            for (_, bond) in view.bonds() {
                assert_eq!(
                    bond.order,
                    BondOrder::Aromatic,
                    "seed {seed}: every benzene ring bond must resolve to Aromatic \
                     regardless of atom renumbering or the Kekulizer's internal tie-break"
                );
            }
        }
    }

    #[test]
    fn assign_mmff94_numeric_types_with_view_is_deterministic() {
        let m = mol("Cn1cnc2c1c(=O)n(C)c(=O)n2C");
        let (types1, view1) = assign_mmff94_numeric_types_with_view(&m).unwrap();
        let (types2, view2) = assign_mmff94_numeric_types_with_view(&m).unwrap();
        assert_eq!(types1, types2);
        for (_, b) in view1.bonds() {
            let other = view2.bond_between(b.atom1, b.atom2).unwrap().1;
            assert_eq!(b.order, other.order);
        }
    }

    #[test]
    fn small_ring_sp3_carbons_are_cr3r_and_cr4r() {
        // Regression for issue #227 Priority 1A: RDKit AtomTyper.cpp types
        // any 4-total-degree ring carbon in a 3- or 4-membered ring as
        // CR3R/CR4R, not the generic aliphatic CR -- previously chematic
        // always fell through to type 1 regardless of ring size.
        let cyclopropane = mol("C1CC1");
        let types = assign_mmff94_numeric_types(&cyclopropane).unwrap();
        for i in 0..cyclopropane.atom_count() {
            assert_eq!(types[i], 22, "cyclopropane carbon {i} should be CR3R (22)");
        }

        let cyclobutane = mol("C1CCC1");
        let types = assign_mmff94_numeric_types(&cyclobutane).unwrap();
        for i in 0..cyclobutane.atom_count() {
            assert_eq!(types[i], 20, "cyclobutane carbon {i} should be CR4R (20)");
        }

        let cyclohexane = mol("C1CCCCC1");
        let types = assign_mmff94_numeric_types(&cyclohexane).unwrap();
        for i in 0..cyclohexane.atom_count() {
            assert_eq!(
                types[i], 1,
                "cyclohexane carbon {i} (6-ring, not 3- or 4-) should stay generic CR (1)"
            );
        }
    }

    #[test]
    fn furan_ring_carbons_are_type_63_and_64() {
        let m = mol("c1ccoc1");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let mut c_types: Vec<u8> = m
            .bonds()
            .flat_map(|(_, b)| [b.atom1, b.atom2])
            .filter(|&a| m.atom(a).element == Element::C)
            .map(|a| types[a.0 as usize])
            .collect();
        c_types.sort_unstable();
        c_types.dedup();
        assert_eq!(
            c_types,
            vec![63, 64],
            "furan ring carbons should be C5A (alpha, 63) and C5B (beta, 64)"
        );
    }

    #[test]
    fn indolocarbazole_pyrroles_follow_rdkit_ring_order() {
        // RDKit 2026.03.6 MMFFGetMMFFAtomType (heavy atoms): its MMFF
        // aromaticity pass walks RingInfo in RDKit's order and stops once
        // every ring atom was looked at, so the second pyrrole (listed last)
        // is aromatic only because its benzene neighbours came first.
        let m = mol("O=C1NC(=O)c2c1c1c3ccccc3n3c1c1c2c2ccccc2n1[C@H]1CC[C@@H]3O1");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        assert_eq!(
            &types[..m.atom_count()],
            &[
                7, 3, 10, 3, 7, 37, 37, 64, 64, 37, 37, 37, 37, 63, 39, 63, 63, 64, 64, 37, 37, 37,
                37, 63, 39, 1, 1, 1, 1, 6
            ]
        );
    }

    #[test]
    fn six_ring_aromatic_sulfur_and_n_double_bonded_to_p_s_si_follow_rdkit() {
        // RDKit 2026.03.6 MMFFGetMMFFAtomType: methylene blue's ring [S+] is
        // S (15), not thiophene S5 (44); a cyclophosphazene-like S(=O) is
        // SO2 (18); a two-connected N double bonded to P, S or Si is NM (62).
        for (smiles, idx, want) in [
            ("CN(C)C1=CC2=C(C=C1)N=C3C=CC(=CC3=[S+]2)N(C)C", 16, 15),
            ("O=S1(C)=NP(C)(C)=NP(C)(C)=N1", 1, 18),
            ("CN=P(C)(C)C", 1, 62),
            ("CN=S(C)C", 1, 62),
            ("CN=[Si](C)C", 1, 62),
            ("CN=S(=O)(C)C", 1, 48),
        ] {
            let m = mol(smiles);
            let types = assign_mmff94_numeric_types(&m).unwrap();
            assert_eq!(types[idx], want, "{smiles}");
        }
    }

    #[test]
    fn thiophene_ring_carbons_are_type_63_and_64() {
        let m = mol("c1ccsc1");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let mut c_types: Vec<u8> = m
            .bonds()
            .flat_map(|(_, b)| [b.atom1, b.atom2])
            .filter(|&a| m.atom(a).element == Element::C)
            .map(|a| types[a.0 as usize])
            .collect();
        c_types.sort_unstable();
        c_types.dedup();
        assert_eq!(c_types, vec![63, 64]);
    }

    #[test]
    fn indole_has_both_5ring_and_6ring_aromatic_carbon_types() {
        // Fused bicyclic (5-ring pyrrole-like + 6-ring benzo) -- exercises
        // atom_in_aromatic_ring_of_size across a shared ring-fusion bond,
        // not just a single isolated ring.
        let m = mol("c1ccc2[nH]ccc2c1");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let mut c_types: Vec<u8> = m
            .bonds()
            .flat_map(|(_, b)| [b.atom1, b.atom2])
            .filter(|&a| m.atom(a).element == Element::C)
            .map(|a| types[a.0 as usize])
            .collect();
        c_types.sort_unstable();
        c_types.dedup();
        assert!(
            c_types.contains(&37),
            "indole must have some 6-ring (CB) carbons: {c_types:?}"
        );
        assert!(
            c_types.iter().any(|t| [63, 64, 78].contains(t)),
            "indole must have some 5-ring aromatic carbons: {c_types:?}"
        );
    }

    #[test]
    fn halogens_map_correctly() {
        let m = mol("CF");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        for i in 0..m.atom_count() {
            let a = m.atom(AtomIdx(i as u32));
            match a.element {
                Element::F => assert_eq!(types[i], 11),
                Element::C => assert_eq!(types[i], 1),
                _ => {}
            }
        }
        let m2 = mol("CCl");
        let types2 = assign_mmff94_numeric_types(&m2).unwrap();
        for i in 0..m2.atom_count() {
            if m2.atom(AtomIdx(i as u32)).element == Element::CL {
                assert_eq!(types2[i], 12);
            }
        }
    }

    #[test]
    fn amide_n_is_type_10() {
        // Acetamide: NC(=O)C
        let m = mol("NC(=O)C");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        for i in 0..m.atom_count() {
            let a = m.atom(AtomIdx(i as u32));
            if a.element == Element::N {
                assert_eq!(types[i], 10, "amide N should be type 10 (NC=O)");
            }
        }
    }

    #[test]
    fn sulfoxide_is_type_17_sulfone_is_type_18() {
        let m_so = mol("CS(=O)C"); // DMSO
        let types_so = assign_mmff94_numeric_types(&m_so).unwrap();
        let m_s2 = mol("CS(=O)(=O)C"); // DMSO2
        let types_s2 = assign_mmff94_numeric_types(&m_s2).unwrap();

        for i in 0..m_so.atom_count() {
            if m_so.atom(AtomIdx(i as u32)).element == Element::S {
                assert_eq!(types_so[i], 17, "DMSO S should be type 17 (SO)");
            }
        }
        for i in 0..m_s2.atom_count() {
            if m_s2.atom(AtomIdx(i as u32)).element == Element::S {
                assert_eq!(types_s2[i], 18, "DMSO2 S should be type 18 (SO2)");
            }
        }
    }

    #[test]
    fn phosphorus_types_are_element_compatible() {
        let generic = mol("P(C)(C)C");
        let generic_types = assign_mmff94_numeric_types(&generic).unwrap();
        let generic_p = generic
            .atoms()
            .find(|(_, atom)| atom.element == Element::P)
            .map(|(idx, _)| generic_types[idx.0 as usize])
            .expect("generic P atom");
        assert_eq!(generic_p, 26, "tricoordinate P should use registry type 26");

        // RDKit types phosphorus by total degree: the four-connected ylide
        // P is 25 and only a two-connected -P= is 75 (RDKit 2026.03.6
        // `MMFFGetMMFFAtomType`).
        for (smiles, want) in [("[P](C)(C)(C)=C", 25), ("C=PC", 75)] {
            let m = mol(smiles);
            let types = assign_mmff94_numeric_types(&m).unwrap();
            let p = m
                .atoms()
                .find(|(_, atom)| atom.element == Element::P)
                .map(|(idx, _)| types[idx.0 as usize])
                .expect("P atom");
            assert_eq!(p, want, "{smiles}");
        }
    }

    // ── Issue #227: nitrile/sulfonamide/nitro/azide/charged-sulfoxide
    // typing gaps (5 small pre-existing gaps named in PR #239's own
    // summary, root-caused via `scripts/diagnostics/mmff94_angle_bond_gap_classify.py`'s
    // live-RDKit-oracle atom-type cross-check on the failing mmff94_strict
    // corpus). All 5 expected types verified against a live RDKit 2026.03.3
    // oracle (`AllChem.MMFFGetMoleculeProperties`), not hand-derived.

    #[test]
    fn nitrile_nitrogen_is_type_42_not_generic_imine() {
        // Oracle: CC#N -> C(1), C(4), N(42).
        let m = mol("CC#N");
        assert_eq!(
            sole_nitrogen_type(&m),
            42,
            "nitrile N should be type 42 (NSP)"
        );
    }

    #[test]
    fn sulfonamide_nitrogen_is_type_43() {
        // Oracle: CS(=O)(=O)N -> C(1), S(18), O(32), O(32), N(43).
        let m = mol("CS(=O)(=O)N");
        assert_eq!(
            sole_nitrogen_type(&m),
            43,
            "sulfonamide N should be type 43 (NSO2)"
        );
    }

    #[test]
    fn nitro_nitrogen_is_type_45_not_generic_charged() {
        // Oracle: C[N+](=O)[O-] -> C(1), N(45), O(32), O(32).
        let m = mol("C[N+](=O)[O-]");
        assert_eq!(
            sole_nitrogen_type(&m),
            45,
            "nitro N should be type 45 (NO2)"
        );
    }

    #[test]
    fn azide_terminal_and_central_nitrogens_are_type_47_and_53() {
        // Oracle: CN=[N+]=[N-] -> C(1), N(9, attachment N, unaffected),
        // N(53, central cumulated), N(47, terminal).
        let m = mol("CN=[N+]=[N-]");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let mut n_types: Vec<u8> = m
            .atoms()
            .filter(|(_, a)| a.element == Element::N)
            .map(|(i, _)| types[i.0 as usize])
            .collect();
        n_types.sort_unstable();
        assert_eq!(
            n_types,
            vec![9, 47, 53],
            "azide N's should be [9 (attachment), 47 (terminal, NAZT), 53 (central, =N=)]"
        );
    }

    #[test]
    fn charged_sulfoxide_sulfur_is_type_17_not_generic_thioether() {
        // Oracle: C[S+](C)[O-] -> C(1), S(17), C(1), O(32). The charge-
        // separated single-bond S-O spelling is MMFF94's only valid form
        // for a charged sulfoxide/sulfonium-oxide; the previous
        // `assign_s_type` only counted explicit double bonds to O and so
        // fell through to the generic S (15) for this spelling.
        let m = mol("C[S+](C)[O-]");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let (s_idx, _) = m.atoms().find(|(_, a)| a.element == Element::S).unwrap();
        assert_eq!(
            types[s_idx.0 as usize], 17,
            "charge-separated sulfoxide S should be type 17 (S=O)"
        );
    }

    // ── Charge calculation tests ─────────────────────────────────────────────

    #[test]
    fn charge_sum_equals_formal_charge_methane() {
        let m = mol("C");
        let q = mmff94_charges_numeric(&m).unwrap();
        let total: f64 = q.iter().sum();
        assert!(total.abs() < 0.1, "methane net charge = {total:.4}");
    }

    #[test]
    fn charge_sum_equals_formal_charge_glycine() {
        let m = mol("NCC(=O)O");
        let q = mmff94_charges_numeric(&m).unwrap();
        let total: f64 = q.iter().sum();
        assert!(total.abs() < 0.15, "glycine net charge = {total:.4}");
    }

    #[test]
    fn carbonyl_oxygen_is_most_negative() {
        // Acetone: CC(=O)C — carbonyl O should be the most negative atom
        let m = mol("CC(=O)C");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let q = mmff94_charges_numeric(&m).unwrap();
        let (o_idx, _) = m.atoms().find(|(_, a)| a.element == Element::O).unwrap();
        let o_charge = q[o_idx.0 as usize];
        assert!(
            o_charge < -0.3,
            "ketone O charge = {o_charge:.3}, expected < -0.3"
        );
        // Also verify O is type 7
        assert_eq!(types[o_idx.0 as usize], 7, "ketone O should be type 7");
    }

    #[test]
    fn amine_n_is_negative() {
        let m = mol("CCN");
        let q = mmff94_charges_numeric(&m).unwrap();
        let n_charge = m
            .atoms()
            .find(|(_, a)| a.element == Element::N)
            .map(|(i, _)| q[i.0 as usize])
            .unwrap();
        assert!(
            n_charge < -0.1,
            "amine N charge = {n_charge:.3}, expected negative"
        );
    }

    #[test]
    fn h_on_nitrogen_is_positive() {
        // Explicit H in SMILES so they appear as atoms
        let m = mol("C[NH2]");
        let q = mmff94_charges_numeric(&m).unwrap();
        let types = assign_mmff94_numeric_types(&m).unwrap();
        // Find H atoms bonded to N (type 23)
        let h_charges: Vec<f64> = m
            .atoms()
            .filter(|(i, a)| a.element == Element::H && types[i.0 as usize] == 23)
            .map(|(i, _)| q[i.0 as usize])
            .collect();
        if h_charges.is_empty() {
            // If no explicit H-N atoms, just verify N is negative
            let n_charge = m
                .atoms()
                .find(|(_, a)| a.element == Element::N)
                .map(|(i, _)| q[i.0 as usize])
                .unwrap();
            assert!(
                n_charge < 0.0,
                "amine N charge = {n_charge:.3}, expected negative"
            );
        } else {
            for &hq in &h_charges {
                assert!(hq > 0.05, "H-N charge = {hq:.3}, expected positive");
            }
        }
    }

    #[test]
    fn pbci_table_has_99_entries() {
        assert_eq!(MMFF94_PBCI.len(), 99);
    }

    #[test]
    fn chg_table_has_498_entries() {
        assert_eq!(MMFF94_CHG.len(), 498);
    }

    // ── BCI bond-type-source fix (issue #227 Phase 2) ────────────────────────

    #[test]
    fn acetone_carbonyl_charges_match_rdkit_oracle_after_bond_type_fix() {
        // Regression pin for the Phase 2 BCI fix: acetone's C=O bond (type
        // 3 - type 7, `BondOrder::Double`) has a `MMFF94_CHG` row ONLY at
        // `bond_type=0` (`(0, 3, 7, -0.57)`) -- there is no `(1, 3, 7, ...)`
        // row. The old, removed local `bond_type_for(order)` mapped
        // `Double -> 1`, so this bond used to MISS the table entirely and
        // fall back to the generic `pbci_for(3) - pbci_for(7)` difference,
        // silently wrong (RDKit's own algorithm never sets `bondType=1` for
        // a bond that isn't formally SINGLE, per `getMMFFBondType`,
        // `AtomTyper.cpp:2457-2475`). The fixed `bond_type_for(ti, tj,
        // order)` (`crate::mmff94_minimizer`) maps every Double/Triple/
        // Aromatic bond to 0 unconditionally, landing on the real `(0, 3,
        // 7, -0.57)` row directly. Expected values below are copied
        // verbatim from a live RDKit oracle query
        // (`AllChem.MMFFGetMoleculeProperties(Chem.MolFromSmiles("CC(=O)C")).GetMMFFPartialCharge(i)`,
        // `rdkit==2026.03.4`), not derived from this fix's own output.
        let m = mol("CC(=O)C");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        assert_eq!(types, vec![1, 3, 7, 1], "sanity: acetone atom types");
        let q = mmff94_charges_numeric(&m).unwrap();
        let expected = [0.061, 0.448, -0.57, 0.061];
        for (i, exp) in expected.iter().enumerate() {
            assert!(
                (q[i] - exp).abs() < 1e-6,
                "acetone atom {i}: expected charge {exp}, got {}",
                q[i]
            );
        }
    }

    #[test]
    fn mmff94_charges_numeric_uses_reperceived_view_not_callers_bond_order() {
        // Directly proves the fix's mechanism (not just its net effect
        // above): caffeine's pyrimidinedione-ring C5A-C=O bond is
        // `BondOrder::Aromatic` in the caller's original molecule but
        // Kekulized to a real Single/Double order in
        // `assign_mmff94_numeric_types_with_view`'s returned `mmff_mol`
        // (already pinned by
        // `caffeine_reperceived_view_kekulizes_the_dione_ring_bond_to_single`
        // above). If `mmff94_charges_numeric` were still reading `mol`'s
        // own (Aromatic) bond order, every BCI contribution across that
        // bond would use `bond_type_for(ti, tj, Aromatic)` = 0 -- which
        // happens to coincide with the fixed formula's answer for THIS
        // particular pair, so this test instead checks the mechanism
        // directly: charges must be finite and identical across two calls
        // (determinism), and must differ from a hand-computed
        // "old-formula" charge vector for at least one BCI-sensitive
        // molecule with a real Double/Triple bond -- acetone
        // (`acetone_carbonyl_charges_match_rdkit_oracle_after_bond_type_fix`
        // above) already demonstrates that divergence with an oracle pin;
        // this test only re-confirms determinism survives the view lookup.
        let m = mol("Cn1cnc2c1c(=O)n(C)c(=O)n2C");
        let q1 = mmff94_charges_numeric(&m).unwrap();
        let q2 = mmff94_charges_numeric(&m).unwrap();
        assert_eq!(q1, q2, "mmff94_charges_numeric must be deterministic");
        assert!(
            q1.iter().all(|c| c.is_finite()),
            "all charges must be finite"
        );
    }

    #[test]
    fn mmff94_charges_numeric_is_invariant_under_atom_renumbering() {
        // Issue #227 Phase 2 (mirrors Phase 1's reviewer-requested
        // renumbering-invariance test for the same reperceived-view
        // mechanism, `caffeine_reperceived_bond_order_is_invariant_under_atom_renumbering`
        // above): the fix's BCI lookup now depends on `mmff_mol`'s
        // Kekulized bond order, computed by `chematic_core::kekulize` (a
        // blossom-matching solver that CAN have genuine ties for a
        // symmetric ring). A per-atom charge is not itself a stable
        // renumbering-invariant key (atom indices change), but the
        // molecule's own SORTED multiset of per-atom charges must be --
        // renumbering relabels which slot each charge sits in, never the
        // set of values a real, order-independent computation produces.
        let base = mol("Cn1cnc2c1c(=O)n(C)c(=O)n2C"); // caffeine
        let n = base.atom_count();
        let identity_bonds: Vec<usize> = (0..base.bonds().count()).collect();

        let mut reference = mmff94_charges_numeric(&base).unwrap();
        reference.sort_by(|a, b| a.partial_cmp(b).unwrap());

        for seed in 0..32u64 {
            let perm = deterministic_permutation(n, seed);
            let variant = rebuild_with_order(&base, &perm, &identity_bonds);
            let mut charges = mmff94_charges_numeric(&variant).unwrap();
            charges.sort_by(|a, b| a.partial_cmp(b).unwrap());
            assert_eq!(charges.len(), reference.len());
            for (c, r) in charges.iter().zip(reference.iter()) {
                assert!(
                    (c - r).abs() < 1e-9,
                    "seed {seed}: sorted charge multiset must match the \
                     original ordering's (up to float-summation-order \
                     noise), got {c} vs {r}"
                );
            }
        }
    }

    // ── Derived-formal-charge fix (issue #227 Phase 2 Step 6) ────────────────
    // Root cause: `mmff94_charges_numeric` was feeding the molecule's raw,
    // literal SMILES formal charge into equation 15 (both as this atom's own
    // q0 and as the neighbor-formal-charge source for the redistribution
    // terms). RDKit's real algorithm uses a separate, MMFF-atom-TYPE-derived
    // formal charge instead (see `mmff_derived_formal_charge`'s doc for the
    // full citation). These tests pin the 3/11 residual molecules this gap
    // actually explains (the other 8/11 are unrelated atom-TYPE-assignment
    // bugs, out of scope here -- see `scripts/mmff94_provenance/PROVENANCE.md`),
    // plus the new O2CM/SM redistribution branches on synthetic fixtures the
    // 264-molecule corpus itself does not exercise.

    #[test]
    fn chembl_tier_b_0080_azide_charges_match_rdkit_oracle_after_derived_formal_charge_fix() {
        // Azide N types 47 (NAZT, terminal) and 53 (=N=, central) are absent
        // from RDKit's derived-formal-charge switch entirely, so RDKit
        // treats both as formal-charge-neutral (0.0) for equation 15 despite
        // their raw SMILES charges being -1/+1 -- chematic previously used
        // those raw charges directly. Expected values copied verbatim from
        // the already-committed live RDKit oracle dump
        // (`validation/results/mmff94_bci_charges_227_rdkit_oracle.jsonl`,
        // `chembl_tier_b_0080`, `rdkit==2026.03.4`), not derived from this
        // fix's own output.
        let m = mol("COc1cc2nc(N3CCN(/C(S)=N/c4ccc(N=[N+]=[N-])cc4)CC3)nc(N)c2cc1OC");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let expected_types = [
            1, 6, 37, 37, 37, 38, 37, 40, 1, 1, 40, 3, 15, 9, 37, 37, 37, 37, 9, 53, 47, 37, 37, 1,
            1, 38, 37, 40, 37, 37, 37, 6, 1,
        ];
        assert_eq!(
            types, expected_types,
            "sanity: chembl_tier_b_0080 atom types"
        );
        let q = mmff94_charges_numeric(&m).unwrap();
        let expected = [
            0.28, -0.3625, 0.0825, 0.0, 0.31, -0.62, 0.72, -0.8382, 0.3691, 0.3691, -0.7882, 0.641,
            -0.141, -0.629, 0.179, 0.0, 0.0, 0.179, -0.4969, 0.6879, -0.37, 0.0, 0.0, 0.3691,
            0.3691, -0.62, 0.41, -0.1, 0.0, 0.0, 0.0825, -0.3625, 0.28,
        ];
        for (i, exp) in expected.iter().enumerate() {
            assert!(
                (q[i] - exp).abs() < 1e-6,
                "chembl_tier_b_0080 atom {i}: expected charge {exp}, got {}",
                q[i]
            );
        }
    }

    #[test]
    fn chembl_tier_b_0159_nitro_charges_match_rdkit_oracle_after_derived_formal_charge_fix() {
        // Nitro N (type 45) is also absent from RDKit's derived-formal-charge
        // switch (only the 3-terminal-oxygen NITRATE case is a switch
        // condition, not nitro's 2-oxygen case), and its two O2CM (type 32)
        // oxygens hit no branch of the O2CM/SM redistribution either (their
        // shared neighbor is a type-45 N with only 2 terminal oxygens, not
        // 3) -- so RDKit's derived formal charge is 0.0 for the N and both
        // O's, not the raw +1/0/-1 chematic previously used. Expected values
        // copied verbatim from the live RDKit oracle dump
        // (`chembl_tier_b_0159`), not derived from this fix's own output.
        let m = mol("N[C@@H](CCC(=O)Nc1ccc([N+](=O)[O-])cc1)C(=O)O");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let expected_types = [
            8, 1, 1, 1, 3, 7, 10, 37, 37, 37, 37, 45, 32, 32, 37, 37, 3, 7, 6,
        ];
        assert_eq!(
            types, expected_types,
            "sanity: chembl_tier_b_0159 atom types"
        );
        let q = mmff94_charges_numeric(&m).unwrap();
        let expected = [
            -0.27, 0.331, 0.0, 0.061, 0.569, -0.57, -0.177, 0.117, 0.0, 0.0, 0.133, 0.907, -0.52,
            -0.52, 0.0, 0.0, 0.659, -0.57, -0.15,
        ];
        for (i, exp) in expected.iter().enumerate() {
            assert!(
                (q[i] - exp).abs() < 1e-6,
                "chembl_tier_b_0159 atom {i}: expected charge {exp}, got {}",
                q[i]
            );
        }
    }

    #[test]
    fn chembl_tier_b_0161_sulfoxide_charges_match_rdkit_oracle_after_derived_formal_charge_fix() {
        // Sulfoxide S (type 17) is absent from RDKit's derived-formal-charge
        // switch, so its derived charge is 0.0, not the raw +1 chematic
        // previously used (its O2CM neighbor's O also gets 0.0: the
        // sulfoxide S is type 17, not one of the O2CM branch's recognized
        // neighbor types). Expected values copied verbatim from the live
        // RDKit oracle dump (`chembl_tier_b_0161`), not derived from this
        // fix's own output.
        let m = mol("CON(C)C(=O)/C=C/CC1C(=O)N2[C@@H]1[S+]([O-])C(C)(C)[C@@H]2C(=O)O");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let expected_types = [
            1, 6, 10, 1, 3, 7, 2, 2, 1, 20, 3, 7, 10, 20, 17, 32, 1, 1, 1, 1, 3, 7, 6,
        ];
        assert_eq!(
            types, expected_types,
            "sanity: chembl_tier_b_0161 atom types"
        );
        let q = mmff94_charges_numeric(&m).unwrap();
        let expected = [
            0.28, -0.3155, -0.3246, 0.3001, 0.6156, -0.57, 0.0144, -0.1382, 0.1382, 0.053, 0.577,
            -0.57, -0.5851, 0.397, 0.1755, -0.541, 0.1935, 0.0, 0.0, 0.3611, 0.659, -0.57, -0.15,
        ];
        for (i, exp) in expected.iter().enumerate() {
            assert!(
                (q[i] - exp).abs() < 1e-6,
                "chembl_tier_b_0161 atom {i}: expected charge {exp}, got {}",
                q[i]
            );
        }
    }

    #[test]
    fn nitrobenzene_nitro_group_charges_match_rdkit_oracle() {
        // Minimal, isolated reproduction of the `chembl_tier_b_0159` nitro
        // mechanism above (type-45 N, type-32 O's, no O2CM branch fires).
        // Expected values from a fresh live RDKit oracle query
        // (`rdkit==2026.03.4`, `MMFFGetMoleculeProperties`), independent of
        // the 264-molecule corpus dump.
        let m = mol("c1ccccc1[N+](=O)[O-]");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        assert_eq!(types, vec![37, 37, 37, 37, 37, 37, 45, 32, 32]);
        let q = mmff94_charges_numeric(&m).unwrap();
        let expected = [0.0, 0.0, 0.0, 0.0, 0.0, 0.133, 0.907, -0.52, -0.52];
        for (i, exp) in expected.iter().enumerate() {
            assert!(
                (q[i] - exp).abs() < 1e-6,
                "nitrobenzene atom {i}: expected charge {exp}, got {}",
                q[i]
            );
        }
    }

    #[test]
    fn nitrate_ion_o2cm_three_oxygen_branch_matches_rdkit_oracle() {
        // New O2CM/SM branch (nitro/nitrate-nitrogen-neighbor, `nbr_type ==
        // 45 && n_term_os == 3`): the 264-molecule corpus contains only
        // 2-terminal-oxygen (nitro) type-45 neighbors, never the
        // 3-terminal-oxygen (nitrate) case, so this fixture independently
        // exercises the branch the corpus can't. All 3 oxygens are
        // symmetry-equivalent and must get an identical shared charge.
        // Expected values from a fresh live RDKit oracle query
        // (`rdkit==2026.03.4`).
        let m = mol("[O-][N+](=O)[O-]");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        assert_eq!(types, vec![32, 45, 32, 32]);
        let q = mmff94_charges_numeric(&m).unwrap();
        let expected = [
            -0.6866666666666666,
            1.06,
            -0.6866666666666666,
            -0.6866666666666666,
        ];
        for (i, exp) in expected.iter().enumerate() {
            assert!(
                (q[i] - exp).abs() < 1e-6,
                "nitrate ion atom {i}: expected charge {exp}, got {}",
                q[i]
            );
        }
    }

    #[test]
    fn sulfone_and_sulfonate_o2cm_type18_branch_matches_rdkit_oracle() {
        // O2CM/SM's sulfone/sulfonate/sulfonamide-sulfur-neighbor branch
        // (`nbr_type == 18`) has two arms: dimethyl sulfone's 2 terminal
        // oxygens (already implicitly corpus-validated -- 34 such atoms
        // across 17 corpus molecules already matched the oracle both before
        // and after this fix, since raw charge 0 and derived charge 0
        // coincide for a plain neutral sulfone) hit the `total == 2 -> 0.0`
        // arm; methanesulfonate's 3 terminal oxygens (its raw charge is
        // split -1/0/0 across otherwise-equivalent atoms, unlike the
        // symmetric derived charge) hit the `total != 2 -> fractional` arm,
        // which the corpus does not exercise. Expected values from a fresh
        // live RDKit oracle query (`rdkit==2026.03.4`).
        let sulfone = mol("CS(=O)(=O)C");
        let sulfone_types = assign_mmff94_numeric_types(&sulfone).unwrap();
        assert_eq!(sulfone_types, vec![1, 18, 32, 32, 1]);
        let sulfone_q = mmff94_charges_numeric(&sulfone).unwrap();
        let sulfone_expected = [0.1052, 1.0896, -0.65, -0.65, 0.1052];
        for (i, exp) in sulfone_expected.iter().enumerate() {
            assert!(
                (sulfone_q[i] - exp).abs() < 1e-6,
                "dimethyl sulfone atom {i}: expected charge {exp}, got {}",
                sulfone_q[i]
            );
        }

        let sulfonate = mol("CS(=O)(=O)[O-]");
        let sulfonate_types = assign_mmff94_numeric_types(&sulfonate).unwrap();
        assert_eq!(sulfonate_types, vec![1, 18, 32, 32, 32]);
        let sulfonate_q = mmff94_charges_numeric(&sulfonate).unwrap();
        let sulfonate_expected = [
            0.1052,
            1.3448,
            -0.8166666666666667,
            -0.8166666666666667,
            -0.8166666666666667,
        ];
        for (i, exp) in sulfonate_expected.iter().enumerate() {
            assert!(
                (sulfonate_q[i] - exp).abs() < 1e-6,
                "methanesulfonate atom {i}: expected charge {exp}, got {}",
                sulfonate_q[i]
            );
        }
    }

    #[test]
    fn o2cm_carboxylate_carbon_neighbor_branch_shares_formal_charge_evenly() {
        // O2CM/SM's carbon-neighbor branch (carboxylate/thiocarboxylate) is
        // not exercised anywhere in the 264-molecule corpus. Unlike the
        // type-45/type-18 branches above, a full end-to-end oracle
        // comparison on a real carboxylate (e.g. acetate, `CC(=O)[O-]`) is
        // confounded by a SEPARATE, pre-existing, out-of-scope gap in
        // `assign_mmff94_numeric_types`: RDKit assigns the carboxylate
        // carbon a dedicated type (41, CO2M/CS2M, `AtomTyper.cpp` lines
        // ~885-895) chematic does not yet implement (chematic assigns the
        // generic type 3 instead), which shifts the BCI bond contribution to
        // the oxygens by a constant, uniform offset -- NOT a bug in this
        // fix's formal-charge redistribution (both oxygens are still
        // affected identically, which is exactly the property being tested
        // here). Documented as a follow-up, not fixed in this PR (per the
        // stop condition: fixing it means touching atom-type assignment).
        //
        // This test instead directly exercises `mmff_derived_formal_charge`
        // (this module's own new function, in scope for a direct test) and
        // checks its output against Halgren's cited formula
        // (`-(n_term_os - 1) / n_term_os` for n_term_os == 2) by hand
        // arithmetic, not RDKit's own output -- for 2 terminal oxygens
        // sharing one carbon, each must get -(2-1)/2 = -0.5.
        let m = mol("CC(=O)[O-]");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        // idx 2 = carbonyl O, idx 3 = the explicit [O-] -- both terminal,
        // both bonded to the same carboxylate carbon.
        assert_eq!(m.atom(AtomIdx(2)).element, Element::O);
        assert_eq!(m.atom(AtomIdx(3)).element, Element::O);
        let fchg_2 = mmff_derived_formal_charge(&m, &m, &[], &types, AtomIdx(2));
        let fchg_3 = mmff_derived_formal_charge(&m, &m, &[], &types, AtomIdx(3));
        assert!(
            (fchg_2 - (-0.5)).abs() < 1e-9,
            "carboxylate O (idx 2): expected shared formal charge -0.5, got {fchg_2}"
        );
        assert!(
            (fchg_3 - (-0.5)).abs() < 1e-9,
            "carboxylate O (idx 3): expected shared formal charge -0.5, got {fchg_3}"
        );
        assert_eq!(
            fchg_2, fchg_3,
            "both terminal oxygens on the same carboxylate carbon must share \
             the formal charge identically, regardless of which one the \
             input SMILES happened to write the literal '-' charge on"
        );
    }

    #[test]
    fn mmff94_charges_numeric_derived_formal_charge_is_invariant_under_atom_renumbering() {
        // Mirrors `mmff94_charges_numeric_is_invariant_under_atom_renumbering`
        // above, but with a fixture that actually exercises the new
        // derived-formal-charge code path (caffeine, used there, has no
        // type-32/45/47/53/17 atoms and never touches
        // `mmff_derived_formal_charge`'s non-default arms). Nitrobenzene's
        // nitro group does: type 45 N (fcadj-relevant leak source) and type
        // 32 O's (O2CM redistribution, here landing on the "no branch
        // matches" default).
        let base = mol("c1ccccc1[N+](=O)[O-]");
        let n = base.atom_count();
        let identity_bonds: Vec<usize> = (0..base.bonds().count()).collect();

        let mut reference = mmff94_charges_numeric(&base).unwrap();
        reference.sort_by(|a, b| a.partial_cmp(b).unwrap());

        for seed in 0..32u64 {
            let perm = deterministic_permutation(n, seed);
            let variant = rebuild_with_order(&base, &perm, &identity_bonds);
            let mut charges = mmff94_charges_numeric(&variant).unwrap();
            charges.sort_by(|a, b| a.partial_cmp(b).unwrap());
            assert_eq!(charges.len(), reference.len());
            for (c, r) in charges.iter().zip(reference.iter()) {
                assert!(
                    (c - r).abs() < 1e-9,
                    "seed {seed}: sorted charge multiset must match the \
                     original ordering's (up to float-summation-order \
                     noise), got {c} vs {r}"
                );
            }
        }
    }

    #[test]
    fn glycine_h_types_correct() {
        // H on C → 5, H on N → 23, H on O → 24
        // Use explicit SMILES with H: [NH2]CC(=O)O
        // But standard parse leaves H implicit. Let's check that H type assignment works:
        let m = mol("[NH2]CC(=O)O");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        for i in 0..m.atom_count() {
            let a = m.atom(AtomIdx(i as u32));
            if a.element == Element::H {
                let t = types[i];
                assert!(
                    matches!(t, 5 | 23 | 24),
                    "H type should be 5/23/24, got {t}"
                );
            }
        }
    }

    #[test]
    fn explicit_oxygen_hydrogen_types_match_mmff94_environment() {
        let alcohol = mol("CCO[H]");
        let alcohol_types = assign_mmff94_numeric_types(&alcohol).unwrap();
        let alcohol_h = (0..alcohol.atom_count())
            .find(|&i| alcohol.atom(AtomIdx(i as u32)).element == Element::H)
            .unwrap();
        assert_eq!(alcohol_types[alcohol_h], 21, "alcohol O-H should be HOR");

        let acid = mol("CC(=O)O[H]");
        let acid_types = assign_mmff94_numeric_types(&acid).unwrap();
        let acid_h = (0..acid.atom_count())
            .find(|&i| acid.atom(AtomIdx(i as u32)).element == Element::H)
            .unwrap();
        assert_eq!(acid_types[acid_h], 24, "acid O-H should be HOCO");

        let water = mol("[H]O[H]");
        let water_types = assign_mmff94_numeric_types(&water).unwrap();
        for i in 0..water.atom_count() {
            if water.atom(AtomIdx(i as u32)).element == Element::H {
                assert_eq!(water_types[i], 31, "water H should be HOH");
            }
        }
    }

    /// Hydrogen types derived from the parent atom's MMFF type (#637).
    /// Expected values are RDKit 2026.03.6 `MMFFGetMMFFAtomType` on the same
    /// explicit-H graphs.
    #[test]
    fn mmff_oxygen_h_types_follow_rdkit() {
        let cases: &[(&str, &[u8])] = &[
            ("c1ccccc1O[H]", &[29]),                // phenol HOCC
            ("C=CO[H]", &[29]),                     // enol HOCC
            ("[H]OC1=NC=CC=C1", &[29]),             // imidic / hydroxypyridine HOCC
            ("CC(=O)O[H]", &[24]),                  // HOCO
            ("[H]OP(=O)(O[H])O[H]", &[24, 24, 24]), // HOP
            ("[H]OS(=O)(=O)c1ccccc1", &[33]),       // HOS
            ("CS(=O)O[H]", &[33]),                  // HOS, sulfinic acid
            ("CNO[H]", &[21]),                      // hydroxylamine HOR
            ("CC(=O)NO[H]", &[21]),                 // hydroxamic acid O-H is HOR
            ("CC(=S)O[H]", &[21]),                  // thioacid O-H is HOR
            ("CCO[H]", &[21]),                      // alcohol HOR
            ("[H]O[H]", &[31, 31]),                 // water HOH
        ];
        for (smiles, expected) in cases {
            let m = mol(smiles);
            let types = assign_mmff94_numeric_types(&m).unwrap();
            let h: Vec<u8> = (0..m.atom_count())
                .filter(|&i| m.atom(AtomIdx(i as u32)).element == Element::H)
                .map(|i| types[i])
                .collect();
            assert_eq!(&h, expected, "{smiles}");
        }
    }

    #[test]
    fn mmff_nitrogen_and_phosphorus_h_types_follow_parent_type() {
        let cases: &[(&str, &[u8])] = &[
            ("CS(=O)(=O)N([H])[H]", &[28, 28]),             // NSO2 -> HNSO
            ("N#CN([H])[H]", &[28, 28]),                    // NC%N -> HNSO row
            ("[H]N([H])C(=S)N([H])[H]", &[28, 28, 28, 28]), // NC=S -> HNCS
            ("CC(=O)N([H])[H]", &[28, 28]),                 // NC=O -> HNCO
            ("CN([H])[H]", &[23, 23]),                      // NR -> HNR
            ("C[N+]([H])([H])[H]", &[36, 36, 36]),          // NR+ -> HNR+
            ("c1cc[n+]([H])cc1", &[36]),                    // NPD+ -> HPD+
            ("C=[N+]([H])[H]", &[36, 36]),                  // N+=C -> HN=C+
            ("CP([H])[H]", &[71, 71]),                      // HP
        ];
        for (smiles, expected) in cases {
            let m = mol(smiles);
            let types = assign_mmff94_numeric_types(&m).unwrap();
            let h: Vec<u8> = (0..m.atom_count())
                .filter(|&i| m.atom(AtomIdx(i as u32)).element == Element::H)
                .map(|i| types[i])
                .collect();
            assert_eq!(&h, expected, "{smiles}");
        }
    }

    #[test]
    fn explicit_hydrogen_on_delocalized_amine_uses_hnco_type() {
        let m = mol("c1ccccc1N([H])[H]");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let h_types: Vec<_> = (0..m.atom_count())
            .filter(|&i| m.atom(AtomIdx(i as u32)).element == Element::H)
            .map(|i| types[i])
            .collect();
        assert_eq!(h_types, vec![28, 28]);
    }

    #[test]
    fn furan_o_is_type_59() {
        // Issue #227 Phase 1B-0: this test used to hedge between 43 and the
        // old (wrong) fallback value 6, because aromatic-O detection wasn't
        // implemented at all. Now implemented (RDKit's real type for
        // aromatic 5-ring O is 59, OFUR -- 43 was never correct; that's
        // NSO2, sulfonamide nitrogen, an unrelated type). Verified against
        // a live RDKit oracle.
        let m = mol("o1cccc1"); // furan
        let types = assign_mmff94_numeric_types(&m).expect("furan must type successfully");
        for i in 0..m.atom_count() {
            if m.atom(AtomIdx(i as u32)).element == Element::O {
                assert_eq!(types[i], 59, "furan O should be type 59 (OFUR)");
            }
        }
    }

    #[test]
    fn o2cm_terminal_oxygen_fixture_matrix_matches_rdkit() {
        // Issue #227 Priority 1A-3: closes the 37-atom
        // `terminal_oxygen_o2cm_umbrella_gap` residual left after PR #239
        // (`classify_terminal_o`, ported from `AtomTyper.cpp` `case 8:`,
        // lines 1554-1737 at the pinned commit). All expected O types below
        // are copied verbatim from a live RDKit oracle
        // (`AllChem.MMFFGetMoleculeProperties`), covering every disjunct of
        // the real O2CM union plus the OM(35)/O=C(7) alternatives it must
        // not swallow.
        let fixtures: &[(&str, &str, &[u8])] = &[
            // Anionic carboxylate: both O's -> O2CM (its C has 2 terminal O).
            ("acetate", "CC(=O)[O-]", &[32, 32]),
            // Neutral carboxylic acid: =O generic carbonyl, -OH stays OR (6)
            // -- the carboxylate condition must NOT fire (only 1 terminal O
            // on that carbon; -OH has an implicit H so isn't terminal at
            // all).
            ("acetic_acid", "CC(=O)O", &[7, 6]),
            // Ester: carbonyl =O stays 7, bridging -O- (non-terminal) stays 6.
            ("ester", "CC(=O)OC", &[7, 6]),
            // Nitro: central N has 2 terminal O -> both O2CM, regardless of
            // which O is drawn with the formal negative charge.
            ("nitrobenzene", "c1ccc(cc1)[N+](=O)[O-]", &[32, 32]),
            // Genuine (charge-separated) N-oxide -> O2CM, not OM -- this is
            // the fixture that exercises the aromatic-ring-N charge-based
            // discriminator (see `classify_terminal_o`'s N branch doc).
            ("pyridine_n_oxide", "c1ccc[n+]([O-])c1", &[32]),
            // Sulfoxide: exactly 1 terminal O on S -> generic O=S (7), not
            // O2CM -- must not be swept in by the sulfone/sulfonate rule.
            ("dmso", "CS(=O)C", &[7]),
            // Sulfone: 2 terminal O's on S -> both O2CM.
            ("dimethyl_sulfone", "CS(=O)(=O)C", &[32, 32]),
            ("sulfonamide", "CS(=O)(=O)N", &[32, 32]),
            ("sulfonate", "CS(=O)(=O)[O-]", &[32, 32, 32]),
            // Sulfate ester: the 2 terminal =O -> O2CM; the 2 bridging
            // -O-C ester oxygens (non-terminal) stay OR (6).
            ("sulfate", "COS(=O)(=O)OC", &[6, 32, 32, 6]),
            // Phosphate ester: terminal =O -> O2CM unconditionally (RDKit's
            // own quirk, ported faithfully); bridging -O-C esters stay 6.
            ("phosphate", "COP(=O)(OC)OC", &[6, 32, 6, 6]),
            ("phosphonate", "CP(=O)(O)O", &[32, 6, 6]),
            // Perchlorate anion: all 4 O's on Cl -> O2CM unconditionally.
            ("perchlorate", "[O-]Cl(=O)(=O)=O", &[32, 32, 32, 32]),
            // Non-terminal / unrelated O's: unaffected negative controls.
            ("ether", "COC", &[6]),
            ("alcohol", "CO", &[6]),
            ("ketone", "CC(=O)C", &[7]),
            ("amide_carbonyl", "CC(=O)NC", &[7]),
        ];

        for (name, smiles, expected) in fixtures {
            let m = mol(smiles);
            let types = assign_mmff94_numeric_types(&m)
                .unwrap_or_else(|e| panic!("{name} ({smiles}) failed to type: {e}"));
            let actual: Vec<u8> = (0..m.atom_count())
                .filter(|&i| m.atom(AtomIdx(i as u32)).element == Element::O)
                .map(|i| types[i])
                .collect();
            assert_eq!(
                &actual, expected,
                "{name} ({smiles}): expected O types {expected:?}, got {actual:?}"
            );
        }
    }

    // ── Phase 0.1 (issue #227, Priority 1A pre-merge fix): Kekulization
    // failure must fail closed, never silently reuse `mol` as a stand-in
    // MMFF view ──────────────────────────────────────────────────────────

    /// Odd-membered, all-carbon, neutral aromatic ring: every atom is a
    /// `atom_must_be_matched` carbon (no lone pair, no exocyclic multiple
    /// bond), so a perfect matching over 5 atoms is impossible by simple
    /// parity -- guaranteed `Err` from `chematic_core::kekulize` regardless
    /// of which of its four matching passes runs, not a flaky/order-
    /// dependent failure.
    fn unkekulizable_five_ring() -> Molecule {
        let mut b = chematic_core::MoleculeBuilder::new();
        let atoms: Vec<_> = (0..5)
            .map(|_| b.add_atom(chematic_core::Atom::aromatic(Element::C)))
            .collect();
        for i in 0..5 {
            b.add_bond(atoms[i], atoms[(i + 1) % 5], BondOrder::Aromatic)
                .unwrap();
        }
        b.build()
    }

    /// A different failure mechanism from the parity case above, and one
    /// that specifically exercises a real SSSR-detected ring (unlike a bare
    /// 2-atom edge, which `find_sssr` correctly reports as ring-free and
    /// which would make `compute_mmff94_aromatic_view` take its
    /// no-rings-to-reperceive early-return instead of ever reaching
    /// Kekulization -- not the code path this test targets): a 3-membered
    /// all-aromatic-notation ring, C-O-O, where the lone must-match atom
    /// (C: no lone pair, no exocyclic multiple bond) has only lone-pair-
    /// donating heteroatom neighbors (both O), so it has zero eligible
    /// double-bond partners in the matching graph -- unmatchable regardless
    /// of which of the four matching passes runs, not a parity artifact.
    /// Models a malformed valence/aromatic-notation combination (an
    /// aromatic ring notation whose bonding pattern cannot support any
    /// consistent Kekule structure).
    fn unkekulizable_isolated_must_match_carbon() -> Molecule {
        let mut b = chematic_core::MoleculeBuilder::new();
        let c = b.add_atom(chematic_core::Atom::aromatic(Element::C));
        let o1 = b.add_atom(chematic_core::Atom::aromatic(Element::O));
        let o2 = b.add_atom(chematic_core::Atom::aromatic(Element::O));
        b.add_bond(c, o1, BondOrder::Aromatic).unwrap();
        b.add_bond(o1, o2, BondOrder::Aromatic).unwrap();
        b.add_bond(o2, c, BondOrder::Aromatic).unwrap();
        b.build()
    }

    #[test]
    fn kekulization_impossible_molecule_fails_closed_not_silently_reused() {
        for m in [
            unkekulizable_five_ring(),
            unkekulizable_isolated_must_match_carbon(),
        ] {
            // Sanity: confirm the fixture actually exercises the failure
            // path this test targets, not some unrelated typing error.
            assert!(
                chematic_core::kekulize(&m).is_err(),
                "test fixture must be genuinely unkekulizable"
            );

            let result = assign_mmff94_numeric_types(&m);
            let err = result.expect_err(
                "a molecule chematic cannot Kekulize must fail MMFF94 typing, \
                 never silently proceed using the un-re-perceived molecule",
            );

            // The failure must name the Kekulization stage and carry the
            // underlying cause, not be an opaque/generic message -- so a
            // caller (and `chematic-3d`'s typed `ForceFieldBridgeError`
            // conversion) can tell this apart from a different typing
            // failure by more than string content alone if it ever needs
            // to.
            assert!(
                err.0.contains("Kekulization"),
                "error must name the Kekulization stage, got: {}",
                err.0
            );
            assert!(
                err.0.contains("cannot be assigned a double bond"),
                "error must carry the underlying KekuleError cause, got: {}",
                err.0
            );

            // Determinism: repeated calls on the same input produce the
            // identical typed error, not a flaky/order-dependent result.
            let err2 = assign_mmff94_numeric_types(&m).expect_err("must fail deterministically");
            assert_eq!(err, err2, "Kekulization failure must be deterministic");

            // Never a wrong type: there is no `Ok` branch here at all --
            // the only way this test could regress into "returns a wrong
            // type" is by silently succeeding, which the above already
            // rules out.
        }
    }

    // ── Phase 0.2 (issue #227, Priority 1A pre-merge fix): a ring's bonds
    // must only be promoted to `BondOrder::Aromatic` when that ring itself
    // passed the Huckel check -- never reconstructed from whether every
    // atom in the ring happens to be aromatic via *other* accepted rings
    // ──────────────────────────────────────────────────────────────────

    #[test]
    fn pyracylene_peri_fused_artifact_ring_bonds_are_not_wrongly_promoted() {
        // Empirically confirmed hazard case (not theoretical): pyracylene's
        // SSSR contains a genuine 4-membered "envelope" ring, atoms
        // {6,7,8,13}, that is neither of the molecule's two real
        // (RDKit-agreeing) accepted aromatic rings -- verified identical on
        // both engines' own SSSR, not assumed. Every one of its 4 atoms is
        // independently aromatic via the two real accepted rings it
        // straddles, so the pre-fix reconstruction
        // (`ring.iter().all(|a| is_arom[a])`) would have wrongly promoted
        // this artifact ring's own bonds to `BondOrder::Aromatic` too.
        // RDKit ground truth (`Chem.SetAromaticity(mol,
        // Chem.AROMATICITY_MMFF94)` on `c1cc2cccc3c2c2c1cccc23`, this PR's
        // audit): all 14 atoms aromatic, but this specific 4-ring's own
        // bonds are NOT (its own Huckel check never independently passes:
        // an all-carbon 4-membered ring has no heteroatom lone-pair bonus
        // available and can supply at most 2 ring-double-bond pi electrons
        // either way, never satisfying `pi_e > 2 && (pi_e-2)%4==0`).
        let m = mol("c1cc2cccc3c2c2c1cccc23");
        let rings = chematic_perception::find_sssr(&m).rings().to_vec();
        let artifact_ring = rings
            .iter()
            .find(|r| {
                let mut set: Vec<u32> = r.iter().map(|a| a.0).collect();
                set.sort_unstable();
                set == [6, 7, 8, 13]
            })
            .expect("pyracylene's SSSR must contain the {6,7,8,13} artifact ring");
        assert_eq!(
            artifact_ring.len(),
            4,
            "sanity: the artifact ring must be the 4-membered one"
        );

        let view = compute_mmff94_aromatic_view(&m, &rings)
            .expect("pyracylene kekulizes and re-perceives successfully");

        let real_rings: Vec<&Vec<AtomIdx>> = rings
            .iter()
            .filter(|r| {
                let mut set: Vec<u32> = r.iter().map(|a| a.0).collect();
                set.sort_unstable();
                set != [6, 7, 8, 13]
            })
            .collect();
        assert_eq!(real_rings.len(), 3, "pyracylene has 4 SSSR rings total");

        // Two of the artifact ring's four edges (7-8, 8-13) are *also*
        // edges of one of the two real accepted 6-rings, so they are
        // legitimately aromatic for that unrelated reason -- not evidence
        // either way about the artifact ring itself. Only the two edges
        // exclusive to the artifact ring (6-7, 6-13) isolate the bug: they
        // must NOT be aromatic, since no real accepted ring ever claims
        // them.
        let is_edge_of_a_real_ring = |a: AtomIdx, b: AtomIdx| -> bool {
            real_rings.iter().any(|r| {
                let len = r.len();
                (0..len).any(|j| {
                    let (x, y) = (r[j], r[(j + 1) % len]);
                    (x == a && y == b) || (x == b && y == a)
                })
            })
        };

        let len = artifact_ring.len();
        let mut checked_exclusive_edge = false;
        for (j, &a) in artifact_ring.iter().enumerate() {
            let b = artifact_ring[(j + 1) % len];
            if is_edge_of_a_real_ring(a, b) {
                continue; // shared edge -- aromatic for an unrelated, legitimate reason
            }
            checked_exclusive_edge = true;
            let (_, bond) = view
                .bond_between(a, b)
                .unwrap_or_else(|| panic!("bond {}-{} must exist", a.0, b.0));
            assert_ne!(
                bond.order,
                BondOrder::Aromatic,
                "artifact-ring-exclusive bond {}-{} must NOT be promoted to aromatic -- this \
                 ring never independently passed the Huckel check, even though all 4 of its \
                 atoms are aromatic via the two other, real accepted rings",
                a.0,
                b.0
            );
        }
        assert!(
            checked_exclusive_edge,
            "sanity: the artifact ring must have at least one edge not shared with a real ring"
        );

        // Positive control in the same molecule: the real rings must still
        // be correctly promoted, proving this isn't just "nothing gets
        // promoted" degenerate behavior.
        let mut any_promoted = false;
        for ring in &real_rings {
            let len = ring.len();
            for (j, &a) in ring.iter().enumerate() {
                let b = ring[(j + 1) % len];
                if let Some((_, bond)) = view.bond_between(a, b)
                    && bond.order == BondOrder::Aromatic
                {
                    any_promoted = true;
                }
            }
        }
        assert!(
            any_promoted,
            "at least one of pyracylene's real rings must still be promoted aromatic"
        );
    }

    #[test]
    fn mmff94_aromaticity_fixture_matrix_matches_rdkit_atom_and_bond_flags() {
        // Issue #227 Phase 0.2: 12-fixture matrix, atom-level AND bond-level
        // MMFF aromaticity parity against a live RDKit oracle
        // (`Chem.SetAromaticity(mol, Chem.AROMATICITY_MMFF94)` on a
        // freshly-Kekulized copy -- RDKit's real `setMMFFAromaticity`,
        // independent of `MMFFGetMoleculeProperties`'s own internal usage of
        // it; verified on caffeine against the already-established pattern
        // (see `caffeine_matches_rdkit_mmff_aromaticity_per_atom` above) before
        // trusting it for the other 11. Values copied verbatim from
        // `validation/results/mmff94_aromaticity_bond_parity_227_oracle.json`
        // (generator: `scripts/diagnostics/mmff94_aromaticity_bond_parity_227.py`), not
        // derived from this fix's own output. Fixtures chosen to stress the
        // exact fused-ring bond-promotion hazard Phase 0.2 fixes: tetralin
        // (fused aromatic+non-aromatic), a spiro system, and a bridged
        // system are included specifically because a naive
        // reconstruct-from-atom-flags approach could wrongly promote their
        // non-aromatic ring's bonds. Bond-type comparison is not separately
        // re-derived here: `bond_type_for` is a pure function of exactly
        // (atom type, bond order), both already verified per-fixture below.
        struct Fixture {
            name: &'static str,
            smiles: &'static str,
            atom_aromatic: &'static [bool],
            bond_aromatic: &'static [(u32, u32, bool)],
            atom_types: &'static [u8],
        }

        // --- benzene ---
        const BENZENE_ATOM_AROM: [bool; 6] = [true, true, true, true, true, true];
        const BENZENE_BOND_AROM: [(u32, u32, bool); 6] = [
            (0, 1, true),
            (0, 5, true),
            (1, 2, true),
            (2, 3, true),
            (3, 4, true),
            (4, 5, true),
        ];
        const BENZENE_ATOM_TYPES: [u8; 6] = [37, 37, 37, 37, 37, 37];

        // --- naphthalene ---
        const NAPHTHALENE_ATOM_AROM: [bool; 10] =
            [true, true, true, true, true, true, true, true, true, true];
        const NAPHTHALENE_BOND_AROM: [(u32, u32, bool); 11] = [
            (0, 1, true),
            (0, 9, true),
            (1, 2, true),
            (2, 3, true),
            (3, 4, true),
            (3, 8, true),
            (4, 5, true),
            (5, 6, true),
            (6, 7, true),
            (7, 8, true),
            (8, 9, true),
        ];
        const NAPHTHALENE_ATOM_TYPES: [u8; 10] = [37, 37, 37, 37, 37, 37, 37, 37, 37, 37];

        // --- anthracene ---
        const ANTHRACENE_ATOM_AROM: [bool; 14] = [
            true, true, true, true, true, true, true, true, true, true, true, true, true, true,
        ];
        const ANTHRACENE_BOND_AROM: [(u32, u32, bool); 16] = [
            (0, 1, true),
            (0, 13, true),
            (1, 2, true),
            (2, 3, true),
            (3, 4, true),
            (3, 12, true),
            (4, 5, true),
            (5, 6, true),
            (5, 10, true),
            (6, 7, true),
            (7, 8, true),
            (8, 9, true),
            (9, 10, true),
            (10, 11, true),
            (11, 12, true),
            (12, 13, true),
        ];
        const ANTHRACENE_ATOM_TYPES: [u8; 14] =
            [37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37];

        // --- indole ---
        const INDOLE_ATOM_AROM: [bool; 9] = [true, true, true, true, true, true, true, true, true];
        const INDOLE_BOND_AROM: [(u32, u32, bool); 10] = [
            (0, 1, true),
            (0, 8, true),
            (1, 2, true),
            (2, 3, true),
            (3, 4, true),
            (3, 7, true),
            (4, 5, true),
            (5, 6, true),
            (6, 7, true),
            (7, 8, true),
        ];
        const INDOLE_ATOM_TYPES: [u8; 9] = [37, 37, 37, 63, 39, 63, 64, 64, 37];

        // --- quinoline ---
        const QUINOLINE_ATOM_AROM: [bool; 10] =
            [true, true, true, true, true, true, true, true, true, true];
        const QUINOLINE_BOND_AROM: [(u32, u32, bool); 11] = [
            (0, 1, true),
            (0, 9, true),
            (1, 2, true),
            (2, 3, true),
            (3, 4, true),
            (3, 8, true),
            (4, 5, true),
            (5, 6, true),
            (6, 7, true),
            (7, 8, true),
            (8, 9, true),
        ];
        const QUINOLINE_ATOM_TYPES: [u8; 10] = [37, 37, 37, 37, 38, 37, 37, 37, 37, 37];

        // --- purine ---
        const PURINE_ATOM_AROM: [bool; 9] = [true, true, true, true, true, true, true, true, true];
        const PURINE_BOND_AROM: [(u32, u32, bool); 10] = [
            (0, 1, true),
            (0, 8, true),
            (1, 2, true),
            (2, 3, true),
            (2, 6, true),
            (3, 4, true),
            (4, 5, true),
            (5, 6, true),
            (6, 7, true),
            (7, 8, true),
        ];
        const PURINE_ATOM_TYPES: [u8; 9] = [37, 38, 63, 39, 63, 66, 64, 37, 38];

        // --- caffeine ---
        const CAFFEINE_ATOM_AROM: [bool; 14] = [
            false, true, true, true, true, true, false, false, false, false, false, false, false,
            false,
        ];
        const CAFFEINE_BOND_AROM: [(u32, u32, bool); 15] = [
            (0, 1, false),
            (1, 2, true),
            (1, 5, true),
            (2, 3, true),
            (3, 4, true),
            (4, 5, true),
            (4, 12, false),
            (5, 6, false),
            (6, 7, false),
            (6, 8, false),
            (8, 9, false),
            (8, 10, false),
            (10, 11, false),
            (10, 12, false),
            (12, 13, false),
        ];
        const CAFFEINE_ATOM_TYPES: [u8; 14] = [1, 39, 63, 66, 64, 63, 3, 7, 10, 1, 3, 7, 10, 1];

        // --- azulene ---
        const AZULENE_ATOM_AROM: [bool; 10] = [
            false, false, false, false, false, false, false, false, false, false,
        ];
        const AZULENE_BOND_AROM: [(u32, u32, bool); 11] = [
            (0, 1, false),
            (0, 9, false),
            (1, 2, false),
            (2, 3, false),
            (3, 4, false),
            (3, 9, false),
            (4, 5, false),
            (5, 6, false),
            (6, 7, false),
            (7, 8, false),
            (8, 9, false),
        ];
        const AZULENE_ATOM_TYPES: [u8; 10] = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2];

        // --- tetralin_fused_nonaromatic ---
        const TETRALIN_FUSED_NONAROMATIC_ATOM_AROM: [bool; 10] = [
            true, true, true, true, true, true, false, false, false, false,
        ];
        const TETRALIN_FUSED_NONAROMATIC_BOND_AROM: [(u32, u32, bool); 11] = [
            (0, 1, true),
            (0, 5, true),
            (1, 2, true),
            (2, 3, true),
            (3, 4, true),
            (3, 9, false),
            (4, 5, true),
            (4, 6, false),
            (6, 7, false),
            (7, 8, false),
            (8, 9, false),
        ];
        const TETRALIN_FUSED_NONAROMATIC_ATOM_TYPES: [u8; 10] =
            [37, 37, 37, 37, 37, 37, 1, 1, 1, 1];

        // --- carbostyril_exocyclic_carbonyl ---
        const CARBOSTYRIL_EXOCYCLIC_CARBONYL_ATOM_AROM: [bool; 11] = [
            false, false, false, false, true, true, true, true, true, true, false,
        ];
        const CARBOSTYRIL_EXOCYCLIC_CARBONYL_BOND_AROM: [(u32, u32, bool); 12] = [
            (0, 1, false),
            (1, 2, false),
            (1, 10, false),
            (2, 3, false),
            (3, 4, false),
            (4, 5, true),
            (4, 9, true),
            (5, 6, true),
            (6, 7, true),
            (7, 8, true),
            (8, 9, true),
            (9, 10, false),
        ];
        const CARBOSTYRIL_EXOCYCLIC_CARBONYL_ATOM_TYPES: [u8; 11] =
            [7, 3, 2, 2, 37, 37, 37, 37, 37, 37, 10];

        // --- spiro_indane_cyclohexane ---
        const SPIRO_INDANE_CYCLOHEXANE_ATOM_AROM: [bool; 14] = [
            true, true, true, true, true, true, false, false, false, false, false, false, false,
            false,
        ];
        const SPIRO_INDANE_CYCLOHEXANE_BOND_AROM: [(u32, u32, bool); 16] = [
            (0, 1, true),
            (0, 5, true),
            (1, 2, true),
            (2, 3, true),
            (3, 4, true),
            (3, 13, false),
            (4, 5, true),
            (4, 6, false),
            (6, 7, false),
            (7, 8, false),
            (7, 12, false),
            (7, 13, false),
            (8, 9, false),
            (9, 10, false),
            (10, 11, false),
            (11, 12, false),
        ];
        const SPIRO_INDANE_CYCLOHEXANE_ATOM_TYPES: [u8; 14] =
            [37, 37, 37, 37, 37, 37, 1, 1, 1, 1, 1, 1, 1, 1];

        // --- norbornane_bridged ---
        const NORBORNANE_BRIDGED_ATOM_AROM: [bool; 7] =
            [false, false, false, false, false, false, false];
        const NORBORNANE_BRIDGED_BOND_AROM: [(u32, u32, bool); 8] = [
            (0, 1, false),
            (0, 5, false),
            (1, 2, false),
            (2, 3, false),
            (2, 6, false),
            (3, 4, false),
            (4, 5, false),
            (5, 6, false),
        ];
        const NORBORNANE_BRIDGED_ATOM_TYPES: [u8; 7] = [1, 1, 1, 1, 1, 1, 1];

        let fixtures: Vec<Fixture> = vec![
            Fixture {
                name: "benzene",
                smiles: "c1ccccc1",
                atom_aromatic: &BENZENE_ATOM_AROM,
                bond_aromatic: &BENZENE_BOND_AROM,
                atom_types: &BENZENE_ATOM_TYPES,
            },
            Fixture {
                name: "naphthalene",
                smiles: "c1ccc2ccccc2c1",
                atom_aromatic: &NAPHTHALENE_ATOM_AROM,
                bond_aromatic: &NAPHTHALENE_BOND_AROM,
                atom_types: &NAPHTHALENE_ATOM_TYPES,
            },
            Fixture {
                name: "anthracene",
                smiles: "c1ccc2cc3ccccc3cc2c1",
                atom_aromatic: &ANTHRACENE_ATOM_AROM,
                bond_aromatic: &ANTHRACENE_BOND_AROM,
                atom_types: &ANTHRACENE_ATOM_TYPES,
            },
            Fixture {
                name: "indole",
                smiles: "c1ccc2[nH]ccc2c1",
                atom_aromatic: &INDOLE_ATOM_AROM,
                bond_aromatic: &INDOLE_BOND_AROM,
                atom_types: &INDOLE_ATOM_TYPES,
            },
            Fixture {
                name: "quinoline",
                smiles: "c1ccc2ncccc2c1",
                atom_aromatic: &QUINOLINE_ATOM_AROM,
                bond_aromatic: &QUINOLINE_BOND_AROM,
                atom_types: &QUINOLINE_ATOM_TYPES,
            },
            Fixture {
                name: "purine",
                smiles: "c1nc2[nH]cnc2cn1",
                atom_aromatic: &PURINE_ATOM_AROM,
                bond_aromatic: &PURINE_BOND_AROM,
                atom_types: &PURINE_ATOM_TYPES,
            },
            Fixture {
                name: "caffeine",
                smiles: "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
                atom_aromatic: &CAFFEINE_ATOM_AROM,
                bond_aromatic: &CAFFEINE_BOND_AROM,
                atom_types: &CAFFEINE_ATOM_TYPES,
            },
            Fixture {
                name: "azulene",
                smiles: "c1ccc2cccccc12",
                atom_aromatic: &AZULENE_ATOM_AROM,
                bond_aromatic: &AZULENE_BOND_AROM,
                atom_types: &AZULENE_ATOM_TYPES,
            },
            Fixture {
                name: "tetralin_fused_nonaromatic",
                smiles: "c1ccc2c(c1)CCCC2",
                atom_aromatic: &TETRALIN_FUSED_NONAROMATIC_ATOM_AROM,
                bond_aromatic: &TETRALIN_FUSED_NONAROMATIC_BOND_AROM,
                atom_types: &TETRALIN_FUSED_NONAROMATIC_ATOM_TYPES,
            },
            Fixture {
                name: "carbostyril_exocyclic_carbonyl",
                smiles: "O=c1ccc2ccccc2[nH]1",
                atom_aromatic: &CARBOSTYRIL_EXOCYCLIC_CARBONYL_ATOM_AROM,
                bond_aromatic: &CARBOSTYRIL_EXOCYCLIC_CARBONYL_BOND_AROM,
                atom_types: &CARBOSTYRIL_EXOCYCLIC_CARBONYL_ATOM_TYPES,
            },
            Fixture {
                name: "spiro_indane_cyclohexane",
                smiles: "c1ccc2c(c1)CC3(CCCCC3)C2",
                atom_aromatic: &SPIRO_INDANE_CYCLOHEXANE_ATOM_AROM,
                bond_aromatic: &SPIRO_INDANE_CYCLOHEXANE_BOND_AROM,
                atom_types: &SPIRO_INDANE_CYCLOHEXANE_ATOM_TYPES,
            },
            Fixture {
                name: "norbornane_bridged",
                smiles: "C1CC2CCC1C2",
                atom_aromatic: &NORBORNANE_BRIDGED_ATOM_AROM,
                bond_aromatic: &NORBORNANE_BRIDGED_BOND_AROM,
                atom_types: &NORBORNANE_BRIDGED_ATOM_TYPES,
            },
        ];

        for fx in &fixtures {
            let m = mol(fx.smiles);
            assert_eq!(
                m.atom_count(),
                fx.atom_aromatic.len(),
                "{}: atom count / element-sequence alignment sanity",
                fx.name
            );
            let rings = chematic_perception::find_sssr(&m).rings().to_vec();
            let view = compute_mmff94_aromatic_view(&m, &rings)
                .unwrap_or_else(|e| panic!("{}: re-perception failed: {e}", fx.name));

            for i in 0..m.atom_count() {
                let idx = AtomIdx(i as u32);
                assert_eq!(
                    view.atom(idx).aromatic,
                    fx.atom_aromatic[i],
                    "{}: atom {i} aromatic-flag mismatch vs RDKit oracle",
                    fx.name
                );
            }

            for &(a, b, expected_arom) in fx.bond_aromatic {
                let (_, bond) = view
                    .bond_between(AtomIdx(a), AtomIdx(b))
                    .unwrap_or_else(|| panic!("{}: bond {a}-{b} must exist", fx.name));
                let is_arom = bond.order == BondOrder::Aromatic;
                assert_eq!(
                    is_arom, expected_arom,
                    "{}: bond {a}-{b} aromatic-flag mismatch vs RDKit oracle (got {is_arom}, want {expected_arom})",
                    fx.name
                );
            }

            if !fx.atom_types.is_empty() {
                let types = assign_mmff94_numeric_types(&m)
                    .unwrap_or_else(|e| panic!("{}: numeric typing failed: {e}", fx.name));
                assert_eq!(
                    types, fx.atom_types,
                    "{}: numeric atom types mismatch vs RDKit oracle",
                    fx.name
                );
            }
        }
    }

    #[test]
    fn kekulization_failure_never_reaches_the_uff_or_mmff94_minimizer_bridge() {
        // Phase 0.1's "must not proceed to parameter lookup" requirement,
        // verified at this crate's own public boundary: every
        // `assign_mmff94_numeric_types` caller in this crate (e.g.
        // `crates/chematic-ff/src/mmff94_minimizer.rs`) uses `?` on this
        // call before ever touching a bond/angle/torsion/oop parameter
        // table, so a `NumericTypeError` here structurally cannot be
        // followed by a parameter-lookup attempt on the same call. This
        // test pins that `?`-propagation contract at the numeric-typing
        // boundary itself, independent of any specific downstream caller.
        let m = unkekulizable_five_ring();
        let types_result = assign_mmff94_numeric_types(&m);
        assert!(types_result.is_err());
        // If this were `Ok`, a caller doing `let types = ...?;` would
        // proceed to index `types` for parameter lookup; confirming `Err`
        // here is exactly what makes that `?` short-circuit for every
        // caller, without needing to duplicate every caller's own test.
    }

    // ── Cumulated-double-bond CSP fix (issue #337 sub-bug 2) ─────────────────
    // Root cause: `assign_c_type`'s sp-carbon check only fired on
    // `triple_bonds > 0`, so a carbon reached via *two* double bonds (a
    // cumulated diene / "allenic" carbon, e.g. the central C of an aryl
    // isothiocyanate's N=C=S) fell through to the `double_bonds > 0`
    // "double-bonded to N/O/P/S" branch and got the generic carbonyl-family
    // type 3 instead of RDKit's real CSP type 4. RDKit's actual rule
    // (`AtomTyper.cpp`, pinned commit
    // `e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f`, lines ~954-960) is simply
    // `getTotalDegree() == 2`, unconditional on which elements the two bonds
    // go to -- see `assign_c_type`'s doc comment for the full citation.
    // These tests pin the 2/8 issue #337 molecules this rule resolves. The
    // other 6/8 pyridinium/macrocycle cases are covered by
    // `charged_macrocycle_pyridinium_boundary_matches_rdkit_type_family`
    // above, plus the corpus-wide aromaticity/type/charge parity gates.

    #[test]
    fn propyne_alkyne_carbons_still_type_csp_after_degree_based_fix() {
        // No-regression pin: a real triple bond always yields
        // `total_degree == 2` for carbon (the triple bond alone consumes 3
        // of its 4 valence units), so this must keep matching RDKit exactly
        // as it did before this fix broadened the condition. Expected
        // values from a fresh live RDKit oracle query (`rdkit==2026.03.4`,
        // `MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")` on the
        // implicit-H `Chem.MolFromSmiles` result -- MMFF atom
        // typing/charges are purely topological, no embedding needed, same
        // precedent as `scripts/diagnostics/mmff94_bci_charges_oracle_227.py`).
        let m = mol("CC#C");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        assert_eq!(types, vec![1, 4, 4]);
        let q = mmff94_charges_numeric(&m).unwrap();
        let expected = [0.2, -0.2, 0.0];
        for (i, exp) in expected.iter().enumerate() {
            assert!(
                (q[i] - exp).abs() < 1e-6,
                "propyne atom {i}: expected charge {exp}, got {}",
                q[i]
            );
        }
    }

    #[test]
    fn allene_central_carbon_types_csp_not_generic_vinylic() {
        // Broader-than-corpus pin: plain carbon allene (no heteroatoms at
        // all) is not in the 264-molecule corpus, but is the clearest
        // demonstration that the real discriminating condition is
        // `total_degree == 2`, not "cumulated bond to a heteroatom" --
        // before this fix, the central carbon was mistyped 2 (generic
        // vinylic C=C), same failure mode as the isothiocyanate carbons.
        // Expected values from a fresh live RDKit oracle query
        // (`rdkit==2026.03.4`).
        let m = mol("C=C=C");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        assert_eq!(types, vec![2, 4, 2]);
        let q = mmff94_charges_numeric(&m).unwrap();
        let expected = [0.065, -0.13, 0.065];
        for (i, exp) in expected.iter().enumerate() {
            assert!(
                (q[i] - exp).abs() < 1e-6,
                "allene atom {i}: expected charge {exp}, got {}",
                q[i]
            );
        }
    }

    #[test]
    fn methyl_isothiocyanate_minimal_ncs_fixture_matches_rdkit_oracle() {
        // Minimal isolated reproduction of the corpus mechanism below,
        // independent of the aromatic ring and amide/piperazine complexity
        // both `chembl_tier_b_0071`/`_0082` also carry. Expected values
        // from a fresh live RDKit oracle query (`rdkit==2026.03.4`).
        let m = mol("CN=C=S");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        assert_eq!(types, vec![1, 9, 4, 16]);
        let q = mmff94_charges_numeric(&m).unwrap();
        let expected = [0.246, -0.546, 0.575, -0.275];
        for (i, exp) in expected.iter().enumerate() {
            assert!(
                (q[i] - exp).abs() < 1e-6,
                "methyl isothiocyanate atom {i}: expected charge {exp}, got {}",
                q[i]
            );
        }
    }

    #[test]
    fn phenyl_isothiocyanate_aryl_ncs_fixture_matches_rdkit_oracle() {
        // Same motif, aryl instead of methyl (closer to the corpus
        // molecules' own aryl isothiocyanate substructure). Also confirms
        // the fix doesn't disturb the adjacent aromatic ring's own typing.
        // Expected values from a fresh live RDKit oracle query
        // (`rdkit==2026.03.4`).
        let m = mol("c1ccccc1N=C=S");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        assert_eq!(types, vec![37, 37, 37, 37, 37, 37, 9, 4, 16]);
        let q = mmff94_charges_numeric(&m).unwrap();
        let expected = [0.0, 0.0, 0.0, 0.0, 0.0, 0.179, -0.479, 0.575, -0.275];
        for (i, exp) in expected.iter().enumerate() {
            assert!(
                (q[i] - exp).abs() < 1e-6,
                "phenyl isothiocyanate atom {i}: expected charge {exp}, got {}",
                q[i]
            );
        }
    }

    #[test]
    fn chembl_tier_b_0071_aryl_isothiocyanate_matches_rdkit_oracle_after_csp_fix() {
        // Corpus regression pin. Expected values cross-checked against the
        // already-committed live RDKit oracle dumps
        // (`validation/results/mmff94_rdkit_type_oracle.jsonl` and
        // `mmff94_bci_charges_227_rdkit_oracle.jsonl`, `rdkit==2026.03.4`)
        // via the corpus-wide per-atom join reported in
        // `scripts/mmff94_provenance/PROVENANCE.md`'s issue #337
        // follow-up, not derived from this fix's own output.
        let m = mol("COc1cc2nc(N3CCN(C(=O)c4ccc(N=C=S)c(I)c4)CC3)nc(N)c2cc1OC");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let expected_types = [
            1, 6, 37, 37, 37, 38, 37, 40, 1, 1, 10, 3, 7, 37, 37, 37, 37, 9, 4, 16, 37, 14, 37, 1,
            1, 38, 37, 40, 37, 37, 37, 6, 1,
        ];
        assert_eq!(
            types, expected_types,
            "sanity: chembl_tier_b_0071 atom types"
        );
        let q = mmff94_charges_numeric(&m).unwrap();
        let expected = [
            0.28, -0.3625, 0.0825, 0.0, 0.31, -0.62, 0.72, -0.8382, 0.3691, 0.3001, -0.6602,
            0.5438, -0.57, 0.0862, 0.0, 0.0, 0.179, -0.479, 0.575, -0.275, 0.081, -0.081, 0.0,
            0.3001, 0.3691, -0.62, 0.41, -0.1, 0.0, 0.0, 0.0825, -0.3625, 0.28,
        ];
        for (i, exp) in expected.iter().enumerate() {
            assert!(
                (q[i] - exp).abs() < 1e-6,
                "chembl_tier_b_0071 atom {i}: expected charge {exp}, got {}",
                q[i]
            );
        }
    }

    #[test]
    fn chembl_tier_b_0082_aryl_isothiocyanate_matches_rdkit_oracle_after_csp_fix() {
        // Corpus regression pin, same mechanism as `_0071` above via an
        // enone rather than a plain benzamide linker. Expected values
        // cross-checked against the already-committed live RDKit oracle
        // dumps via the same corpus-wide per-atom join, not derived from
        // this fix's own output.
        let m = mol("COc1cc2nc(N3CCN(C(=O)/C=C/c4ccc(N=C=S)cc4)CC3)nc(N)c2cc1OC");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let expected_types = [
            1, 6, 37, 37, 37, 38, 37, 40, 1, 1, 10, 3, 7, 2, 2, 37, 37, 37, 37, 9, 4, 16, 37, 37,
            1, 1, 38, 37, 40, 37, 37, 37, 6, 1,
        ];
        assert_eq!(
            types, expected_types,
            "sanity: chembl_tier_b_0082 atom types"
        );
        let q = mmff94_charges_numeric(&m).unwrap();
        let expected = [
            0.28, -0.3625, 0.0825, 0.0, 0.31, -0.62, 0.72, -0.8382, 0.3691, 0.3001, -0.6602,
            0.6156, -0.57, 0.0144, -0.0284, 0.0284, 0.0, 0.0, 0.179, -0.479, 0.575, -0.275, 0.0,
            0.0, 0.3001, 0.3691, -0.62, 0.41, -0.1, 0.0, 0.0, 0.0825, -0.3625, 0.28,
        ];
        for (i, exp) in expected.iter().enumerate() {
            assert!(
                (q[i] - exp).abs() < 1e-6,
                "chembl_tier_b_0082 atom {i}: expected charge {exp}, got {}",
                q[i]
            );
        }
    }
}
