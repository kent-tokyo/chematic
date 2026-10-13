#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn re_perception_shortcut_only_skips_non_extending_verdicts() {
        for smi in [
            "c1ccc2c(c1)CCC2",
            "O=C1NCCc2ccccc21",
            "c1ccc2c(c1)-c1ccccc1-2",
            "O=c1[nH]c2ccccc2c2ccccc12",
            "C1=Cc2ccccc2C1",
            "O=C1C=CC(=O)c2ccccc21",
            "c1ccc2c(c1)[nH]c1ccccc12",
            "O=C1c2ccccc2C(=O)N1C",
        ] {
            let mol = chematic_smiles::parse(smi).unwrap();
            if !explicit_aromaticity_is_consistent(&mol) {
                continue;
            }
            let cleared = clear_aromatic_flags(&mol);
            let Ok(k) = chematic_core::kekulize(&cleared) else {
                continue;
            };
            let kekulized = chematic_core::apply_kekule(&cleared, &k);
            if re_perception_can_extend(&mol, &kekulized) {
                continue;
            }
            let (atoms, bonds) = rdkit_parity_aromaticity(&kekulized);
            for a in &atoms {
                assert!(mol.atom(*a).aromatic, "{smi}: new aromatic atom {a:?}");
            }
            for b in &bonds {
                assert_eq!(
                    mol.bond(*b).order,
                    BondOrder::Aromatic,
                    "{smi}: new aromatic bond {b:?}"
                );
            }
        }
    }

    fn mol_kekulized(smiles: &str) -> Molecule {
        let mol = chematic_smiles::parse(smiles).expect("valid SMILES");
        let k = chematic_core::kekulize(&mol).expect("kekulizable");
        chematic_core::apply_kekule(&mol, &k)
    }

    // Calibration battery, RDKit-atom-index-verified (not guessed): every
    // entry here was checked against a live `rdkit.Chem.MolFromSmiles(...)`
    // atom-aromaticity dump before being pinned. Covers the exact cases that
    // motivated this module: simple monocyclics (benzene/pyrrole/furan/
    // thiophene), the exocyclic-carbonyl-in-ring rule (tropone/2-pyridone/
    // 4-pyranone), a genuine bridgehead spanning two valid rings
    // (indolizine), a non-alternant fused bicyclic needing the whole-perimeter
    // candidate (azulene), plain fused benzenoids (naphthalene/anthracene),
    // fused heteroaromatics (indole/quinoline/purine), and both open findings
    // from Aromaticity-A1-1a (the false-positive reproducer, purine).
    #[test]
    fn calibration_battery_matches_rdkit() {
        let cases: &[(&str, &str, &[u32])] = &[
            ("benzene", "c1ccccc1", &[0, 1, 2, 3, 4, 5]),
            ("pyrrole", "c1cc[nH]c1", &[0, 1, 2, 3, 4]),
            ("furan", "c1ccoc1", &[0, 1, 2, 3, 4]),
            ("thiophene", "c1ccsc1", &[0, 1, 2, 3, 4]),
            (
                "selenophene (Se analog control, pre-Kekulized input)",
                "C1=C[Se]C=C1",
                &[0, 1, 2, 3, 4],
            ),
            (
                "tellurophene (pre-Kekulized input; regression test for the Te \
                 normal_valences() gap fixed in chematic-core's element.rs — Te previously \
                 had no valence-table entry, so default_valence/count_atom_pi_electrons/\
                 get_atom_electron_donor_type all returned None and \
                 is_atom_candidate_for_aromaticity rejected it outright)",
                "C1=C[Te]C=C1",
                &[0, 1, 2, 3, 4],
            ),
            ("tropone", "O=c1cccccc1", &[1, 2, 3, 4, 5, 6, 7]),
            ("2-pyridone", "O=c1cccc[nH]1", &[1, 2, 3, 4, 5, 6]),
            ("4-pyranone", "O=c1ccocc1", &[1, 2, 3, 4, 5, 6]),
            (
                "indolizine (true bridgehead, both rings valid)",
                "c1ccn2ccccc12",
                &[0, 1, 2, 3, 4, 5, 6, 7, 8],
            ),
            (
                "azulene (non-alternant, needs whole-perimeter candidate)",
                "C1=CC2=CC=CC=CC2=C1",
                &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            ),
            (
                "naphthalene",
                "c1ccc2ccccc2c1",
                &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            ),
            (
                "anthracene",
                "c1ccc2cc3ccccc3cc2c1",
                &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
            ),
            ("indole", "c1ccc2[nH]ccc2c1", &[0, 1, 2, 3, 4, 5, 6, 7, 8]),
            (
                "quinoline",
                "c1ccc2ncccc2c1",
                &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            ),
            (
                "purine (Aromaticity-A1-1a open finding, fixed here)",
                "c1cnc2[nH]cnc2n1",
                &[0, 1, 2, 3, 4, 5, 6, 7, 8],
            ),
            (
                "PR #86 false-positive reproducer (Aromaticity-A1-1a open finding, fixed here)",
                "C1=Cc2ccccc2C2=NCCCN12",
                &[2, 3, 4, 5, 6, 7],
            ),
        ];

        for (name, smi, expected) in cases {
            let mol = mol_kekulized(smi);
            let (atoms, _bonds) = rdkit_parity_aromaticity(&mol);
            let mut got: Vec<u32> = atoms.iter().map(|a| a.0).collect();
            got.sort();
            assert_eq!(&got, expected, "{name} ({smi}): should match RDKit exactly");
        }
    }

    #[test]
    fn issue_767_neutral_hypercoordinate_ring_chalcogens_are_not_aromatic() {
        for (name, smi) in [
            (
                "neutral sulfur with two exocyclic single bonds",
                "C1=CS(O)(O)C=C1",
            ),
            (
                "neutral sulfur with one exocyclic single bond",
                "C1=CS(O)C=C1",
            ),
            (
                "explicit aromatic spelling of hypercoordinate sulfur",
                "c1s(O)(O)ccc1",
            ),
            ("explicit aromatic thiophene 1-oxide", "O=s1cccc1"),
            ("explicit aromatic thiophene 1,1-dioxide", "O=s1(=O)cccc1"),
        ] {
            let mol = chematic_smiles::parse(smi).expect("valid issue #767 SMILES");
            let applied = apply_aromaticity_rdkit_parity_experimental(&mol)
                .expect("issue #767 input must be re-perceived");
            assert_eq!(
                applied.atoms().filter(|(_, atom)| atom.aromatic).count(),
                0,
                "{name}: {smi}"
            );
            assert_eq!(
                applied
                    .bonds()
                    .filter(|(_, bond)| bond.order == BondOrder::Aromatic)
                    .count(),
                0,
                "{name}: {smi}"
            );
        }
    }

    #[test]
    fn issue_767_aromaticity_controls_remain_unchanged() {
        let cases = [
            ("thiophene", "C1=CC=CS1", 5),
            ("thiophene 1-oxide", "C1=CS(C=C1)=O", 0),
            ("thiophene 1,1-dioxide", "C1=CS(=O)(=O)C=C1", 0),
            ("furan", "C1=CC=CO1", 5),
            ("N-methylpyrrole", "Cn1cccc1", 5),
            ("N-methylpyrrole N-oxide", "C1=C[N+](C)([O-])C=C1", 0),
            ("pyridine", "C1=CC=NC=C1", 6),
            ("pyridine N-oxide", "C1=C[N+]([O-])=CC=C1", 6),
            ("isoxazole", "Cc1cc(N)no1", 5),
        ];

        for (name, smi, expected) in cases {
            let mol = chematic_smiles::parse(smi).expect("valid control SMILES");
            let applied = apply_aromaticity_rdkit_parity_experimental(&mol)
                .expect("control must remain perceivable");
            assert_eq!(
                applied.atoms().filter(|(_, atom)| atom.aromatic).count(),
                expected,
                "{name}: {smi}"
            );
        }
    }

    #[test]
    fn te_default_valence_and_aromatic_bond_flags() {
        // Direct check of the dependent-path fix: default_valence(52) must now
        // resolve to Some(2) (was None before the chematic-core element.rs
        // Te normal_valences() fix).
        assert_eq!(default_valence(52), Some(2));

        // The calibration battery above only checks the aromatic ATOM set for
        // tellurophene. RDKit also reports tellurophene's ring BONDS as
        // aromatic (bond order 1.5 on every ring bond, oracle-verified via
        // rdkit.Chem.MolFromSmiles("C1=C[Te]C=C1")); confirm chematic's
        // aromatic bond set matches too, not just the atom set.
        let mol = mol_kekulized("C1=C[Te]C=C1");
        let (atoms, bonds) = rdkit_parity_aromaticity(&mol);
        assert_eq!(
            atoms.len(),
            5,
            "all 5 tellurophene ring atoms must be aromatic"
        );
        assert_eq!(
            bonds.len(),
            5,
            "all 5 tellurophene ring bonds must be aromatic"
        );
    }

    // Purine's Aromaticity-A1-0 finding was that production's answer depends
    // on whether the input was Kekulized before `apply_aromaticity` ran.
    // rdkit_parity_aromaticity must NOT reintroduce that: both a raw
    // aromatic-lowercase parse (kekulized here identically to every other
    // corpus entry, so this mostly re-confirms `mol_kekulized`'s own
    // determinism) and chematic's own `kekulize()` choice must agree with
    // each other and with RDKit.
    #[test]
    fn purine_representation_stable() {
        let smi = "c1cnc2[nH]cnc2n1";
        let raw = chematic_smiles::parse(smi).expect("valid SMILES");
        let k = chematic_core::kekulize(&raw).expect("purine should kekulize");
        let via_own_kekulize = chematic_core::apply_kekule(&raw, &k);

        let (atoms_a, _) = rdkit_parity_aromaticity(&mol_kekulized(smi));
        let (atoms_b, _) = rdkit_parity_aromaticity(&via_own_kekulize);

        let mut a: Vec<u32> = atoms_a.iter().map(|x| x.0).collect();
        let mut b: Vec<u32> = atoms_b.iter().map(|x| x.0).collect();
        a.sort();
        b.sort();
        assert_eq!(a, b, "purine: two Kekulization paths disagree");
        assert_eq!(
            a,
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8],
            "purine: should match RDKit (all 9 atoms aromatic)"
        );
    }

    #[test]
    fn porphyrin_like_fused_system_matches_rdkit_aromatic_atom_set() {
        // RDKit 2025.09.3 and the official RDKit.js 2026.03.6 package agree
        // that the two exocyclic vinyl branches are not aromatic. This is one
        // of the retained fixed-10k Morgan residuals; retain the atom-level
        // oracle here so an eventual fingerprint fix cannot mask a wrong
        // aromaticity partition behind folded-bit collisions.
        let smi =
            "CC1=C2NC(=C1CCC(O)=O)C=C3N=C(C=C4NC(=CC5=NC(=C2)C(=C5C)C=C)C(=C4C)C=C)C(=C3CCC(O)=O)C";
        let parsed = chematic_smiles::parse(smi).expect("valid porphyrin-like SMILES");
        let applied =
            apply_aromaticity_rdkit_parity_experimental(&parsed).expect("aromaticity perception");
        let aromatic: Vec<u32> = applied
            .atoms()
            .filter_map(|(idx, atom)| atom.aromatic.then_some(idx.0))
            .collect();
        assert_eq!(
            aromatic,
            vec![
                1, 2, 3, 4, 5, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 29, 30,
            ],
            "must retain RDKit's aromatic partition before Morgan expansion"
        );
    }

    #[test]
    fn fullerene_kekule_copy_uses_rdkit_symmetrized_rings() {
        // RDKit 2026.03.6 perceives over its symmetrized SSSR; on this
        // Kekule-written fullerene adduct one hexagon only becomes aromatic
        // through a ring `symmetrizeSSSR` adds.
        let smi = "COCCOCCOCCN1CC23C4=C5C6=C7C8=C9C%10=C%11C%12=C%13C%14=C%15C%16=C%17C%14=C%14C%12=C%12C%10=C%10C%18=C%19C%12=C%14C%12=C%17C%14C%17=C%12C%19=C%12C%18=C(C6=C9%10)C5C5C%12=C%17C(=C6C9=C2C2C4=C7C4=C8C%11=C%13C7=C4C2C(=C%157)C9=C%16C6%14)C53C1COCCOCCOC";
        let parsed = chematic_smiles::parse(smi).expect("parses");
        let applied = apply_aromaticity_rdkit_parity_experimental(&parsed).expect("perceives");
        let aromatic: Vec<u32> = applied
            .atoms()
            .filter_map(|(idx, atom)| atom.aromatic.then_some(idx.0))
            .collect();
        assert_eq!(
            aromatic,
            vec![
                13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,
                34, 35, 36, 37, 38, 40, 41, 42, 43, 44, 45, 46, 47, 50, 51, 57, 58, 59, 60, 61, 62,
                63, 64, 66, 67, 68, 69
            ]
        );
    }

    #[test]
    fn production_api_assign_matches_engine_on_benzene() {
        let mol = chematic_smiles::parse("c1ccccc1").expect("valid SMILES");
        let model = assign_aromaticity_rdkit_parity_experimental(&mol).expect("benzene kekulizes");
        assert_eq!(model.aromatic_atom_count(), 6);
        for (idx, _) in mol.atoms() {
            assert!(
                model.is_atom_aromatic(idx),
                "atom {idx:?} should be aromatic"
            );
        }
        // This engine only determines atom/bond sets, not per-ring
        // classification or antiaromaticity -- both must be empty.
        assert!(model.ring_classifications().is_empty());
        assert!(model.antiaromatic_rings().is_empty());
    }

    #[test]
    fn production_api_apply_sets_aromatic_flags_and_bond_orders() {
        let mol = chematic_smiles::parse("C1=CC=CC=C1").expect("valid SMILES"); // Kekule benzene
        let applied = apply_aromaticity_rdkit_parity_experimental(&mol).expect("benzene kekulizes");
        assert_eq!(applied.atom_count(), mol.atom_count());
        for (_, atom) in applied.atoms() {
            assert!(atom.aromatic, "every benzene atom should end up aromatic");
        }
        for (_, bond) in applied.bonds() {
            assert_eq!(bond.order, BondOrder::Aromatic);
        }
    }

    #[test]
    fn mixed_aromatic_morphine_promotes_bridge_oxygen_like_rdkit() {
        // RDKit canonicalizes this mixed aromatic/Kekule spelling to
        // `CN1CCc2oc3c(O)ccc4c3c2C1C4`: the degree-two bridge oxygen is
        // aromatic. Preserving only the parser-supplied `c1ccc...` subgraph
        // missed that promotion and shifted both TPSA and Crippen LogP.
        let mol =
            chematic_smiles::parse("Oc1ccc2CC3N(CCC4=C3c2c1O4)C").expect("valid morphine SMILES");
        let applied =
            apply_aromaticity_rdkit_parity_experimental(&mol).expect("morphine kekulizes");
        let aromatic_oxygen_count = applied
            .atoms()
            .filter(|(idx, atom)| {
                atom.element.atomic_number() == 8 && atom.aromatic && applied.degree(*idx) == 2
            })
            .count();
        assert_eq!(
            aromatic_oxygen_count, 1,
            "RDKit marks exactly the bridge oxygen aromatic"
        );
    }

    #[test]
    fn mixed_fused_cyclic_ether_does_not_promote_oxygen() {
        // Negative control for the morphine regression: an aromatic neighbour
        // plus a vinylic neighbour is not sufficient to make a cyclic ether
        // oxygen aromatic.
        let mol = chematic_smiles::parse("CC(=O)c1c(O)c(C)c(O)c2c1OC1=Cc3c(c(C)nn3C)C(=O)C12C")
            .expect("valid fused cyclic ether SMILES");
        let applied = apply_aromaticity_rdkit_parity_experimental(&mol)
            .expect("fused cyclic ether kekulizes");
        assert!(
            applied
                .atoms()
                .filter(|(_, atom)| atom.element.atomic_number() == 8)
                .all(|(_, atom)| !atom.aromatic),
            "RDKit keeps every oxygen in this fused cyclic ether non-aromatic"
        );
    }

    #[test]
    fn production_api_preserves_r_group_sidecar() {
        use chematic_core::{MoleculeBuilder, RGroupLabel};

        let parsed = chematic_smiles::parse("c1ccccc1.[*]").expect("valid disconnected SMILES");
        let mut builder = MoleculeBuilder::from_molecule(&parsed);
        let r1 = AtomIdx(6);
        builder.set_r_group(r1, RGroupLabel::numbered(1).unwrap());
        let mol = builder.build();

        let applied = apply_aromaticity_rdkit_parity_experimental(&mol)
            .expect("aromaticity perception preserves pseudoatom metadata");
        assert_eq!(applied.r_group_label(r1), mol.r_group_label(r1));
    }

    #[test]
    fn production_api_preserves_valid_explicit_aromatic_input() {
        // RDKit accepts this fused purine-like graph even though chematic's
        // matching-based Kekule conversion cannot represent it. The explicit
        // aromatic input is self-consistent, so it must be preserved rather
        // than rejected or silently sent through the ordinary Hückel model.
        let smi = "Cc1cn2c(=O)c3ncn(COCCO)c3nc2n1C";
        let mol = chematic_smiles::parse(smi).expect("valid SMILES");

        let model = assign_aromaticity_rdkit_parity_experimental(&mol).expect("valid input");
        assert!(model.aromatic_atom_count() > 0);
        let applied = apply_aromaticity_rdkit_parity_experimental(&mol).expect("valid input");
        assert_eq!(applied.atom_count(), mol.atom_count());
    }

    #[test]
    fn preperceived_fast_path_keeps_disconnected_saturated_ring_non_aromatic() {
        let mol = chematic_smiles::parse("c1ccccc1.C1CCCCC1").expect("valid SMILES");
        let applied = apply_aromaticity_rdkit_parity_experimental(&mol).expect("valid aromaticity");
        assert_eq!(applied.atoms().filter(|(_, atom)| atom.aromatic).count(), 6);
        assert_eq!(
            applied
                .bonds()
                .filter(|(_, bond)| bond.order == BondOrder::Aromatic)
                .count(),
            6
        );
    }

    #[test]
    fn kekule_ring_bypasses_fast_path_and_is_still_perceived() {
        let mol = chematic_smiles::parse("c1ccccc1.C1=CC=CC=C1").expect("valid SMILES");
        let applied = apply_aromaticity_rdkit_parity_experimental(&mol).expect("valid aromaticity");
        assert_eq!(
            applied.atoms().filter(|(_, atom)| atom.aromatic).count(),
            12,
            "the separate Kekulé benzene must still use full perception"
        );
        assert_eq!(
            applied
                .bonds()
                .filter(|(_, bond)| bond.order == BondOrder::Aromatic)
                .count(),
            12
        );
    }

    /// The full kekulize-and-perceive path, with no shortcut.
    fn full_path(mol: &Molecule) -> Result<Molecule, AromaticityError> {
        let (kekulized, verdict) = kekulize_for_rdkit_parity_with_verdict(mol)?;
        if kekulized
            .bonds()
            .any(|(_, bond)| bond.order == BondOrder::Aromatic)
        {
            return Ok(kekulized);
        }
        let model = match verdict {
            Some(v) => assign_from_verdict(&kekulized, v)?,
            None => assign_from_kekulized(&kekulized)?,
        };
        Ok(crate::aromaticity::build_molecule_from_model(
            &kekulized, &model,
        ))
    }

    fn fingerprint(mol: &Molecule) -> String {
        let atoms: String = mol
            .atoms()
            .map(|(i, a)| {
                format!(
                    "{}{}{}{:?}{},",
                    a.element.symbol(),
                    a.charge,
                    a.aromatic as u8,
                    a.hydrogen_count,
                    chematic_core::implicit_hcount(mol, i)
                )
            })
            .collect();
        let bonds: String = mol
            .bonds()
            .map(|(_, b)| format!("{}-{}:{:?},", b.atom1.0, b.atom2.0, b.order))
            .collect();
        format!("{atoms}|{bonds}")
    }

    #[test]
    fn kekulized_copy_matches_the_rebuilds() {
        for smi in [
            "c1ccccc1",
            "c1ccc2[nH]ccc2c1",
            "C/C=C/c1ccccc1",
            "[C@H](F)(Cl)c1ccccc1",
            "[*]c1ccccc1",
            "[*:1]c1ccc(-c2ccccn2)cc1",
            "CCO",
        ] {
            let mut mol = chematic_smiles::parse(smi).expect("valid SMILES");
            if mol.atom(AtomIdx(0)).wildcard {
                // An R-group label (dropped by apply_kekule's rebuild).
                mol.set_r_group(AtomIdx(0), chematic_core::RGroupLabel::unnumbered());
            }
            let k = chematic_core::kekulize(&mol).expect("kekulizable");
            let a = kekulized_copy(&mol, &k);
            let b = chematic_core::apply_kekule(&clear_aromatic_flags(&mol), &k);
            assert_eq!(fingerprint(&a), fingerprint(&b), "{smi}");
            assert_eq!(
                chematic_smiles::canonical_smiles(&a),
                chematic_smiles::canonical_smiles(&b),
                "{smi}"
            );
            for i in 0..a.atom_count() {
                let idx = AtomIdx(i as u32);
                assert_eq!(
                    a.neighbors(idx).collect::<Vec<_>>(),
                    b.neighbors(idx).collect::<Vec<_>>(),
                    "{smi}"
                );
                assert_eq!(
                    a.stereo_neighbor_order(idx),
                    b.stereo_neighbor_order(idx),
                    "{smi}"
                );
                assert_eq!(a.r_group_label(idx), b.r_group_label(idx), "{smi}");
            }
            for bi in 0..a.bond_count() {
                let bidx = BondIdx(bi as u32);
                assert_eq!(a.bond_direction(bidx), b.bond_direction(bidx), "{smi}");
            }
        }
    }

    #[test]
    fn shortcuts_agree_with_the_full_path() {
        for smi in [
            // explicit aromatic, pure
            "c1ccccc1",
            "c1ccc2[nH]ccc2c1",
            // mixed components (explicit + Kekulé/saturated in one ring system)
            "O=C1c2ccccc2C(=O)N1CCCCN1CCN(c2cccc3ccccc23)CC1",
            "CSc1ccc2c(c1)N(CCC1CCCCN1C)c1ccccc1S2",
            "c1ccc2c(c1)CCC2",
            "O=C1NCCc2ccccc21",
            "C1=Cc2ccccc2C1",
            "c1ccccc1.C1=CC=CC=C1",
            // non-aromatic components next to aromatic ones
            "O=C1NC(=O)/C(=C/c2ccccc2)N1",
            "O=C1C=CNC(=O)N1Cc1ccccc1",
            "C1=CC=CC=C1c1ccccc1",
            "O=C1C=CC(=O)C=C1c1ccccc1",
            "C1=CNC=C1c1ccccc1",
            "O=C1NNC(=O)N1c1ccccc1",
            "S=C1NNC(=S)S1",
            // no aromaticity at all
            "CC12CCC3C(CCC4=CC(=O)CCC34C)C1CCC2=O",
            "C1=CC=CC=C1",
            "C1=CNC=C1",
            "O=C1C=CC(=O)C=C1",
            "O=C1NNC(=O)N1",
            "C1CCNCC1",
            "C1=CCC=C1C1=CCCC1",
            // flags without aromatic bonds / inconsistent explicit forms
            "[cH]1[cH][cH][cH][cH][cH]1",
        ] {
            let mol = chematic_smiles::parse(smi).expect("valid SMILES");
            let full = full_path(&mol);
            let fast = apply_aromaticity_rdkit_parity_uncached(&mol);
            match (&full, &fast) {
                (Ok(a), Ok(b)) => assert_eq!(fingerprint(a), fingerprint(b), "{smi}"),
                (Err(a), Err(b)) => assert_eq!(a, b, "{smi}"),
                _ => panic!(
                    "{smi}: full ok={} vs fast ok={}",
                    full.is_ok(),
                    fast.is_ok()
                ),
            }
            if rdkit_parity_view_is_identity(&mol) {
                assert_eq!(
                    fingerprint(&mol),
                    fingerprint(full.as_ref().unwrap()),
                    "{smi}"
                );
            }
        }
    }

    #[test]
    fn kekule_verdict_identity_is_sound() {
        let mut proven = 0;
        for smi in [
            "c1ccccc1",
            "c1ccc2[nH]ccc2c1",
            "c1ccc2ccccc2c1",
            "c1ccc2c(c1)ccc1ccccc12",
            "c1cc2ccc3cccc4ccc(c1)c2c34",
            "c1ccc2cc3ccccc3cc2c1",
            "c1ccc-2c(c1)-c1ccccc1-2",
            "c1cc2cccc3cccc(c1)c23",
            "c1ccc2c(c1)c1ccccc1n2C",
            "O=c1cccc[nH]1",
            "O=c1ccc2ccccc2o1",
            "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
            "c1cc[n+]cc1",
            "c1ccn2ccnc2c1",
            "c1cnc2nccn2c1",
            "c1ccc2c(c1)[se]c1ccccc12",
            "c1cc2ccc1cc2",
            "c1cccc1",
            "c1ccc1",
            "c1cccccc1",
            "c1ccccccc1",
            "c1ccc2cccc2cc1",
            "C1=Cc2ccccc2C1",
            "c1ccc2c(c1)CCC2",
            "O=C1c2ccccc2C(=O)N1CCCCN1CCN(c2cccc3ccccc23)CC1",
            "c1ccc2c(c1)[nH]c1ccccc12",
            "c1cc2c3c(c1)ccc1cccc(c13)cc2",
            "B1c2ccccc2-c2ccccc12",
            "c1ccc2c(c1)oc1ccccc12",
            "O=c1c2ccccc2c(=O)c2ccccc12",
            "c12c3c4c5c1c1c6c7c2c2c8c3c3c9c4c4c%10c5c5c1c1c6c6c%11c7c2c2c7c8c3c3c8c9c4c4c9c%10c5c5c1c1c6c6c%11c2c2c7c3c3c8c4c4c9c5c1c1c6c2c3c41",
        ] {
            let mol = chematic_smiles::parse(smi).expect("valid SMILES");
            let info = ComponentInfo::new(&mol);
            if !explicit_kept_by_kekule_verdict(&mol, &info) {
                continue;
            }
            proven += 1;
            match full_path(&mol) {
                Ok(full) => assert_eq!(fingerprint(&mol), fingerprint(&full), "{smi}"),
                Err(e) => panic!("{smi}: full path failed: {e:?}"),
            }
        }
        assert!(proven > 10, "only {proven} proven");
    }

    #[test]
    fn ring_with_exocyclic_multiple_bond_bypasses_fast_path() {
        for smiles in ["O=C1NNC(=O)N1", "CN=C1SSC(=O)N1C", "S=C1NNC(=S)S1"] {
            let mol = chematic_smiles::parse(smiles).expect("valid SMILES");
            assert!(
                complete_preperceived_aromaticity(&mol).is_none(),
                "exocyclic unsaturation can participate in RDKit aromaticity: {smiles}"
            );
            apply_aromaticity_rdkit_parity_experimental(&mol)
                .expect("full aromaticity perception succeeds");
        }
    }

    #[test]
    fn production_api_does_not_mutate_input_on_explicit_aromatic_path() {
        // "元の分子を途中まで書き換えてから失敗する経路は作らないでください" --
        // `mol` is only ever taken by `&Molecule` throughout the fallible
        // path (`clear_aromatic_flags`/`kekulize`/`apply_kekule` all build
        // new molecules rather than mutating in place), so this is enforced
        // by the type system. Pin it as an explicit regression: the input's
        // own atom/bond flags and counts are unchanged after a failed call.
        let smi = "Cc1cn2c(=O)c3ncn(COCCO)c3nc2n1C";
        let mol = chematic_smiles::parse(smi).expect("valid SMILES");
        let atom_count_before = mol.atom_count();
        let bond_count_before = mol.bond_count();
        let aromatic_before: Vec<bool> = mol.atoms().map(|(_, a)| a.aromatic).collect();

        let result = assign_aromaticity_rdkit_parity_experimental(&mol);
        assert!(result.is_ok(), "self-consistent explicit aromatic input");

        assert_eq!(mol.atom_count(), atom_count_before);
        assert_eq!(mol.bond_count(), bond_count_before);
        let aromatic_after: Vec<bool> = mol.atoms().map(|(_, a)| a.aromatic).collect();
        assert_eq!(aromatic_before, aromatic_after);
    }

    /// Was a known kekulize-gap case (found during Morgan M4-A0, `chematic-fp`'s
    /// `rdkit_morgan_hash.rs`, full-corpus validation): pyridinium's protonated
    /// `[nH+]` used to make `chematic_core::kekulize()` hard-fail because
    /// `atom_must_be_matched`'s N-with-H lone-pair-donor rule was charge-blind
    /// (docs/rfcs/aromaticity_rdkit_parity_rfc.md §1, root cause A). Fixed by
    /// `fix/kekulize-charge-aware-k1` (see
    /// docs/rfcs/kekulize_charge_aware_rdkit_parity.md): the rule now requires
    /// `atom.charge <= 0`, so a protonated ring N routes back to "must be
    /// matched" -- same as neutral pyridine's bare N -- instead of being
    /// wrongly treated like neutral pyrrole's `[nH]`. Kept as a regression
    /// test (not deleted) so a future re-introduction of the charge-blind rule
    /// is caught here, not just in the 40-fixture diagnosis corpus.
    ///
    /// NOTE for whoever picks up the companion fp-side fix: `chematic-fp`'s
    /// `rdkit_morgan_ecfp4.rs` has two tests
    /// (`kekule_pyridinium_reports_kekulization_failed_not_a_fallback_result`,
    /// `hueckel_fallback_would_be_detectable_if_silently_reintroduced`) that
    /// also used this exact SMILES as their "kekulize fails" positive control
    /// for a *different* invariant (the fallible ECFP4 path must not silently
    /// fall back to Hueckel) -- those are out of scope for K1 (chematic-fp is
    /// off limits for this fix) and still fail as of this commit. They need
    /// the same swap this test got: same invariant, a still-failing molecule
    /// as the example (`Cc1cn2c(=O)c3ncn(COCCO)c3nc2n1C`, see
    /// `production_api_reports_kekulize_failure_not_panic` above), not a
    /// deletion. `validation/README.md`'s Morgan M4-A0 section and several
    /// `validation/*.json` artifacts also cite a now-stale "2 of 5,048"
    /// preprocessing-failure count (only the purine molecule remains).
    #[test]
    fn kekulize_charge_aware_k1_fixes_protonated_pyridinium() {
        let smi = "c1cc[nH+]cc1";
        let mol = chematic_smiles::parse(smi).expect("valid SMILES");
        let result = apply_aromaticity_rdkit_parity_experimental(&mol);
        match result {
            Ok(applied) => {
                for (_, atom) in applied.atoms() {
                    assert!(
                        atom.aromatic,
                        "pyridinium's ring must end up fully aromatic"
                    );
                }
            }
            Err(other) => {
                panic!("expected pyridinium to kekulize successfully post-K1, got {other:?}")
            }
        }
    }

    #[test]
    fn production_api_stale_aromatic_flag_is_overridden_not_leaked() {
        // A non-aromatic atom that happens to carry a stale `aromatic=true`
        // flag on input must not leak that flag into the output -- the
        // engine's own verdict is the sole source of truth for the result.
        use chematic_core::MoleculeBuilder;

        let mut base = chematic_smiles::parse("CC").expect("valid SMILES"); // ethane, acyclic
        let mut builder = MoleculeBuilder::new();
        for (_, atom) in base.atoms() {
            let mut a = atom.clone();
            a.aromatic = true; // stale/bogus annotation
            builder.add_atom(a);
        }
        for (_, bond) in base.bonds() {
            let _ = builder.add_bond(bond.atom1, bond.atom2, bond.order);
        }
        base = builder.build();
        assert!(base.atoms().all(|(_, a)| a.aromatic), "test setup sanity");

        let applied = apply_aromaticity_rdkit_parity_experimental(&base)
            .expect("acyclic molecule kekulizes trivially (no-op)");
        assert!(
            applied.atoms().all(|(_, a)| !a.aromatic),
            "stale aromatic=true on an acyclic atom must not survive"
        );
    }

    #[test]
    fn neutral_oxygen_and_sulfur_are_not_aromatic_macrocycle_donors() {
        // RDKit Aromaticity.cpp suppresses neutral divalent O/S two-electron
        // donation when the containing SSSR ring has at least nine atoms.
        for smiles in ["o1cccccccc1", "s1cccccccc1", "o1cccccccccccc1"] {
            let mol = chematic_smiles::parse(smiles).expect("valid macrocycle");
            let perceived = apply_aromaticity_rdkit_parity_experimental(&mol)
                .expect("macrocycle can be kekulized");
            assert!(
                perceived.atoms().all(|(_, atom)| !atom.aromatic),
                "{smiles} must be non-aromatic in RDKit's model"
            );
            assert!(
                perceived
                    .bonds()
                    .all(|(_, bond)| bond.order != BondOrder::Aromatic),
                "{smiles} must have a Kekule representation"
            );
        }

        for smiles in ["o1cccc1", "s1cccc1"] {
            let mol = chematic_smiles::parse(smiles).expect("valid five-membered ring");
            let perceived = apply_aromaticity_rdkit_parity_experimental(&mol)
                .expect("five-membered aromatic ring is supported");
            assert!(
                perceived.atoms().all(|(_, atom)| atom.aromatic),
                "{smiles} must remain aromatic"
            );
        }
    }
}
