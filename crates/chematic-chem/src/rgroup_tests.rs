use super::*;
use chematic_smiles::{parse, random_smiles};

fn mol(s: &str) -> Molecule {
    parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"))
}

fn decompose<'a>(
    scaffold: &str,
    smiles: impl IntoIterator<Item = &'a str>,
) -> Vec<Option<RGroupResult>> {
    let mols: Vec<Molecule> = smiles.into_iter().map(mol).collect();
    let refs: Vec<&Molecule> = mols.iter().collect();
    rgroup_decompose(scaffold, &refs).expect("decompose failed")
}

#[test]
fn test_rgroup_monosubstituted_benzene() {
    // Three alkylbenzenes — scaffold = benzene ring (6 C, no attachment markers).
    let results = decompose("c1ccccc1", ["Cc1ccccc1", "CCc1ccccc1", "CCCc1ccccc1"]);
    assert_eq!(results.len(), 3);
    for (i, r) in results.iter().enumerate() {
        let r = r.as_ref().unwrap_or_else(|| panic!("mol {i} should match"));
        assert_eq!(r.mol_idx, i);
        // Must have exactly one R-group (one substituted position).
        assert_eq!(r.r_groups.len(), 1, "mol {i}: expected 1 R-group");
        let r1 = r.r_groups.get(&1).expect("R1 missing");
        // R1 must contain a wildcard marker and the appropriate alkyl chain.
        assert!(
            r1.contains('*'),
            "R1 must contain [*] attachment, got {r1:?}"
        );
    }
}

#[test]
fn test_rgroup_mol_idx_preserved() {
    let results = decompose("c1ccccc1", ["Cc1ccccc1", "CCCC", "CCCc1ccccc1"]);
    // CCCC does not match benzene scaffold → None
    assert!(results[0].is_some());
    assert!(results[1].is_none());
    assert!(results[2].is_some());
    assert_eq!(results[0].as_ref().unwrap().mol_idx, 0);
    assert_eq!(results[2].as_ref().unwrap().mol_idx, 2);
}

#[test]
fn test_rgroup_no_substituents() {
    // Benzene itself: scaffold matches, zero R-groups.
    let results = decompose("c1ccccc1", ["c1ccccc1"]);
    let r = results[0]
        .as_ref()
        .expect("benzene should match its own scaffold");
    assert_eq!(r.r_groups.len(), 0, "benzene has no R-groups vs itself");
}

#[test]
fn test_rgroup_two_substituents() {
    // para-xylene: two methyl groups → two R-groups.
    let results = decompose("c1ccccc1", ["Cc1ccc(C)cc1"]);
    let r = results[0].as_ref().expect("should match");
    assert_eq!(r.r_groups.len(), 2, "para-xylene has two R-groups");
}

#[test]
fn test_rgroup_invalid_smarts() {
    let mols = [mol("c1ccccc1")];
    let refs: Vec<&Molecule> = mols.iter().collect();
    let err = rgroup_decompose("((invalid", &refs);
    assert!(matches!(err, Err(RGroupError::InvalidSmarts(_))));
}

#[test]
fn test_rgroup_different_heteroatom_substituents() {
    // Aniline and phenol vs benzene scaffold.
    let results = decompose("c1ccccc1", ["Nc1ccccc1", "Oc1ccccc1"]);
    let aniline = results[0].as_ref().expect("aniline matches");
    let phenol = results[1].as_ref().expect("phenol matches");
    assert_eq!(aniline.r_groups.len(), 1);
    assert_eq!(phenol.r_groups.len(), 1);
    // Both R1 values must start with N or O (amine / hydroxyl).
    let r1_a = aniline.r_groups.get(&1).unwrap();
    let r1_p = phenol.r_groups.get(&1).unwrap();
    assert!(
        r1_a.starts_with('N') || r1_a.contains('N'),
        "aniline R1 should contain N, got {r1_a}"
    );
    assert!(
        r1_p.starts_with('O') || r1_p.contains('O'),
        "phenol R1 should contain O, got {r1_p}"
    );
}

#[test]
fn multiple_bond_rgroups_retain_the_cut_bond_order() {
    // R-group strings and core graphs checked against RDKit 2026.09.1.
    let mols = [mol("O=Cc1ccccc1"), mol("N#Cc1ccccc1")];
    let refs: Vec<&Molecule> = mols.iter().collect();

    let labeled = rgroup_decompose_labeled("c1ccccc1[C:1]", &refs).unwrap();
    let expected_core = canonical_smiles(&mol("c1ccc(C[*:1])cc1"));
    assert_eq!(labeled[0].as_ref().unwrap().core_smiles, expected_core);
    assert_eq!(labeled[1].as_ref().unwrap().core_smiles, expected_core);
    assert_eq!(labeled[0].as_ref().unwrap().r_groups[&1], "O=[*:1]");
    assert_eq!(labeled[1].as_ref().unwrap().r_groups[&1], "N#[*:1]");

    let unlabeled = rgroup_decompose("c1ccccc1C", &refs).unwrap();
    for (row, expected_order) in unlabeled.iter().zip([BondOrder::Double, BondOrder::Triple]) {
        let fragment = mol(&row.as_ref().unwrap().r_groups[&1]);
        assert_eq!(fragment.bond(BondIdx(0)).order, expected_order);
    }
}

#[test]
fn alkene_stereo_survives_a_cut_on_the_rgroup_side() {
    // Source-built RDKit 2026.09.1 emits these same R1 structures. The
    // boundary slash is needed to distinguish the two alkene isomers.
    for (target, expected) in [
        ("c1ccccc1/C=C/Cl", "Cl/C=C/[*:1]"),
        ("c1ccccc1/C=C\\Cl", "Cl/C=C\\[*:1]"),
    ] {
        let target = mol(target);
        let rows = rgroup_decompose_labeled("c1ccccc1", &[&target]).unwrap();
        assert_eq!(
            rows[0].as_ref().unwrap().r_groups[&1],
            expected,
            "R-group alkene geometry must survive the cut"
        );
    }
}

#[test]
fn alkene_stereo_survives_a_cut_on_the_core_side() {
    // Both alkene substituents are outside the C=C query. Their boundary
    // directions are still required to identify the core's E/Z geometry.
    let mut core_smiles = Vec::new();
    for target in ["F/C=C/Cl", "F/C=C\\Cl"] {
        let target = mol(target);
        let rows = rgroup_decompose_labeled("C=C", &[&target]).unwrap();
        core_smiles.push(rows[0].as_ref().unwrap().core_smiles.clone());
    }
    assert_ne!(core_smiles[0], core_smiles[1]);
    assert!(core_smiles.iter().all(|smiles| smiles.contains('/')));
}

#[test]
fn alkene_core_distinguishes_ez_isomers_across_substituent_pairs() {
    // Every pair has the same constitution and two independently specified
    // alkene geometries. R-group extraction must retain that distinction.
    let substituents = ["F", "Cl", "Br", "I", "C", "CO"];
    for (left_index, left) in substituents.iter().enumerate() {
        for right in &substituents[left_index + 1..] {
            let first = format!("{left}/C=C/{right}");
            let second = format!("{left}/C=C\\{right}");
            let molecules = [mol(&first), mol(&second)];
            let refs = molecules.iter().collect::<Vec<_>>();
            let rows = rgroup_decompose_labeled("C=C", &refs).unwrap();
            let first_row = rows[0].as_ref().expect("first geometry must match");
            let second_row = rows[1].as_ref().expect("second geometry must match");
            assert_eq!(
                first_row.r_groups, second_row.r_groups,
                "{first} / {second}"
            );
            assert_ne!(
                first_row.core_smiles, second_row.core_smiles,
                "{first} / {second}"
            );
        }
    }
}

#[test]
fn tetrahedral_stereo_keeps_its_neighbor_frame_after_cutting() {
    // A stereocenter need not touch the boundary to change handedness if
    // its parser-side neighbor frame is discarded during atom renumbering.
    for (target, expected) in [
        ("CC[C@H](F)c1ccccc1", "CC[C@H](F)[*:1]"),
        ("CC[C@@H](F)c1ccccc1", "CC[C@@H](F)[*:1]"),
        ("C[C@H](F)CCc1ccccc1", "C[C@H](F)CC[*:1]"),
    ] {
        let target = mol(target);
        let rows = rgroup_decompose_labeled("c1ccccc1", &[&target]).unwrap();
        let fragment = &rows[0].as_ref().unwrap().r_groups[&1];
        assert_eq!(
            canonical_smiles(&mol(fragment)),
            canonical_smiles(&mol(expected))
        );
    }
}

#[test]
fn stereo_rgroups_are_stable_under_input_atom_order() {
    for (target, expected) in [
        ("c1ccccc1/C=C/Cl", "Cl/C=C/[*:1]"),
        ("c1ccccc1/C=C\\Cl", "Cl/C=C\\[*:1]"),
        ("CC[C@H](F)c1ccccc1", "CC[C@H](F)[*:1]"),
        ("C[C@H](F)CCc1ccccc1", "C[C@H](F)CC[*:1]"),
    ] {
        let target = mol(target);
        let expected = canonical_smiles(&mol(expected));
        for seed in 0..12 {
            let reordered = mol(&random_smiles(&target, seed));
            let rows = rgroup_decompose_labeled("c1ccccc1", &[&reordered]).unwrap();
            let fragment = &rows[0].as_ref().unwrap().r_groups[&1];
            assert_eq!(canonical_smiles(&mol(fragment)), expected, "seed {seed}");
        }
    }
}

#[test]
fn explicit_multiple_bond_placeholder_keeps_core_and_rgroup_orders() {
    for (scaffold, target, order, expected_core, expected_rgroup) in [
        (
            "c1ccccc1C=[*:1]",
            "O=Cc1ccccc1",
            BondOrder::Double,
            "C(c1ccccc1)=[*:1]",
            "O=[*:1]",
        ),
        (
            "c1ccccc1C#[*:1]",
            "N#Cc1ccccc1",
            BondOrder::Triple,
            "C(c1ccccc1)#[*:1]",
            "N#[*:1]",
        ),
    ] {
        let molecule = mol(target);
        let rows = rgroup_decompose_labeled(scaffold, &[&molecule]).unwrap();
        let row = rows[0].as_ref().expect("multiple-bond core should match");
        assert_eq!(row.core_smiles, canonical_smiles(&mol(expected_core)));
        assert_eq!(row.r_groups[&1], expected_rgroup);
        let core = mol(&row.core_smiles);
        let rgroup = mol(&row.r_groups[&1]);
        assert!(core.bonds().any(|(_, bond)| bond.order == order));
        assert_eq!(rgroup.bond(BondIdx(0)).order, order);
    }
}

#[test]
fn bridged_rgroups_keep_both_core_attachments() {
    // RDKit 2026.09.1 emits the same two-label fragment in R1 and R2.
    for (target, expected) in [
        ("c1ccc2c(c1)CCC2", "C(C[*:1])C[*:2]"),
        ("c1ccc2c(c1)OCC2", "C(C[*:1])O[*:2]"),
        ("c1ccc2ccccc2c1", "c(ccc:[*:2]):[*:1]"),
    ] {
        let molecule = mol(target);
        let rows = rgroup_decompose_labeled("c1ccccc1", &[&molecule]).unwrap();
        let row = rows[0].as_ref().expect("fused core should match");
        let expected = canonical_smiles(&mol(expected));
        assert_eq!(row.r_groups.len(), 2);
        assert_eq!(row.r_groups[&1], expected);
        assert_eq!(row.r_groups[&2], expected);
        assert!(row.core_smiles.contains("*:1"));
        assert!(row.core_smiles.contains("*:2"));
    }
}

#[test]
fn unlabeled_bridged_rgroups_keep_two_attachment_stubs() {
    let rows = decompose("c1ccccc1", ["c1ccc2c(c1)CCC2"]);
    let row = rows[0].as_ref().expect("fused core should match");
    assert_eq!(row.r_groups.len(), 2);
    for fragment in row.r_groups.values() {
        assert_eq!(fragment.matches('*').count(), 2);
    }
}

#[test]
fn bridged_rgroup_labels_ignore_input_atom_order() {
    let original = mol("c1ccc2c(c1)OCC2");
    let respelled = mol("c1ccc2OCCc2c1");
    let first = rgroup_decompose_labeled("c1ccccc1", &[&original]).unwrap();
    let second = rgroup_decompose_labeled("c1ccccc1", &[&respelled]).unwrap();
    assert_eq!(first, second);
}

#[test]
fn labelled_dummy_atoms_align_symmetric_columns_like_rdkit() {
    let mols = [mol("Cc1ccc(CC)cc1"), mol("CCc1ccc(N)cc1")];
    let refs: Vec<&Molecule> = mols.iter().collect();
    let rows =
        rgroup_decompose_labeled("c1cc([*:1])ccc1[*:2]", &refs).expect("labelled decomposition");
    let first = rows[0].as_ref().expect("first molecule matches");
    let second = rows[1].as_ref().expect("second molecule matches");

    assert_eq!(first.core_smiles, second.core_smiles);
    assert_eq!(first.r_groups.keys().copied().collect::<Vec<_>>(), [1, 2]);
    assert!(first.r_groups[&1].contains("CC"));
    assert!(first.r_groups[&2].contains('C'));
    assert!(second.r_groups[&1].contains("CC"));
    assert!(second.r_groups[&2].contains('N'));
    assert!(first.core_smiles.contains("*:1"));
    assert!(first.core_smiles.contains("*:2"));
    assert!(first.r_groups[&1].contains("*:1"));
    assert!(first.r_groups[&2].contains("*:2"));
}

#[test]
fn labelled_core_fills_hydrogen_only_for_established_columns() {
    for (core, unsubstituted, substituted) in [
        ("[c:1]1ccccc1", "c1ccccc1", "Cc1ccccc1"),
        ("[C:1]1CCCCC1", "C1CCCCC1", "CC1CCCCC1"),
        ("[O:1]C", "CO", "COC"),
        ("[N:1]C", "CN", "CNC"),
    ] {
        let mols = [mol(unsubstituted), mol(substituted)];
        let refs: Vec<&Molecule> = mols.iter().collect();
        let rows = rgroup_decompose_labeled(core, &refs).unwrap();
        let first = rows[0].as_ref().unwrap();
        let second = rows[1].as_ref().unwrap();
        assert_eq!(
            first.r_groups.get(&1).map(String::as_str),
            Some("[*:1][H]"),
            "{core}"
        );
        assert_eq!(first.core_smiles, second.core_smiles, "{core}");

        let one = rgroup_decompose_labeled(core, &refs[..1]).unwrap();
        assert!(one[0].as_ref().unwrap().r_groups.is_empty(), "{core}");
    }
}

#[test]
fn automatically_labelled_site_fills_hydrogen_in_established_column() {
    let mols = [mol("c1ccncc1"), mol("Cc1ccncc1"), mol("Clc1ccncc1")];
    let refs: Vec<&Molecule> = mols.iter().collect();
    let core = "c1ccncc1[*:1]";
    let rows = rgroup_decompose_labeled(core, &refs).unwrap();
    let first = rows[0].as_ref().unwrap();
    assert_eq!(first.r_groups.get(&2).map(String::as_str), Some("[*:2][H]"));
    assert_eq!(first.core_smiles, rows[1].as_ref().unwrap().core_smiles);
    assert_eq!(first.core_smiles, rows[2].as_ref().unwrap().core_smiles);
    for (row, expected_group) in rows.iter().zip(["[H][*:2]", "C[*:2]", "Cl[*:2]"]) {
        let row = row.as_ref().unwrap();
        assert_eq!(
            canonical_smiles(&mol(&row.core_smiles)),
            canonical_smiles(&mol("c1cc([*:2])ccn1"))
        );
        assert_eq!(
            canonical_smiles(&mol(&row.r_groups[&2])),
            canonical_smiles(&mol(expected_group))
        );
    }
    let alone = rgroup_decompose_labeled(core, &refs[..1]).unwrap();
    assert!(alone[0].as_ref().unwrap().r_groups.is_empty());
}

#[test]
fn terminal_mapped_wildcard_accepts_hydrogen_only_for_established_columns() {
    for (core, unsubstituted, substituted, label) in [
        ("c1ccccc1[*:1]", "c1ccccc1", "Cc1ccccc1", 1),
        ("CC([*:1])O", "CCO", "CC(N)O", 1),
        ("[C:1]([*:2])O", "CO", "CCO", 2),
    ] {
        let mols = [mol(unsubstituted), mol(substituted)];
        let refs: Vec<&Molecule> = mols.iter().collect();
        let rows = rgroup_decompose_labeled(core, &refs).unwrap();
        let first = rows[0].as_ref().expect("hydrogen site matches");
        let second = rows[1].as_ref().expect("substituent matches");
        assert_eq!(
            first.r_groups.get(&label).map(String::as_str),
            Some(format!("[*:{label}][H]").as_str()),
            "{core}"
        );
        assert_eq!(first.core_smiles, second.core_smiles, "{core}");

        let alone = rgroup_decompose_labeled(core, &refs[..1]).unwrap();
        assert!(alone[0].as_ref().unwrap().r_groups.is_empty(), "{core}");
    }
}

#[test]
fn two_terminal_wildcards_allow_one_hydrogen_site() {
    let mols = [mol("Cc1ccccc1"), mol("Cc1ccc(CC)cc1")];
    let refs: Vec<&Molecule> = mols.iter().collect();
    let rows = rgroup_decompose_labeled("c1cc([*:1])ccc1[*:2]", &refs).unwrap();
    let first = rows[0].as_ref().expect("one-substituent row matches");
    let second = rows[1].as_ref().expect("two-substituent row matches");
    assert_eq!(first.r_groups[&1], "[*:1]C");
    assert_eq!(first.r_groups[&2], "[*:2][H]");
    assert_eq!(first.core_smiles, second.core_smiles);
}

#[test]
fn two_terminal_wildcards_on_one_atom_align_substituent_and_hydrogen() {
    let mols = [mol("CO"), mol("CCO"), mol("CC(O)N")];
    let refs: Vec<&Molecule> = mols.iter().collect();
    let rows = rgroup_decompose_labeled("C([*:1])([*:2])O", &refs).unwrap();
    let rows: Vec<_> = rows.iter().map(|row| row.as_ref().unwrap()).collect();
    assert!(
        rows.iter()
            .all(|row| row.core_smiles == rows[0].core_smiles)
    );
    assert_eq!(
        canonical_smiles(&mol(&rows[0].core_smiles)),
        canonical_smiles(&mol("OC([*:1])[*:2]"))
    );
    assert_eq!(rows[0].r_groups[&1], "[*:1][H]");
    assert_eq!(rows[0].r_groups[&2], "[*:2][H]");
    assert_eq!(rows[1].r_groups[&1], "[*:1]C");
    assert_eq!(rows[1].r_groups[&2], "[*:2][H]");
    assert_eq!(rows[2].r_groups[&1], "[*:1]C");
    assert_eq!(rows[2].r_groups[&2], "N[*:2]");

    let alone = rgroup_decompose_labeled("C([*:1])([*:2])O", &refs[..1]).unwrap();
    assert!(alone[0].as_ref().unwrap().r_groups.is_empty());

    let pair = rgroup_decompose_labeled("C([*:1])([*:2])O", &refs[..2]).unwrap();
    for row in pair.iter().flatten() {
        assert_eq!(row.r_groups.keys().copied().collect::<Vec<_>>(), [1]);
        assert_eq!(
            canonical_smiles(&mol(&row.core_smiles)),
            canonical_smiles(&mol("OC[*:1]"))
        );
    }
    assert_eq!(pair[0].as_ref().unwrap().r_groups[&1], "[*:1][H]");
    assert_eq!(pair[1].as_ref().unwrap().r_groups[&1], "[*:1]C");

    for (smiles, heteroatom) in [("CC(O)N", "N"), ("CC(O)Cl", "Cl"), ("CC(O)F", "F")] {
        let single = [mol(smiles)];
        let single_refs: Vec<&Molecule> = single.iter().collect();
        for inputs in [&single_refs[..], &[&mols[0], &single[0]][..]] {
            let rows = rgroup_decompose_labeled("C([*:1])([*:2])O", inputs).unwrap();
            let row = rows.last().unwrap().as_ref().unwrap();
            assert_eq!(
                canonical_smiles(&mol(&row.r_groups[&1])),
                canonical_smiles(&mol(&format!("{heteroatom}[*:1]"))),
                "{smiles} with {} rows",
                inputs.len()
            );
            assert_eq!(
                canonical_smiles(&mol(&row.r_groups[&2])),
                canonical_smiles(&mol("C[*:2]")),
                "{smiles} with {} rows",
                inputs.len()
            );
        }
    }
}

#[test]
fn three_terminal_wildcards_align_two_groups_and_one_hydrogen() {
    let core = "C([*:1])([*:2])([*:3])O";
    let mols = [mol("CO"), mol("CCO"), mol("CC(O)N")];
    let refs: Vec<&Molecule> = mols.iter().collect();

    let single = rgroup_decompose_labeled(core, &refs[2..]).unwrap();
    let row = single[0].as_ref().unwrap();
    assert_eq!(row.r_groups[&1], "N[*:1]");
    assert_eq!(row.r_groups[&2], "[*:2]C");
    assert!(!row.r_groups.contains_key(&3));

    let pair = rgroup_decompose_labeled(core, &[refs[0], refs[2]]).unwrap();
    let first = pair[0].as_ref().unwrap();
    let second = pair[1].as_ref().unwrap();
    assert_eq!(first.r_groups[&1], "[*:1][H]");
    assert_eq!(first.r_groups[&2], "[*:2][H]");
    assert_eq!(second.r_groups[&1], "N[*:1]");
    assert_eq!(second.r_groups[&2], "[*:2]C");
    assert!(!second.r_groups.contains_key(&3));

    let carbon_first = rgroup_decompose_labeled(core, &refs[1..]).unwrap();
    let last = carbon_first[1].as_ref().unwrap();
    assert_eq!(last.r_groups[&1], "[*:1]C");
    assert_eq!(last.r_groups[&2], "N[*:2]");

    let established = rgroup_decompose_labeled(core, &refs).unwrap();
    let last = established[2].as_ref().unwrap();
    assert_eq!(last.r_groups[&1], "[*:1]C");
    assert_eq!(last.r_groups[&2], "N[*:2]");
    assert!(!last.r_groups.contains_key(&3));
}

#[test]
fn optional_wildcard_fallback_does_not_invent_multiple_bonded_hydrogen() {
    let mols = [mol("C")];
    let refs: Vec<&Molecule> = mols.iter().collect();
    let rows = rgroup_decompose_labeled("C=[*:1]", &refs).unwrap();
    assert!(rows[0].is_none());
}

#[test]
fn prepared_labelled_core_matches_free_function() {
    let mols = [
        mol("CCc1cc(O)ccc1"),
        mol("Oc1ccc(CC)cc1"),
        mol("Clc1ccc(CC)cc1"),
    ];
    let refs: Vec<&Molecule> = mols.iter().collect();
    let core = "c1cc([*:1])ccc1[*:2]";

    let expected = rgroup_decompose_labeled(core, &refs).unwrap();
    let prepared = PreparedRGroupCore::new(core).unwrap();
    let actual = prepared.decompose_labeled(&refs).unwrap();

    assert_eq!(actual, expected);
}

#[test]
fn labelled_symmetric_columns_use_fragment_similarity_like_rdkit() {
    // The third row has no exact payload seen in either column.  RDKit's
    // table alignment keeps the carbon chain in R1 and the heteroatom in
    // R2; lexical tie-breaking alone swapped those labels.
    let mols = [
        mol("Cc1cccc(N)c1"),
        mol("COc1cccc(Cl)c1"),
        mol("CCc1cccc(O)c1"),
    ];
    let refs: Vec<&Molecule> = mols.iter().collect();
    let rows =
        rgroup_decompose_labeled("c1c([*:1])cccc1[*:2]", &refs).expect("labelled decomposition");
    let third = rows[2].as_ref().expect("third molecule matches");

    assert!(third.r_groups[&1].contains("CC"), "R1 must keep ethyl");
    assert!(third.r_groups[&2].contains('O'), "R2 must keep hydroxyl");

    let respelled = [
        mol("Nc1cccc(C)c1"),
        mol("Clc1cccc(OC)c1"),
        mol("Oc1cccc(CC)c1"),
    ];
    let respelled_refs: Vec<&Molecule> = respelled.iter().collect();
    let respelled_rows = rgroup_decompose_labeled("c1c([*:1])cccc1[*:2]", &respelled_refs)
        .expect("respelled labelled decomposition");
    for (original, respelled) in rows.iter().zip(&respelled_rows) {
        assert_eq!(
            original, respelled,
            "input atom order changed the table row"
        );
    }
}

#[test]
fn mapped_core_atom_labels_its_substituent() {
    let mols = [mol("Cc1ccccc1")];
    let refs: Vec<&Molecule> = mols.iter().collect();
    let rows = rgroup_decompose_labeled("[c:7]1ccccc1", &refs).expect("mapped-core decomposition");
    let row = rows[0].as_ref().expect("toluene matches");
    assert_eq!(row.r_groups.keys().copied().collect::<Vec<_>>(), [7]);
    assert!(row.core_smiles.contains("*:7"));
    assert!(row.r_groups[&7].contains("*:7"));
}

#[test]
fn duplicate_and_nonterminal_labels_are_typed_errors() {
    let mols = [mol("Cc1ccccc1")];
    let refs: Vec<&Molecule> = mols.iter().collect();

    let duplicate = rgroup_decompose_labeled("[c:1]1cccc[c:1]1", &refs);
    assert!(matches!(duplicate, Err(RGroupError::DuplicateLabel(1))));

    let nonterminal = rgroup_decompose_labeled("c1cc([*:2]C)ccc1", &refs);
    assert!(matches!(
        nonterminal,
        Err(RGroupError::NonTerminalAttachmentLabel(2))
    ));
}

#[test]
fn invalid_and_ambiguous_labelled_cores_fail_closed() {
    let mols = [mol("CC(C)(C)C")];
    let refs: Vec<&Molecule> = mols.iter().collect();

    assert!(matches!(
        rgroup_decompose_labeled("[C:0]", &refs),
        Err(RGroupError::InvalidLabel(0))
    ));
    assert!(matches!(
        rgroup_decompose_labeled("[C;D4:1]", &refs),
        Err(RGroupError::MultipleAttachments {
            mol_idx: 0,
            label: 1
        })
    ));

    let radical = [mol("[C]")];
    let radical_refs: Vec<&Molecule> = radical.iter().collect();
    let missing_hydrogen = rgroup_decompose_labeled("C([*:1])", &radical_refs).unwrap();
    assert!(missing_hydrogen[0].is_none());
}

#[test]
fn labelled_core_refuses_when_match_cap_is_reached() {
    let input = vec!["C"; LABELED_MATCH_LIMIT].join(".");
    let mols = [mol(&input)];
    let refs: Vec<&Molecule> = mols.iter().collect();
    assert!(matches!(
        rgroup_decompose_labeled("C", &refs),
        Err(RGroupError::MatchLimitExceeded {
            mol_idx: 0,
            limit: LABELED_MATCH_LIMIT
        })
    ));
}
