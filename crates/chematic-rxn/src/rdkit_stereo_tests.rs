use super::*;

#[test]
fn tokenizer_orders_ring_closures_last_by_number() {
    let atoms = tokenize_template("C[C@@H]1CC[C@H](O)CC1").unwrap();
    assert_eq!(atoms[1].order, vec![0, 2, 7]);
    assert_eq!(atoms[1].tag, Some(false));
    assert_eq!(atoms[4].order, vec![3, 5, 6]);
    assert_eq!(atoms[4].tag, Some(false));
}

#[test]
fn tokenizer_covers_grouped_components_large_rings_and_invalid_boundaries() {
    assert_eq!(
        split_templates("([C:1].[O:2]).[N:3]"),
        vec!["[C:1].[O:2]", "[N:3]"]
    );
    assert_eq!(
        split_templates("([C:1])tail.[O:2]"),
        vec!["([C:1])tail", "[O:2]"]
    );
    assert!(tokenize_template("[C:1]%12CCCCC%12").is_some());
    assert!(tokenize_template("[C:1]%(123)CCCCC%(123)").is_some());
    assert!(tokenize_template("[C:1]([O:2])[N:3]").is_some());
    assert!(tokenize_template("[C:1].[O:2]").is_some());
    // The reaction parser rejects an unclosed branch before this helper;
    // an unmatched closing branch is nevertheless defensive input here.
    assert!(tokenize_template("[C:1])").is_none());
    assert!(tokenize_template("[C:1]%(").is_none());
    assert!(tokenize_template("[C:1]%x1").is_none());
    assert!(tokenize_template("[C:1]1CC").is_none());
}

#[test]
fn parser_cleanup_covers_charge_bond_and_ez_edits() {
    for text in [
        "N(=O)=O",
        "[O-][N+](=O)O",
        "[NH3+]->[Cu+2]",
        "C/N=C(/N)N",
        "F/C=C/F",
        "C/C=C(/F)\\Cl",
    ] {
        let mol = chematic_smiles::parse(text).unwrap();
        let first = rdkit_parse_cleanup(&mol);
        let second = rdkit_parse_cleanup(&mol);
        assert_eq!(
            first.as_ref().as_ref().map(Molecule::atom_count),
            second.as_ref().as_ref().map(Molecule::atom_count)
        );
    }
}

#[test]
fn tokenizer_h_only_in_simple_spelling() {
    assert_eq!(
        tokenize_template("[C@@H](N)(C)O").unwrap()[0].tag,
        Some(false)
    );
    assert_eq!(
        tokenize_template("[C@@;H1:2](N)(C)O").unwrap()[0].tag,
        Some(true)
    );
    assert_eq!(
        tokenize_template("N[C@@H:2](C)O").unwrap()[1].tag,
        Some(true)
    );
    assert_eq!(
        tokenize_template("[C@@H:2](N)(C)O").unwrap()[0].map,
        Some(2)
    );
}

#[test]
fn count_swaps_matches_rdkit() {
    assert_eq!(count_swaps(&[0, 1, 2, 3], &[1, 0, 2, 3]), Some(1));
    assert_eq!(count_swaps(&[0, 1, 2, 3], &[1, 2, 3, 0]), Some(3));
    assert_eq!(count_swaps(&[0, 1, 2, 3], &[1, 2, 0, 3]), Some(2));
    assert_eq!(count_swaps(&[0, 1], &[0]), None);
    assert_eq!(count_swaps(&[0, 1], &[0, 2]), None);
}

#[test]
fn tokenizer_and_chiral_h_reject_unsupported_syntax_without_panicking() {
    for accepted in [
        "ClBrBNOPSFIbcnopsaA*",
        "[13C@@H+:7](F)(Cl)Br",
        "[#6@H:1](F)(Cl)Br",
        "[C;$([C](F)):1](N)O",
    ] {
        assert!(tokenize_template(accepted).is_some(), "{accepted}");
    }
    for rejected in [
        "[",
        "C(C.C)",
        "C?C",
        "[C@@@:1](F)(Cl)Br",
        "[C@TH1:1](F)(Cl)Br",
        "[C@?:1](F)(Cl)Br",
    ] {
        assert!(tokenize_template(rejected).is_none(), "{rejected}");
    }

    for accepted in ["C@H", "13C@@H+", "#6@H1-", "Cl@H", "*@H"] {
        assert!(simple_chiral_h(accepted), "{accepted}");
    }
    for rejected in ["", "1@H", "#@H", "C", "C@", "C@@", "C@TH1", "C@Hfoo"] {
        assert!(!simple_chiral_h(rejected), "{rejected}");
    }
}

#[test]
fn permutation_and_template_overlap_helpers_cover_bounded_failures() {
    assert_eq!(permutation_parity(&[0, 1], &[1, 0]), Some(false));
    assert_eq!(permutation_parity(&[0], &[0, 1]), None);
    assert_eq!(permutation_parity(&[0, 1], &[0, 2]), None);
    assert_eq!(permutation_parity(&[0, 1], &[0, 0]), None);

    let perms = permutations(&[1, 2, 3], MAX_ORDER_OPTIONS).unwrap();
    assert_eq!(perms.len(), 6);
    assert!(permutations(&[0, 1, 2, 3, 4, 5, 6], MAX_ORDER_OPTIONS).is_none());

    let mut with_slot = [1, -1, 3];
    assert!(order_overlap(&mut with_slot, 1, &[1, 2, 3]));
    assert_eq!(with_slot, [1, 2, 3]);
    let mut no_slot = [1, 3];
    assert!(!order_overlap(&mut no_slot, 0, &[1, 2, 3]));
    assert!(!order_overlap(&mut no_slot, 1, &[1, 2, 3]));
}

#[test]
fn inversion_flags_cover_retain_invert_remove_create_and_unmapped_atoms() {
    let info = ReactionStereo::new(
        "[C@@H:1](F)(Cl)Br.[C@@H:2](F)(Cl)Br.[C:3](F)(Cl)Br",
        "[C@@H:1](F)(Cl)Br.[C@H:2](F)(Cl)Br.[C:4](F)(Cl)Br.[C@@H](F)(Cl)Br",
    )
    .expect("supported reaction stereo templates");
    assert_eq!(info.flags[0][0], Some(2));
    assert_eq!(info.flags[1][0], Some(1));
    assert_eq!(info.flags[2][0], Some(4));
    assert_eq!(info.flags[3][0], Some(4));

    let removed = ReactionStereo::new("[C@@H:1](F)(Cl)Br", "[C:1](F)(Cl)Br").unwrap();
    assert_eq!(removed.flags[0][0], Some(3));
}
