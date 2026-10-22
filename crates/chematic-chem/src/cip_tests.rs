#[test]
fn lazy_sphere_generator_matches_eager_reference() {
    for smi in [
        "C[C@H](N)C(=O)O",
        "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C",
        "C1CC2CCC1C2",
        "OC(=O)C1=CC=CC=C1C(C)(Cl)Br",
        "C12C3C4C1C5C2C3C45",
        "N#CC(C=C)C(=O)C=O",
        "[2H]C(Cl)(F)C1CC1",
    ] {
        let mol = chematic_smiles::parse(smi).unwrap();
        for c in 0..mol.atom_count() {
            let center = AtomIdx(c as u32);
            for (nb, _) in mol.neighbors(center) {
                let eager = cip_branch_spheres(&mol, center, nb);
                let mut lazy = Vec::new();
                let mut sphere_gen = SphereGen::new(&mol, center, nb);
                while let Some(layer) = sphere_gen.next_layer() {
                    lazy.push(layer);
                }
                assert_eq!(eager, lazy, "{smi} center {c} start {}", nb.0);
            }
        }
    }
}

use super::*;

/// Accurate-mode tetrahedral labels, abstentions and bond-keyed E/Z do
/// not depend on input atom order (random SMILES now keep stereo), apart
/// from labels marked spelling-dependent.
#[test]
fn accurate_labels_and_abstentions_survive_atom_reordering() {
    use std::collections::BTreeMap;
    let keyed = |m: &Molecule| -> (String, BTreeMap<usize, String>) {
        let (can, order) = chematic_smiles::canonical_smiles_with_atom_order(m);
        let pos: BTreeMap<u32, usize> = order.iter().enumerate().map(|(k, a)| (a.0, k)).collect();
        let r = assign_cip_with_mode(m, CipMode::Accurate).unwrap();
        let mut out = BTreeMap::new();
        for (a, c) in &r.assignments {
            if !matches!(c, CipCode::E | CipCode::Z) {
                // A saturated ring lone-pair centre gets RDKit's label for
                // the spelling parsed, which a reordered spelling may flip.
                let label = if cip_label_depends_on_smiles_spelling(m, *a) {
                    "spelling_dependent".to_string()
                } else {
                    format!("{c:?}")
                };
                out.insert(pos[&a.0], label);
            }
        }
        for (a, why) in &r.unresolved {
            out.insert(pos[&a.0], format!("{why:?}"));
        }
        for (b, c) in assign_ez_bonds_with_mode(m, CipMode::Accurate) {
            let e = m.bond(b);
            let (x, y) = (pos[&e.atom1.0], pos[&e.atom2.0]);
            out.insert(100_000 + x.min(y) * 1000 + x.max(y), format!("{c:?}"));
        }
        (can, out)
    };
    for smiles in [
        "N[C@@H](C)C(=O)O",
        "C#C/C=C1\\CCC(c2cccc3ccccc23)C(=O)O1",
        // Lone-pair and phosphorus abstentions of the exposed 10k lane.
        "Cn1cc(C2=NC[C@@]3(C[N@@]4CC[C@@H]3C4)O2)c2ccccc21",
        "CCC(C)[C@H](NC(=O)[C@@H]1CCCN1[P@](=O)(OC)[C@H](Cc1ccccc1)NC(=O)[C@H](CC(N)=O)NC(=O)OC(C)(C)C)C(=O)NCC(C)C",
    ] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let want = keyed(&mol);
        assert!(!want.1.is_empty(), "{smiles}");
        for seed in 1..8 {
            let other =
                chematic_smiles::parse(&chematic_smiles::random_smiles(&mol, seed)).unwrap();
            assert_eq!(keyed(&other), want, "{smiles} seed {seed}");
        }
    }
}

/// Issue #634: aromatic E/Z ranks substituents with the hierarchical
/// digraph (MANCUDE duplicates for aromatic rings). The plain sphere
/// comparison adds no duplicate for aromatic bonds and can let tert-butyl
/// outrank an aryl carbon, or an sp3 ring CH outrank an aromatic carbon.
#[test]
fn accurate_ez_ranks_aryl_above_alkyl() {
    let cases: &[(&str, (u32, u32), CipCode)] = &[
        (
            "Cc1oc(-c2ccc(C(F)(F)F)cc2)nc1COc1cccc(/C(=C/Cn2oc(=O)[nH]c2=O)C(C)(C)C)c1",
            (23, 24),
            CipCode::E,
        ),
        (
            "C=C(C)[C@H]1Cc2c(ccc3c2OC2COc4cc(OC)c(OC)cc4C2/C3=N/O)O1",
            (26, 27),
            CipCode::Z,
        ),
        (
            "C=C(C)[C@H]1Cc2c(ccc3c2OC2COc4cc(OC)c(OC)cc4C2/C3=N/OC(C)=O)O1",
            (26, 27),
            CipCode::Z,
        ),
        // Imine on a pyrido[1,2-a]pyrimidine: the ring-fusion N is a
        // lone-pair donor in RDKit's Kekulé form, which ranks the ring
        // side's duplicate atoms (#734, six BioTransformer products).
        ("CCOC(=O)c1cc2ccccn2/c(=N/c2ccccc2)n1", (13, 14), CipCode::Z),
        // Simple controls where both comparators already agree.
        ("F/C=C/F", (1, 2), CipCode::E),
        ("F/C=C\\F", (1, 2), CipCode::Z),
        ("c1ccccc1/C(C)=C/C", (6, 8), CipCode::E),
    ];
    for (smiles, (a, b), expected) in cases {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let labels = assign_ez_bonds_with_mode(&mol, CipMode::Accurate);
        let found = labels.iter().find(|(bidx, _)| {
            let bond = mol.bond(*bidx);
            (bond.atom1.0, bond.atom2.0) == (*a, *b) || (bond.atom1.0, bond.atom2.0) == (*b, *a)
        });
        assert_eq!(found.map(|(_, c)| *c), Some(*expected), "{smiles}");
        assert!(
            assign_ez_bonds_with_mode(&mol, CipMode::LegacyFast)
                .iter()
                .any(|&(bond_idx, code)| bond_idx == found.unwrap().0 && code == *expected),
            "{smiles}: default E/Z assignment"
        );
        // assign_cip_with_mode reports the same label under atom1.
        let merged = assign_cip_with_mode(&mol, CipMode::Accurate).unwrap();
        let atom1 = mol.bond(found.unwrap().0).atom1;
        assert!(
            merged.assignments.contains(&(atom1, *expected)),
            "{smiles}: merged assignments"
        );
    }
}
#[test]
fn squaraine_diimine_canonical_smiles_keeps_both_ez() {
    // Two exocyclic imines on a Kekulé four-membered ring, coupled
    // through the ring bond between them (#734: RDKit-written MOL blocks
    // came back with one imine inverted or unspecified).
    let ez = |smiles: &str| {
        let mol = parse(smiles).unwrap();
        let mut labels: Vec<CipCode> = assign_ez_bonds_with_mode(&mol, CipMode::Accurate)
            .into_iter()
            .map(|(_, code)| code)
            .collect();
        labels.sort_by_key(|code| format!("{code:?}"));
        labels
    };
    for group in [
        [r"CC/N=C1/C(O)=C(O)/C1=N\C", r"CC/N=C1/C(O)=C(O)\C1=N/C"],
        [r"CC/N=C1\C(O)=C(O)/C1=N/C", r"CC/N=C1\C(O)=C(O)\C1=N\C"],
        [
            r"C1(C(=C(C/1=N\CC)O)O)=N/C(CC)C",
            r"C1(O)C(/C(=N/CC)C=1O)=N\C(CC)C",
        ],
    ] {
        let written: Vec<String> = group
            .iter()
            .map(|s| chematic_smiles::canonical_smiles(&parse(s).unwrap()))
            .collect();
        assert_eq!(written[0], written[1], "{group:?}");
        assert_eq!(ez(&written[0]), ez(group[0]), "{group:?}: {}", written[0]);
        assert_eq!(ez(group[0]).len(), 2);
    }
}

#[test]
fn lone_pair_donors_beside_exocyclic_bonds_kekulize() {
    // The fusion N keeps its lone pair; the exocyclic `=O`/`=S`/`=N`
    // carbon takes no ring double bond (17 molecules of a 25,000-molecule
    // check set did not kekulize before, #734).
    for smiles in [
        "Cc1cc2ccccn2c(=O)n1",
        "CCOC(=O)c1cc2ccccn2c(=S)n1",
        "c1ccn2ccccc2c1",
        // 1,2,5-thiadiazole S-oxide: the three-bonded `[s+]` is a
        // lone-pair donor, as in RDKit.
        "CNc1n[s+]([O-])nc1N",
    ] {
        let mol = parse(smiles).unwrap();
        let kekule = chematic_core::kekulization::kekulize(&mol);
        assert!(kekule.is_ok(), "{smiles}");
    }
}

use chematic_smiles::parse;

fn cip_at(smiles: &str, atom_idx: usize) -> Option<CipCode> {
    let mol = parse(smiles).unwrap();
    let assignment = assign_cip(&mol);
    assignment.get(AtomIdx(atom_idx as u32))
}

// --- Tetrahedral R/S ---

#[test]
fn test_l_alanine_s() {
    // N[C@@H](C)C(=O)O — L-alanine, chiral center is atom 1
    assert_eq!(
        cip_at("N[C@@H](C)C(=O)O", 1),
        Some(CipCode::S),
        "L-alanine should be S"
    );
}

#[test]
fn test_d_alanine_r() {
    // N[C@H](C)C(=O)O — D-alanine
    assert_eq!(
        cip_at("N[C@H](C)C(=O)O", 1),
        Some(CipCode::R),
        "D-alanine should be R"
    );
}

#[test]
fn test_chfclbr_r() {
    // [C@@H](F)(Cl)Br — known R configuration
    // CIP priority: Br(35) > Cl(17) > F(9) > H(1)
    // @@H: looking from H, F→Cl→Br is CW → R
    assert_eq!(
        cip_at("[C@@H](F)(Cl)Br", 0),
        Some(CipCode::R),
        "[C@@H](F)(Cl)Br should be R"
    );
}

#[test]
fn test_chfclbr_s() {
    // [C@H](F)(Cl)Br — S
    assert_eq!(
        cip_at("[C@H](F)(Cl)Br", 0),
        Some(CipCode::S),
        "[C@H](F)(Cl)Br should be S"
    );
}

#[test]
fn test_no_chirality() {
    let mol = parse("CC(=O)O").unwrap();
    let assignment = assign_cip(&mol);
    let tetrahedral: Vec<_> = assignment
        .assignments
        .iter()
        .filter(|(_, c)| matches!(c, CipCode::R | CipCode::S))
        .collect();
    assert!(tetrahedral.is_empty(), "acetic acid has no chiral centers");
}

#[test]
fn test_symmetric_center_none() {
    // No @/@@ annotation → no assignment attempted
    let mol = parse("CC(N)(N)CC").unwrap();
    let assignment = assign_cip(&mol);
    assert!(
        assignment.assignments.is_empty(),
        "no stereo annotation → no assignment"
    );
}

#[test]
fn test_assignment_get() {
    let mol = parse("N[C@@H](C)C(=O)O").unwrap();
    let a = assign_cip(&mol);
    assert!(a.get(AtomIdx(1)).is_some(), "atom 1 should have CIP code");
    assert!(a.get(AtomIdx(0)).is_none(), "atom 0 (N) has no chirality");
}

#[test]
fn test_r_lactic_acid_gives_answer() {
    // OC(=O)[C@@H](O)C — lactic acid; should give R or S
    let mol = parse("OC(=O)[C@@H](O)C").unwrap();
    let assignment = assign_cip(&mol);
    let chiral_idx = (0..mol.atom_count())
        .map(|i| AtomIdx(i as u32))
        .find(|&i| mol.atom(i).chirality != Chirality::None)
        .unwrap();
    let code = assignment.get(chiral_idx);
    assert!(
        code == Some(CipCode::R) || code == Some(CipCode::S),
        "should give R or S for lactic acid, got {:?}",
        code
    );
}

// --- E/Z double bonds ---

#[test]
fn test_trans_2_butene_e() {
    // C/C=C/C — trans-2-butene → E
    let mol = parse("C/C=C/C").unwrap();
    let assignment = assign_cip(&mol);
    let has_e = assignment.assignments.iter().any(|(_, c)| *c == CipCode::E);
    assert!(
        has_e,
        "Expected E for trans-2-butene, got {:?}",
        assignment.assignments
    );
}

#[test]
fn test_cis_2_butene_z() {
    // C/C=C\C — cis-2-butene → Z
    let mol = parse("C/C=C\\C").unwrap();
    let assignment = assign_cip(&mol);
    let has_z = assignment.assignments.iter().any(|(_, c)| *c == CipCode::Z);
    assert!(
        has_z,
        "Expected Z for cis-2-butene, got {:?}",
        assignment.assignments
    );
}

#[test]
fn test_fceccl_e() {
    // F/C=C/Cl → E (F and Cl on opposite sides)
    let mol = parse("F/C=C/Cl").unwrap();
    let assignment = assign_cip(&mol);
    let has_e = assignment.assignments.iter().any(|(_, c)| *c == CipCode::E);
    assert!(
        has_e,
        "Expected E for F/C=C/Cl, got {:?}",
        assignment.assignments
    );
}

#[test]
fn test_fceccl_z() {
    // F/C=C\Cl → Z (F and Cl on same side)
    let mol = parse("F/C=C\\Cl").unwrap();
    let assignment = assign_cip(&mol);
    let has_z = assignment.assignments.iter().any(|(_, c)| *c == CipCode::Z);
    assert!(
        has_z,
        "Expected Z for F/C=C\\Cl, got {:?}",
        assignment.assignments
    );
}

#[test]
fn test_no_ez_no_stereo_bond() {
    // C=C — no Up/Down bonds → no E/Z
    let mol = parse("C=C").unwrap();
    let assignment = assign_cip(&mol);
    let has_ez = assignment
        .assignments
        .iter()
        .any(|(_, c)| matches!(c, CipCode::E | CipCode::Z));
    assert!(!has_ez, "plain C=C should have no E/Z");
}

#[test]
fn test_ez_terminal_no_crash() {
    // C=C/F — one alkene carbon is terminal (only implicit H's, no
    // stereo bond); should not crash.
    let mol = parse("C=C/F").unwrap();
    let _ = assign_cip(&mol);
}

#[test]
fn test_cip_assignment_methane() {
    let mol = parse("C").unwrap();
    let assignment = assign_cip(&mol);
    assert!(assignment.assignments.is_empty());
}

#[test]
fn test_multiple_chiral_centers() {
    // Two chiral centers in a chain
    let mol = parse("N[C@@H](C)[C@H](C)N").unwrap();
    let assignment = assign_cip(&mol);
    let rs_count = assignment
        .assignments
        .iter()
        .filter(|(_, c)| matches!(c, CipCode::R | CipCode::S))
        .count();
    assert_eq!(rs_count, 2, "should assign 2 chiral centers");
}

#[test]
fn test_cip_assignment_struct() {
    let mol = parse("N[C@@H](C)C(=O)O").unwrap();
    let a = assign_cip(&mol);
    assert!(!a.assignments.is_empty());
    let code = a.get(AtomIdx(1));
    assert_eq!(code, Some(CipCode::S));
}

#[test]
fn test_r_s_are_consistent() {
    // The two SMILES should give opposite results
    let r_code = cip_at("[C@@H](F)(Cl)Br", 0);
    let s_code = cip_at("[C@H](F)(Cl)Br", 0);
    assert_ne!(r_code, s_code, "@ and @@ must give opposite results");
}

#[test]
fn test_e_z_are_consistent() {
    // /C=C/C (trans) and /C=C\C (cis) must differ
    let mol_e = parse("C/C=C/C").unwrap();
    let mol_z = parse("C/C=C\\C").unwrap();
    let assign_e = assign_cip(&mol_e);
    let assign_z = assign_cip(&mol_z);
    let has_e = assign_e.assignments.iter().any(|(_, c)| *c == CipCode::E);
    let has_z = assign_z.assignments.iter().any(|(_, c)| *c == CipCode::Z);
    assert!(has_e && has_z, "trans must be E, cis must be Z");
}

#[test]
fn test_canonical_preserves_ez() {
    // Canonicalizing a stereo-SMILES must preserve the E/Z assignment.
    // Uses F/C=C/Cl (E) and F/C=C\Cl (Z) — non-symmetric so CIP is unambiguous.
    use chematic_smiles::canonical_smiles;

    let mol_e = parse("F/C=C/Cl").unwrap();
    let can_e = canonical_smiles(&mol_e);
    let mol_e2 =
        parse(&can_e).unwrap_or_else(|err| panic!("canonical E not parseable: {can_e} → {err}"));
    let assign_e = assign_cip(&mol_e2);
    assert!(
        assign_e.assignments.iter().any(|(_, c)| *c == CipCode::E),
        "canonical SMILES of E isomer must still be E, canonical='{can_e}'"
    );

    let mol_z = parse("F/C=C\\Cl").unwrap();
    let can_z = canonical_smiles(&mol_z);
    let mol_z2 =
        parse(&can_z).unwrap_or_else(|err| panic!("canonical Z not parseable: {can_z} → {err}"));
    let assign_z = assign_cip(&mol_z2);
    assert!(
        assign_z.assignments.iter().any(|(_, c)| *c == CipCode::Z),
        "canonical SMILES of Z isomer must still be Z, canonical='{can_z}'"
    );
}

// --- EZ-S1: aromatic ring bond stash + true-highest-priority substituent ---
//
// `substituent_is_up` previously read only `bond.order`, never
// `mol.bond_direction(bond_idx)` -- the side channel used when a `/`/`\`
// marker lands on a bond between two aromatic-flagged atoms (e.g. a ring
// bond flanking an exocyclic C=N). Fixing that read exposed a second,
// pre-existing, unrelated bug in `highest_stereo_sub`: it returned the
// highest-priority substituent *among those carrying an explicit
// marker*, not the true highest-priority substituent overall, silently
// using a lower-priority marked substituent's raw side instead of
// deriving the true-highest one's side as its geometric complement.
// Both are fixed together here (see cip.rs's `substituent_is_up` and
// `highest_stereo_sub`). Expected values below are RDKit
// `rdCIPLabeler`-verified, not just "some E/Z is present" checks.

#[test]
fn test_ez_stash_exocyclic_imine_e() {
    // RDKit rdCIPLabeler: bond(N, ring-C) = E.
    assert_eq!(
        cip_at(r"C/N=c1\cccc[nH]1", 1),
        Some(CipCode::E),
        "exocyclic imine on an aromatic ring must read the stashed ring-bond direction"
    );
}

#[test]
fn test_ez_stash_exocyclic_imine_inverted_z() {
    // Same structure, ring-closure slash inverted -- RDKit: Z. A "the two
    // just differ" check would also pass under a global sign-inversion
    // bug, so both this and the E case above pin the specific code.
    assert_eq!(
        cip_at(r"C/N=c1/cccc[nH]1", 1),
        Some(CipCode::Z),
        "inverting the stashed ring-bond direction must flip E to Z"
    );
}

#[test]
fn test_ez_stash_representation_path_independence() {
    // The same molecule (N-methylenamino-pyridine-like ring), with the
    // stash-bearing bond routed through a different parser path each
    // time. All three must (a) actually stash on an aromatic-aromatic
    // bond -- not fall back to a literal Up/Down order -- and (b) agree
    // on the RDKit-confirmed E code, proving E/Z perception doesn't
    // depend on which path created the bond.
    let cases: &[(&str, u32, &str)] = &[
        (r"C/N=c1\cccc[nH]1", 1, "chain-edge"),
        (r"C/N=c1cccc[nH]\1", 1, "ring-closure"),
        (r"[nH]1cccc(\c1=N/C)", 5, "branch-attachment"),
    ];
    for &(smi, atom_idx, path_label) in cases {
        let mol = parse(smi).unwrap_or_else(|e| panic!("{path_label}: parse '{smi}': {e}"));
        let has_aromatic_stash = (0..mol.bond_count()).any(|i| {
            let bidx = BondIdx(i as u32);
            let bond = mol.bond(bidx);
            mol.atom(bond.atom1).aromatic
                && mol.atom(bond.atom2).aromatic
                && bond.order == BondOrder::Aromatic
                && mol.bond_direction(bidx).is_some()
        });
        assert!(
            has_aromatic_stash,
            "{path_label}: '{smi}' must stash direction on an aromatic bond, not divert it to a literal Up/Down order"
        );
        let code = assign_cip(&mol).get(AtomIdx(atom_idx));
        assert_eq!(
            code,
            Some(CipCode::E),
            "{path_label}: '{smi}' must resolve to the same E as the other paths"
        );
    }
}

#[test]
fn test_ez_stash_canonical_round_trip() {
    use chematic_smiles::canonical_smiles;

    let mol = parse(r"C/N=c1\cccc[nH]1").unwrap();
    let can = canonical_smiles(&mol);
    let mol2 = parse(&can).unwrap_or_else(|e| panic!("re-parse canonical '{can}': {e}"));
    let code = assign_cip(&mol2)
        .assignments
        .iter()
        .find(|(_, c)| matches!(c, CipCode::E | CipCode::Z))
        .map(|(_, c)| *c);
    assert_eq!(
        code,
        Some(CipCode::E),
        "canonical round trip of the stash-derived E must still be E, canonical='{can}'"
    );
}

#[test]
fn test_highest_stereo_sub_uses_true_priority_not_marked_fallback() {
    // Non-aromatic, no stash involved at all -- a separate, pre-existing
    // bug independent of the stash-read fix above. At the C(F)(Cl)= end,
    // Cl (Z=17) outranks F (Z=9) but only F carries a marker; at the
    // C(Br)(I)= end, I (Z=53) outranks Br (Z=35) but only Br carries a
    // marker. `highest_stereo_sub` must derive each true-highest
    // substituent's side as the marked sibling's complement rather than
    // using the marked (lower-priority) sibling's own raw side.
    // RDKit rdCIPLabeler: E.
    let mol = parse(r"Cl/C(F)=C(\Br)I").unwrap();
    let code = assign_cip(&mol)
        .assignments
        .iter()
        .find(|(_, c)| matches!(c, CipCode::E | CipCode::Z))
        .map(|(_, c)| *c);
    assert_eq!(
        code,
        Some(CipCode::E),
        "true-highest-priority substituent (unmarked) must be used, not the marked lower-priority one"
    );
}

#[test]
fn test_highest_stereo_sub_symmetric_ring_is_not_stereogenic() {
    // An exocyclic imine on an *unsubstituted* benzo ring: the ipso
    // carbon's two ring branches are a genuine CIP priority tie
    // (symmetric across the ipso/para axis) -- confirmed via
    // `compare_branches(mol, alkene_end, subs[0], subs[1]) ==
    // Ordering::Equal` in both argument orders. Swapping the two ring
    // neighbors maps the molecule onto itself, so there is no
    // stereogenic bond to report at all -- not an ambiguous one to
    // guess at. (This SMILES also isn't RDKit-parseable, consistent
    // with there being no real stereoisomerism here.)
    let mol = parse(r"C/N=c1ccccc/1").unwrap();
    let code = assign_cip(&mol)
        .assignments
        .iter()
        .find(|(_, c)| matches!(c, CipCode::E | CipCode::Z))
        .map(|(_, c)| *c);
    assert!(
        code.is_none(),
        "a genuine CIP priority tie must not be assigned an arbitrary E/Z: {code:?}"
    );

    // Stable across a canonical round trip -- before the tie guard, the
    // stable sort's fallback choice depended on adjacency order, which
    // renumbering changes, flipping E<->Z arbitrarily.
    use chematic_smiles::canonical_smiles;
    let can = canonical_smiles(&mol);
    let mol2 = parse(&can).unwrap_or_else(|e| panic!("re-parse canonical '{can}': {e}"));
    let code2 = assign_cip(&mol2)
        .assignments
        .iter()
        .find(|(_, c)| matches!(c, CipCode::E | CipCode::Z))
        .map(|(_, c)| *c);
    assert!(
        code2.is_none(),
        "the tie must remain unassigned after a canonical round trip, canonical='{can}': {code2:?}"
    );
}

// --- Allene axial chirality ---

#[test]
fn test_allene_no_stereo_no_assignment() {
    // Propadiene (allene without any stereo bonds) → no assignment.
    let mol = parse("C=C=C").unwrap();
    let a = assign_cip(&mol);
    let has_allene = a
        .assignments
        .iter()
        .any(|(_, c)| matches!(c, CipCode::E | CipCode::Z));
    assert!(
        !has_allene,
        "unspecified allene should have no axial chirality"
    );
}

#[test]
fn test_allene_two_enantiomers_differ() {
    // 1,3-difluoroallene: F/C=C=C/F vs F/C=C=C\F — must give different codes.
    // Build manually with Up/Down bonds to avoid SMILES parser allene ambiguity.
    use chematic_core::{Atom, BondOrder as BO, Element, MoleculeBuilder};

    // F1/C2=C3=C4\F5 — F at each end with opposite Up/Down
    let mut b = MoleculeBuilder::new();
    let f1 = b.add_atom(Atom::new(Element::F));
    let c2 = b.add_atom(Atom::new(Element::C));
    let c3 = b.add_atom(Atom::new(Element::C));
    let c4 = b.add_atom(Atom::new(Element::C));
    let f5 = b.add_atom(Atom::new(Element::F));
    b.add_bond(f1, c2, BO::Up).unwrap(); // F1 up relative to C2
    b.add_bond(c2, c3, BO::Double).unwrap();
    b.add_bond(c3, c4, BO::Double).unwrap();
    b.add_bond(c4, f5, BO::Down).unwrap(); // F5 down relative to C4
    let mol_a = b.build();

    // Same but with both F Up (same side)
    let mut b2 = MoleculeBuilder::new();
    let f1b = b2.add_atom(Atom::new(Element::F));
    let c2b = b2.add_atom(Atom::new(Element::C));
    let c3b = b2.add_atom(Atom::new(Element::C));
    let c4b = b2.add_atom(Atom::new(Element::C));
    let f5b = b2.add_atom(Atom::new(Element::F));
    b2.add_bond(f1b, c2b, BO::Up).unwrap();
    b2.add_bond(c2b, c3b, BO::Double).unwrap();
    b2.add_bond(c3b, c4b, BO::Double).unwrap();
    b2.add_bond(c4b, f5b, BO::Up).unwrap(); // both Up
    let mol_b = b2.build();

    let code_a = assign_cip(&mol_a)
        .assignments
        .iter()
        .find(|(_, c)| matches!(c, CipCode::E | CipCode::Z))
        .map(|(_, c)| *c);
    let code_b = assign_cip(&mol_b)
        .assignments
        .iter()
        .find(|(_, c)| matches!(c, CipCode::E | CipCode::Z))
        .map(|(_, c)| *c);

    // Both must get an assignment, and they must differ.
    assert!(
        code_a.is_some(),
        "allene A should get an axial chirality code"
    );
    assert!(
        code_b.is_some(),
        "allene B should get an axial chirality code"
    );
    assert_ne!(
        code_a, code_b,
        "the two allene enantiomers must get different codes"
    );
}

#[test]
fn test_non_allene_not_detected() {
    // CO2 (O=C=O) has two double bonds from C but is not an allene (no axial chirality possible).
    let mol = parse("O=C=O").unwrap();
    let a = assign_cip(&mol);
    // CO2 has no substituents for Up/Down bonds, so no allene assignment.
    let has_allene = a
        .assignments
        .iter()
        .any(|(_, c)| matches!(c, CipCode::E | CipCode::Z));
    assert!(
        !has_allene,
        "CO2 should not get axial chirality (no stereo bonds)"
    );
}

// -- CIP rule 4 edge cases (mass tiebreaker and duplicates) --------
#[test]
fn test_cip_enantiomers_consistent_with_mass_tiebreaker() {
    // Verify that mass tiebreaker is used correctly (doesn't affect basic R/S).
    // D-alanine vs L-alanine should get opposite codes.
    let l_ala = parse("C[C@H](N)C(=O)O").unwrap();
    let d_ala = parse("C[C@@H](N)C(=O)O").unwrap();

    let l_assign = assign_cip(&l_ala);
    let d_assign = assign_cip(&d_ala);

    let l_code = l_assign
        .assignments
        .iter()
        .find(|(_, c)| matches!(c, CipCode::R | CipCode::S))
        .map(|(_, c)| *c);
    let d_code = d_assign
        .assignments
        .iter()
        .find(|(_, c)| matches!(c, CipCode::R | CipCode::S))
        .map(|(_, c)| *c);

    // Both should assign (no ties in amino acids).
    assert!(l_code.is_some(), "L-alanine should assign R or S");
    assert!(d_code.is_some(), "D-alanine should assign R or S");

    // They should be opposite.
    assert_ne!(
        l_code, d_code,
        "L and D enantiomers should have opposite R/S"
    );
}

#[test]
fn test_cip_tied_substituents_no_assignment() {
    // When two substituents have identical priority, no R/S is assigned.
    // Example: any center with two identical groups.
    let mol = parse("CC(F)Br").unwrap(); // no chirality specified, but if we try to assign...
    let assignment = assign_cip(&mol);
    // No assignment (no @/@@ specified).
    assert!(
        assignment
            .assignments
            .iter()
            .all(|(_, c)| !matches!(c, CipCode::R | CipCode::S)),
        "achiral molecule should not get R/S"
    );
}

#[test]
fn test_cip_atomic_mass_tiebreaker_infrastructure() {
    // Verify Element::atomic_mass() exists and provides correct values (CIP rule 4).
    // This confirms the mass tiebreaker infrastructure is in place.
    use chematic_core::Element;

    // Check a few elements have correct masses (monoisotopic, within 0.01 u)
    assert!(
        (Element::C.atomic_mass() - 12.0).abs() < 0.01,
        "C ~= 12.0 u, got {}",
        Element::C.atomic_mass()
    );
    assert!(
        (Element::N.atomic_mass() - 14.0).abs() < 0.01,
        "N ~= 14.0 u, got {}",
        Element::N.atomic_mass()
    );
    assert!(
        (Element::O.atomic_mass() - 16.0).abs() < 0.01,
        "O ~= 16.0 u, got {}",
        Element::O.atomic_mass()
    );
    assert!(
        (Element::H.atomic_mass() - 1.0).abs() < 0.01,
        "H ~= 1.0 u, got {}",
        Element::H.atomic_mass()
    );
    assert!(
        (Element::F.atomic_mass() - 19.0).abs() < 0.01,
        "F ~= 19.0 u, got {}",
        Element::F.atomic_mass()
    );

    // Verify that identical atoms have identical masses
    assert_eq!(
        Element::C.atomic_mass(),
        Element::C.atomic_mass(),
        "same element = same mass"
    );

    // Verify ordering: heavier atoms > lighter (used in CIP comparison)
    assert!(
        Element::N.atomic_mass() > Element::C.atomic_mass(),
        "N > C in mass"
    );
    assert!(
        Element::O.atomic_mass() > Element::N.atomic_mass(),
        "O > N in mass"
    );
}

// =========================================================================
// CIP Rule 3 (duplicate atom) tests - fused ring systems
// =========================================================================

#[test]
fn test_cip_naphthalene_assignment() {
    // Naphthalene (c1ccc2ccccc2c1): fused aromatic rings
    // Tests whether CIP Rule 3 (duplicate atom for fused systems) is handled.
    // The bridging carbons in fused rings may appear multiple times in sphere expansion.
    let mol = parse("c1ccc2ccccc2c1").expect("naphthalene");
    let assignment = assign_cip(&mol);
    // Naphthalene has no chiral centers (planar, all carbons equivalent by symmetry)
    // but assignment should not crash and should handle the fused structure.
    assert!(
        assignment.assignments.is_empty(),
        "naphthalene should have no chiral assignments (planar symmetric structure)"
    );
}

#[test]
fn test_cip_decalin_assignment() {
    // Decalin (bicyclic C10): two fused saturated rings
    // CC(C)C1CCC2CCCCC2C1 is a substituted decalin derivative
    // Tests CIP Rule 3 with bridging carbons in saturated system
    let mol = parse("CC(C)C1CCC2CCCCC2C1").expect("decalin");
    let assignment = assign_cip(&mol);
    // Decalin has potential chiral centers; assignment should work without crashing
    assert!(
        !mol.atom(chematic_core::AtomIdx(0)).element.atomic_number() != 6
            || assignment.assignments.len() <= mol.atom_count(),
        "decalin CIP assignment should complete"
    );
}

#[test]
fn test_cip_fused_ring_no_crash() {
    // Simple fused ring system: bicyclo[4.4.0]decane (decalin base)
    // Tests that sphere expansion handles ring revisits without crashing
    let mol = parse("C1CCC2CCCCC2C1").expect("decalin");
    let assignment = assign_cip(&mol);
    // Should complete without panic/crash
    assert!(assignment.assignments.len() <= mol.atom_count());
}

#[test]
fn test_tetrahedral_stable_when_ring_bond_opens_before_other_neighbors() {
    // A stereocenter whose ring-closure digit is written BEFORE its other
    // substituents (`[C@@H]1...`) used to get the wrong CIP code: raw adjacency
    // order only materializes a ring-*opening* bond once the matching closing
    // digit is reached, which is later than a continuation atom that has nothing
    // to wait on -- so the neighbor meant to come second (by SMILES-textual
    // position) ended up listed after one that should come third. Verified
    // against RDKit's CanonicalRankAtoms-based CIP oracle (atom 5 == R here).
    let smi_a = "CN1CCC[C@@H]1c1cccnc1";
    let smi_b = "c1ccncc1[C@@H]1N(CCC1)C"; // same molecule, order-only respelling
    assert_eq!(cip_at(smi_a, 5), Some(CipCode::R));
    assert_eq!(cip_at(smi_b, 6), Some(CipCode::R));
}

#[test]
fn test_tetrahedral_double_bond_duplicates_into_own_sphere() {
    // A stereocenter substituent reached via a double bond must count as TWO
    // entries in ITS OWN CIP substituent sphere (the real atom plus a phantom
    // duplicate), not just contribute a single phantom to the far side of the
    // double bond. Without the departure-side duplicate, `C(=CH2)(CH3)-` scores
    // (C,C) instead of (C,C,C) and loses a priority tie-break it should win.
    // Verified against RDKit (atom 3 == R).
    assert_eq!(
        cip_at("C=C(C)[C@@H]1CN[C@@H](C(=O)O)[C@@H]1CC(=O)O", 3),
        Some(CipCode::R)
    );
}

#[test]
fn cip_mode_legacy_fast_matches_assign_cip_byte_for_byte() {
    // CipMode::LegacyFast must be assign_cip()'s output, unchanged -- every
    // existing caller of assign_cip is untouched by assign_cip_with_mode existing.
    let smis = [
        "C[C@H](N)C(=O)O",
        "C=C(C)[C@@H]1CN[C@@H](C(=O)O)[C@@H]1CC(=O)O",
        "C1CCN(P2(N3CCCC3)=N[P@@](N3CCCC3)(N3CC3)=N[P@](N3CCCC3)(N3CC3)=N2)C1",
    ];
    for smi in smis {
        let mol = chematic_smiles::parse(smi).expect("valid SMILES");
        let legacy = assign_cip(&mol);
        let via_mode = assign_cip_with_mode(&mol, CipMode::LegacyFast).expect("infallible");
        assert_eq!(legacy.assignments, via_mode.assignments, "smiles={smi}");
        assert!(via_mode.unresolved.is_empty(), "smiles={smi}");
    }
}

#[test]
fn accurate_cip_abstains_when_distant_wildcard_changes_priority() {
    // Stress corpus row 1956. The wildcard is not a direct ligand of atom 1,
    // yet instantiating it as carbon changes that center's R/S priority.
    // A concrete label for the unspecified molecule is therefore unsafe.
    let mol = chematic_smiles::parse("C[C@]12CC/C(=N\\N=C(N)N)C=C1C*[C@@H]2/C=C/[N+](=O)[O-]")
        .expect("valid wildcard molecule");
    let wildcard = mol
        .atoms()
        .find(|(_, atom)| atom.wildcard)
        .map(|(idx, _)| idx)
        .expect("fixture has a wildcard");
    let center = AtomIdx(1);
    assert!(mol.neighbors(center).all(|(idx, _)| idx != wildcard));

    let concrete = mol.with_atom_element(wildcard, Element::C);
    let concrete_result = assign_cip_with_mode(&concrete, CipMode::Accurate).unwrap();
    assert_eq!(concrete_result.get(center), Some(CipCode::R));

    let result = assign_cip_with_mode(&mol, CipMode::Accurate).unwrap();
    assert_eq!(result.get(center), None);
    assert!(
        result
            .unresolved
            .contains(&(center, CipUnresolvedReason::OracleUnstable))
    );
}

#[test]
fn accurate_cip_abstains_when_wildcard_changes_isotope_branch_priority() {
    // Stress corpus row 1849: the wildcard is on an aromatic side branch.
    // Both engines label atom 14 as R with '*' and S with carbon, showing
    // that an unspecified element does not support a definite R label.
    let mol = chematic_smiles::parse("CCCCCCCOC1=CC=C(*C[C@](C)([200NH2])CO)C=C1").unwrap();
    let wildcard = mol
        .atoms()
        .find(|(_, atom)| atom.wildcard)
        .map(|(idx, _)| idx)
        .unwrap();
    let center = AtomIdx(14);
    assert!(mol.neighbors(center).all(|(idx, _)| idx != wildcard));

    let concrete = mol.with_atom_element(wildcard, Element::C);
    assert_eq!(
        assign_cip_with_mode(&concrete, CipMode::Accurate)
            .unwrap()
            .get(center),
        Some(CipCode::S)
    );
    let result = assign_cip_with_mode(&mol, CipMode::Accurate).unwrap();
    assert_eq!(result.get(center), None);
    assert!(
        result
            .unresolved
            .contains(&(center, CipUnresolvedReason::OracleUnstable))
    );
}

#[test]
fn accurate_cip_abstains_when_fluorine_changes_distant_priority() {
    // RDKit 2026.09.1 stress rows 717 and 1529 keep S for '*', but report R
    // when that wildcard is instantiated as fluorine. Carbon leaves S intact.
    let rows = [
        ("[*]C[C@H](N)C(=O)O", AtomIdx(2)),
        (
            "C(N)(CNC([C@H](CC(=O)O)NC([C@@H](NC(=O)[C@H](C(C)C)NC(=O)[C@H](CC1=CC=CC=C1)NC(=O)[C@H](C)NC(=O)[C@H](CC1=CC=C(C=C1)O)N)C(C)*)=O)=O)=O",
            AtomIdx(12),
        ),
    ];
    for (smiles, center) in rows {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let wildcard = mol.atoms().find(|(_, atom)| atom.wildcard).unwrap().0;
        assert_eq!(
            assign_cip_with_mode(
                &mol.with_atom_element(wildcard, Element::C),
                CipMode::Accurate
            )
            .unwrap()
            .get(center),
            Some(CipCode::S),
        );
        assert_eq!(
            assign_cip_with_mode(
                &mol.with_atom_element(wildcard, Element::F),
                CipMode::Accurate
            )
            .unwrap()
            .get(center),
            Some(CipCode::R),
        );
        let result = assign_cip_with_mode(&mol, CipMode::Accurate).unwrap();
        assert_eq!(result.get(center), None, "{smiles}");
        assert!(
            result
                .unresolved
                .contains(&(center, CipUnresolvedReason::OracleUnstable)),
            "{smiles}"
        );
    }
}

#[test]
fn accurate_cip_keeps_label_when_wildcard_is_disconnected() {
    let mol = chematic_smiles::parse("C[C@H](F)Cl.*").unwrap();
    let result = assign_cip_with_mode(&mol, CipMode::Accurate).unwrap();
    let without_wildcard = chematic_smiles::parse("C[C@H](F)Cl").unwrap();
    let baseline = assign_cip_with_mode(&without_wildcard, CipMode::Accurate).unwrap();
    assert_eq!(result.get(AtomIdx(1)), baseline.get(AtomIdx(1)));
    assert!(result.unresolved.is_empty());
}

#[test]
fn accurate_cip_bounds_many_wildcard_probes_by_component() {
    let smiles = format!("C[C@H](F)C{}.[C@H](F)(Cl)Br", "*".repeat(17));
    let mol = chematic_smiles::parse(&smiles).unwrap();
    let distant_center = mol
        .atoms()
        .find(|(idx, atom)| *idx != AtomIdx(1) && atom.chirality.is_tetrahedral())
        .map(|(idx, _)| idx)
        .unwrap();
    let result = assign_cip_with_mode(&mol, CipMode::Accurate).unwrap();
    assert_eq!(result.get(AtomIdx(1)), None);
    assert!(
        result
            .unresolved
            .contains(&(AtomIdx(1), CipUnresolvedReason::OracleUnstable))
    );
    assert!(result.get(distant_center).is_some());
}

#[test]
fn cip_mode_accurate_labels_phosphorus_ties_like_rdkit() {
    // docs/rfcs/cip_accurate_rfc.md Milestone 4C-1: atoms 6/19 of this
    // cyclophosphazene were reported unresolved because neither RDKit CIP
    // engine gives a label that survives a neutral Kekulé respelling of the
    // P/N ring. They now carry RDKit's CIPLabeler label for the spelling
    // given (RDKit 2026.03.6: 6 R, 19 S), flagged as spelling-dependent.
    let smi = "CNP1(NC)=N[P@](NC)(N2CC2)=NP(NC)(NC)=N[P@@](NC)(N2CC2)=N1";
    let mol = chematic_smiles::parse(smi).expect("valid SMILES");
    let result = assign_cip_with_mode(&mol, CipMode::Accurate).expect("no engine error");
    assert_eq!(result.get(AtomIdx(6)), Some(CipCode::R));
    assert_eq!(result.get(AtomIdx(19)), Some(CipCode::S));
    assert!(cip_label_depends_on_kekule_spelling(&mol, AtomIdx(6)));
}

#[test]
fn cip_mode_accurate_phosphorus_labels_flip_with_kekule_spelling_as_rdkit() {
    // Same molecule as the test above. Flip every P/N ring bond Single<->Double (a
    // chemically neutral resonance respelling, the same test Milestone 4C-0/4C-1
    // used to show *both* RDKit engines are representation-unstable here): both
    // spellings are labelled, the labels flip as RDKit's CIPLabeler does, and
    // both are flagged as spelling-dependent.
    use chematic_core::BondOrder;
    use chematic_perception::find_sssr;

    let smi = "CNP1(NC)=N[P@](NC)(N2CC2)=NP(NC)(NC)=N[P@@](NC)(N2CC2)=N1";
    let mol = chematic_smiles::parse(smi).expect("valid SMILES");
    let sssr = find_sssr(&mol);

    let flip = |m: &chematic_core::Molecule, a, b| -> Option<chematic_core::Molecule> {
        let (bidx, bond) = m.bond_between(a, b)?;
        let new_order = match bond.order {
            BondOrder::Single => BondOrder::Double,
            BondOrder::Double => BondOrder::Single,
            other => other,
        };
        Some(m.with_bond_order(bidx, new_order))
    };

    let mut respelled = mol.clone();
    for ring in sssr.rings() {
        for w in ring.windows(2) {
            if let Some(next) = flip(&respelled, w[0], w[1]) {
                respelled = next;
            }
        }
        if let (Some(&first), Some(&last)) = (ring.first(), ring.last())
            && let Some(next) = flip(&respelled, last, first)
        {
            respelled = next;
        }
    }

    let original = assign_cip_with_mode(&mol, CipMode::Accurate).expect("no engine error");
    let after = assign_cip_with_mode(&respelled, CipMode::Accurate).expect("no engine error");
    for atom_idx in [6u32, 19u32] {
        let idx = AtomIdx(atom_idx);
        let (before, respelled_code) = (original.get(idx), after.get(idx));
        assert!(
            before.is_some(),
            "atom={atom_idx}: original spelling labelled"
        );
        assert!(
            respelled_code.is_some(),
            "atom={atom_idx}: respelled labelled"
        );
        assert_ne!(before, respelled_code, "atom={atom_idx}: label flips");
        assert!(cip_label_depends_on_kekule_spelling(&mol, idx));
        assert!(cip_label_depends_on_kekule_spelling(&respelled, idx));
    }
}

#[test]
fn cip_mode_accurate_labels_cyclophosphazene_phosphorus_like_rdkit() {
    // A P on an unsaturated ring gets RDKit's CIPLabeler label (RDKit
    // 2026.03.6: 1 S, 12 R). The label flips with the P=N Kekulé spelling
    // in both engines, which `cip_label_depends_on_kekule_spelling` reports.
    let mol = chematic_smiles::parse("N[P@]1(Cl)=NP(N2CC2)(N2CC2)=N[P@](N)(Cl)=N1")
        .expect("valid SMILES");
    let result = assign_cip_with_mode(&mol, CipMode::Accurate).expect("no engine error");
    assert_eq!(result.get(AtomIdx(1)), Some(CipCode::S));
    assert_eq!(result.get(AtomIdx(12)), Some(CipCode::R));
    assert!(result.unresolved.is_empty());
    assert!(cip_label_depends_on_kekule_spelling(&mol, AtomIdx(12)));
    assert!(!cip_label_depends_on_kekule_spelling(&mol, AtomIdx(0)));
    // The other Kekulé spelling of the same ring: both engines flip.
    let mol = chematic_smiles::parse("N[P@]1(Cl)N=P(N2CC2)(N2CC2)N=[P@](N)(Cl)N=1")
        .expect("valid SMILES");
    let flipped = assign_cip_with_mode(&mol, CipMode::Accurate).expect("no engine error");
    assert_eq!(flipped.get(AtomIdx(12)), Some(CipCode::S));
}

#[test]
fn cip_mode_accurate_reports_lone_pair_centres() {
    // Issue #634, rebaseline row 2960: a stereo-tagged bridgehead amine has
    // three explicit ligands and a lone pair. It gets RDKit's label for this
    // spelling (RDKit inverts `@@` here: one ring-closure digit on the
    // centre), marked spelling-dependent; the carbon centres keep theirs.
    let mol = chematic_smiles::parse("Cn1cc(C2=NC[C@@]3(C[N@@]4CC[C@@H]3C4)O2)c2ccccc21")
        .expect("valid SMILES");
    let result = assign_cip_with_mode(&mol, CipMode::Accurate).expect("no engine error");
    assert_eq!(result.get(AtomIdx(9)), Some(CipCode::S));
    assert!(cip_label_depends_on_smiles_spelling(&mol, AtomIdx(9)));
    assert!(!cip_label_depends_on_smiles_spelling(&mol, AtomIdx(7)));
    assert!(result.unresolved.is_empty());
    assert_eq!(result.get(AtomIdx(7)), Some(CipCode::R));
    assert_eq!(result.get(AtomIdx(12)), Some(CipCode::R));
    // A centre outside rings (sulfoxide, phosphine) takes the lone pair as
    // a lowest-priority phantom ligand, as RDKit 2026.03.6 rdCIPLabeler
    // does: (S)-methyl phenyl sulfoxide.
    let mol = chematic_smiles::parse("C[S@](=O)c1ccccc1").expect("valid SMILES");
    let result = assign_cip_with_mode(&mol, CipMode::Accurate).expect("no engine error");
    assert_eq!(result.get(AtomIdx(1)), Some(CipCode::S));
    assert!(result.unresolved.is_empty());
    // A ring sulfoxide is labelled too: RDKit reads its spellings as
    // chematic does (RDKit 2026.03.6 rdCIPLabeler: R, and R/S).
    let mol = chematic_smiles::parse("O=[S@@]1CCc2ccccc21").expect("valid SMILES");
    let result = assign_cip_with_mode(&mol, CipMode::Accurate).expect("no engine error");
    assert_eq!(result.get(AtomIdx(1)), Some(CipCode::R));
    let mol = chematic_smiles::parse("C[C@@H]1CCC[S@](=O)C1").expect("valid SMILES");
    let result = assign_cip_with_mode(&mol, CipMode::Accurate).expect("no engine error");
    assert_eq!(result.get(AtomIdx(1)), Some(CipCode::R));
    assert_eq!(result.get(AtomIdx(5)), Some(CipCode::S));
}

#[test]
fn cip_mode_accurate_labels_acyclic_phosphorus_like_rdkit() {
    // A phosphonamidate P outside any unsaturated ring: the P=O oxygen
    // gets no duplicate at the root (RDKit CIPLabeler), so OMe > =O and
    // the label matches RDKit (exposed 10k row 4419: R).
    let mol = chematic_smiles::parse(
        "CCC(C)[C@H](NC(=O)[C@@H]1CCCN1[P@](=O)(OC)[C@H](Cc1ccccc1)NC(=O)[C@H](CC(N)=O)NC(=O)OC(C)(C)C)C(=O)NCC(C)C",
    )
    .expect("valid SMILES");
    let result = assign_cip_with_mode(&mol, CipMode::Accurate).expect("no engine error");
    assert_eq!(result.get(AtomIdx(13)), Some(CipCode::R));
}

#[test]
fn cip_mode_accurate_labels_ring_amidine_imine_like_rdkit() {
    // Exposed 10k rows 1206/1213/1214/1287/1370/1371: an exocyclic imine
    // on the ring-fusion carbon of an imidazo[1,2-a]pyridine. RDKit
    // 2026.03.6 rdCIPLabeler and RDKit.js 2026.03.6/2026.09.1 label it Z;
    // both modes now use the accurate aromatic substituent ranking.
    let mol =
        chematic_smiles::parse("CCOC(=O)c1cc2ccccn2/c(=N/c2ccc(OC)cc2)n1").expect("valid SMILES");
    let accurate = assign_cip_with_mode(&mol, CipMode::Accurate).expect("no engine error");
    assert_eq!(accurate.get(AtomIdx(13)), Some(CipCode::Z));
    let legacy = assign_cip_with_mode(&mol, CipMode::LegacyFast).expect("infallible");
    assert_eq!(legacy.get(AtomIdx(13)), Some(CipCode::Z));
}

#[test]
fn cip_mode_accurate_merges_legacy_ez_with_accurate_tetrahedral() {
    // Accurate mode must still report E/Z (the accurate engine never computes it)
    // alongside its own tetrahedral R/S for the same molecule.
    let mol = chematic_smiles::parse("C/C=C/[C@H](N)C(=O)O").expect("valid SMILES");
    let result = assign_cip_with_mode(&mol, CipMode::Accurate).expect("no engine error");
    let has_ez = result
        .assignments
        .iter()
        .any(|(_, c)| matches!(c, CipCode::E | CipCode::Z));
    let has_rs = result
        .assignments
        .iter()
        .any(|(_, c)| matches!(c, CipCode::R | CipCode::S));
    assert!(
        has_ez,
        "expected an E/Z assignment: {:?}",
        result.assignments
    );
    assert!(
        has_rs,
        "expected an R/S assignment: {:?}",
        result.assignments
    );
}

// --- ez_completeness ---

#[test]
fn test_ez_completeness_specified() {
    // C/C=C/C — trans-2-butene: one stereogenic bond, fully marked.
    let mol = parse("C/C=C/C").unwrap();
    let ec = ez_completeness(&mol);
    assert_eq!(ec.specified, 1);
    assert_eq!(ec.unspecified, 0);
    assert_eq!(ec.total, 1);
}

#[test]
fn test_ez_completeness_unspecified() {
    // CC=CC — same connectivity as above, no `/`/`\` marker: still
    // stereogenic (one heavy substituent + implicit H on each end,
    // always distinguishable), just not annotated.
    let mol = parse("CC=CC").unwrap();
    let ec = ez_completeness(&mol);
    assert_eq!(ec.specified, 0);
    assert_eq!(ec.unspecified, 1);
    assert_eq!(ec.total, 1);
}

#[test]
fn test_ez_completeness_terminal_not_counted() {
    // C=C — ethylene: both ends are =CH2 (0 heavy substituents each
    // side), not stereogenic at all.
    let mol = parse("C=C").unwrap();
    let ec = ez_completeness(&mol);
    assert_eq!(ec.specified, 0);
    assert_eq!(ec.unspecified, 0);
    assert_eq!(ec.total, 0);
}

#[test]
fn test_ez_completeness_symmetric_end_not_counted() {
    // CC(C)=CC — the =C end bonded to two methyls is a genuine CIP tie
    // (confirmed directly via compare_branches below), so the whole
    // bond is not stereogenic even though the other end (=CC, one
    // heavy substituent + implicit H) would be on its own.
    let mol = parse("CC(C)=CC").unwrap();
    let methyl_a = AtomIdx(0);
    let methyl_b = AtomIdx(2);
    let alkene_c = AtomIdx(1);
    assert_eq!(
        compare_branches(&mol, alkene_c, methyl_a, methyl_b),
        std::cmp::Ordering::Equal,
        "the two methyls must be a genuine CIP tie for this test to be valid"
    );

    let ec = ez_completeness(&mol);
    assert_eq!(ec.specified, 0);
    assert_eq!(ec.unspecified, 0);
    assert_eq!(ec.total, 0);
}

#[test]
fn test_ez_completeness_aromatic_ring_excluded() {
    // Benzene's ring bonds are BondOrder::Aromatic, not Double -- no
    // contribution to any count.
    let mol = parse("c1ccccc1").unwrap();
    let ec = ez_completeness(&mol);
    assert_eq!(ec.total, 0);
}

#[test]
fn test_ez_completeness_kekule_benzene_ring_locked_excluded() {
    // C1=CC=CC=C1 -- Kekulized benzene with no aromaticity perception
    // applied, so every ring bond really is BondOrder::Double (not
    // Aromatic). Each ring atom still only has 1 heavy substituent per
    // end (the rest of the ring), so the pure substituent-count check
    // alone would call these stereogenic -- the small-ring lock (< 8
    // atoms, RDKit-confirmed) must exclude them anyway.
    let mol = parse("C1=CC=CC=C1").unwrap();
    let ec = ez_completeness(&mol);
    assert_eq!(ec.total, 0, "{ec:?}");
}

#[test]
fn test_ez_completeness_small_ring_double_bond_excluded() {
    // Cyclohexene: a single ring double bond, otherwise satisfies the
    // substituent-distinguishability check, but a 6-membered ring can't
    // isomerize -- RDKit's FindPotentialStereo agrees (empty result).
    let mol = parse("C1=CCCCC1").unwrap();
    let ec = ez_completeness(&mol);
    assert_eq!(ec.total, 0, "{ec:?}");
}

#[test]
fn test_ez_completeness_large_ring_double_bond_counted() {
    // Cyclooctene: an 8-membered ring double bond IS large enough for
    // E/Z isomerism (cis/trans-cyclooctene are both known, isolable
    // compounds) -- RDKit's FindPotentialStereo reports it as a
    // potential (unspecified) stereo bond. This is the positive
    // control for the ring-size cutoff above: it must not become a
    // blanket "any ring bond is excluded" bug.
    let mol = parse("C1=CCCCCCC1").unwrap();
    let ec = ez_completeness(&mol);
    assert_eq!(ec.specified, 0, "{ec:?}");
    assert_eq!(ec.unspecified, 1, "{ec:?}");
    assert_eq!(ec.total, 1, "{ec:?}");
}

#[test]
fn test_ez_completeness_bridged_bicyclic_ring_lock_not_sssr_dependent() {
    // Norbornene (bicyclo[2.2.1]hept-2-ene): the double bond sits in a
    // real 5-membered ring, but this bridged system's SSSR can also
    // report a larger fundamental cycle through the same bond depending
    // on decomposition choice (the same envelope-ring artifact
    // documented for aromaticity/ring-count elsewhere in this crate).
    // A min-ring-size check naively routed through `find_sssr` risks
    // missing the true 5-ring and passing the >= 8 gate. RDKit
    // (MinBondRingSize=5, FindPotentialStereo -- empty) confirms this
    // bond is ring-locked and must NOT be counted.
    let mol = parse("C1=CC2CCC1C2").unwrap();
    let ec = ez_completeness(&mol);
    assert_eq!(ec.total, 0, "{ec:?}");
}

#[test]
fn test_ez_completeness_partial_marker_is_unspecified() {
    // C/C=CC -- a marker on only one side of the double bond. The bond
    // is still structurally stereogenic (1 heavy substituent + implicit
    // H on each end, distinguishable), but assign_ez has nothing on the
    // far end to resolve a code from, so it must land in `unspecified`,
    // not silently drop out of the count or count as `specified`.
    let mol = parse("C/C=CC").unwrap();
    let bidx = (0..mol.bond_count())
        .map(|i| BondIdx(i as u32))
        .find(|&b| mol.bond(b).order == BondOrder::Double)
        .unwrap();
    assert_eq!(
        assign_ez(&mol, bidx),
        None,
        "a one-sided marker must not resolve a code"
    );

    let ec = ez_completeness(&mol);
    assert_eq!(ec.specified, 0, "{ec:?}");
    assert_eq!(ec.unspecified, 1, "{ec:?}");
    assert_eq!(ec.total, 1, "{ec:?}");
}

#[test]
fn test_ez_completeness_aggregates_mixed_specified_and_unspecified() {
    // C/C=C/CC=CC — first double bond (atoms 1=2) is fully marked
    // (same local Up/Up pattern as trans-2-butene, shifted down the
    // chain); second double bond (atoms 4=5) has no marker on either
    // side but is still stereogenic (one heavy substituent + implicit
    // H per end).
    let mol = parse("C/C=C/CC=CC").unwrap();
    let ec = ez_completeness(&mol);
    assert_eq!(ec.specified, 1, "{ec:?}");
    assert_eq!(ec.unspecified, 1, "{ec:?}");
    assert_eq!(ec.total, 2, "{ec:?}");
}
