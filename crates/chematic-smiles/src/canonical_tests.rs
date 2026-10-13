#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse;

    /// Positive control (section 3/19): directly verify -- by calling the
    /// legacy exhaustive enumeration itself, not by trusting a doc comment
    /// -- that budget exhaustion in the *old* design silently returns a
    /// partial/wrong result rather than failing closed. With `budget = 0`,
    /// `enumerate_discrete_ranks`'s very first non-singleton cell hits the
    /// `*budget == 0 { break; }` guard on its first loop iteration and
    /// returns an empty `Vec`; the caller's `.min_by(..).unwrap_or_default()`
    /// then silently produces `("", vec![])` -- an EMPTY canonical SMILES,
    /// not an error. This is exactly the anti-pattern the new
    /// `canonical_smiles_with_limits`/`CanonicalizationError::
    /// SearchBudgetExceeded` API (in `canonical_search.rs`) exists to
    /// retire for callers who opt into a bounded search.
    #[test]
    fn old_design_zero_budget_silently_returns_empty_string_not_an_error() {
        let mol = parse("c1ccccc1").unwrap(); // benzene: 1 non-singleton cell of all 6 atoms
        let plateaued = morgan_ranks(&mol);
        let mut budget = 0usize;
        let branches = enumerate_discrete_ranks(&mol, plateaued, &mut budget);
        assert!(
            branches.is_empty(),
            "zero budget should yield zero explored branches from the old design"
        );

        // Reproduce exactly what `legacy_winning_individualized_ranks` (and,
        // before this PR, `winning_individualized_ranks`) does with that
        // empty branch list.
        let (ranks, s) = branches
            .into_iter()
            .map(|r: Vec<u64>| {
                let s = CanonicalWriter::new(&mol, &r).write_all();
                (r, s)
            })
            .min_by(|(_, a), (_, b)| a.cmp(b))
            .unwrap_or_default();
        assert_eq!(
            s, "",
            "old design: zero budget silently yields an EMPTY canonical string"
        );
        assert!(ranks.is_empty());

        // Confirm the NEW engine does NOT reproduce this failure mode: the
        // same molecule under an equivalently tiny bounded search fails
        // closed (a typed error), never an empty string.
        let bounded = crate::canonical_search::canonical_smiles_with_limits(
            &mol,
            &crate::canonical_search::CanonicalizationLimits {
                max_search_nodes: Some(0),
                max_automorphism_tests: None,
            },
        );
        assert!(
            matches!(
                bounded,
                Err(crate::canonical_search::CanonicalizationError::SearchBudgetExceeded { .. })
            ),
            "new engine must fail closed, got {bounded:?}"
        );
    }

    /// Build a copy of `mol` with atoms reordered by `perm` (perm[new_idx] = old_idx).
    /// Bonds are remapped to the new indices; stereo/direction metadata is
    /// intentionally dropped since this helper only exists to test whether the
    /// *skeleton* rank partition is invariant under atom relabeling.
    fn permute_molecule(mol: &Molecule, perm: &[usize]) -> Molecule {
        let mut old_to_new = vec![0u32; perm.len()];
        for (new_idx, &old_idx) in perm.iter().enumerate() {
            old_to_new[old_idx] = new_idx as u32;
        }
        let mut builder = chematic_core::MoleculeBuilder::new();
        for &old_idx in perm {
            builder.add_atom(mol.atom(AtomIdx(old_idx as u32)).clone());
        }
        for (_, bond) in mol.bonds() {
            let a = AtomIdx(old_to_new[bond.atom1.0 as usize]);
            let b = AtomIdx(old_to_new[bond.atom2.0 as usize]);
            let _ = builder.add_bond(a, b, bond.order);
        }
        builder.build()
    }

    /// Like [`permute_molecule`], but preserves bond order (including
    /// literal `Up`/`Down` E/Z markers) and the `bond_direction` stash,
    /// remapped onto the new atom labeling -- needed to test canonical
    /// *convergence* across genuinely different atom-labelings of the same
    /// molecule (unlike `permute_molecule`, which deliberately drops this
    /// metadata since it only tests skeleton rank-partition invariance).
    /// Bonds are added in the same relative order as `mol.bonds()` yields
    /// them, so each new bond's index equals its original one -- the
    /// `bond_direction` stash (keyed by `BondIdx`) can therefore be
    /// reapplied directly after `build()`.
    fn relabel_molecule_preserving_ez(mol: &Molecule, perm: &[usize]) -> Molecule {
        let mut old_to_new = vec![0u32; perm.len()];
        for (new_idx, &old_idx) in perm.iter().enumerate() {
            old_to_new[old_idx] = new_idx as u32;
        }
        let mut builder = chematic_core::MoleculeBuilder::new();
        for &old_idx in perm {
            builder.add_atom(mol.atom(AtomIdx(old_idx as u32)).clone());
        }
        let mut direction_stash = Vec::new();
        for (bidx, bond) in mol.bonds() {
            let a = AtomIdx(old_to_new[bond.atom1.0 as usize]);
            let b = AtomIdx(old_to_new[bond.atom2.0 as usize]);
            let new_bidx = builder
                .add_bond(a, b, bond.order)
                .expect("relabeling a valid molecule's own bonds cannot fail");
            if let Some(dir) = mol.bond_direction(bidx) {
                let anchor = mol
                    .bond_direction_anchor(bidx)
                    .map(|old| AtomIdx(old_to_new[old.0 as usize]));
                direction_stash.push((new_bidx, dir, anchor));
            }
        }
        // Tetrahedral chirality (`Atom.chirality`) is meaningless without
        // its accompanying `stereo_neighbor_order` (the original neighbor
        // listing @/@@ is interpreted against) -- remapped here the same
        // way `Molecule::remove_atom` remaps it on atom removal (see that
        // method's own doc comment), just with no removed atom to drop.
        for &old_idx in perm {
            let old_atom = AtomIdx(old_idx as u32);
            if let Some(order) = mol.stereo_neighbor_order(old_atom) {
                let new_order: Vec<u32> = order
                    .iter()
                    .map(|&v| {
                        if v == chematic_core::STEREO_H_SENTINEL {
                            v
                        } else {
                            old_to_new[v as usize]
                        }
                    })
                    .collect();
                builder.set_stereo_neighbor_order(AtomIdx(old_to_new[old_idx]), new_order);
            }
        }
        let mut relabeled = builder.build();
        for (bidx, dir, anchor) in direction_stash {
            relabeled.set_bond_direction(bidx, dir);
            if let Some(anchor) = anchor {
                relabeled.set_bond_direction_anchor(bidx, anchor);
            }
        }
        relabeled
    }

    /// Relabel a rank vector into "first-seen order" group ids, so partitions
    /// can be compared structurally (which atoms are grouped together)
    /// independent of the actual numeric rank values assigned.
    fn partition_key(ranks: &[u64]) -> Vec<usize> {
        let mut seen: Vec<u64> = Vec::new();
        ranks
            .iter()
            .map(|&r| match seen.iter().position(|&s| s == r) {
                Some(pos) => pos,
                None => {
                    seen.push(r);
                    seen.len() - 1
                }
            })
            .collect()
    }

    /// Sanity check demanded before any individualize-refine rewrite: the
    /// refinement-to-plateau partition (which atoms share a rank, not the raw
    /// numeric values) must be invariant under input atom permutation. If this
    /// fails, the root cause is a non-invariant initial invariant / refinement
    /// step itself (over-splitting orbits), not a missing individualize step
    /// (under-splitting ties) -- a completely different bug to fix.
    #[test]
    fn morgan_ranks_partition_is_permutation_invariant() {
        let corpus = [
            "O=C(NCc1cccnc1)NC[C@H]1CCC[C@H](OCc2cc(C(F)(F)F)cc(C(F)(F)F)c2)[C@@H]1c1ccccc1",
            "c1ccccc1",
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "O=C1CCC(=O)N1",
            "c1ccc2ccc3ccccc3c2c1",
            "O=C(NCc1cccnc1)NCC1CCCC(OCc2cc(C(F)(F)F)cc(C(F)(F)F)c2)C1c1ccccc1",
        ];

        for smi in corpus {
            let mol = parse(smi).unwrap_or_else(|e| panic!("{smi}: {e}"));
            let n = mol.atom_count();
            let part_orig = partition_key(&morgan_ranks(&mol));

            // A few deterministic, non-identity permutations (no RNG dependency).
            let perms: Vec<Vec<usize>> = vec![
                (0..n).rev().collect(),
                {
                    let mut p: Vec<usize> = (0..n).collect();
                    if n > 2 {
                        p.rotate_left(n / 3 + 1);
                    }
                    p
                },
                {
                    let mut p: Vec<usize> = (0..n).rev().collect();
                    if n > 3 {
                        p.swap(1, n - 2);
                        p.rotate_right(2);
                    }
                    p
                },
            ];

            for perm in perms {
                let permuted = permute_molecule(&mol, &perm);
                let part_perm = partition_key(&morgan_ranks(&permuted));

                // inverse: where did old atom `old_idx` land in the permuted molecule?
                let mut new_of_old = vec![0usize; n];
                for (new_idx, &old_idx) in perm.iter().enumerate() {
                    new_of_old[old_idx] = new_idx;
                }

                for i in 0..n {
                    for j in (i + 1)..n {
                        let same_orig = part_orig[i] == part_orig[j];
                        let same_perm = part_perm[new_of_old[i]] == part_perm[new_of_old[j]];
                        assert_eq!(
                            same_orig, same_perm,
                            "partition not permutation-invariant for '{smi}': \
                             atoms {i},{j} (perm {perm:?})"
                        );
                    }
                }
            }
        }
    }

    /// `canonical_atom_order` must be permutation-invariant: relabeling the
    /// same molecule's atoms (different parse order) must not change WHICH
    /// symmetry class of atom appears 1st, 2nd, 3rd, ... in the returned
    /// order. Unlike `canonical_smiles`, `canonical_atom_order` does not run
    /// individualize-refine -- it sorts raw `morgan_ranks` with no tie-break,
    /// so this is expected to fail on any molecule with a genuine
    /// (non-singleton) rank tie. This is a diagnostic probe for that gap, not
    /// an already-passing invariant.
    #[test]
    fn canonical_atom_order_permutation_invariance_probe() {
        let corpus = [
            "O=C(NCc1cccnc1)NC[C@H]1CCC[C@H](OCc2cc(C(F)(F)F)cc(C(F)(F)F)c2)[C@@H]1c1ccccc1",
            "c1ccccc1",
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "O=C1CCC(=O)N1",
            "c1ccc2ccc3ccccc3c2c1",
            "O=C(NCc1cccnc1)NCC1CCCC(OCc2cc(C(F)(F)F)cc(C(F)(F)F)c2)C1c1ccccc1",
        ];

        let mut bad = 0;
        let mut total = 0;
        for smi in corpus {
            let mol = parse(smi).unwrap_or_else(|e| panic!("{smi}: {e}"));
            let n = mol.atom_count();
            // Ground-truth atom-class labels, in the ORIGINAL molecule's index
            // space only -- never relabeled independently for the permuted
            // copy, or two arbitrary first-seen-order numberings would be
            // compared against each other and any mismatch would be a
            // methodology artifact, not a real instability (see the
            // project's canonicalization-tie-break-theory / measurement-
            // harness-controls lessons).
            let part_orig = partition_key(&morgan_ranks(&mol));
            let order_orig = canonical_atom_order(&mol);
            let profile_orig: Vec<usize> = order_orig.iter().map(|&i| part_orig[i]).collect();

            let perms: Vec<Vec<usize>> = vec![(0..n).rev().collect(), {
                let mut p: Vec<usize> = (0..n).collect();
                if n > 2 {
                    p.rotate_left(n / 3 + 1);
                }
                p
            }];

            for perm in perms {
                total += 1;
                // perm[new_idx] = old_idx (see permute_molecule's contract).
                let permuted = permute_molecule(&mol, &perm);
                let order_perm = canonical_atom_order(&permuted);
                // Map each returned NEW index back to the class of the
                // corresponding OLD atom, via `part_orig` -- the same
                // ground-truth labeling used for profile_orig.
                let profile_perm: Vec<usize> = order_perm
                    .iter()
                    .map(|&new_i| part_orig[perm[new_i]])
                    .collect();

                if profile_orig != profile_perm {
                    bad += 1;
                    eprintln!(
                        "canonical_atom_order NOT permutation-invariant for '{smi}' (perm {perm:?}): \
                         {profile_orig:?} != {profile_perm:?}"
                    );
                }
            }
        }
        eprintln!("canonical_atom_order instability: {bad}/{total} permutation trials");
        assert_eq!(
            bad, 0,
            "{bad}/{total} permutation trials were unstable -- see stderr"
        );
    }

    /// Direct probe (no permutation needed): does `canonical_atom_order`'s
    /// naive `morgan_ranks`-only sort ever disagree with the FULLY
    /// individualized/resolved rank order that `canonical_smiles` actually
    /// verified-correct output is built from? Comparing against
    /// same-partition-class labels (as the permutation probe above does) is
    /// blind to intra-class reordering among genuinely symmetric
    /// (automorphism-equivalent) atoms, where any order is harmless -- this
    /// test instead reconstructs the winning individualized branch directly,
    /// so it also catches disagreement WITHIN a Morgan-rank tie that is not
    /// a true automorphism (exactly the class of bug individualize-refine
    /// was built to fix for `canonical_smiles`, see Round 10-12 history).
    #[test]
    fn canonical_atom_order_matches_individualized_ranks_probe() {
        let corpus = [
            "O=C(NCc1cccnc1)NC[C@H]1CCC[C@H](OCc2cc(C(F)(F)F)cc(C(F)(F)F)c2)[C@@H]1c1ccccc1",
            "c1ccccc1",
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "O=C1CCC(=O)N1",
            "c1ccc2ccc3ccccc3c2c1",
            "O=C(NCc1cccnc1)NCC1CCCC(OCc2cc(C(F)(F)F)cc(C(F)(F)F)c2)C1c1ccccc1",
            // Extra cases picked for symmetry that Weisfeiler-Leman-style
            // refinement is known to struggle with (fused/bridged systems).
            "C1CC2CCC1CC2",
            "C1CC2CC1CC2",
            "c1ccc(-c2ccccc2)cc1",
            "OC1CCC(O)CC1",
            "C12CC3CC(CC(C3)C1)C2",
        ];

        let mut needed_individualization = 0;
        let mut mismatched = 0;
        for smi in corpus {
            let mol = parse(smi).unwrap_or_else(|e| panic!("{smi}: {e}"));
            let plateaued = morgan_ranks(&mol);
            let mut budget = MAX_INDIVIDUALIZE_BRANCHES;
            let branches = enumerate_discrete_ranks(&mol, plateaued, &mut budget);
            if branches.len() > 1 {
                needed_individualization += 1;
            }
            let winning_ranks = branches
                .into_iter()
                .min_by_key(|ranks| CanonicalWriter::new(&mol, ranks).write_all())
                .expect("at least one branch");

            let n = mol.atom_count();
            let mut winning_order: Vec<usize> = (0..n).collect();
            winning_order.sort_by(|&a, &b| winning_ranks[b].cmp(&winning_ranks[a]));

            let naive_order = canonical_atom_order(&mol);
            if naive_order != winning_order {
                mismatched += 1;
                eprintln!(
                    "canonical_atom_order disagrees with resolved canonical order for '{smi}': \
                     naive={naive_order:?} resolved={winning_order:?}"
                );
            }
        }
        eprintln!(
            "{needed_individualization}/{} molecules needed individualization; \
             {mismatched}/{} disagreed with canonical_atom_order",
            corpus.len(),
            corpus.len()
        );
        assert_eq!(
            mismatched, 0,
            "canonical_atom_order must match the individualized/resolved order -- see stderr"
        );
    }

    /// Canonical SMILES must be stable: applying it twice gives the same result.
    fn is_stable(smiles: &str) -> bool {
        let mol1 = parse(smiles).expect(smiles);
        let c1 = canonical_smiles(&mol1);
        assert!(
            !c1.is_empty(),
            "canonical_smiles returned empty for '{smiles}'"
        );
        let mol2 =
            parse(&c1).unwrap_or_else(|e| panic!("canonical SMILES '{c1}' is not parseable: {e}"));
        let c2 = canonical_smiles(&mol2);
        c1 == c2
    }

    /// Two SMILES representing the same molecule must give the same canonical form.
    fn same_canonical(a: &str, b: &str) -> bool {
        let mol_a = parse(a).expect(a);
        let mol_b = parse(b).expect(b);
        canonical_smiles(&mol_a) == canonical_smiles(&mol_b)
    }

    fn high_cycle_fan() -> Molecule {
        // 103 atoms and 204 bonds give a cycle rank of 102, requiring extended
        // ring labels. Unique isotopes keep this a deterministic writer test
        // instead of a symmetry-search stress test.
        let mut builder = chematic_core::MoleculeBuilder::with_capacity(103, 204);
        let center = builder.add_atom(chematic_core::Atom::new(chematic_core::Element::C));
        let mut rim = Vec::with_capacity(102);
        for isotope in 1..=102 {
            let mut atom = chematic_core::Atom::new(chematic_core::Element::C);
            atom.isotope = Some(isotope);
            rim.push(builder.add_atom(atom));
        }
        for &atom in &rim {
            builder.add_bond(center, atom, BondOrder::Single).unwrap();
        }
        for edge in rim.windows(2) {
            builder
                .add_bond(edge[0], edge[1], BondOrder::Single)
                .unwrap();
        }
        builder.build()
    }

    #[test]
    fn canonical_writer_round_trips_extended_ring_labels() {
        let mol = high_cycle_fan();
        let canonical = canonical_smiles(&mol);
        assert!(canonical.contains("%(100)"), "{canonical}");
        let reparsed = parse(&canonical).expect("canonical extended labels must parse");
        assert_eq!(reparsed.atom_count(), mol.atom_count());
        assert_eq!(reparsed.bond_count(), mol.bond_count());
        assert_eq!(canonical_smiles(&reparsed), canonical);
    }

    #[test]
    fn test_methane_stable() {
        assert!(is_stable("C"));
    }
    #[test]
    fn test_ethane_stable() {
        assert!(is_stable("CC"));
    }
    #[test]
    fn test_ethanol_stable() {
        assert!(is_stable("CCO"));
    }
    #[test]
    fn test_acetic_acid_stable() {
        assert!(is_stable("CC(=O)O"));
    }
    #[test]
    fn test_benzene_stable() {
        assert!(is_stable("c1ccccc1"));
    }
    #[test]
    fn test_pyridine_stable() {
        assert!(is_stable("c1ccncc1"));
    }
    #[test]
    fn test_naphthalene_stable() {
        assert!(is_stable("c1ccc2ccccc2c1"));
    }
    #[test]
    fn test_aspirin_stable() {
        assert!(is_stable("CC(=O)Oc1ccccc1C(=O)O"));
    }
    #[test]
    fn test_caffeine_stable() {
        assert!(is_stable("Cn1cnc2c1c(=O)n(c(=O)n2C)C"));
    }

    #[test]
    fn test_ethanol_same_from_different_starts() {
        assert!(same_canonical("CCO", "OCC"));
    }

    #[test]
    fn test_isobutane_same_canonical() {
        // CC(C)C and C(C)(C)C are the same molecule.
        assert!(same_canonical("CC(C)C", "C(C)(C)C"));
    }

    #[test]
    fn test_wildcard_roundtrip() {
        let mol = parse("[*]CC").unwrap();
        let c = canonical_smiles(&mol);
        assert!(!c.is_empty());
        let mol2 = parse(&c).unwrap();
        assert_eq!(mol.atom_count(), mol2.atom_count());
        assert!(is_stable("[*]CC"));
    }

    #[test]
    fn canonical_wildcard_preserves_atom_map() {
        let mol = parse("[*:17]C").unwrap();
        let canonical = canonical_smiles(&mol);
        assert!(
            canonical.contains("[*:17]"),
            "canonical wildcard output must retain its atom map: {canonical}"
        );
        let reparsed = parse(&canonical).unwrap();
        let wildcard_maps = reparsed
            .atoms()
            .filter_map(|(_, atom)| atom.wildcard.then_some(atom.atom_map))
            .collect::<Vec<_>>();
        assert_eq!(wildcard_maps, vec![Some(17)]);
        assert_eq!(canonical_smiles(&reparsed), canonical);
    }

    #[test]
    fn test_disconnected_stable() {
        assert!(is_stable("[Na+].[Cl-]"));
    }

    // E/Z stereo bond direction tests.
    #[test]
    fn test_ez_e_stable() {
        assert!(is_stable("C/C=C/C"));
    }
    #[test]
    fn test_ez_z_stable() {
        assert!(is_stable("C/C=C\\C"));
    }
    #[test]
    fn test_ez_fluoro_e_stable() {
        assert!(is_stable("F/C=C/Cl"));
    }
    #[test]
    fn test_ez_fluoro_z_stable() {
        assert!(is_stable("F/C=C\\Cl"));
    }
    #[test]
    fn test_ez_e_ne_z() {
        // E and Z isomers of 1-fluoro-2-chloroethylene must yield different canonical forms.
        let mol_e = parse("F/C=C/Cl").unwrap();
        let mol_z = parse("F/C=C\\Cl").unwrap();
        assert_ne!(canonical_smiles(&mol_e), canonical_smiles(&mol_z));
    }

    // ── Tetrahedral stereo parity tests ─────────────────────────────────────

    #[test]
    fn test_tetrahedral_stable_no_from_atom() {
        // Bracket-H form at start of fragment — no from-atom.
        assert!(is_stable("[C@@H](F)(Cl)Br"));
        assert!(is_stable("[C@H](F)(Cl)Br"));
    }

    #[test]
    fn test_tetrahedral_stable_with_from_atom() {
        // L-alanine: chiral atom has a from-atom (N).
        assert!(is_stable("N[C@@H](C)C(=O)O"));
        assert!(is_stable("N[C@H](C)C(=O)O"));
    }

    #[test]
    fn test_enantiomers_differ() {
        // R and S configurations must give distinct canonical SMILES.
        assert!(!same_canonical("N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O"));
        assert!(!same_canonical("[C@@H](F)(Cl)Br", "[C@H](F)(Cl)Br"));
    }

    #[test]
    fn test_tetrahedral_same_from_different_starts() {
        // L-alanine from N vs from methyl — odd permutation, parity correction required.
        // RDKit: N[C@@H](C)C(=O)O and C[C@H](N)C(=O)O both → C[C@H](N)C(=O)O.
        assert!(same_canonical("N[C@@H](C)C(=O)O", "C[C@H](N)C(=O)O"));
        // D-alanine must differ from L-alanine.
        assert!(!same_canonical("N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O"));
    }

    #[test]
    fn test_rdkit_agreement_alanine() {
        // Pairs where the Morgan ranks distinguish all atoms unambiguously.
        // N[C@@H](C)C(=O)O and C[C@H](N)C(=O)O: same L-alanine (RDKit agrees).
        assert!(same_canonical("N[C@@H](C)C(=O)O", "C[C@H](N)C(=O)O"));
        // Enantiomers must differ (RDKit: C[C@@H](N)C(=O)O for D-alanine).
        assert!(!same_canonical("N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O"));
        // Stability: L-alanine canonical is self-stable.
        assert!(is_stable("N[C@@H](C)C(=O)O"));
        assert!(is_stable("C[C@H](N)C(=O)O"));
    }

    #[test]
    fn test_tetrahedral_all_heavy_substituents_stable() {
        // Chiral centre with no bracket H (all four heavy substituents).
        assert!(is_stable("[C@](F)(Cl)(Br)I"));
        assert!(is_stable("[C@@](F)(Cl)(Br)I"));
    }

    #[test]
    fn test_tetrahedral_all_heavy_enantiomers_differ() {
        assert!(!same_canonical("[C@](F)(Cl)(Br)I", "[C@@](F)(Cl)(Br)I"));
    }

    #[test]
    fn test_ring_stereocentre_stable() {
        // Chiral atom inside a ring — tests ring-closure partner resolution.
        assert!(is_stable("[C@@H]1CCCC1F"));
        assert!(is_stable("[C@H]1CCCC1F"));
    }

    #[test]
    fn test_ring_stereocentre_enantiomers_differ() {
        assert!(!same_canonical("[C@@H]1CCCC1F", "[C@H]1CCCC1F"));
    }

    #[test]
    fn test_chirality_from_different_entry_points() {
        // Same chiral molecule, two SMILES with different traversal order.
        // F[C@@H](Cl)Br  ≡  Cl[C@H](F)Br  (same S-configuration, just written
        // from different entry atoms — verified by signed-tetrahedral-volume).
        // Their canonical SMILES must be identical.
        let c1 = canonical_smiles(&parse("F[C@@H](Cl)Br").unwrap());
        let c2 = canonical_smiles(&parse("Cl[C@H](F)Br").unwrap());
        assert_eq!(c1, c2, "same molecule from different starts should match");

        // Cross-check: the enantiomer gives a different canonical form.
        let c3 = canonical_smiles(&parse("F[C@H](Cl)Br").unwrap());
        assert_ne!(c1, c3, "enantiomers must differ");
    }

    // ── Bond-order canonicality tests (#14 fix) ──────────────────────────

    #[test]
    fn test_acetic_acid_canonical_same_from_different_starts() {
        // Bug #14: both oxygens in acetic acid had the same Morgan rank because
        // the refinement loop omitted bond orders.  After the fix, O= (double)
        // and O-H (single) get distinct ranks regardless of atom insertion order.
        assert!(same_canonical("CC(=O)O", "OC(C)=O"));
        assert!(same_canonical("CC(=O)O", "O=C(O)C"));
        assert!(same_canonical("CC(=O)O", "C(C)(=O)O"));
    }

    #[test]
    fn test_oxygens_in_acetic_acid_not_equivalent() {
        // The two oxygens (O= vs O-H) are chemically distinct and must receive
        // different Morgan symmetry classes.
        let mol = parse("CC(=O)O").unwrap();
        let classes = equivalent_atom_classes(&mol);
        let o_classes: Vec<usize> = mol
            .atoms()
            .filter(|(_, a)| a.element.atomic_number() == 8)
            .map(|(i, _)| classes[i.0 as usize])
            .collect();
        assert_eq!(o_classes.len(), 2);
        assert_ne!(
            o_classes[0], o_classes[1],
            "O= and O-H must be in different symmetry classes"
        );
    }

    #[test]
    fn test_formic_acid_canonical_consistent() {
        // OC=O and O=CO — same formic acid, should canonicalize identically.
        assert!(same_canonical("OC=O", "O=CO"));
    }

    // ── RDKit PR #9066: conjugated E/Z round-trip ────────────────────────────

    #[test]
    fn conjugated_double_bond_ez_round_trip() {
        // RDKit PR #9066: removeRedundantBondDirSpecs() could strip bond directions
        // on conjugated double bonds, losing E/Z stereo.  Chematic does not apply
        // aggressive direction removal, but this test guards against regressions.
        for smi in &[
            r"F/C=C/C=C/Cl", // all-E conjugated diene
            r"F/C=C\C=C\Cl", // E then Z
            r"F/C=C/C=C\Cl", // E then inverted-Z
        ] {
            let mol = parse(smi).unwrap_or_else(|e| panic!("parse {smi}: {e:?}"));
            let out = canonical_smiles(&mol);
            let mol2 = parse(&out).unwrap_or_else(|e| panic!("re-parse {out}: {e:?}"));
            let out2 = canonical_smiles(&mol2);
            assert_eq!(
                out, out2,
                "conjugated E/Z must be stable after two rounds: {smi} → {out} → {out2}"
            );
        }
    }

    // ── Round 10: ring-closure directional bond flip ─────────────────────────
    //
    // A directional marker (`/`, `\`) is read "toward" the ring digit from
    // wherever it's written. At the ring-OPENING occurrence that's already
    // the open->close direction; at the CLOSING occurrence it's close->open
    // (the opposite traversal direction over the same physical bond) and must
    // be flipped before use (parser.rs `close_or_open_ring`). Before the fix,
    // the closing-side marker was stored raw/unflipped, which silently
    // produced a *different stereoisomer* whenever a random SMILES spelling
    // routed a conjugated system's connecting single bond through a
    // ring-closure digit instead of a plain adjacent chain bond -- confirmed
    // via a corpus-wide worst-of-10 sweep (RDKit-checked structural
    // correctness, not just self-stability/idempotency, which this bug class
    // passed trivially since it was deterministic-but-wrong on each input).

    #[test]
    fn ring_closure_direction_flip_real_world_repro() {
        // Real molecule found via corpus sweep. `variant` is an RDKit
        // doRandom=True re-spelling of the exact same molecule as `orig`,
        // routing the diene's connecting single bond through ring-closure
        // digit "1" instead of a plain chain bond. Before the parser fix,
        // chematic silently emitted a different (RDKit-confirmed
        // non-equivalent) stereoisomer for `variant`.
        let orig = r"CC1CCOC(=O)/C=C/C=C\C(=O)O[C@@H]2C[C@H]3O[C@@H]4C[C@@H](C)C(=O)C[C@]4(COC(=O)C1O)[C@]2(C)C31CO1";
        let variant = r"C1=C\C(=O)O[C@@H]2C[C@H]3O[C@H]4[C@@]([C@@]2(C32OC2)C)(CC(=O)[C@H](C)C4)COC(=O)C(O)C(C)CCOC(=O)/C=C/1";
        assert!(
            same_canonical(orig, variant),
            "ring-closure-routed diene must canonicalize identically to the \
             chain-form spelling of the same molecule"
        );
    }

    #[test]
    fn ring_closure_direction_minimal_ez_agreement() {
        // Minimal case isolating the same mechanism: a ring-closure bond
        // (distinct from the exocyclic C=C double bond itself) whose
        // directional markers are specified at BOTH the opening and closing
        // occurrences of the ring digit. Per the flip rule, opposite raw
        // symbols (one `/`, one `\`) describe one consistent bond and must
        // parse successfully; same-symbol at both ends is the conflicting
        // case (unchanged by this fix -- only Up/Down are flipped, so a
        // same-vs-different Double/Single conflict, e.g. "C=1CC-1", is
        // unaffected).
        let mol = parse(r"F/C=C/1CCCC\1").unwrap_or_else(|e| panic!("{e:?}"));
        let out = canonical_smiles(&mol);
        let mol2 = parse(&out).unwrap_or_else(|e| panic!("re-parse {out}: {e:?}"));
        assert_eq!(
            canonical_smiles(&mol2),
            out,
            "ring-closure E/Z with opposite-symbol agreement must round-trip stably"
        );

        // Same-symbol at both ends of a ring-closure directional bond is now
        // (correctly) the conflicting combination.
        assert!(matches!(
            parse(r"F/C=C/1CCCC/1"),
            Err(crate::error::SmilesError::ConflictingRingBond { ring_num: 1, .. })
        ));
    }

    // ── Round 10: ring-digit reuse racing PendingRing resolution ─────────────
    //
    // A stereocenter that OPENS a ring whose partner closes INSIDE the
    // stereocenter's own branch subtree (e.g. `[C@]1(...[C@H]1...)`) has its
    // own stereo record still unfinalized at the moment of that first
    // closure -- the immediate-resolution fast path in `close_or_open_ring`
    // only patches already-finalized records, so this case falls through to
    // the end-of-parse fallback. Before the fix, that fallback resolved by
    // raw ring DIGIT via `ring_close_partners: HashMap<u8, AtomIdx>` -- if the
    // same digit was reused later for an unrelated ring (e.g. a trailing
    // phenyl `c1ccccc1`), the later reuse's closer silently overwrote the
    // earlier, still-pending resolution, corrupting the stereocenter's
    // neighbor order with a foreign atom index (confirmed: the wrong index
    // pointed at an aromatic carbon in the unrelated trailing ring, not
    // anywhere near the stereocenter). Fixed by keying resolution on a
    // per-occurrence slot id (`next_ring_slot`) that is never reused,
    // regardless of how many times the same ring digit is.

    #[test]
    fn ring_digit_reuse_inside_stereocenter_branch_real_world_repro() {
        // Real molecule found via corpus sweep. `variant` is an RDKit
        // doRandom=True re-spelling of the exact same molecule as `orig`,
        // where the stereocenter's ring-1 partner closes inside its own
        // branch AND ring digit 1 is reused later for a trailing phenyl.
        let orig = r"COc1ccc2c3c1OC1[C@H](O)[C@](CO)(CCCCCc4ccccc4)CC4C(C2)N(C)CCC341";
        let variant = r"C([C@@]1(CC2C34CCN(C2Cc2ccc(c(c24)OC3[C@@H]1O)OC)C)CO)CCCCc1ccccc1";
        assert!(
            same_canonical(orig, variant),
            "ring-digit reuse must not corrupt a stereocenter whose own ring \
             partner closes inside its branch"
        );
    }

    #[test]
    fn ring_digit_reuse_inside_stereocenter_branch_minimal() {
        // Minimal case matching the real repro's precondition exactly:
        // `[C@H]1` opens ring 1 and its partner closes INSIDE its own first
        // branch `(CC1)` -- i.e. before the parser ever advances to a new
        // *chain* atom for atom0, so atom0's stereo record is still
        // unfinalized at the moment of that closure (the immediate-resolution
        // fast path in `close_or_open_ring` cannot catch it; only the
        // end-of-parse fallback does). Ring digit 1 is then reused by an
        // unrelated, disconnected fragment. Before the fix, the fallback
        // resolved by raw digit and the benzene's ring closure silently
        // stole atom0's still-pending resolution.
        let smi = r"[C@H]1(CC1)Cl.c1ccccc1";
        let mol = parse(smi).unwrap_or_else(|e| panic!("{smi}: {e:?}"));
        let out = canonical_smiles(&mol);
        let mol2 = parse(&out).unwrap_or_else(|e| panic!("re-parse {out}: {e:?}"));
        assert_eq!(
            canonical_smiles(&mol2),
            out,
            "stereocenter with in-branch ring closure + later digit reuse must be stable"
        );
    }

    // ── Allene cumulated double bond stereo ──────────────────────────────────

    #[test]
    fn allene_stereo_two_enantiomers_differ() {
        // F[C@@H]=[C]=[C@H]Cl and F[C@H]=[C]=[C@@H]Cl must produce different canonical SMILES.
        let mol_r = parse("F[C@@H]=[C]=[C@H]Cl").unwrap();
        let mol_s = parse("F[C@H]=[C]=[C@@H]Cl").unwrap();
        let smi_r = canonical_smiles(&mol_r);
        let smi_s = canonical_smiles(&mol_s);
        assert_ne!(
            smi_r, smi_s,
            "allene enantiomers must produce different canonical SMILES: {smi_r}"
        );
    }

    #[test]
    fn allene_stereo_round_trip_stable() {
        for smi in &["F[C@@H]=[C]=[C@H]Cl", "F[C@H]=[C]=[C@@H]Cl"] {
            let mol = parse(smi).unwrap();
            let out = canonical_smiles(&mol);
            let mol2 = parse(&out).unwrap();
            let out2 = canonical_smiles(&mol2);
            assert_eq!(
                out, out2,
                "allene stereo must be stable: {smi} -> {out} -> {out2}"
            );
        }
    }

    #[test]
    fn issue_9678_central_allene_tag_round_trips_and_distinguishes_enantiomers() {
        let first = canonical_smiles(&parse("CC=[C@AL1]=CCC").unwrap());
        let second = canonical_smiles(&parse("CC=[C@AL2]=CCC").unwrap());
        assert!(first.contains("@AL"), "axial tag was dropped: {first}");
        assert!(second.contains("@AL"), "axial tag was dropped: {second}");
        assert_ne!(first, second, "opposite allene stereoisomers collapsed");
        for canonical in [first, second] {
            let reparsed = parse(&canonical).unwrap();
            assert_eq!(canonical_smiles(&reparsed), canonical);
        }
    }

    #[test]
    fn issue_9678_shorthand_and_atom_renumbering_preserve_allene_stereo() {
        let explicit = parse("CC=[C@AL1]=CCC").unwrap();
        let shorthand = parse("CC=[C@]=CCC").unwrap();
        let expected = canonical_smiles(&explicit);
        assert_eq!(canonical_smiles(&shorthand), expected);

        let n = explicit.atom_count();
        for permutation in [
            (0..n).rev().collect::<Vec<_>>(),
            vec![2, 0, 5, 3, 1, 4],
        ] {
            let mut relabeled = relabel_molecule_preserving_ez(&explicit, &permutation);
            let mut old_to_new = vec![0u32; n];
            for (new_idx, &old_idx) in permutation.iter().enumerate() {
                old_to_new[old_idx] = new_idx as u32;
            }
            for old_idx in 0..n {
                let old = AtomIdx(old_idx as u32);
                if let Some(order) = explicit.stereo_neighbor_order(old) {
                    let mapped = order
                        .iter()
                        .map(|&id| {
                            if id == STEREO_H_SENTINEL {
                                id
                            } else {
                                old_to_new[id as usize]
                            }
                        })
                        .collect();
                    relabeled.set_stereo_neighbor_order(
                        AtomIdx(old_to_new[old_idx]),
                        mapped,
                    );
                }
            }
            assert_eq!(canonical_smiles(&relabeled), expected);
        }
    }

    /// Pinned golden-value regression for the allene-end-carbon bug this PR
    /// found and fixed: an allene end carbon (sp2, one real double-bond
    /// partner standing in for the 4th tetrahedral-like position) has a
    /// **3-element** `stereo_neighbor_order`, not 4 -- routing it through
    /// `chematic_core::remap_tetrahedral_parity` (fixed at `[u32; 4]`) would
    /// silently fall through to the "unchanged" fallback for the wrong
    /// reason (array-conversion failure, not "no verifiable order"), which
    /// produces a DIFFERENT valid-looking tag than the length-generic
    /// `permutation_is_odd` correctly computes -- this was caught only by a
    /// byte-identical before/after diff during development
    /// (`allene_stereo_two_enantiomers_differ`/`allene_stereo_round_trip_stable`
    /// above kept passing throughout, since neither checks an exact golden
    /// value, only relative invariants -- see the RFC's §8/§12.1 for the
    /// full incident writeup). `corrected_chirality`'s fallback dispatch
    /// (`canonical.rs`, the `original_arr`/`canonical_arr` `<[u32; 4]>::try_from`
    /// pair) explicitly falls back to `permutation_is_odd` whenever either
    /// side isn't exactly 4 elements -- exercised here, not just described.
    /// Pinned exact values so a future regression fails a normal
    /// `cargo test`, not only a throwaway diff tool.
    #[test]
    fn allene_stereo_exact_canonical_value_is_pinned() {
        assert_eq!(
            canonical_smiles(&parse("F[C@@H]=[C]=[C@H]Cl").unwrap()),
            "C(=[C@H]Cl)=[C@H]F"
        );
        assert_eq!(
            canonical_smiles(&parse("F[C@H]=[C]=[C@@H]Cl").unwrap()),
            "C(=[C@@H]Cl)=[C@@H]F"
        );
    }

    /// Mirrors `square_planar_stereo.rs`'s own
    /// `untagged_four_coordinate_pt_is_never_auto_promoted` for the
    /// tetrahedral case: a plain, untagged 4-distinct-substituent carbon
    /// must stay `Chirality::None` through parsing and canonicalization --
    /// never auto-assigned a `@`/`@@` tag just because it happens to have 4
    /// distinct substituents that *could* form a stereocenter. Guaranteed
    /// structurally by `corrected_chirality`'s very first line (`if stored
    /// == Chirality::None { return Chirality::None; }`, unchanged by this
    /// PR): the new `stereo_geometry` module's `canonicalize_configuration`/
    /// `remap_*` functions are never even called for an atom with no
    /// declared chirality, so there is no code path in this PR that could
    /// invent one.
    #[test]
    fn untagged_tetrahedral_center_is_never_auto_promoted() {
        let mol = parse("C(F)(Cl)(Br)I").unwrap();
        for (idx, atom) in mol.atoms() {
            assert_eq!(
                atom.chirality,
                Chirality::None,
                "atom {idx:?} unexpectedly carries stereo in an untagged fixture"
            );
        }
        let smi = canonical_smiles(&mol);
        assert!(
            !smi.contains('@'),
            "untagged input must never gain a stereo tag on write: {smi}"
        );
    }

    // ── RDKit PR #8957: fused-ring stereo round-trip ────────────────────────

    #[test]
    fn ring_stereo_stable_in_fused_system() {
        // RDKit PR #8957: "modern stereo" perception inverted R/S in fused polycyclic
        // systems with multiple stereocenters.  The canonical SMILES form may use @
        // instead of @@ (both are valid encodings of the same stereoisomer depending
        // on traversal order), so the invariant is round-trip stability, not literal
        // @@ count.
        let smi = r"CC[C@@]1(C)C[C@@](CC)(c2ccccc2)CCO1";
        let mol = parse(smi).expect("fused ring stereo mol");
        let out = canonical_smiles(&mol);
        // Round-trip must be stable: canonical(canonical(x)) == canonical(x).
        let mol2 = parse(&out).expect("canonical re-parse");
        let out2 = canonical_smiles(&mol2);
        assert_eq!(
            out, out2,
            "fused ring stereo must be stable after canonical round-trip"
        );
        // The canonical SMILES must still contain at least 2 stereocenters (@/@@ count ≥ 2).
        let stereo_count = out.matches('@').count();
        assert!(
            stereo_count >= 2,
            "both stereocenters must be encoded (got {stereo_count}): {out}"
        );
    }

    // ── E/Z directional-bond canonical stability (issue: Sprint 8) ──────────
    //
    // The canonical writer emits `/`,`\` directional bonds with traversal-direction
    // correction but no separate "normalization" pass. These tests lock in that the
    // direction choice is already deterministic and idempotent for stable skeletons,
    // so a future writer change cannot silently regress E/Z output. (The residual
    // canonical_diff idempotency failures are large fused-polycyclic atom-ranking
    // non-convergence, not a `/`,`\` direction bug — see docs/rdkit_compat.md.)

    /// E/Z parity of the first stereo double bond: `Some(true)` = E (the two
    /// *reference* substituents are on opposite sides), `Some(false)` = Z,
    /// `None` = no specified geometry.
    ///
    /// At each end, the reference substituent is whichever of its (at most
    /// two) non-double-bond neighbors has the lower `morgan_ranks` value —
    /// a deterministic, structure-derived choice (not CIP priority, but
    /// consistent across any re-parse of the same molecule, which is all
    /// this self-consistency check needs) — falling back to the sibling's
    /// marker via the trigonal-carbon complement identity when the
    /// reference substituent's own bond isn't the one marked.
    ///
    /// A prior version of this helper picked "whichever substituent bond
    /// happens to be marked, in adjacency order" instead of a fixed
    /// structural reference. That is parse-order-dependent: at a
    /// trisubstituted end, adjacency order is just first-encountered-in-the-
    /// input order, not tied to which substituent is higher priority, so
    /// comparing (say) the low-priority substituent at one end against the
    /// high-priority one at the other can silently invert the apparent
    /// parity. It happened to still pass before `resolve_ez_markers`
    /// existed only because canonicalization never relocated a marker to a
    /// *different* substituent, so both sides of the before/after
    /// comparison picked the same (arbitrary) atom by construction. Once a
    /// marker can legitimately move to a different, equally-valid carrier,
    /// that construction no longer holds, and the fixed rank-based
    /// reference is what actually verifies "the encoded geometry didn't
    /// change" rather than "the same accidental artifact survived".
    fn double_bond_is_e(smiles: &str) -> Option<bool> {
        let mol = parse(smiles).unwrap();
        double_bond_is_e_mol(&mol)
    }

    fn raw_dir(mol: &Molecule, bidx: BondIdx) -> Option<BondOrder> {
        let order = mol.bond(bidx).order;
        if matches!(order, BondOrder::Up | BondOrder::Down) {
            return Some(order);
        }
        mol.bond_direction(bidx)
    }

    /// Whether `end`'s substituent side of a double bond points "up",
    /// reading whichever substituent actually carries a direction marker
    /// via the fixed, rank-based reference (lower-`ranks` substituent) with
    /// sibling-complement fallback -- the same reference `resolve_ez_marker_
    /// for_end` uses to pick a canonical carrier. Shared by
    /// `double_bond_is_e_mol` (single representative bond) and
    /// `geometry_fingerprint` (every stereogenic bond) so there is exactly
    /// one implementation of "which substituent to trust" for tests, not
    /// two that could silently disagree.
    fn up_of_reference(mol: &Molecule, ranks: &[u64], end: AtomIdx) -> Option<bool> {
        let subs: Vec<(AtomIdx, BondIdx)> = mol
            .neighbors(end)
            .filter(|&(_, b)| mol.bond(b).order != BondOrder::Double)
            .collect();
        match subs.len() {
            1 => {
                let (_, bidx) = subs[0];
                let dir = raw_dir(mol, bidx)?;
                Some(CanonicalWriter::direction_is_up(
                    dir,
                    mol.bond(bidx).atom1,
                    end,
                ))
            }
            2 => {
                let reference = *subs
                    .iter()
                    .min_by_key(|&&(a, _)| ranks[a.0 as usize])
                    .expect("subs has 2 elements");
                let sibling = if reference.0 == subs[0].0 {
                    subs[1]
                } else {
                    subs[0]
                };
                if let Some(dir) = raw_dir(mol, reference.1) {
                    Some(CanonicalWriter::direction_is_up(
                        dir,
                        mol.bond(reference.1).atom1,
                        end,
                    ))
                } else {
                    let dir = raw_dir(mol, sibling.1)?;
                    Some(!CanonicalWriter::direction_is_up(
                        dir,
                        mol.bond(sibling.1).atom1,
                        end,
                    ))
                }
            }
            _ => None,
        }
    }

    /// Like `double_bond_is_e`, but tries every double bond in the molecule
    /// (not just the first one found) and returns the first that yields a
    /// defined answer -- a molecule can have an earlier, non-stereogenic
    /// double bond (e.g. a ketone's `C=O`, whose oxygen side has no
    /// substituent at all) ahead of its actual stereo alkene in bond order.
    fn double_bond_is_e_mol(mol: &Molecule) -> Option<bool> {
        let ranks = morgan_ranks(mol);
        mol.bonds()
            .filter(|(_, b)| b.order == BondOrder::Double)
            .find_map(|(_, b)| {
                let ua = up_of_reference(mol, &ranks, b.atom1)?;
                let ub = up_of_reference(mol, &ranks, b.atom2)?;
                Some(ua != ub)
            })
    }

    /// A geometry fingerprint suitable for comparing two *different* parses
    /// of the same molecule (e.g. original input vs. a canonical-output
    /// reparse) -- one E/Z fact per stereogenic double bond (both ends have
    /// at least one substituent), via `up_of_reference`. Ordered by each
    /// bond's lower-ranked endpoint: `ranks` is a molecule-intrinsic, permutation-
    /// invariant key (unlike `BondIdx`/bond-parse-order, which differs
    /// between two independently-atom-ordered spellings of the same
    /// molecule), so this ordering lines up positionally across the two
    /// parses being compared as long as neither parse gained or lost a
    /// stereogenic double bond.
    fn geometry_fingerprint(mol: &Molecule) -> Vec<Option<bool>> {
        let ranks = morgan_ranks(mol);
        let mut doubles: Vec<(u64, BondIdx)> = mol
            .bonds()
            .filter(|(_, b)| b.order == BondOrder::Double)
            .filter(|(_, b)| {
                CanonicalWriter::end_has_substituent(mol, b.atom1)
                    && CanonicalWriter::end_has_substituent(mol, b.atom2)
            })
            .map(|(bidx, b)| {
                let key = ranks[b.atom1.0 as usize].min(ranks[b.atom2.0 as usize]);
                (key, bidx)
            })
            .collect();
        doubles.sort_by_key(|&(k, _)| k);

        doubles
            .into_iter()
            .map(|(_, bidx)| {
                let bond = mol.bond(bidx);
                let ua = up_of_reference(mol, &ranks, bond.atom1)?;
                let ub = up_of_reference(mol, &ranks, bond.atom2)?;
                Some(ua != ub)
            })
            .collect()
    }

    /// A spelling-independent view of the production E/Z geometry extractor.
    /// The raw atom/bond indices deliberately do not escape this helper;
    /// rank keys and each double bond's same-side relation are the facts a
    /// later canonical plan is allowed to compare.
    type EzGeometrySignature = Vec<((u64, u64), (u64, u64), bool)>;

    fn ez_geometry_signature(mol: &Molecule) -> EzGeometrySignature {
        let (ranks, _) = winning_individualized_ranks(mol);
        let writer = CanonicalWriter::new(mol, &ranks);
        writer
            .extract_ez_geometry_facts()
            .into_iter()
            .map(|fact| {
                let mut references = [fact.ends[0].reference_rank, fact.ends[1].reference_rank];
                references.sort_unstable();
                (
                    fact.double_key,
                    (references[0], references[1]),
                    fact.same_side,
                )
            })
            .collect()
    }

    /// A canonical output occurrence for one candidate E/Z carrier bond.
    /// `kind` is tree edge, ring opening, or ring closing respectively.  The
    /// signature contains only canonical ranks, never parse-order indices.
    type EzOutputSlotSignature = Vec<((u64, u64), u64, u8)>;

    fn ez_output_slot_signature(mol: &Molecule) -> EzOutputSlotSignature {
        let (ranks, _) = winning_individualized_ranks(mol);
        let ends = CanonicalWriter::compute_stereo_alkene_ends(mol);
        let candidates: HashSet<BondIdx> = CanonicalWriter::coupling_components(mol, &ends)
            .into_iter()
            .flatten()
            .flat_map(|end| CanonicalWriter::substituents(mol, end))
            .map(|(_, bond)| bond)
            .collect();
        let mut writer = CanonicalWriter::new(mol, &ranks);
        let starts = writer.canonical_atom_list();
        writer.find_ring_closures(&starts);

        let mut slots = Vec::new();
        for bond in candidates {
            let bond_ref = mol.bond(bond);
            let mut edge = (
                ranks[bond_ref.atom1.0 as usize],
                ranks[bond_ref.atom2.0 as usize],
            );
            if edge.1 < edge.0 {
                edge = (edge.1, edge.0);
            }
            if writer.ring_bonds[bond.0 as usize] {
                for (atom, rings) in writer.atom_ring_nums.iter().enumerate() {
                    for &(_, is_open, _, ring_bond) in rings {
                        if ring_bond == bond {
                            slots.push((edge, ranks[atom], if is_open { 1 } else { 2 }));
                        }
                    }
                }
            } else if let Some((parent, _)) = writer
                .canonical_tree_parent
                .iter()
                .flatten()
                .find(|(_, tree_bond)| *tree_bond == bond)
            {
                slots.push((edge, ranks[parent.0 as usize], 0));
            }
        }
        slots.sort_unstable();
        slots
    }

    /// Full rank-keyed canonical DFS shape, excluding only directional token
    /// polarity.  Unlike `ez_output_slot_signature`, this records every tree
    /// parent and every ring-token occurrence so #503 diagnostics can tell a
    /// skeleton divergence from a stereochemical-token divergence without
    /// inspecting parse-order atom or bond indices.
    type CanonicalSkeletonSignature = (Vec<(u64, u64)>, Vec<(u64, u64, u32, bool)>);

    fn canonical_skeleton_signature(mol: &Molecule) -> CanonicalSkeletonSignature {
        let (ranks, _) = winning_individualized_ranks(mol);
        let mut writer = CanonicalWriter::new(mol, &ranks);
        let starts = writer.canonical_atom_list();
        writer.find_ring_closures(&starts);

        let mut tree = writer
            .canonical_tree_parent
            .iter()
            .enumerate()
            .filter_map(|(child, parent)| {
                parent.map(|(parent, _)| (ranks[parent.0 as usize], ranks[child]))
            })
            .collect::<Vec<_>>();
        tree.sort_unstable();

        let mut rings = writer
            .atom_ring_nums
            .iter()
            .enumerate()
            .flat_map(|(atom, entries)| {
                let ranks = &ranks;
                entries.iter().map(move |&(number, is_open, partner, _)| {
                    (ranks[atom], ranks[partner.0 as usize], number, is_open)
                })
            })
            .collect::<Vec<_>>();
        rings.sort_unstable();
        (tree, rings)
    }

    const EZ_STABLE_CORPUS: &[&str] = &[
        "C/C=C/C",     // (E)-2-butene
        "C/C=C\\C",    // (Z)-2-butene
        "F/C=C/F",     // (E)-1,2-difluoroethene
        "F/C=C\\F",    // (Z)
        "CC/C=C/CC",   // (E)-3-hexene
        "CC/C=C\\CC",  // (Z)-3-hexene
        "C/C=C/C=C/C", // (2E,4E)-hexadiene
        "Cl/C=C/Br",
        "C/C=C/c1ccccc1", // (E)-propenylbenzene
        "C/C(F)=C(\\F)C",
    ];

    #[test]
    fn stashed_directions_read_from_atom2_keep_their_geometry() {
        // The same markers moved into stashes read from the bond's atom2
        // (Up from atom2 is Down from atom1): carrier resolution wrote
        // atom1-relative markers that emission then read against the
        // anchor, mirroring the E/Z (#734, 92 of 272 E/Z molecules).
        for smiles in [
            r"COc1cc2nc(N3CCN(C(=O)/C(F)=C/c4ccccc4)CC3)nc(N)c2cc1OC",
            r"COc1cc2nc(N3CCN(/C(S)=N/C4CCCCC4)CC3)nc(N)c2cc1OC",
            r"O=C1OC(=C/I)/CCC1c1cccc2ccccc12",
            r"C/C=C/C(C)=C/C=C\C",
        ] {
            let literal = crate::parse(smiles).unwrap();
            let mut stashed = literal.clone();
            let marked: Vec<(BondIdx, BondOrder, AtomIdx)> = literal
                .bonds()
                .filter(|(_, b)| matches!(b.order, BondOrder::Up | BondOrder::Down))
                .map(|(i, b)| (i, b.order, b.atom2))
                .collect();
            for (bidx, order, atom2) in marked {
                let from_atom2 = if order == BondOrder::Up {
                    BondOrder::Down
                } else {
                    BondOrder::Up
                };
                stashed.set_bond_order(bidx, BondOrder::Single);
                stashed.set_bond_direction(bidx, from_atom2);
                stashed.set_bond_direction_anchor(bidx, atom2);
            }
            assert_eq!(
                canonical_smiles(&stashed),
                canonical_smiles(&literal),
                "{smiles}"
            );
        }
    }

    #[test]
    fn ez_canonical_smiles_is_idempotent() {
        for s in EZ_STABLE_CORPUS {
            assert!(
                is_stable(s),
                "E/Z canonical SMILES must be idempotent for {s}"
            );
        }
    }

    #[test]
    fn ez_geometry_preserved_through_canonicalization() {
        for s in EZ_STABLE_CORPUS {
            let want = double_bond_is_e(s)
                .unwrap_or_else(|| panic!("input {s} must have specified geometry"));
            let canon = canonical_smiles(&parse(s).unwrap());
            let got = double_bond_is_e(&canon)
                .unwrap_or_else(|| panic!("canonical {canon} dropped geometry from {s}"));
            assert_eq!(got, want, "E/Z geometry changed: {s} -> {canon}");
        }
    }

    #[test]
    fn ez_e_and_z_differ_for_each_skeleton() {
        // Each E form must canonicalize differently from its Z form.
        for (e, z) in [
            ("C/C=C/C", "C/C=C\\C"),
            ("F/C=C/F", "F/C=C\\F"),
            ("CC/C=C/CC", "CC/C=C\\CC"),
        ] {
            assert_ne!(
                canonical_smiles(&parse(e).unwrap()),
                canonical_smiles(&parse(z).unwrap()),
                "E and Z must produce different canonical SMILES ({e} vs {z})"
            );
        }
    }

    // ── Fused-aromatic canonical idempotency (Sprint 9) ─────────────────────
    //
    // Lock in the fused aromatics that DO round-trip consistently. The residual
    // ~1.6% canonical idempotency failures on large fused polycyclics are caused
    // by aromaticity-perception round-trip inconsistency (a molecule vs the
    // re-parse of its own canonical SMILES can disagree on which bonds are
    // aromatic — e.g. 16 vs 17 on a fluorene-type linkage — which shifts Morgan
    // ranks). That is an aromaticity/parser-core issue, not a canonical-ranking
    // bug; see docs/rdkit_compat.md. These cases are stable and guarded here.

    #[test]
    fn fused_aromatic_canonical_is_idempotent() {
        for s in [
            "c1ccc2ccccc2c1",         // naphthalene
            "c1ccc2ncccc2c1",         // quinoline
            "c1ccc2c(c1)cc[nH]2",     // indole
            "c1ccc2cc3ccccc3cc2c1",   // anthracene
            "c1ccc2[nH]c3ccccc3c2c1", // carbazole
            "c1ccc2c(c1)oc1ccccc12",  // dibenzofuran
        ] {
            assert!(
                is_stable(s),
                "fused-aromatic canonical SMILES must be idempotent for {s}"
            );
        }
    }

    // ── Round 12: simple E/Z direction normalization ────────────────────────
    //
    // `/N=N/` and `\N=N\` are two equally valid SMILES spellings of the same
    // geometry (flipping every directional bond of one connected E/Z system
    // preserves meaning). Before this fix, the writer just propagated
    // whichever direction the parser happened to read, so two spellings of
    // the same molecule could canonicalize to two different strings. The fix
    // normalizes each connected E/Z system so its first directional bond (in
    // canonical write order) is always `/`.

    #[test]
    fn ez_simple_bond_direction_normalized_azo() {
        assert!(
            same_canonical("CN(C)/N=N/c1ccccc1", r"CN(C)\N=N\c1ccccc1",),
            "isolated E/Z double bond must canonicalize identically regardless \
             of which of the two equally-valid slash spellings was parsed"
        );
    }

    #[test]
    fn ez_simple_bond_direction_normalized_symmetric() {
        assert!(same_canonical("F/C=C/F", r"F\C=C\F"));
        assert!(same_canonical("C(/F)=C/F", r"C(\F)=C\F"));
    }

    // Bracket atoms with `hydrogen_count: None` must still get their implicit
    // hydrogens written by the canonical writer too — same bug/fixtures as
    // writer::tests::test_bracket_implicit_h_*.
    use chematic_core::{Atom, Element};

    #[test]
    fn canonical_bracket_implicit_h_ammonium_charge_only() {
        let mut b = chematic_core::MoleculeBuilder::new();
        let mut n = Atom::new(Element::N);
        n.charge = 1;
        b.add_atom(n);
        assert_eq!(canonical_smiles(&b.build()), "[NH4+]");
    }

    #[test]
    fn canonical_bracket_implicit_h_isotope_only() {
        let mut b = chematic_core::MoleculeBuilder::new();
        let mut c = Atom::new(Element::C);
        c.isotope = Some(13);
        b.add_atom(c);
        assert_eq!(canonical_smiles(&b.build()), "[13CH4]");
    }

    #[test]
    fn canonical_bracket_implicit_h_atom_map_only() {
        let mut b = chematic_core::MoleculeBuilder::new();
        let mut c0 = Atom::new(Element::C);
        c0.atom_map = Some(7);
        let c0 = b.add_atom(c0);
        let c1 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c0, c1, BondOrder::Single).unwrap();
        assert_eq!(canonical_smiles(&b.build()), "C[CH3:7]");
    }

    #[test]
    fn canonical_bracket_implicit_h_isotope_and_atom_map() {
        let mut b = chematic_core::MoleculeBuilder::new();
        let mut c0 = Atom::new(Element::C);
        c0.isotope = Some(13);
        c0.atom_map = Some(7);
        let c0 = b.add_atom(c0);
        let c1 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c0, c1, BondOrder::Single).unwrap();
        assert_eq!(canonical_smiles(&b.build()), "C[13CH3:7]");
    }

    #[test]
    fn canonical_bracket_implicit_h_hydroxide_charge_only() {
        let mut b = chematic_core::MoleculeBuilder::new();
        let mut o = Atom::new(Element::O);
        o.charge = -1;
        b.add_atom(o);
        assert_eq!(canonical_smiles(&b.build()), "[OH-]");
    }

    // ── Standalone wedge/hash bond must not emit a meaningless SMILES
    // directional token in canonical_smiles either (same fix as
    // `writer::tests`, applied to the canonical writer's own emission
    // sites) ──────────────────────────────────────────────────────────────

    /// Minimal repro from docs/rfcs/stereo2d_reader_integration_rfc.md §3: a
    /// wedge bond with no adjacent double bond must not become `/`/`\`.
    #[test]
    fn canonical_standalone_solid_wedge_not_written_as_slash() {
        let mut b = chematic_core::MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        b.add_bond(c, f, BondOrder::Single).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        let wedge_bond = b.add_bond(c, br, BondOrder::Up).unwrap();
        let mol = b.build();

        let out = canonical_smiles(&mol);
        assert!(
            !out.contains('/') && !out.contains('\\'),
            "standalone solid wedge (no adjacent double bond) must not be \
             written as a directional token: got '{out}'"
        );
        assert_eq!(mol.bond(wedge_bond).order, BondOrder::Up);

        let mol2 = crate::parser::parse(&out).unwrap_or_else(|e| panic!("re-parse '{out}': {e}"));
        assert_eq!(mol2.atom_count(), mol.atom_count());
        assert_eq!(mol2.bond_count(), mol.bond_count());
    }

    /// Same shape but with a hash wedge (`BondOrder::Down`).
    #[test]
    fn canonical_standalone_hash_wedge_not_written_as_backslash() {
        let mut b = chematic_core::MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        b.add_bond(c, f, BondOrder::Single).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        let wedge_bond = b.add_bond(c, br, BondOrder::Down).unwrap();
        let mol = b.build();

        let out = canonical_smiles(&mol);
        assert!(
            !out.contains('/') && !out.contains('\\'),
            "standalone hash wedge (no adjacent double bond) must not be \
             written as a directional token: got '{out}'"
        );
        assert_eq!(mol.bond(wedge_bond).order, BondOrder::Down);

        let mol2 = crate::parser::parse(&out).unwrap_or_else(|e| panic!("re-parse '{out}': {e}"));
        assert_eq!(mol2.atom_count(), mol.atom_count());
        assert_eq!(mol2.bond_count(), mol.bond_count());
    }

    /// A wedge bond landing on a ring-closure edge must be suppressed
    /// identically in the canonical writer's ring-closure emission site.
    #[test]
    fn canonical_standalone_wedge_on_ring_closure_bond() {
        let mut b = chematic_core::MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let c3 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        b.add_bond(c2, c3, BondOrder::Single).unwrap();
        b.add_bond(c1, c3, BondOrder::Up).unwrap();
        let mol = b.build();

        let out = canonical_smiles(&mol);
        assert!(
            !out.contains('/') && !out.contains('\\'),
            "standalone wedge on a ring-closure bond must not be written as \
             a directional token: got '{out}'"
        );
        let mol2 = crate::parser::parse(&out).unwrap_or_else(|e| panic!("re-parse '{out}': {e}"));
        assert_eq!(mol2.atom_count(), 3);
        assert_eq!(mol2.bond_count(), 3);
    }

    /// The fix must not disturb `Atom.chirality`/`Molecule::stereo_neighbor_order`
    /// while suppressing a standalone wedge bond's spurious token on the
    /// same stereocenter -- mirrors `writer::tests::
    /// test_standalone_wedge_does_not_disturb_stereocenter_chirality` for
    /// the canonical writer's own (rank-dependent) chirality-correction path.
    #[test]
    fn canonical_standalone_wedge_does_not_disturb_stereocenter_chirality() {
        let mut center = Atom::new(Element::C);
        center.chirality = Chirality::Clockwise;
        // Force bracket notation (h=0 still forces `[..]` without printing
        // an `H` token) so the chirality symbol is actually emitted.
        center.hydrogen_count = Some(0);
        let mut b = chematic_core::MoleculeBuilder::new();
        let c = b.add_atom(center);
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let iodine = b.add_atom(Atom::new(Element::I));
        b.add_bond(c, f, BondOrder::Single).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        b.add_bond(c, br, BondOrder::Single).unwrap();
        let wedge_bond = b.add_bond(c, iodine, BondOrder::Up).unwrap();
        b.set_stereo_neighbor_order(c, vec![f.0, cl.0, br.0, iodine.0]);
        let mol = b.build();

        let out = canonical_smiles(&mol);
        assert!(
            !out.contains('/') && !out.contains('\\'),
            "standalone wedge on a stereocenter's substituent must still be \
             suppressed: got '{out}'"
        );
        assert!(
            out.contains('@'),
            "the stereocenter's own chirality symbol must still be emitted \
             alongside the (now-suppressed) wedge bond: got '{out}'"
        );

        assert_eq!(mol.atom(c).chirality, Chirality::Clockwise);
        assert_eq!(mol.bond(wedge_bond).order, BondOrder::Up);
        assert_eq!(
            mol.stereo_neighbor_order(c),
            Some([f.0, cl.0, br.0, iodine.0].as_slice())
        );
    }

    /// Regression guard: genuine E/Z directional markers must still be
    /// emitted by canonical_smiles (already covered by `ez_e_stable` etc.
    /// above; this pins the specific double-bond-adjacency mechanism the
    /// standalone-wedge fix now gates on).
    #[test]
    fn canonical_genuine_ez_directional_bond_still_written() {
        let mol = parse("F/C=C/F").unwrap();
        let out = canonical_smiles(&mol);
        assert!(
            out.contains('/') || out.contains('\\'),
            "genuine double-bond-adjacent directional marker must survive \
             canonical_smiles(): got '{out}'"
        );
    }

    /// The pre-existing "E/Z direction stashed on an aromatic bond" fixtures
    /// (`parser::tests::direction_stash_*`) all involve a real exocyclic
    /// double bond and must keep round-tripping unaffected by this fix --
    /// exercised end-to-end here via `canonical_smiles` directly (those
    /// tests live in parser.rs and already re-run on every `cargo test`).
    #[test]
    fn canonical_aromatic_stash_with_real_double_bond_still_emits_direction() {
        let mol = parse(r"C/N=c1\c(O)c(O)\c1=N/C").unwrap();
        let out = canonical_smiles(&mol);
        assert!(
            out.contains('/') || out.contains('\\'),
            "genuine exocyclic-double-bond-adjacent aromatic stash must \
             still emit a directional token: got '{out}'"
        );
    }

    #[test]
    fn canonical_rejects_same_side_ring_carriers_on_fused_imine() {
        let input = r"CCOC(=O)c1cc2ccccn2/c(=N/c2ccc(Br)cc2)n1";
        let old_conflicting = r"N(/c3ccc(Br)cc3)=c/1/n2ccccc2cc(n1)C(=O)OCC";
        assert!(
            !ez_direction_assignments_consistent(&parse(old_conflicting).unwrap()),
            "the historical spelling has two same-side carriers at one alkene end"
        );

        let output = canonical_smiles(&parse(input).unwrap());
        assert!(
            ez_direction_assignments_consistent(&parse(&output).unwrap()),
            "canonical output must not contain a conflicting E/Z specification: {output}"
        );
    }

    #[test]
    fn canonical_rejects_same_side_carriers_on_coupled_bis_imine() {
        let input = r"CC/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O";
        let old_conflicting =
            r"c/3(c(/c(/c3=N\CC)=N\[C@@H](Cc1ccc(NC(=O)c2c(cncc2Cl)Cl)cc1)C(O)=O)O)O";
        assert!(
            !ez_direction_assignments_consistent(&parse(old_conflicting).unwrap()),
            "the historical spelling has two same-side carriers at one alkene end"
        );

        let output = canonical_smiles(&parse(input).unwrap());
        assert!(
            ez_direction_assignments_consistent(&parse(&output).unwrap()),
            "canonical output must not contain a conflicting E/Z specification: {output}"
        );
    }

    #[test]
    fn cyclic_imine_carrier_normalization_is_input_order_invariant() {
        let fixtures = [
            (
                r"C[C@]12CC/C(=N\N=C(N)N)C=C1CC[C@@H]2CO",
                r"C1(=N/N=C(N)N)\C=C2CC[C@H](CO)[C@@]2(C)CC1",
            ),
            (
                r"C[C@]12CC/C(=N\N=C(N)N)C=C1CC[C@@H]2/C=C/[N+](=O)[O-]",
                r"C[C@@]12[C@@H](/C=C/[N+](=O)[O-])CCC1=C/C(=N/N=C(N)N)CC2",
            ),
        ];
        for (a, b) in fixtures {
            assert_eq!(
                canonical_smiles(&parse(a).unwrap()),
                canonical_smiles(&parse(b).unwrap()),
                "equivalent cyclic-imine carrier spellings must converge"
            );
        }
    }

    // ── E/Z marker-carrier normalization (fix/canonical-ez-carrier-
    // normalization) ─────────────────────────────────────────────────────
    //
    // At a tri-/tetra-substituted stereo alkene end, SMILES only requires
    // ONE of its ≥2 substituent bonds to carry the `/`/`\` marker -- the
    // other's position is implied (a trigonal alkene carbon has exactly two
    // sides). *Which* substituent gets the mark used to be whatever the
    // parser happened to read, so two RDKit-valid respellings of the same
    // molecule that mark different substituents produced two different
    // canonical outputs; the current evidence is recorded in
    // docs/canonical-ez-plan.md. `resolve_ez_markers` picks the marker carrier deterministically
    // from canonical rank instead. Every case below is a real molecule from
    // the residual corpus (`validation/results/
    // canonical_residual_diagnosis_summary.json`'s `permutation_invariance_
    // failures_sample`), pinned as a regression guard, not a synthetic
    // approximation of the bug.

    /// Assert `a` and `b` -- two real, already-observed-divergent canonical
    /// outputs of the SAME molecule -- now canonicalize identically, AND
    /// that the E/Z geometry each encodes is unchanged by that
    /// canonicalization (checked via `double_bond_is_e`, not string
    /// comparison, per the "zero new semantic changes" requirement).
    #[track_caller]
    fn assert_ez_carrier_pair_resolved(a: &str, b: &str) {
        let want_a = double_bond_is_e(a);
        let want_b = double_bond_is_e(b);
        let canon_a = canonical_smiles(&parse(a).unwrap());
        let canon_b = canonical_smiles(&parse(b).unwrap());
        assert_eq!(
            canon_a, canon_b,
            "two divergent real-corpus spellings of the same molecule must \
             now canonicalize identically: '{a}' -> '{canon_a}' vs '{b}' -> '{canon_b}'"
        );
        assert_eq!(
            double_bond_is_e(&canon_a),
            want_a,
            "canonicalizing '{a}' -> '{canon_a}' changed its E/Z geometry"
        );
        assert_eq!(
            double_bond_is_e(&canon_b),
            want_b,
            "canonicalizing '{b}' -> '{canon_b}' changed its E/Z geometry"
        );
    }

    /// Simplest case: a disubstituted-vs-trisubstituted amidine carbon with
    /// exactly two non-H substituents (methylamino / 4-iodobenzylamino), no
    /// ring at all. The two spellings mark different substituents.
    #[test]
    fn ez_carrier_trisub_no_ring() {
        assert_ez_carrier_pair_resolved("Ic1ccc(cc1)CN/C(=N/C)NC", "Ic1ccc(cc1)CNC(=N/C)/NC");
    }

    /// Tetrasubstituted alkene: both ends have two real substituents (no
    /// implicit H on either alkene carbon), so both ends independently pick
    /// a canonical carrier.
    #[test]
    fn ez_carrier_tetrasub_no_ring() {
        assert_ez_carrier_pair_resolved(
            "COc1ccc(/C=C(\\C)C(=O)c2cc(OC)c(OC)c(OC)c2)cc1O",
            "COc1ccc(/C=C(C(=O)c2cc(OC)c(OC)c(OC)c2)\\C)cc1O",
        );
    }

    /// One candidate substituent is reached via a plain tree edge, the
    /// other via the SAME alkene end's ring-closure bond (an aliphatic
    /// cyclopentylidene hydrazone) -- covers "E/Z bond that is also part of
    /// a ring closure" directly, with no aromaticity involved at all.
    #[test]
    fn ez_carrier_ring_closure_candidate() {
        assert_ez_carrier_pair_resolved(
            "OC(=O)CCCCCCC/1CCCC1=N/NCCCCCC",
            "OC(=O)CCCCCCC1CCC/C1=N\\NCCCCCC",
        );
    }

    /// The aromatic-bond-direction stash (`resolve_aromatic_direction_stash`
    /// in parser.rs) combined with a ring-closure candidate: a four-membered
    /// ring carbon carries an exocyclic C=N imine, and the mark can validly
    /// land on either its plain ring tree-edge or its ring-closure bond
    /// (both aromatic-aromatic, so both route through the stash side
    /// channel, not a literal `Up`/`Down` bond order). This is the exact
    /// intersection the RFC's dominant root cause calls out: tri-substituted
    /// + aromatic stash + ring closure, all at once.
    ///
    /// This same ring also has a ketone (`c(=O)`) adjacent to the imine,
    /// sharing a ring bond with it -- a ketone's carbon has two ring-bond
    /// "substituents" of its own even though a ketone has no E/Z stereo at
    /// all (the oxygen side has none). `resolve_ez_markers` must recognize
    /// the ketone end is not stereogenic (both ends of a double bond need a
    /// substituent) and never touch that shared bond on the ketone's
    /// account, or it can move/erase the imine's own marker -- this pin
    /// exercises that guard too, not just the ring-closure/stash mechanism.
    #[test]
    fn canonical_reorient_for_write_uses_stash_anchor() {
        use chematic_core::{Atom, Element, MoleculeBuilder};

        let mut builder = MoleculeBuilder::new();
        let c1 = builder.add_atom(Atom::aromatic(Element::C));
        let c2 = builder.add_atom(Atom::aromatic(Element::C));
        let bidx = builder.add_bond(c2, c1, BondOrder::Aromatic).unwrap();
        let mut mol = builder.build();
        mol.set_bond_direction(bidx, BondOrder::Up);
        mol.set_bond_direction_anchor(bidx, c1);

        let writer = CanonicalWriter::new(&mol, &[]);
        assert_eq!(
            writer.raw_direction_anchor(bidx),
            c1,
            "the parser-side stash anchor must survive reversed atom storage"
        );
        assert_eq!(
            CanonicalWriter::reorient_for_write(
                writer.raw_direction_anchor(bidx),
                c1,
                BondOrder::Up
            ),
            BondOrder::Up
        );
        assert_eq!(
            CanonicalWriter::reorient_for_write(
                writer.raw_direction_anchor(bidx),
                c2,
                BondOrder::Up
            ),
            BondOrder::Down
        );
    }

    #[test]
    fn ez_carrier_aromatic_stash_ring_closure() {
        assert_ez_carrier_pair_resolved(
            "c4(c(cncc4Cl)Cl)C(=O)Nc3ccc(cc3)C[C@H](/N=c/2c(=O)c(c2N1CCSCC1)O)C(=O)O",
            "c4(c(cncc4Cl)Cl)C(=O)Nc3ccc(cc3)C[C@H](/N=c2\\c(=O)c(c2N1CCSCC1)O)C(=O)O",
        );
    }

    /// A deterministic (seeded, reproducible -- never true randomness)
    /// Fisher-Yates permutation of `0..n`, for exercising many atom
    /// relabelings without external randomness.
    ///
    /// Deliberately NOT built on `chematic_smiles::random_smiles`: probing
    /// it during this investigation found that `random_smiles(&mol, 2)` on
    /// `CC1=C2CC[C@H](/C=N/N=C(N)N)[C@@]2(C)CC/C1=N\N=C(N)N` (one of the
    /// [`EZ_SHARED_CARRIER_FULLY_RESOLVED`] fixtures) produces a
    /// canonical key differing at a SECOND, unrelated ring stereocenter's
    /// `@`/`@@` tag from the original -- independently confirmed via RDKit
    /// (`Chem.MolToSmiles`/`compare_molecules(..., StandardInchiString)` ->
    /// `Distinct`, a genuinely different diastereomer, not just a
    /// respelling). Whether that is a real `random_smiles` bug or a
    /// legitimate edge case was NOT investigated further here -- it is
    /// out of scope for this PR -- but using it as this gate's relabeling
    /// source would risk silently conflating that separate, unconfirmed
    /// question with the E/Z-carrier claim this test exists to check.
    fn deterministic_permutation(n: usize, seed: u64) -> Vec<usize> {
        let mut state = seed ^ 0x9E37_79B9_7F4A_7C15;
        let mut next_u32 = move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 33) as u32
        };
        let mut perm: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            let j = (next_u32() as usize) % (i + 1);
            perm.swap(i, j);
        }
        perm
    }

    /// Every alternate encoding of `mol` reachable by moving ONE raw
    /// `/`/`\` mark from its current candidate bond to the SIBLING
    /// candidate at the same stereo-alkene end, via the same sibling-
    /// complement identity `end_votes` itself relies on -- i.e. a
    /// different, but geometrically identical, valid SMILES-level
    /// spelling of the exact same real molecule.
    ///
    /// This is a DIFFERENT degree of freedom than atom-index relabeling
    /// (`relabel_molecule_preserving_ez`/`deterministic_permutation`),
    /// which always preserves which physical bond carries the mark.
    /// Measured directly: pure index relabeling alone (even 16+ deterministic
    /// permutations) never reproduces the divergence the corpus diagnosis
    /// found for several [`EZ_SHARED_CARRIER_FULLY_RESOLVED`] fixtures,
    /// because that divergence is specifically about the solver's decision
    /// changing based on WHICH candidate bond happens to carry the input's
    /// raw mark -- this helper is what actually exercises that axis.
    ///
    /// Restricted to the plain (non-aromatic-stash) case: both candidates'
    /// real chemical order must be `Single` (one currently `Up`/`Down`, the
    /// other genuinely unmarked) -- sufficient for every fixture in this
    /// file (none of the 18 issue #149 fixtures route an E/Z-relevant
    /// candidate through the aromatic-bond-direction stash), not a claim
    /// this covers every possible molecule.
    ///
    /// Self-validating: a candidate bond can be shared between TWO
    /// different ends (the exact coupling this whole module resolves), so
    /// moving its mark away for ONE end's sake, considered in isolation,
    /// can silently strip the OTHER end's only source of geometry --
    /// measured directly on one of the
    /// [`EZ_SHARED_CARRIER_FULLY_RESOLVED`] fixtures (the
    /// enol-ether double bond's own mark vanishing as a side effect of
    /// relocating the SEPARATE, merely-bond-adjacent coupled imine's own
    /// mark). Every candidate alternate is checked against
    /// `geometry_fingerprint` before being returned -- only alternates that
    /// provably encode the exact same real-world geometry as `mol` are
    /// ever handed to a caller, so this helper can never itself be the
    /// source of a "lost stereo" false positive in the tests that use it.
    fn alternate_ez_markings(mol: &Molecule) -> Vec<Molecule> {
        let baseline_geo = geometry_fingerprint(mol);
        let ends = CanonicalWriter::compute_stereo_alkene_ends(mol);
        let mut alternates = Vec::new();
        for &end in &ends {
            let subs = CanonicalWriter::substituents(mol, end);
            if subs.len() != 2 {
                continue;
            }
            for i in 0..2 {
                let chosen = subs[i];
                let other = subs[1 - i];
                let chosen_order = mol.bond(chosen.1).order;
                let other_order = mol.bond(other.1).order;
                let chosen_dir = match chosen_order {
                    BondOrder::Up | BondOrder::Down => chosen_order,
                    _ => continue,
                };
                if other_order != BondOrder::Single {
                    continue; // other already marked, or non-plain -- not this helper's case
                }

                let up =
                    CanonicalWriter::direction_is_up(chosen_dir, mol.bond(chosen.1).atom1, end);
                let new_other_dir =
                    CanonicalWriter::direction_for_up(mol.bond(other.1).atom1, end, !up);

                let alt = mol
                    .with_bond_order(chosen.1, BondOrder::Single)
                    .with_bond_order(other.1, new_other_dir);
                if geometry_fingerprint(&alt) != baseline_geo {
                    continue; // this end's own move stripped a coupled partner's geometry -- discard
                }
                alternates.push(alt);
            }
        }
        alternates
    }

    /// The full set of relabelings + alternate spellings this file's
    /// permutation-invariance gates check: the original, its reversed-
    /// atom-order relabeling, 16 deterministic Fisher-Yates relabelings
    /// (all preserving which bond carries each mark), PLUS every
    /// mark-relocated alternate from [`alternate_ez_markings`] (which
    /// instead moves a mark to its sibling candidate, atom labeling
    /// unchanged) -- the two axes together are what the corpus-measured
    /// residual divergence actually needs to be reproduced locally.
    fn ez_carrier_test_variants(mol: &Molecule) -> Vec<Molecule> {
        let n = mol.atom_count();
        let mut variants: Vec<Molecule> = vec![
            relabel_molecule_preserving_ez(mol, &(0..n).collect::<Vec<_>>()),
            relabel_molecule_preserving_ez(mol, &(0..n).rev().collect::<Vec<_>>()),
        ];
        for seed in 0..16u64 {
            variants.push(relabel_molecule_preserving_ez(
                mol,
                &deterministic_permutation(n, seed),
            ));
        }
        variants.extend(alternate_ez_markings(mol));
        variants
    }

    /// A permutation-invariant tetrahedral-stereo fact per stereocenter,
    /// ordered by the CENTER atom's own ascending canonical rank so it
    /// lines up positionally across two different atom-labelings of the
    /// same molecule -- the tetrahedral counterpart to
    /// `geometry_fingerprint` (E/Z). Reuses [`permutation_is_odd`], the
    /// same primitive `corrected_chirality` uses to compute the write-time
    /// parity correction, but against a FIXED rank-sorted reference order
    /// instead of "whichever order the writer happened to traverse in," so
    /// relabeling can never change the reported polarity. Implicit H
    /// (`STEREO_H_SENTINEL`) gets a fixed extreme rank purely as an
    /// internally-consistent tie-break -- this does not need to match real
    /// CIP H-priority, only be the same choice on both sides of a
    /// comparison. `None` for a stereocenter whose 4 neighbor ranks are
    /// not fully discriminated (would make the by-value permutation lookup
    /// ambiguous) -- rare (only possible for locally-automorphic
    /// substituents) and never silently miscompared.
    fn tetrahedral_fingerprint(mol: &Molecule) -> Vec<Option<bool>> {
        let ranks = morgan_ranks(mol);
        let mut centers: Vec<(u64, AtomIdx)> = mol
            .atoms()
            .filter(|(_, a)| a.chirality.is_tetrahedral())
            .map(|(idx, _)| (ranks[idx.0 as usize], idx))
            .collect();
        centers.sort_by_key(|&(r, _)| r);

        centers
            .into_iter()
            .map(|(_, idx)| {
                let stored = mol.atom(idx).chirality;
                let original = mol.stereo_neighbor_order(idx)?;
                let rank_of = |v: u32| -> u64 {
                    if v == STEREO_H_SENTINEL {
                        u64::MAX
                    } else {
                        ranks[v as usize]
                    }
                };
                let original_ranks: Vec<u64> = original.iter().map(|&v| rank_of(v)).collect();
                let mut sorted_ranks = original_ranks.clone();
                sorted_ranks.sort_unstable();
                let mut deduped = sorted_ranks.clone();
                deduped.dedup();
                if deduped.len() != sorted_ranks.len() {
                    return None; // neighbor ranks not fully discriminated
                }
                let is_odd = permutation_is_odd(&original_ranks, &sorted_ranks);
                Some(match (stored, is_odd) {
                    (Chirality::Clockwise, false) => true,
                    (Chirality::Clockwise, true) => false,
                    (Chirality::CounterClockwise, false) => false,
                    (Chirality::CounterClockwise, true) => true,
                    (Chirality::None, _) => unreachable!("filtered out above"),
                    (_, _) => {
                        unreachable!("centers is filtered to is_tetrahedral() above")
                    }
                })
            })
            .collect()
    }

    /// For every shared-candidate coupling component in `mol` recorded as
    /// abstained (production's own `ez_shared_bond_abstains`), NONE of its
    /// ends' candidate bonds carry a resolved marker -- proving "abstain"
    /// really means untouched, never a partially-applied plan left behind
    /// by an early-return path that ran after some `ez_marker` inserts
    /// already happened.
    ///
    /// Deliberately NOT asserting the stronger "0 or all of an end's own
    /// 2 candidates are marked" for every end unconditionally: that is
    /// FALSE in a real, legitimate case measured directly on one of the
    /// [`EZ_SHARED_CARRIER_FULLY_RESOLVED`] fixtures -- an end
    /// with zero raw geometry of its own (`reference_up` returns `None`,
    /// so `end_votes` contributes nothing) can still be a bystander whose
    /// shared bond gets marked purely by its *coupled partner's* own vote.
    /// That end's OTHER (private) candidate correctly stays unmarked. This
    /// is not corruption -- it is what `geometry_fingerprint`/
    /// `tetrahedral_fingerprint` (checked separately) confirm produces the
    /// right molecule either way -- so it must not be flagged here as if
    /// it were a bug.
    fn assert_no_partial_marker_application(mol: &Molecule) {
        let ends = CanonicalWriter::compute_stereo_alkene_ends(mol);
        let components = CanonicalWriter::coupling_components(mol, &ends);
        let (ranks, _) = winning_individualized_ranks(mol);
        let mut writer = CanonicalWriter::new(mol, &ranks);
        writer.resolve_ez_markers();

        for component in &components {
            let component_abstained = component
                .iter()
                .any(|end| writer.ez_shared_bond_abstains.contains(end));
            if !component_abstained {
                continue;
            }
            for &end in component {
                let subs = CanonicalWriter::substituents(mol, end);
                let marked = subs
                    .iter()
                    .filter(|(_, b)| writer.ez_marker.contains_key(b))
                    .count();
                assert_eq!(
                    marked, 0,
                    "end {end:?} belongs to a component recorded as \
                     abstained, but still has a marked candidate bond"
                );
            }
        }
    }

    /// 18 of these 19 are real-corpus molecules from the 282-molecule
    /// `has_ez_marker` diagnosis subset (issue #149) where two independently
    /// stereogenic double bonds share one physical candidate carrier bond --
    /// the joint component solver (`resolve_component_jointly`) resolves
    /// **all** of them fully: canonical output is invariant under every
    /// relabeling tested below, not just one. The 19th (see this list's own
    /// doc comment) is a separately-found real corpus molecule with the same
    /// shared-carrier shape, added once measured to resolve with the same
    /// rigor.
    ///
    /// Originally split 10/8: the last 8 (`CC1=C2CC[C@H](/C=N/N=C(N)N)...`
    /// through `CCO/C(O)=C(\C1=NCCN1)...`) were a genuine, still-open
    /// residual until `compute_stereo_alkene_ends` gained a ring-size gate
    /// (`double_bond_endocyclic_in_small_ring`, issue #149's Wave 2C audit,
    /// `docs/rfcs/ez_ring_constrained_residual_audit.md`) that excludes an
    /// end whose own double bond is endocyclic in a ring smaller than 8
    /// atoms -- exactly the shape all 8 shared (an endocyclic culprit with
    /// no real stereochemical freedom, coupled via a shared candidate bond
    /// to a genuinely stereogenic partner). Measured directly: all 8 now
    /// converge to one canonical output, and the original 10 (re-verified,
    /// not assumed -- 5 of them share the identical endocyclic shape per
    /// the audit's own finding) remain fully resolved. Merged into one list
    /// since the split no longer describes anything real.
    const EZ_SHARED_CARRIER_FULLY_RESOLVED: &[&str] = &[
        r"CCCCC/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
        r"O=C(Nc1ccc(C[C@H](/N=c2\c(O)c(O)\c2=N/Cc2ccccc2)C(=O)O)cc1)c1c(Cl)cncc1Cl",
        r"CCC/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
        r"O=C(Nc1ccc(C[C@H](/N=c2\c(O)c(O)\c2=N/c2ccccc2)C(=O)O)cc1)c1c(Cl)cncc1Cl",
        r"CC(C)(C)/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
        r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1cccc(I)c1",
        r"CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1ccccc1C(F)(F)F",
        r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1ccccc1OC",
        r"CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1ccc([N+](=O)[O-])cc1",
        r"CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1cccc(C(F)(F)F)c1",
        r"CC1=C2CC[C@H](/C=N/N=C(N)N)[C@@]2(C)CC/C1=N\N=C(N)N",
        r"CC1=C2CC[C@@H](/C=N/N=C(N)N)[C@@]2(C)CC/C1=N\N=C(N)N",
        r"COC(=O)/C=C/[C@H]1CCC2=C(C)/C(=N/N=C(N)N)CC[C@@]21C",
        r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1ccc(I)cc1",
        r"CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1ccccc1[N+](=O)[O-]",
        r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1cccc(C)c1",
        r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1cccc(OC)c1",
        r"CCO/C(O)=C(\C1=NCCN1)c1nnc(N)s1",
        r"OC(=O)[C@H](Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)/N=c1/c(c(c1O)O)=N/CCCCC",
    ];

    /// Three held-out issue #149 residuals measured by the 2026-09-06
    /// 64-relabeling audit. Each is a genuine coupled E/Z system whose
    /// aromatic direction-stash carriers can still produce two valid,
    /// semantically identical canonical spellings. Keep these separate from
    /// [`EZ_SHARED_CARRIER_FULLY_RESOLVED`]: they are deliberately not a
    /// claim that the general ranking problem is solved.
    const EZ_SHARED_CARRIER_HELD_OUT_RESIDUALS: &[&str] = &[
        r"CC/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
        r"CCCC(C)/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
        r"COCC/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
    ];

    /// Enumerate every DFS priority mask for the candidate bonds of each
    /// coupled E/Z component. A future production search must not mistake
    /// this local traversal space for a complete solution: the observed
    /// spellings expose no common geometry-preserving output here.
    #[test]
    fn issue503_component_traversal_only_has_no_common_geometry_preserving_output() {
        let observed_pairs = [
            (
                r"c/3(c(/c(c3=N\CC)=N\[C@@H](Cc1ccc(NC(=O)c2c(cncc2Cl)Cl)cc1)C(O)=O)O)O",
                r"c3(c(c(/c3=N/CC)=N\[C@@H](Cc1ccc(NC(c2c(Cl)cncc2Cl)=O)cc1)C(O)=O)O)O",
            ),
            (
                r"c/1(c(/c(c1=N\[C@H](C(O)=O)Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)=N\C(C)CCC)O)O",
                r"c1(c(c(/c1=N/[C@H](C(O)=O)Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)=N\C(C)CCC)O)O",
            ),
            (
                r"c/1(O)c(O)/c(=N\[C@@H](Cc3ccc(cc3)NC(=O)c2c(Cl)cncc2Cl)C(=O)O)c1=N\CCOC",
                r"c1(O)c(O)c(=N/[C@@H](Cc3ccc(cc3)NC(=O)c2c(Cl)cncc2Cl)C(=O)O)\c1=N\CCOC",
            ),
        ];

        for (&input, &(observed_a, observed_b)) in EZ_SHARED_CARRIER_HELD_OUT_RESIDUALS
            .iter()
            .zip(observed_pairs.iter())
        {
            let mut per_spelling = Vec::new();
            for spelling in [input, observed_a, observed_b] {
                let mol = parse(spelling).unwrap();
                let geometry = geometry_fingerprint(&mol);
                let ends = CanonicalWriter::compute_stereo_alkene_ends(&mol);
                let mut component_bonds: Vec<_> = CanonicalWriter::coupling_components(&mol, &ends)
                    .into_iter()
                    .filter(|component| component.len() > 1)
                    .flat_map(|component| {
                        component.into_iter().flat_map(|end| {
                            CanonicalWriter::substituents(&mol, end)
                                .into_iter()
                                .map(|(_, bond)| bond)
                        })
                    })
                    .collect();
                component_bonds.sort_unstable_by_key(|bond| bond.0);
                component_bonds.dedup();
                assert!(component_bonds.len() <= 8, "diagnostic bound exceeded");

                let (ranks, _) = winning_individualized_ranks(&mol);
                let mut outputs = HashSet::new();
                for mask in 0usize..(1usize << component_bonds.len()) {
                    let preferences: HashMap<BondIdx, bool> = component_bonds
                        .iter()
                        .enumerate()
                        .map(|(index, &bond)| (bond, (mask & (1 << index)) != 0))
                        .collect();
                    let mut writer = CanonicalWriter::new(&mol, &ranks);
                    writer.traversal_bond_preference = Some(preferences);
                    let output = writer.write_all();
                    let reparsed = parse(&output).unwrap_or_else(|e| {
                        panic!("{spelling}: traversal candidate did not parse: {e}: {output}")
                    });
                    if geometry_fingerprint(&reparsed) == geometry {
                        outputs.insert(output);
                    }
                }
                per_spelling.push(outputs);
            }
            let common = per_spelling
                .into_iter()
                .reduce(|left, right| left.intersection(&right).cloned().collect())
                .unwrap();
            assert!(
                common.is_empty(),
                "{input}: candidate-bond DFS priorities unexpectedly found a common output; \
                 update the complete-plan boundary before relying on it"
            );
        }
    }

    /// The input spelling variants behind #503 must first reduce to exactly
    /// the same chemical E/Z facts.  This pins the phase boundary for the
    /// future full component solver: it may change output token locations,
    /// but it must preserve these rank-keyed relations.
    #[test]
    fn issue503_ez_geometry_facts_are_spelling_invariant() {
        let observed_pairs = [
            (
                r"c/3(c(/c(c3=N\CC)=N\[C@@H](Cc1ccc(NC(=O)c2c(cncc2Cl)Cl)cc1)C(O)=O)O)O",
                r"c3(c(c(/c3=N/CC)=N\[C@@H](Cc1ccc(NC(c2c(Cl)cncc2Cl)=O)cc1)C(O)=O)O)O",
            ),
            (
                r"c/1(c(/c(c1=N\[C@H](C(O)=O)Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)=N\C(C)CCC)O)O",
                r"c1(c(c(/c1=N/[C@H](C(O)=O)Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)=N\C(C)CCC)O)O",
            ),
            (
                r"c/1(O)c(O)/c(=N\[C@@H](Cc3ccc(cc3)NC(=O)c2c(Cl)cncc2Cl)C(=O)O)c1=N\CCOC",
                r"c1(O)c(O)c(=N/[C@@H](Cc3ccc(cc3)NC(=O)c2c(Cl)cncc2Cl)C(=O)O)\c1=N\CCOC",
            ),
        ];

        for (&input, &(observed_a, observed_b)) in EZ_SHARED_CARRIER_HELD_OUT_RESIDUALS
            .iter()
            .zip(observed_pairs.iter())
        {
            let expected = ez_geometry_signature(&parse(input).unwrap());
            assert!(
                expected.len() >= 2,
                "{input}: setup must expose both coupled E/Z facts"
            );
            for spelling in [observed_a, observed_b] {
                assert_eq!(
                    ez_geometry_signature(&parse(spelling).unwrap()),
                    expected,
                    "{input}: spelling '{spelling}' changed the extracted E/Z facts"
                );
            }
        }
    }

    /// Before a whole-component E/Z plan can choose token polarity, every
    /// equivalent spelling must expose the same canonical output occurrences
    /// for the coupled carrier bonds.  Ring bonds deliberately expose both
    /// legal digit occurrences; the future solver, not parse order, decides
    /// which one carries a marker.
    #[test]
    fn issue503_ez_output_slots_are_spelling_invariant() {
        let observed_pairs = [
            (
                r"c/3(c(/c(c3=N\CC)=N\[C@@H](Cc1ccc(NC(=O)c2c(cncc2Cl)Cl)cc1)C(O)=O)O)O",
                r"c3(c(c(/c3=N/CC)=N\[C@@H](Cc1ccc(NC(c2c(Cl)cncc2Cl)=O)cc1)C(O)=O)O)O",
            ),
            (
                r"c/1(c(/c(c1=N\[C@H](C(O)=O)Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)=N\C(C)CCC)O)O",
                r"c1(c(c(/c1=N/[C@H](C(O)=O)Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)=N\C(C)CCC)O)O",
            ),
            (
                r"c/1(O)c(O)/c(=N\[C@@H](Cc3ccc(cc3)NC(=O)c2c(Cl)cncc2Cl)C(=O)O)c1=N\CCOC",
                r"c1(O)c(O)c(=N/[C@@H](Cc3ccc(cc3)NC(=O)c2c(Cl)cncc2Cl)C(=O)O)\c1=N\CCOC",
            ),
        ];

        for (&input, &(observed_a, observed_b)) in EZ_SHARED_CARRIER_HELD_OUT_RESIDUALS
            .iter()
            .zip(observed_pairs.iter())
        {
            let expected = ez_output_slot_signature(&parse(input).unwrap());
            assert!(
                !expected.is_empty(),
                "{input}: setup must expose candidate E/Z output slots"
            );
            for spelling in [observed_a, observed_b] {
                assert_eq!(
                    ez_output_slot_signature(&parse(spelling).unwrap()),
                    expected,
                    "{input}: spelling '{spelling}' changed canonical output slots"
                );
            }
        }
    }

    /// The output-slot gate is intentionally component-local.  Before a new
    /// #503 design changes canonical traversal, establish whether the *full*
    /// rank-keyed DFS tree and ring-token layout already agree across the
    /// equivalent spellings.  A mismatch would put the next investigation in
    /// traversal/rank construction; agreement leaves only the writer's
    /// stereochemical token and parity handling.
    #[test]
    fn issue503_full_canonical_skeleton_is_spelling_invariant() {
        let observed_pairs = [
            (
                r"c/3(c(/c(c3=N\CC)=N\[C@@H](Cc1ccc(NC(=O)c2c(cncc2Cl)Cl)cc1)C(O)=O)O)O",
                r"c3(c(c(/c3=N/CC)=N\[C@@H](Cc1ccc(NC(c2c(Cl)cncc2Cl)=O)cc1)C(O)=O)O)O",
            ),
            (
                r"c/1(c(/c(c1=N\[C@H](C(O)=O)Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)=N\C(C)CCC)O)O",
                r"c1(c(c(/c1=N/[C@H](C(O)=O)Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)=N\C(C)CCC)O)O",
            ),
            (
                r"c/1(O)c(O)/c(=N\[C@@H](Cc3ccc(cc3)NC(=O)c2c(Cl)cncc2Cl)C(=O)O)c1=N\CCOC",
                r"c1(O)c(O)c(=N/[C@@H](Cc3ccc(cc3)NC(=O)c2c(Cl)cncc2Cl)C(=O)O)\c1=N\CCOC",
            ),
        ];

        for (&input, &(observed_a, observed_b)) in EZ_SHARED_CARRIER_HELD_OUT_RESIDUALS
            .iter()
            .zip(observed_pairs.iter())
        {
            let expected = canonical_skeleton_signature(&parse(input).unwrap());
            for spelling in [observed_a, observed_b] {
                assert_eq!(
                    canonical_skeleton_signature(&parse(spelling).unwrap()),
                    expected,
                    "{input}: spelling '{spelling}' changed full canonical DFS skeleton"
                );
            }
        }
    }

    /// Once all tree and ring occurrences agree, compare the rendered
    /// canonical strings with only E/Z direction characters removed.  This
    /// does not validate chemistry; it isolates whether the remaining #503
    /// divergence is exclusively the placement/polarity of those characters
    /// rather than atom, branch, ring, or tetrahedral-token serialization.
    #[test]
    fn issue503_directionless_canonical_text_is_spelling_invariant() {
        let observed_pairs = [
            (
                r"c/3(c(/c(c3=N\CC)=N\[C@@H](Cc1ccc(NC(=O)c2c(cncc2Cl)Cl)cc1)C(O)=O)O)O",
                r"c3(c(c(/c3=N/CC)=N\[C@@H](Cc1ccc(NC(c2c(Cl)cncc2Cl)=O)cc1)C(O)=O)O)O",
            ),
            (
                r"c/1(c(/c(c1=N\[C@H](C(O)=O)Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)=N\C(C)CCC)O)O",
                r"c1(c(c(/c1=N/[C@H](C(O)=O)Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)=N\C(C)CCC)O)O",
            ),
            (
                r"c/1(O)c(O)/c(=N\[C@@H](Cc3ccc(cc3)NC(=O)c2c(Cl)cncc2Cl)C(=O)O)c1=N\CCOC",
                r"c1(O)c(O)c(=N/[C@@H](Cc3ccc(cc3)NC(=O)c2c(Cl)cncc2Cl)C(=O)O)\c1=N\CCOC",
            ),
        ];

        for (&input, &(observed_a, observed_b)) in EZ_SHARED_CARRIER_HELD_OUT_RESIDUALS
            .iter()
            .zip(observed_pairs.iter())
        {
            let directionless = |smiles: &str| {
                canonical_smiles(&parse(smiles).unwrap())
                    .chars()
                    .filter(|ch| !matches!(ch, '/' | '\\'))
                    .collect::<String>()
            };
            let expected = directionless(input);
            for spelling in [observed_a, observed_b] {
                assert_eq!(
                    directionless(spelling),
                    expected,
                    "{input}: spelling '{spelling}' changed canonical text beyond E/Z directions"
                );
            }
        }
    }

    /// The slot-space probes retain candidates using `geometry_fingerprint`.
    /// Establish that this older diagnostic has the same spelling-invariant
    /// precondition as the rank-keyed geometry-fact extractor before using an
    /// empty output intersection as evidence about the writer.
    #[test]
    fn issue503_geometry_fingerprint_is_spelling_invariant() {
        let observed_pairs = [
            (
                r"c/3(c(/c(c3=N\CC)=N\[C@@H](Cc1ccc(NC(=O)c2c(cncc2Cl)Cl)cc1)C(O)=O)O)O",
                r"c3(c(c(/c3=N/CC)=N\[C@@H](Cc1ccc(NC(c2c(Cl)cncc2Cl)=O)cc1)C(O)=O)O)O",
            ),
            (
                r"c/1(c(/c(c1=N\[C@H](C(O)=O)Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)=N\C(C)CCC)O)O",
                r"c1(c(c(/c1=N/[C@H](C(O)=O)Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)=N\C(C)CCC)O)O",
            ),
            (
                r"c/1(O)c(O)/c(=N\[C@@H](Cc3ccc(cc3)NC(=O)c2c(Cl)cncc2Cl)C(=O)O)c1=N\CCOC",
                r"c1(O)c(O)c(=N/[C@@H](Cc3ccc(cc3)NC(=O)c2c(Cl)cncc2Cl)C(=O)O)\c1=N\CCOC",
            ),
        ];

        for (&input, &(observed_a, observed_b)) in EZ_SHARED_CARRIER_HELD_OUT_RESIDUALS
            .iter()
            .zip(observed_pairs.iter())
        {
            let expected = geometry_fingerprint(&parse(input).unwrap());
            for spelling in [observed_a, observed_b] {
                assert_eq!(
                    geometry_fingerprint(&parse(spelling).unwrap()),
                    expected,
                    "{input}: spelling '{spelling}' changed the slot-probe geometry fingerprint"
                );
            }
        }
    }

    /// The complete slot probe includes every actual raw directional carrier,
    /// not only substituent bonds reached through a coupled component.  The
    /// latter misses aromatic-stash positions in the #503 canonical output
    /// and falsely makes the semantic candidate intersection appear empty.
    #[test]
    fn issue503_complete_slot_probe_covers_every_canonical_direction_carrier() {
        let observed_pairs = [
            (
                r"c/3(c(/c(c3=N\CC)=N\[C@@H](Cc1ccc(NC(=O)c2c(cncc2Cl)Cl)cc1)C(O)=O)O)O",
                r"c3(c(c(/c3=N/CC)=N\[C@@H](Cc1ccc(NC(c2c(Cl)cncc2Cl)=O)cc1)C(O)=O)O)O",
            ),
            (
                r"c/1(c(/c(c1=N\[C@H](C(O)=O)Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)=N\C(C)CCC)O)O",
                r"c1(c(c(/c1=N/[C@H](C(O)=O)Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)=N\C(C)CCC)O)O",
            ),
            (
                r"c/1(O)c(O)/c(=N\[C@@H](Cc3ccc(cc3)NC(=O)c2c(Cl)cncc2Cl)C(=O)O)c1=N\CCOC",
                r"c1(O)c(O)c(=N/[C@@H](Cc3ccc(cc3)NC(=O)c2c(Cl)cncc2Cl)C(=O)O)\c1=N\CCOC",
            ),
        ];

        for (&input, &(observed_a, observed_b)) in EZ_SHARED_CARRIER_HELD_OUT_RESIDUALS
            .iter()
            .zip(observed_pairs.iter())
        {
            for spelling in [input, observed_a, observed_b] {
                let output = canonical_smiles(&parse(spelling).unwrap());
                let mol = parse(&output).unwrap();
                let ends = CanonicalWriter::compute_stereo_alkene_ends(&mol);
                let candidates: HashSet<BondIdx> =
                    CanonicalWriter::coupling_components(&mol, &ends)
                        .into_iter()
                        .flatten()
                        .flat_map(|end| CanonicalWriter::substituents(&mol, end))
                        .map(|(_, bond)| bond)
                        .collect();
                let ranks = morgan_ranks(&mol);
                let writer = CanonicalWriter::new(&mol, &ranks);
                let raw_directional: HashSet<BondIdx> = mol
                    .bonds()
                    .filter_map(|(bond, _)| {
                        matches!(
                            writer.raw_input_direction(bond),
                            Some(BondOrder::Up | BondOrder::Down)
                        )
                        .then_some(bond)
                    })
                    .collect();
                assert!(
                    raw_directional
                        .iter()
                        .any(|bond| !candidates.contains(bond)),
                    "{spelling}: canonical output '{output}' no longer exercises the known \
                     component-only slot-universe omission"
                );
                let complete_candidates: HashSet<_> =
                    candidates.union(&raw_directional).copied().collect();
                assert!(
                    raw_directional
                        .iter()
                        .all(|bond| complete_candidates.contains(bond)),
                    "{spelling}: complete slot universe omitted a canonical direction carrier"
                );
            }
        }
    }

    /// The closing side of a ring digit is syntactically capable of carrying
    /// its directional token (the parser applies the required orientation
    /// flip).  Exhaust the two remaining *local* writer freedoms together:
    /// every candidate-bond DFS priority and every candidate ring-marker
    /// side.  A positive result would identify an admissible bounded plan;
    /// a negative result proves that neither independent local preference nor
    /// their combination can close #503.
    #[test]
    fn issue503_local_traversal_and_ring_marker_side_have_no_common_output() {
        let observed_pairs = [
            (
                r"c/3(c(/c(c3=N\CC)=N\[C@@H](Cc1ccc(NC(=O)c2c(cncc2Cl)Cl)cc1)C(O)=O)O)O",
                r"c3(c(c(/c3=N/CC)=N\[C@@H](Cc1ccc(NC(c2c(Cl)cncc2Cl)=O)cc1)C(O)=O)O)O",
            ),
            (
                r"c/1(c(/c(c1=N\[C@H](C(O)=O)Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)=N\C(C)CCC)O)O",
                r"c1(c(c(/c1=N/[C@H](C(O)=O)Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)=N\C(C)CCC)O)O",
            ),
            (
                r"c/1(O)c(O)/c(=N\[C@@H](Cc3ccc(cc3)NC(=O)c2c(Cl)cncc2Cl)C(=O)O)c1=N\CCOC",
                r"c1(O)c(O)c(=N/[C@@H](Cc3ccc(cc3)NC(=O)c2c(Cl)cncc2Cl)C(=O)O)\c1=N\CCOC",
            ),
        ];

        for (&input, &(observed_a, observed_b)) in EZ_SHARED_CARRIER_HELD_OUT_RESIDUALS
            .iter()
            .zip(observed_pairs.iter())
        {
            let mut per_spelling = Vec::new();
            for spelling in [input, observed_a, observed_b] {
                let mol = parse(spelling).unwrap();
                let geometry = geometry_fingerprint(&mol);
                let ends = CanonicalWriter::compute_stereo_alkene_ends(&mol);
                let mut component_bonds: Vec<_> = CanonicalWriter::coupling_components(&mol, &ends)
                    .into_iter()
                    .flatten()
                    .flat_map(|end| CanonicalWriter::substituents(&mol, end))
                    .map(|(_, bond)| bond)
                    .collect();
                component_bonds.sort_unstable_by_key(|bond| bond.0);
                component_bonds.dedup();
                assert!(component_bonds.len() <= 8, "diagnostic bound exceeded");
                let (ranks, _) = winning_individualized_ranks(&mol);
                let mut outputs = HashSet::new();
                for traversal_mask in 0usize..(1usize << component_bonds.len()) {
                    let preferences: HashMap<BondIdx, bool> = component_bonds
                        .iter()
                        .enumerate()
                        .map(|(index, &bond)| (bond, (traversal_mask & (1 << index)) != 0))
                        .collect();
                    for close_mask in 0usize..(1usize << component_bonds.len()) {
                        let mut writer = CanonicalWriter::new(&mol, &ranks);
                        writer.traversal_bond_preference = Some(preferences.clone());
                        writer.ring_marker_on_close = component_bonds
                            .iter()
                            .enumerate()
                            .filter_map(|(index, &bond)| {
                                ((close_mask & (1 << index)) != 0).then_some(bond)
                            })
                            .collect();
                        let output = writer.write_all();
                        let reparsed = parse(&output).unwrap_or_else(|e| {
                            panic!("{spelling}: local-plan candidate did not parse: {e}: {output}")
                        });
                        if geometry_fingerprint(&reparsed) == geometry {
                            outputs.insert(output);
                        }
                    }
                }
                per_spelling.push(outputs);
            }
            let common = per_spelling
                .into_iter()
                .reduce(|left, right| left.intersection(&right).cloned().collect())
                .unwrap();
            assert!(
                common.is_empty(),
                "{input}: local writer plan unexpectedly found {common:?}; \
                 promote it from a probe only after its full contract is verified"
            );
        }
    }

    /// Local traversal choices cannot solve the #503 residuals.  This probe
    /// holds the canonical DFS skeleton fixed, then exhausts the remaining
    /// grammar-level freedom directly: every component candidate bond can be
    /// plain, `/`, or `\\`, and every ring candidate may print its token at
    /// either legal digit occurrence.  Candidates are reparsed and retained
    /// only when their complete E/Z geometry fingerprint agrees with the
    /// source.  It is intentionally test-only: a future production solver
    /// must derive a rank-keyed component plan, not adopt a parse-indexed
    /// enumeration.
    #[test]
    fn issue503_full_slot_polarity_space_has_common_output() {
        let observed_pairs = [
            (
                r"c/3(c(/c(c3=N\CC)=N\[C@@H](Cc1ccc(NC(=O)c2c(cncc2Cl)Cl)cc1)C(O)=O)O)O",
                r"c3(c(c(/c3=N/CC)=N\[C@@H](Cc1ccc(NC(c2c(Cl)cncc2Cl)=O)cc1)C(O)=O)O)O",
            ),
            (
                r"c/1(c(/c(c1=N\[C@H](C(O)=O)Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)=N\C(C)CCC)O)O",
                r"c1(c(c(/c1=N/[C@H](C(O)=O)Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)=N\C(C)CCC)O)O",
            ),
            (
                r"c/1(O)c(O)/c(=N\[C@@H](Cc3ccc(cc3)NC(=O)c2c(Cl)cncc2Cl)C(=O)O)c1=N\CCOC",
                r"c1(O)c(O)c(=N/[C@@H](Cc3ccc(cc3)NC(=O)c2c(Cl)cncc2Cl)C(=O)O)\c1=N\CCOC",
            ),
        ];

        for (&input, &(observed_a, observed_b)) in EZ_SHARED_CARRIER_HELD_OUT_RESIDUALS
            .iter()
            .zip(observed_pairs.iter())
        {
            let mut per_spelling = Vec::new();
            let mut reranked_per_spelling = Vec::new();
            for spelling in [input, observed_a, observed_b] {
                let mol = parse(spelling).unwrap();
                let geometry = geometry_fingerprint(&mol);
                let ends = CanonicalWriter::compute_stereo_alkene_ends(&mol);
                let mut component_bonds: Vec<_> = CanonicalWriter::coupling_components(&mol, &ends)
                    .into_iter()
                    .flatten()
                    .flat_map(|end| CanonicalWriter::substituents(&mol, end))
                    .map(|(_, bond)| bond)
                    .collect();
                let (ranks, _) = winning_individualized_ranks(&mol);
                let raw_directional: Vec<_> = (0..mol.bond_count())
                    .map(|index| BondIdx(index as u32))
                    .filter(|&bond| {
                        matches!(
                            CanonicalWriter::new(&mol, &ranks).raw_input_direction(bond),
                            Some(BondOrder::Up | BondOrder::Down)
                        )
                    })
                    .collect();
                component_bonds.extend(raw_directional);
                component_bonds.sort_unstable_by_key(|bond| bond.0);
                component_bonds.dedup();
                assert!(component_bonds.len() <= 8, "diagnostic bound exceeded");

                let mut skeleton = CanonicalWriter::new(&mol, &ranks);
                skeleton.find_ring_closures(&skeleton.canonical_atom_list());
                let ring_component_bonds: Vec<_> = component_bonds
                    .iter()
                    .copied()
                    .filter(|bond| skeleton.ring_bonds[bond.0 as usize])
                    .collect();
                assert!(
                    ring_component_bonds.len() <= 8,
                    "ring-side diagnostic bound exceeded"
                );

                let mut outputs = HashSet::new();
                let mut reranked_outputs = HashSet::new();
                let polarity_plans = 3usize.pow(component_bonds.len() as u32);
                for polarity_code in 0..polarity_plans {
                    let mut code = polarity_code;
                    let orders: HashMap<BondIdx, BondOrder> = component_bonds
                        .iter()
                        .map(|&bond| {
                            let choice = code % 3;
                            code /= 3;
                            let order = match choice {
                                0 => CanonicalWriter::plain_order(mol.bond(bond).order),
                                1 => BondOrder::Up,
                                2 => BondOrder::Down,
                                _ => unreachable!(),
                            };
                            (bond, order)
                        })
                        .collect();
                    for close_mask in 0usize..(1usize << ring_component_bonds.len()) {
                        let mut writer = CanonicalWriter::new(&mol, &ranks);
                        writer.direct_ez_orders = Some(orders.clone());
                        writer.ring_marker_on_close = ring_component_bonds
                            .iter()
                            .enumerate()
                            .filter_map(|(index, &bond)| {
                                ((close_mask & (1 << index)) != 0).then_some(bond)
                            })
                            .collect();
                        let output = writer.write_all();
                        let reparsed = parse(&output).unwrap_or_else(|e| {
                            panic!("{spelling}: full-slot candidate did not parse: {e}: {output}")
                        });
                        if geometry_fingerprint(&reparsed) == geometry {
                            outputs.insert(output.clone());
                            let (_, reranked) = winning_individualized_ranks_single(&reparsed);
                            let reranked_mol = parse(&reranked).unwrap_or_else(|e| {
                                panic!(
                                    "{spelling}: reranked carrier candidate did not parse: {e}: {reranked}"
                                )
                            });
                            if geometry_fingerprint(&reranked_mol) == geometry {
                                reranked_outputs.insert(reranked);
                            }
                        }
                    }
                }
                per_spelling.push(outputs);
                reranked_per_spelling.push(reranked_outputs);
            }
            let selected: Vec<_> = per_spelling
                .iter()
                .map(|outputs| {
                    outputs
                        .iter()
                        .min()
                        .cloned()
                        .expect("each spelling must have a semantic candidate")
                })
                .collect();
            assert!(
                selected.windows(2).all(|pair| pair[0] == pair[1]),
                "{input}: independent full-slot minimization diverged: {selected:?}"
            );
            let common = per_spelling
                .into_iter()
                .reduce(|left, right| left.intersection(&right).cloned().collect())
                .unwrap();
            assert!(
                !common.is_empty(),
                "{input}: the complete carrier space unexpectedly has no common output"
            );
            let reranked_selected: Vec<_> = reranked_per_spelling
                .iter()
                .map(|outputs| {
                    outputs
                        .iter()
                        .min()
                        .cloned()
                        .expect("each spelling must have a reranked semantic candidate")
                })
                .collect();
            assert!(
                reranked_selected.windows(2).all(|pair| pair[0] == pair[1]),
                "{input}: independent reranked full-slot minimization diverged: \
                 {reranked_selected:?}"
            );
            let reranked_common = reranked_per_spelling
                .into_iter()
                .reduce(|left, right| left.intersection(&right).cloned().collect())
                .unwrap();
            assert!(
                !reranked_common.is_empty(),
                "{input}: reranking the complete carrier space lost every common output"
            );
        }
    }

    /// The bounded complete-slot planner must converge the formerly
    /// fail-closed residual spellings and admit their self-stable keys.
    #[test]
    fn ez_shared_carrier_held_out_residuals_converge_to_stable_keys() {
        let observed_pairs = [
            (
                r"c/3(c(/c(c3=N\CC)=N\[C@@H](Cc1ccc(NC(=O)c2c(cncc2Cl)Cl)cc1)C(O)=O)O)O",
                r"c3(c(c(/c3=N/CC)=N\[C@@H](Cc1ccc(NC(=O)c2c(cncc2Cl)Cl)cc1)C(O)=O)O)O",
            ),
            (
                r"c/1(c(/c(c1=N\[C@H](C(O)=O)Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)=N\C(C)CCC)O)O",
                r"c1(c(c(/c1=N/[C@H](C(O)=O)Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)=N\C(C)CCC)O)O",
            ),
            (
                r"c/1(O)c(O)/c(=N\[C@@H](Cc3ccc(cc3)NC(=O)c2c(Cl)cncc2Cl)C(=O)O)c1=N\CCOC",
                r"c1(O)c(O)c(=N/[C@@H](Cc3ccc(cc3)NC(=O)c2c(Cl)cncc2Cl)C(=O)O)\c1=N\CCOC",
            ),
        ];

        for (&s, &(a, b)) in EZ_SHARED_CARRIER_HELD_OUT_RESIDUALS
            .iter()
            .zip(observed_pairs.iter())
        {
            let mol = parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"));
            let outputs: HashSet<String> = [a, b]
                .into_iter()
                .map(|variant| canonical_smiles(&parse(variant).unwrap()))
                .collect();
            assert_eq!(outputs.len(), 1, "'{s}': held-out spellings must converge");
            assert!(
                canonical_smiles_stable_key(&mol).is_some(),
                "'{s}': converged aromatic E/Z residual must become a cache key"
            );
        }
    }

    /// Re-run the held-out Wave 3 residuals through the same deterministic
    /// relabeling axis used by the corpus audit. The complete planner must
    /// converge every relabeling and observed equivalent spelling to one
    /// semantic, self-stable canonical string.
    #[test]
    fn ez_shared_carrier_held_out_residuals_converge_under_relabeling() {
        let observed_pairs = [
            (
                r"c/3(c(/c(c3=N\CC)=N\[C@@H](Cc1ccc(NC(=O)c2c(cncc2Cl)Cl)cc1)C(O)=O)O)O",
                r"c3(c(c(/c3=N/CC)=N\[C@@H](Cc1ccc(NC(=O)c2c(cncc2Cl)Cl)cc1)C(O)=O)O)O",
            ),
            (
                r"c/1(c(/c(c1=N\[C@H](C(O)=O)Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)=N\C(C)CCC)O)O",
                r"c1(c(c(/c1=N/[C@H](C(O)=O)Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)=N\C(C)CCC)O)O",
            ),
            (
                r"c/1(O)c(O)/c(=N\[C@@H](Cc3ccc(cc3)NC(=O)c2c(Cl)cncc2Cl)C(=O)O)c1=N\CCOC",
                r"c1(O)c(O)c(=N/[C@@H](Cc3ccc(cc3)NC(=O)c2c(Cl)cncc2Cl)C(=O)O)\c1=N\CCOC",
            ),
        ];

        for (&input, &(observed_a, observed_b)) in EZ_SHARED_CARRIER_HELD_OUT_RESIDUALS
            .iter()
            .zip(observed_pairs.iter())
        {
            let mol = parse(input).unwrap_or_else(|e| panic!("parse '{input}': {e}"));
            let before_geo = geometry_fingerprint(&mol);
            let n = mol.atom_count();
            let mut variants = vec![
                relabel_molecule_preserving_ez(&mol, &(0..n).collect::<Vec<_>>()),
                relabel_molecule_preserving_ez(&mol, &(0..n).rev().collect::<Vec<_>>()),
            ];
            for seed in 0..256u64 {
                variants.push(relabel_molecule_preserving_ez(
                    &mol,
                    &deterministic_permutation(n, seed),
                ));
            }

            let relabeled_outputs: HashSet<String> = variants
                .iter()
                .map(|variant| {
                    assert_eq!(
                        geometry_fingerprint(variant),
                        before_geo,
                        "'{input}': relabeling changed the encoded E/Z geometry"
                    );
                    canonical_smiles(variant)
                })
                .collect();
            assert_eq!(
                relabeled_outputs.len(),
                1,
                "'{input}': relabeling-only audit must remain deterministic"
            );

            let mut outputs = relabeled_outputs;
            for observed in [observed_a, observed_b] {
                let observed_mol = parse(observed)
                    .unwrap_or_else(|e| panic!("parse observed residual '{observed}': {e}"));
                assert_eq!(
                    geometry_fingerprint(&observed_mol),
                    before_geo,
                    "'{input}': observed residual spelling changed E/Z geometry"
                );
                outputs.insert(canonical_smiles(&observed_mol));
            }
            assert_eq!(outputs.len(), 1, "'{input}': 256-seed audit must converge");
            assert!(
                canonical_smiles_stable_key(&mol).is_some(),
                "'{input}': converged residual must have a stable key"
            );
        }
    }

    /// Proves all 19 [`EZ_SHARED_CARRIER_FULLY_RESOLVED`] fixtures are
    /// genuinely, fully permutation-invariant -- not just under the one
    /// relabeling a weaker test might check. Per fixture: the original
    /// parse, its reversed-atom-order relabeling, and 16 deterministic
    /// (seeded) Fisher-Yates relabelings (see [`deterministic_permutation`]
    /// for why `random_smiles` is not used) -- 18 spellings total -- must
    /// all canonicalize to exactly ONE string. That string must be
    /// idempotent (`canonical(canonical(x)) == canonical(x)`), must
    /// re-parse without error, and a reparse of it must preserve both the
    /// E/Z ([`geometry_fingerprint`]) and tetrahedral
    /// ([`tetrahedral_fingerprint`]) stereo facts read from the original
    /// parse. Also asserts the joint solver never abstains for these 19
    /// (production's own `ez_shared_bond_abstains` record, not re-derived
    /// topology) and never applies a partial marker plan.
    ///
    /// The 19th fixture (two independently stereogenic exocyclic imines on
    /// the same four-membered ring, sharing the ring-closure bond each
    /// would otherwise use as a candidate carrier) was previously believed
    /// to be a genuine, documented residual that a fully general carrier
    /// choice could not resolve to one canonical string -- tracked by a
    /// separate, weaker test that only checked E/Z-preservation across two
    /// specific hand-picked respellings, not full permutation-invariance.
    /// Re-measured directly (not assumed) while auditing a stale-doc-comment
    /// flag in ROADMAP.md: it passes this test's full 18-spelling rigor,
    /// same as the other 18 -- the residual is resolved, and the split has
    /// been merged into this one list for the same reason the original 8/10
    /// split above was.
    #[test]
    fn ez_shared_carrier_fully_resolved_are_permutation_invariant() {
        for &s in EZ_SHARED_CARRIER_FULLY_RESOLVED {
            let mol = parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"));

            let (ranks, _) = winning_individualized_ranks(&mol);
            let mut writer = CanonicalWriter::new(&mol, &ranks);
            writer.resolve_ez_markers();
            assert!(
                writer.ez_shared_bond_abstains.is_empty(),
                "'{s}': expected the joint component solver to resolve \
                 without abstaining -- this fixture may need re-diagnosis"
            );
            assert_no_partial_marker_application(&mol);

            let variants = ez_carrier_test_variants(&mol);
            assert!(
                variants.len() >= 18,
                "sanity: at least 16 deterministic + 2 fixed relabelings"
            );

            let canonical_outputs: Vec<String> = variants.iter().map(canonical_smiles).collect();
            let unique: HashSet<&String> = canonical_outputs.iter().collect();
            assert_eq!(
                unique.len(),
                1,
                "'{s}': expected ONE canonical output across {} relabelings/\
                 markings (fully resolved), got {} distinct: {:?}",
                variants.len(),
                unique.len(),
                unique
            );
            let canon = canonical_outputs[0].clone();

            // Idempotence.
            let reparsed_once =
                parse(&canon).unwrap_or_else(|e| panic!("'{s}': reparse of '{canon}': {e}"));
            let canon_twice = canonical_smiles(&reparsed_once);
            assert_eq!(
                canon, canon_twice,
                "'{s}': canonical(canonical(x)) must equal canonical(x)"
            );

            // Zero corruption: E/Z and tetrahedral facts survive the
            // canonicalize -> reparse round trip.
            let before_geo = geometry_fingerprint(&mol);
            let after_geo = geometry_fingerprint(&reparsed_once);
            assert_eq!(
                before_geo, after_geo,
                "'{s}': E/Z geometry must be preserved by canonicalization"
            );
            let before_tet = tetrahedral_fingerprint(&mol);
            let after_tet = tetrahedral_fingerprint(&reparsed_once);
            assert_eq!(
                before_tet, after_tet,
                "'{s}': tetrahedral stereo must be preserved by canonicalization"
            );
            assert!(
                before_geo.iter().any(|f| f.is_some()),
                "test setup sanity: '{s}' must have at least one defined \
                 E/Z geometry fact"
            );
        }
    }

    /// Issue #390 witness: `O/N=C/C(C=N/O)=N\NC`. A prior version of
    /// `normalize_ez` silently changed this molecule's atom3=atom7 double
    /// bond from Z to E -- confirmed genuinely wrong (not just a different,
    /// equally valid spelling) via an independent RDKit `MolToInchi`/
    /// `GetStereo()` check, before any fix existed.
    ///
    /// Root cause, in two parts:
    /// 1. `resolve_ez_markers`'s carrier election for the ambiguous end
    ///    (2 candidates: the bond toward the atom1=atom2 double bond's own
    ///    side, raw-marked and load-bearing for atom1=atom2's geometry; the
    ///    bond toward the atom4=atom5 double bond's side, unmarked --
    ///    atom4=atom5 is itself genuinely undefined in this molecule, only
    ///    one of its two flanking bonds ever carries a mark, confirmed via
    ///    InChI's own `?` stereo descriptor for it) could elect the
    ///    unmarked candidate. That demotes the marked one -- silently
    ///    under-specifying atom1=atom2 -- while handing atom4=atom5 a
    ///    geometry it never had. Neither move is neutral. Fixed by
    ///    [`Self::is_load_bearing_elsewhere`]: an election must not demote a
    ///    raw-marked candidate that is some other, non-ambiguous double
    ///    bond's only geometric anchor.
    /// 2. Independently, `build_ez_groups`/`normalize_ez` chained the
    ///    elected carrier's group together with atom4=atom5's own group
    ///    (both touch the same connecting bond), and the group's sign was
    ///    seeded from a value that had already been re-oriented for one
    ///    specific DFS write direction -- which bond within a group gets
    ///    visited "forward" vs "backward" varies across candidate canonical
    ///    numberings for reasons unrelated to this bond's own geometry, so
    ///    the seeded sign could vary too. Fixed by splitting `normalize_ez`
    ///    into a mol-relative propagation step (always operates on
    ///    [`Self::effective_order`], never an already-write-oriented value)
    ///    and a write-perspective anchor-seeding step (uses the write atom
    ///    ONLY to decide the group's shared sign once, never to decide what
    ///    gets propagated).
    #[test]
    fn issue390_witness_geometry_preserved_and_stable() {
        // `C=N/O` has a marker on one side only, which specifies nothing
        // and is dropped (as RDKit does); the two defined bonds stay.
        let smi = "O/N=C/C(C=N/O)=N\\NC";
        let mol = parse(smi).unwrap();
        let canon = canonical_smiles(&mol);
        assert_eq!(canon, "O/N=C/C(C=NO)=N\\NC");

        let reparsed = parse(&canon).unwrap();
        let canon_twice = canonical_smiles(&reparsed);
        assert_eq!(
            canon, canon_twice,
            "canonical(canonical(x)) must equal canonical(x)"
        );

        let before_geo = geometry_fingerprint(&mol);
        let after_geo = geometry_fingerprint(&reparsed);
        assert_eq!(
            before_geo, after_geo,
            "E/Z geometry must be preserved by canonicalization"
        );
        assert!(
            before_geo.iter().any(|f| f.is_some()),
            "test setup sanity: witness must have at least one defined E/Z fact"
        );
    }

    /// The issue #390 witness's ambiguous end has exactly one *safe*
    /// candidate: the bond toward atom1=atom2, which is load-bearing for
    /// that unrelated double bond's own geometry. Its sibling (toward
    /// atom4=atom5, itself undefined) is therefore NOT a valid
    /// geometry-preserving alternate carrier -- moving the mark there would
    /// silently strip atom1=atom2's geometry, exactly the corruption
    /// [`Self::is_load_bearing_elsewhere`] exists to forbid.
    /// [`alternate_ez_markings`] independently verifies geometry
    /// preservation before offering an alternate, so it returning none here
    /// is itself a direct, positive confirmation the fix is in effect --
    /// not a vacuous negative.
    #[test]
    fn issue390_witness_has_no_valid_alternate_marking() {
        let mol = parse("O/N=C/C(C=N/O)=N\\NC").unwrap();
        let alternates = alternate_ez_markings(&mol);
        assert!(
            alternates.is_empty(),
            "the only other candidate is load-bearing for atom1=atom2's own geometry; \
             moving the mark there must not be offered as a valid alternate"
        );
    }

    /// Mirror/stereoisomer distinctness: the issue #390 witness's E and Z
    /// forms must remain genuinely distinguishable after canonicalization,
    /// never collapse into the same string. Direct regression test for a
    /// mid-investigation defect (never shipped) where an earlier attempted
    /// fix made canonicalization informationally lossy for this coupled
    /// shape -- both the true-Z and true-E spellings canonicalized to
    /// identical output regardless of which geometry was actually input.
    #[test]
    fn issue390_mirror_stereoisomers_stay_distinct() {
        let z = canonical_smiles(&parse("O/N=C/C(C=N/O)=N\\NC").unwrap());
        let e = canonical_smiles(&parse("O/N=C/C(C=N/O)=N/NC").unwrap());
        assert_ne!(
            z, e,
            "genuinely different E/Z witnesses must not canonicalize to the same string"
        );
    }

    /// Atom-order permutation invariance for the issue #390 witness: every
    /// deterministic relabeling [`ez_carrier_test_variants`] produces
    /// (identity, reversed, 16 seeded Fisher-Yates shuffles, plus any
    /// geometry-preserving alternate markings) must canonicalize to the
    /// identical string, and that string must be idempotent and preserve
    /// the input's own E/Z geometry.
    #[test]
    fn issue390_witness_permutation_invariant() {
        let mol = parse("O/N=C/C(C=N/O)=N\\NC").unwrap();
        let variants = ez_carrier_test_variants(&mol);
        assert!(
            variants.len() >= 18,
            "sanity: at least 16 deterministic + 2 fixed relabelings"
        );

        let outputs: Vec<String> = variants.iter().map(canonical_smiles).collect();
        let unique: HashSet<&String> = outputs.iter().collect();
        assert_eq!(
            unique.len(),
            1,
            "expected ONE canonical output across {} relabelings/markings, got {}: {:?}",
            variants.len(),
            unique.len(),
            unique
        );

        let canon = outputs[0].clone();
        let reparsed = parse(&canon).unwrap_or_else(|e| panic!("reparse of '{canon}': {e}"));
        assert_eq!(
            canonical_smiles(&reparsed),
            canon,
            "canonical(canonical(x)) must equal canonical(x)"
        );

        let before_geo = geometry_fingerprint(&mol);
        let after_geo = geometry_fingerprint(&reparsed);
        assert_eq!(
            before_geo, after_geo,
            "E/Z geometry must be preserved across permutation + canonicalization"
        );
    }

    /// Direct empirical check of the invariant [`Self::resolve_component_
    /// jointly`]'s tie-abstain path relies on: `winning_individualized_
    /// ranks` is a strict total order over EVERY atom, with no ties, even
    /// for a genuinely symmetric molecule (swapping the two coupled ends
    /// is an actual graph automorphism here, not just a superficial
    /// resemblance) -- built via `MoleculeBuilder` directly: methyl-
    /// hydrazone-hydrazone-methyl, `P_A-A(=N_A-N_A2)-B(=N_B-N_B2)-P_B` with
    /// `A`/`B` directly bonded (the shared/coupled candidate) and each end
    /// independently wedge-marked so there is real E/Z content at stake.
    /// The extra terminal `N_A2`/`N_B2` atoms are required, not decorative
    /// -- `compute_stereo_alkene_ends` only counts a double bond as
    /// stereogenic when BOTH ends have >=1 substituent (matching
    /// `chematic_chem::cip::assign_ez`'s own guard), and a bare `C=N` imine
    /// nitrogen has none. If this ever starts reporting equal ranks for
    /// the two automorphic ends, the tie-abstain path documented on
    /// `resolve_component_jointly` becomes reachable in practice, not just
    /// defensively coded -- this test exists to catch that change, not to
    /// assert it can never happen.
    #[test]
    fn individualized_ranks_never_tie_even_for_symmetric_molecules() {
        use chematic_core::{Atom, Element, MoleculeBuilder};
        let mut b = MoleculeBuilder::new();
        let pa = b.add_atom(Atom::new(Element::C)); // P_A methyl
        let a = b.add_atom(Atom::new(Element::C)); // A (coupled alkene end)
        let na = b.add_atom(Atom::new(Element::N)); // N_A
        let na2 = b.add_atom(Atom::new(Element::N)); // N_A's own substituent
        let bb = b.add_atom(Atom::new(Element::C)); // B (coupled alkene end)
        let nb = b.add_atom(Atom::new(Element::N)); // N_B
        let nb2 = b.add_atom(Atom::new(Element::N)); // N_B's own substituent
        let pb = b.add_atom(Atom::new(Element::C)); // P_B methyl
        b.add_bond(pa, a, BondOrder::Up).unwrap();
        b.add_bond(a, na, BondOrder::Double).unwrap();
        b.add_bond(na, na2, BondOrder::Single).unwrap();
        b.add_bond(a, bb, BondOrder::Single).unwrap(); // shared/coupled candidate
        b.add_bond(bb, nb, BondOrder::Double).unwrap();
        b.add_bond(nb, nb2, BondOrder::Single).unwrap();
        b.add_bond(bb, pb, BondOrder::Up).unwrap();
        let mol = b.build();

        let (ranks, _) = winning_individualized_ranks(&mol);
        assert_ne!(
            ranks[a.0 as usize], ranks[bb.0 as usize],
            "individualized ranks are documented as a strict total order \
             even for automorphic atoms -- if this now ties, the \
             tie-abstain path this test's sibling exercises via a forced \
             rank array has become reachable through real input too"
        );

        // Sanity: confirm the two ends really are coupled (share the A-B
        // bond as a mutual candidate) and the solver resolves them (no
        // rank tie means no abstain), rather than this being a vacuous
        // check on a topology that never reaches the tie-break logic.
        let ends = CanonicalWriter::compute_stereo_alkene_ends(&mol);
        let components = CanonicalWriter::coupling_components(&mol, &ends);
        assert!(
            components.iter().any(|c| c.len() == 2),
            "test setup sanity: A and B must form a genuine size-2 coupling \
             component"
        );
        let mut writer = CanonicalWriter::new(&mol, &ranks);
        writer.resolve_ez_markers();
        assert!(
            writer.ez_shared_bond_abstains.is_empty(),
            "test setup sanity: with distinct ranks, the solver should \
             resolve rather than abstain"
        );
    }

    #[test]
    fn issue149_aromatic_stash_matches_exhaustive_oracle() {
        let fixtures = [
            (
                r"CC/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
                [
                    r"c/3(c(/c(c3=N\CC)=N\[C@@H](Cc1ccc(NC(=O)c2c(cncc2Cl)Cl)cc1)C(O)=O)O)O",
                    r"c3(c(c(/c3=N/CC)=N\[C@@H](Cc1ccc(NC(=O)c2c(cncc2Cl)Cl)cc1)C(O)=O)O)O",
                ],
            ),
            (
                r"CCCC(C)/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
                [
                    r"c/1(c(/c(c1=N\[C@H](C(O)=O)Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)=N\C(C)CCC)O)O",
                    r"c1(c(c(/c1=N/[C@H](C(O)=O)Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)=N\C(C)CCC)O)O",
                ],
            ),
            (
                r"COCC/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
                [
                    r"c/1(O)c(O)/c(=N\[C@@H](Cc3ccc(cc3)NC(=O)c2c(Cl)cncc2Cl)C(=O)O)c1=N\CCOC",
                    r"c1(O)c(O)c(=N/[C@@H](Cc3ccc(cc3)NC(=O)c2c(Cl)cncc2Cl)C(=O)O)\c1=N\CCOC",
                ],
            ),
        ];

        // The exhaustive oracle must agree not only for the corpus spelling,
        // but also for both independently observed aromatic-stash spellings.
        // This proves that orbit pruning has not introduced a second source
        // of disagreement; it deliberately does not choose between the two
        // representation-dependent winners.
        for (input, observed_variants) in fixtures {
            for spelling in std::iter::once(input).chain(observed_variants) {
                let mol = crate::parser::parse(spelling).expect("residual parses");
                let current = canonical_smiles(&mol);
                let oracle = canonical_smiles_exhaustive_oracle(&mol);
                assert_eq!(current, oracle, "search/oracle mismatch for {spelling}");
            }
        }
    }

    /// White-box probe of the tie-abstain path, and a genuine (not assumed)
    /// finding about why it cannot be exercised through a size-2 component
    /// -- the ONLY size ever observed in real molecules (measured, see
    /// [`Self::MAX_JOINT_COMPONENT_SIZE`]'s doc comment).
    ///
    /// Built the SAME symmetric molecule as
    /// [`individualized_ranks_never_tie_even_for_symmetric_molecules`],
    /// then FORCED a rank tie between the two coupled ends via a
    /// hand-crafted `ranks` array (bypassing `winning_individualized_ranks`
    /// entirely, since that function never produces one on its own) --
    /// tried under two different mark configurations (`P_A`/`P_B` both
    /// `Up`, and `Up`/`Down`). In BOTH cases the solver still resolved
    /// uniquely, never abstaining, even with the forced tie. Traced why:
    /// for a size-2 component, a "mixed" combination (one end's chosen
    /// candidate is the shared bond, the other's is its own private one)
    /// ALWAYS conflicts by construction -- the end treating the shared
    /// bond as "chosen" votes a real direction on it, while the end
    /// treating it as "other" always votes to demote it to plain
    /// (`Self::plain_order`), and a direction can never equal "plain". So
    /// only two combinations are ever valid at size 2: BOTH ends use the
    /// shared bond, or BOTH use their own private one -- never a case
    /// where the two ends' deviation bits differ, which is the specific
    /// precondition `is_tie` checks for (`winner_bits[i] != winner_bits[j]`).
    /// A rank tie between the two ends is therefore insufficient on its
    /// own to create genuine ambiguity at size 2, REGARDLESS of which
    /// marks are present -- there is no "swap which one deviates" case to
    /// be ambiguous about, only "both agree" (resolved) or "both would
    /// conflict, so both fall back" (also resolved, uniquely).
    ///
    /// This does NOT prove `is_tie` can never fire -- a genuinely
    /// unreachable-in-practice size-3+ component (never observed in any
    /// corpus scanned) could still expose a real swap ambiguity between
    /// two of its rank-tied members while a third breaks the symmetry
    /// differently; that combination was not constructed here (todo,
    /// tracked as a documented gap rather than silently assumed safe).
    /// What IS established: the defensive `is_tie` check is unreachable
    /// for the only topology this crate's own corpus scans have ever
    /// produced, for two independent, now-verified reasons (ranks don't
    /// tie; and even a forced tie can't create a swap-style ambiguity at
    /// size 2) -- not merely "happens not to trigger today".
    #[test]
    fn size_two_component_forced_rank_tie_still_resolves_uniquely() {
        use chematic_core::{Atom, Element, MoleculeBuilder};
        fn build(pb_order: BondOrder) -> (Molecule, AtomIdx, AtomIdx) {
            let mut b = MoleculeBuilder::new();
            let pa = b.add_atom(Atom::new(Element::C));
            let a = b.add_atom(Atom::new(Element::C));
            let na = b.add_atom(Atom::new(Element::N));
            let na2 = b.add_atom(Atom::new(Element::N));
            let bb = b.add_atom(Atom::new(Element::C));
            let nb = b.add_atom(Atom::new(Element::N));
            let nb2 = b.add_atom(Atom::new(Element::N));
            let pb = b.add_atom(Atom::new(Element::C));
            b.add_bond(pa, a, BondOrder::Up).unwrap();
            b.add_bond(a, na, BondOrder::Double).unwrap();
            b.add_bond(na, na2, BondOrder::Single).unwrap();
            b.add_bond(a, bb, BondOrder::Single).unwrap();
            b.add_bond(bb, nb, BondOrder::Double).unwrap();
            b.add_bond(nb, nb2, BondOrder::Single).unwrap();
            b.add_bond(bb, pb, pb_order).unwrap();
            (b.build(), a, bb)
        }

        for pb_order in [BondOrder::Up, BondOrder::Down] {
            let (mol, a, bb) = build(pb_order);
            let (mut ranks, _) = winning_individualized_ranks(&mol);
            // Force a tie between the two coupled ends, preserving every
            // other atom's real rank -- the minimal, targeted way to
            // inject the one precondition `is_tie` checks for.
            ranks[a.0 as usize] = ranks[bb.0 as usize];

            let mut writer = CanonicalWriter::new(&mol, &ranks);
            writer.resolve_ez_markers();
            assert!(
                writer.ez_shared_bond_abstains.is_empty(),
                "pb_order={pb_order:?}: a forced rank tie at size 2 should \
                 still resolve uniquely (see this test's own doc comment \
                 for why no swap-style ambiguity is possible here) -- an \
                 abstain here would mean the structural argument above no \
                 longer holds and needs re-examining, not just re-asserting"
            );
            assert_no_partial_marker_application(&mol);
        }
    }

    /// [`Self::coupling_components`] deliberately sorts its `HashSet`
    /// iteration on raw `AtomIdx` purely to avoid process-random traversal
    /// order (see that method's own doc comment) -- never as a canonical
    /// tie-break. This proves that choice is sufficient: rebuilding the
    /// SAME set of stereo-alkene ends via several different insertion
    /// orders always yields the SAME partition into components (compared
    /// as sets of atoms, never as the returned `Vec`'s order), across
    /// every fixture in [`EZ_SHARED_CARRIER_FULLY_RESOLVED`] -- i.e.
    /// component MEMBERSHIP is a pure function of molecule structure, never
    /// of `HashSet` build/iteration order.
    #[test]
    fn coupling_components_membership_is_insertion_order_independent() {
        fn normalized_components(
            mol: &Molecule,
            ends: &HashSet<AtomIdx>,
        ) -> Vec<std::collections::BTreeSet<u32>> {
            let mut components: Vec<std::collections::BTreeSet<u32>> =
                CanonicalWriter::coupling_components(mol, ends)
                    .into_iter()
                    .map(|c| c.into_iter().map(|a| a.0).collect())
                    .collect();
            components.sort();
            components
        }

        for &s in EZ_SHARED_CARRIER_FULLY_RESOLVED {
            let mol = parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"));
            let ends = CanonicalWriter::compute_stereo_alkene_ends(&mol);
            let baseline = normalized_components(&mol, &ends);

            let mut as_vec: Vec<AtomIdx> = ends.iter().copied().collect();
            as_vec.sort_by_key(|a| a.0);
            let mut orders: Vec<Vec<AtomIdx>> = vec![as_vec.clone()];
            let mut reversed = as_vec.clone();
            reversed.reverse();
            orders.push(reversed);
            for seed in 0..8u64 {
                let perm = deterministic_permutation(as_vec.len(), seed);
                orders.push(perm.into_iter().map(|i| as_vec[i]).collect());
            }

            for order in orders {
                let rebuilt: HashSet<AtomIdx> = order.into_iter().collect();
                let variant = normalized_components(&mol, &rebuilt);
                assert_eq!(
                    baseline, variant,
                    "'{s}': coupling_components membership must not depend \
                     on HashSet insertion order"
                );
            }
        }
    }

    /// Positive control for [`tetrahedral_fingerprint`]: without this, the
    /// preservation checks above would pass vacuously if the fingerprint
    /// function were broken and always returned e.g. all-`None` or some
    /// other input-insensitive constant. Flipping one real `[C@H]` ->
    /// `[C@@H]` tag in one of the fully-resolved fixtures must change the
    /// fingerprint, proving the function actually reads the input.
    #[test]
    fn tetrahedral_fingerprint_is_sensitive_to_a_real_flip() {
        let original = "CCCCC/N=c1\\c(O)c(O)\\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O";
        let flipped = "CCCCC/N=c1\\c(O)c(O)\\c1=N/[C@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O";
        let fp_original = tetrahedral_fingerprint(&parse(original).unwrap());
        let fp_flipped = tetrahedral_fingerprint(&parse(flipped).unwrap());
        assert_ne!(
            fp_original, fp_flipped,
            "flipping a real [C@H]->[C@@H] tag must change the fingerprint"
        );
    }

    /// Positive control for `geometry_fingerprint`: without this, the
    /// "before == after" check above would pass vacuously if the fingerprint
    /// function were broken and always returned e.g. all-`None` or some
    /// other input-insensitive constant. Flipping one real mark in one of
    /// the residual fixtures (`/N=c1\...` -> `/N=c1/...`) must change the
    /// fingerprint, proving the function actually reads the input.
    #[test]
    fn geometry_fingerprint_is_sensitive_to_a_real_flip() {
        let original = "CCCCC/N=c1\\c(O)c(O)\\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O";
        let flipped = "CCCCC/N=c1/c(O)c(O)\\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O";
        let fp_original = geometry_fingerprint(&parse(original).unwrap());
        let fp_flipped = geometry_fingerprint(&parse(flipped).unwrap());
        assert_ne!(
            fp_original, fp_flipped,
            "geometry_fingerprint must be sensitive to a real E/Z flip, or \
             the no-corruption check above is vacuous"
        );
    }
}

// ── Explicit-vs-implicit hydrogen-count canonicalization invariance ────────
//
// Regression suite for issue #205 (kent-tokyo/renkin PR #65 finding): an
// atom's `hydrogen_count` field is `Some(n)` when it was written with
// bracket notation and `None` when written organic-subset, but these can
// represent the *same* chemical state (n happens to equal what valence
// inference would give anyway). Before this fix, `initial_invariant`
// treated `None` and `Some(n)` as different Morgan-rank seeds, and
// `emit_atom` treated any `Some(_)` as forcing bracket notation regardless
// of whether it was redundant -- so two representations of one molecule
// (e.g. a bracket `[Cl]` substituent vs. organic-subset `Cl`) could
// canonicalize to different strings. Fixed by routing both decisions
// through `implicit_hcount`/`valence_inferred_hcount` instead of the raw
// field. These tests prove the invariant mechanically (structural/data
// comparisons), not by assuming two SMILES "look like" the same molecule.
#[cfg(test)]
mod explicit_implicit_h_invariance {
    use super::*;
    use crate::parser::parse;
    use chematic_core::{Atom, Element, MoleculeBuilder};

    fn atom_multiset(mol: &Molecule) -> Vec<(u8, i8, Option<u16>, bool, u8)> {
        // (atomic_number, charge, isotope, aromatic, effective_h_count) --
        // effective (not raw) H count, since that's the structural fact that
        // must be invariant; two representations legitimately differ in
        // whether hydrogen_count is None or Some(same value).
        let mut v: Vec<_> = mol
            .atoms()
            .map(|(idx, a)| {
                (
                    a.element.atomic_number(),
                    a.charge,
                    a.isotope,
                    a.aromatic,
                    chematic_core::implicit_hcount(mol, idx),
                )
            })
            .collect();
        v.sort();
        v
    }

    fn bond_multiset(mol: &Molecule) -> Vec<(u8, u8, u8)> {
        // (min_atomic_number, max_atomic_number, bond_order) at each edge,
        // order-independent identification of the bond multiset.
        let mut v: Vec<_> = mol
            .bonds()
            .map(|(_, b)| {
                let e1 = mol.atom(b.atom1).element.atomic_number();
                let e2 = mol.atom(b.atom2).element.atomic_number();
                (e1.min(e2), e1.max(e2), bond_order_value(b.order) as u8)
            })
            .collect();
        v.sort();
        v
    }

    /// Asserts `a` and `b` are the same chemical graph (same atom multiset,
    /// same bond multiset) *before* checking canonical string equality --
    /// this is what makes the string-equality assertions below a proof of
    /// canonicalization invariance rather than an assumption that two
    /// differently-written SMILES denote the same molecule.
    fn assert_same_graph_then_same_canonical(a: &Molecule, b: &Molecule, label: &str) {
        assert_eq!(
            a.atom_count(),
            b.atom_count(),
            "{label}: atom count must match for this to be a meaningful test"
        );
        assert_eq!(
            atom_multiset(a),
            atom_multiset(b),
            "{label}: atom multisets (element/charge/isotope/aromatic/effective-H) differ -- \
             these are not actually the same molecule, so canonical equality would prove nothing"
        );
        assert_eq!(
            bond_multiset(a),
            bond_multiset(b),
            "{label}: bond multisets differ -- not the same graph"
        );
        assert_eq!(
            canonical_smiles(a),
            canonical_smiles(b),
            "{label}: same graph must canonicalize identically"
        );
    }

    /// Isolated version of the kent-tokyo/renkin PR #65 fixture with zero
    /// reaction machinery involved -- pure SMILES parsing, bracket vs
    /// organic notation. `chematic-rxn`'s own test suite
    /// (`transform.rs::reaction_derived_matches_direct_parse_chlorobenzene`)
    /// covers the actual `run_reactants`-derived case, since chematic-smiles
    /// cannot depend on chematic-rxn.
    #[test]
    fn bracket_vs_organic_notation_isolated_no_reaction() {
        let cases = [
            ("[Cl]c1ccccc1", "Clc1ccccc1"), // aromatic, symmetric ring
            ("[Cl]CCC", "ClCCC"),           // aliphatic, asymmetric
            ("[OH]CC", "OCC"),              // hydroxyl group
            ("C[NH2]", "CN"),               // amine
            ("[CH3]C", "CC"),               // methyl written explicit
        ];
        for (bracket, organic) in cases {
            let a = parse(bracket).unwrap();
            let b = parse(organic).unwrap();
            assert_same_graph_then_same_canonical(&a, &b, &format!("{bracket} vs {organic}"));
        }
    }

    /// Randomized atom relabeling: build the same graph via `MoleculeBuilder`
    /// with atoms inserted in several different orders (including a
    /// deterministic pseudo-random shuffle), holding chemistry fixed. Proves
    /// insertion order alone -- independent of bracket/organic notation --
    /// cannot change `canonical_smiles`'s output.
    #[test]
    fn randomized_atom_relabeling_does_not_change_canonical_form() {
        // A 5-membered ring with a substituent breaking full symmetry:
        // cyclopentane with one CH2 replaced conceptually by using a
        // pendant Cl so every construction order is exercised meaningfully.
        fn ring_with_substituent(insertion_order: &[usize]) -> Molecule {
            // Logical atoms: 0..5 ring carbons, 5 = Cl substituent on ring atom 0.
            let atoms: Vec<Atom> = (0..5)
                .map(|_| Atom {
                    element: Element::C,
                    isotope: None,
                    charge: 0,
                    hydrogen_count: None,
                    aromatic: false,
                    chirality: Chirality::None,
                    wildcard: false,
                    atom_map: None,
                    cip_code: None,
                })
                .chain(std::iter::once(Atom {
                    element: Element::CL,
                    isotope: None,
                    charge: 0,
                    hydrogen_count: None,
                    aromatic: false,
                    chirality: Chirality::None,
                    wildcard: false,
                    atom_map: None,
                    cip_code: None,
                }))
                .collect();
            // Logical bonds: ring 0-1-2-3-4-0, plus 0-5 (Cl).
            let bonds = [
                (0usize, 1usize, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Single),
                (3, 4, BondOrder::Single),
                (4, 0, BondOrder::Single),
                (0, 5, BondOrder::Single),
            ];

            let mut logical_to_new = [0u32; 6];
            let mut b = MoleculeBuilder::new();
            for (new_idx, &logical_idx) in insertion_order.iter().enumerate() {
                logical_to_new[logical_idx] = new_idx as u32;
                b.add_atom(atoms[logical_idx].clone());
            }
            for (l1, l2, order) in bonds {
                let a1 = AtomIdx(logical_to_new[l1]);
                let a2 = AtomIdx(logical_to_new[l2]);
                let _ = b.add_bond(a1, a2, order);
            }
            b.build()
        }

        let orders: Vec<Vec<usize>> = vec![
            vec![0, 1, 2, 3, 4, 5],
            vec![5, 4, 3, 2, 1, 0],
            vec![2, 0, 4, 1, 5, 3],
            vec![5, 0, 1, 2, 3, 4],
            vec![3, 1, 5, 0, 4, 2],
        ];
        let reference = canonical_smiles(&ring_with_substituent(&orders[0]));
        for order in &orders[1..] {
            let mol = ring_with_substituent(order);
            assert_eq!(
                canonical_smiles(&mol),
                reference,
                "insertion order {order:?} must not change canonical_smiles"
            );
        }
    }

    /// Disconnected structures: a salt/mixture where one fragment uses
    /// bracket notation and the other uses organic-subset notation for
    /// what is, on the relevant atom, the same effective H count.
    #[test]
    fn disconnected_structures_bracket_organic_mix() {
        let a = parse("[Cl]CC.[Na+]").unwrap();
        let b = parse("ClCC.[Na+]").unwrap();
        assert_same_graph_then_same_canonical(&a, &b, "disconnected salt, bracket vs organic Cl");
    }

    /// A Kekulized ring (explicit alternating single/double bonds, no
    /// aromaticity perception applied), combined with a bracket/organic
    /// notation difference on its substituent, must still canonicalize
    /// identically. Kept separate from the aromatic-flagged case in
    /// `bracket_vs_organic_notation_isolated_no_reaction` -- raw `parse()`
    /// does not itself perceive aromaticity from Kekulized input (that is a
    /// separate, opt-in `chematic-perception` step), so an aromatic-flagged
    /// molecule and a Kekulized one are genuinely different `Atom.aromatic`
    /// states, not merely different spellings of the same data; comparing
    /// them here would test aromaticity perception, not this fix.
    #[test]
    fn kekulized_ring_bracket_vs_organic_substituent() {
        let bracket = parse("[Cl]C1=CC=CC=C1").unwrap();
        let organic = parse("ClC1=CC=CC=C1").unwrap();
        assert_same_graph_then_same_canonical(
            &bracket,
            &organic,
            "Kekulized ring, bracket vs organic Cl substituent",
        );
    }

    /// Isotope and charge must NOT be unified away -- only redundant
    /// explicit H-count information is. An explicit isotope or nonzero
    /// charge is always genuinely distinguishing (there is no "implicit
    /// isotope"/"implicit charge" to fall back to), so these must still
    /// force bracket notation and distinct canonical forms from their
    /// natural-abundance/neutral counterparts.
    #[test]
    fn isotope_and_charge_remain_distinguishing_not_unified() {
        let normal = parse("CCl").unwrap();
        let isotopic = parse("C[37Cl]").unwrap();
        assert_ne!(
            canonical_smiles(&normal),
            canonical_smiles(&isotopic),
            "an explicit isotope must remain distinguishing"
        );

        let neutral = parse("CC(=O)[O-]").unwrap();
        let anion_str = canonical_smiles(&neutral);
        assert!(
            anion_str.contains('-'),
            "explicit negative charge must remain visible: {anion_str}"
        );
    }

    /// A real stereocenter (explicit chirality) is unaffected by the
    /// explicit/implicit-H unification: canonicalizing a chiral molecule
    /// from two different atom orderings must still agree, and the
    /// canonical form must still be stable under repeated canonicalization
    /// (idempotence).
    #[test]
    fn stereocenter_canonicalization_idempotent_and_order_independent() {
        let a = parse("N[C@@H](C)C(=O)O").unwrap(); // L-alanine
        let b = parse("OC(=O)[C@H](C)N").unwrap(); // same molecule, other end first
        assert_same_graph_then_same_canonical(&a, &b, "L-alanine, two parse orders");

        let canon = canonical_smiles(&a);
        let reparsed = parse(&canon).unwrap();
        assert_eq!(
            canonical_smiles(&reparsed),
            canon,
            "canonical_smiles must be idempotent: canonicalize(parse(canonicalize(x))) == canonicalize(x)"
        );
    }

    /// General idempotence check across the bracket/organic pairs above:
    /// canonicalizing an already-canonical string must reproduce it
    /// byte-for-byte, regardless of which side of a bracket/organic pair it
    /// started from.
    #[test]
    fn canonicalization_is_idempotent_for_all_bracket_organic_pairs() {
        for (bracket, organic) in [
            ("[Cl]c1ccccc1", "Clc1ccccc1"),
            ("[Cl]CCC", "ClCCC"),
            ("[OH]CC", "OCC"),
        ] {
            for smi in [bracket, organic] {
                let mol = parse(smi).unwrap();
                let canon = canonical_smiles(&mol);
                let reparsed = parse(&canon).unwrap();
                assert_eq!(
                    canonical_smiles(&reparsed),
                    canon,
                    "not idempotent starting from {smi}: {canon}"
                );
            }
        }
    }

    // ── Round 11: ring-closure explicit-bond-order aromaticity check (#395) ──
    //
    // `write_chain`'s ring-closure marker decision checked only the
    // currently-written atom's own aromaticity, never its ring-closure
    // partner's -- unlike the equivalent tree-edge decision, which correctly
    // checks both endpoints. A bare ring-closure digit between two aromatic
    // atoms is read back as an *aromatic* bond, so a genuinely `Single`
    // ring-closure bond joining two atoms that each individually happen to
    // be aromatic (e.g. a non-aromatic fusion bond between two separately-
    // aromatic ring systems, `c1-2`) silently became aromatic on re-parse
    // whenever the writer omitted the `-` marker. Confirmed via corpus sweep
    // (130/10,000 molecules) and RDKit InChI cross-check.

    #[test]
    fn ring_closure_explicit_single_bond_between_aromatic_atoms_real_world_repro() {
        // Smallest real-corpus repro from issue #395: the final atom closes
        // both ring 1 (implicit aromatic, correct) and ring 2 (explicit
        // `-2`, a non-aromatic fusion bond) at once.
        let smi = "Oc1[nH]c(Br)nc2nnc(Br)c1-2";
        let mol = parse(smi).unwrap();
        let once = canonical_smiles(&mol);
        let reparsed = parse(&once).unwrap();
        let twice = canonical_smiles(&reparsed);
        assert_eq!(
            once, twice,
            "explicit-bond-order ring closure not idempotent: {smi} -> once={once} twice={twice}"
        );
    }

    #[test]
    fn ring_closure_explicit_bond_order_survives_reparse_for_various_orders() {
        // A minimal two-ring system where the fusion bond (ring digit `2`)
        // is deliberately non-aromatic, exercised at each explicit bond
        // order the writer can emit on a ring closure.
        for smi in [
            "c1ccc2c1-c1ccccc-21", // single fusion bond
            "c1ccc2c1=CC=CC2",     // double fusion bond into a non-aromatic ring
        ] {
            let mol = match parse(smi) {
                Ok(m) => m,
                Err(_) => continue, // not every hand-written combination is valid; skip malformed ones
            };
            let once = canonical_smiles(&mol);
            let reparsed = parse(&once).unwrap_or_else(|e| {
                panic!("re-parse of canonical output '{once}' (from {smi}) failed: {e}")
            });
            let twice = canonical_smiles(&reparsed);
            assert_eq!(
                once, twice,
                "not idempotent: {smi} -> once={once} twice={twice}"
            );
        }
    }
}

#[cfg(test)]
mod legacy_boundary_contract_tests {
    use super::*;
    use chematic_core::MoleculeBuilder;

    #[test]
    fn equivalent_atom_api_agrees_with_known_chain_and_ring_symmetries() {
        for (text, a, b, expected) in [
            ("CC", 0, 1, true),
            ("CCC", 0, 2, true),
            ("CCC", 0, 1, false),
            ("c1ccccc1", 0, 3, true),
            ("CC(=O)O", 2, 3, false),
            ("[13CH3]C", 0, 1, false),
            ("[NH3+]CCN", 0, 3, false),
        ] {
            let mol = crate::parse(text).unwrap();
            assert_eq!(
                are_atoms_equivalent(&mol, AtomIdx(a), AtomIdx(b)),
                expected,
                "{text}"
            );
            assert_eq!(
                are_atoms_equivalent(&mol, AtomIdx(b), AtomIdx(a)),
                expected,
                "{text}"
            );
            assert!(!are_atoms_equivalent(&mol, AtomIdx(a), AtomIdx(u32::MAX)));
        }
        let empty = MoleculeBuilder::new().build();
        assert!(!are_atoms_equivalent(&empty, AtomIdx(0), AtomIdx(0)));
        assert_eq!(legacy_branch_count_for_benchmark(&empty), 0);
        assert_eq!(legacy_canonical_smiles_for_benchmark(&empty), "");
        assert_eq!(canonical_smiles_exhaustive_oracle(&empty), "");
        assert_eq!(
            canonical_smiles_exhaustive_oracle_with_ranks(&empty),
            (vec![], String::new())
        );
    }

    #[test]
    fn capped_legacy_and_tied_rank_writers_preserve_graph_and_stereo() {
        for text in [
            "CCO",
            "CC(C)C",
            "[13CH3]C[O-]",
            "C[NH2+]CC[O-]",
            "c1ccccc1",
            "C1CC1",
            "F[C@H](Cl)Br",
            "F/C=C/F",
            "[Mo]$[Mo]",
            "[Co@OH1](F)(Cl)(Br)(I)(N)O",
        ] {
            let mol = crate::parse(text).unwrap_or_else(|e| panic!("{text}: {e}"));
            let expected = canonical_smiles(&mol);
            assert!(legacy_branch_count_for_benchmark(&mol) > 0, "{text}");
            let ties = vec![0; mol.atom_count()];
            let oversized = vec![mol.atom_count() as u64; mol.atom_count()];
            for written in [
                legacy_canonical_smiles_for_benchmark(&mol),
                CanonicalWriter::new(&mol, &ties).write_all(),
                CanonicalWriter::new(&mol, &oversized).write_all(),
            ] {
                let restored =
                    crate::parse(&written).unwrap_or_else(|e| panic!("{text} -> {written}: {e}"));
                assert_eq!(restored.atom_count(), mol.atom_count(), "{text}");
                assert_eq!(restored.bond_count(), mol.bond_count(), "{text}");
                assert_eq!(canonical_smiles(&restored), expected, "{text} -> {written}");
            }
        }
    }
}
