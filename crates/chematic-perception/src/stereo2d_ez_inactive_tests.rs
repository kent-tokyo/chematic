//! Existing pinned inputs and exact pre-change API/classifier oracle.
use super::*;

fn reference_ez(
    mol: &mut Molecule,
    coords: &[(f64, f64)],
    explicitly_unspecified: &HashSet<BondIdx>,
) -> Vec<EzDirectionDiagnostic> {
    let double_bonds: Vec<BondIdx> = mol
        .bonds()
        .filter(|(_, b)| b.order == BondOrder::Double)
        .map(|(bidx, _)| bidx)
        .collect();

    // With no literal double bond, every later phase has no outcomes and
    // cannot write a direction or diagnostic. Avoid the unused ring work.
    if double_bonds.is_empty() {
        return Vec::new();
    }
    // Exclude Kekule aromatic ring bonds with the same non-mutating query.
    let aromaticity = assign_aromaticity(mol);
    // Classify independently; no double bond reads another bond's outcome.
    let mut outcomes = reference_classify_double_bonds(
        mol,
        coords,
        &double_bonds,
        explicitly_unspecified,
        &aromaticity,
    );

    // Phase 2: a physical substituent bond can be the chosen carrier for TWO
    // different double bonds at once (the ordinary conjugated-diene shape,
    // see module docs). Group every Assigned outcome's claims by physical
    // bond and check they agree; any bond with disagreeing claims takes
    // every double bond relying on it down with it, symmetrically -- never
    // "first-claimed wins" (that would make the result depend on
    // `double_bonds` iteration order, i.e. bond-index/parse order, which the
    // spec explicitly rules out).
    let mut claims: HashMap<BondIdx, Vec<(BondIdx, BondOrder)>> = HashMap::new();
    for (&db, outcome) in &outcomes {
        if let EzOutcome::Assigned {
            carrier_a1,
            carrier_a2,
        } = outcome
        {
            claims
                .entry(carrier_a1.0)
                .or_default()
                .push((db, carrier_a1.1));
            claims
                .entry(carrier_a2.0)
                .or_default()
                .push((db, carrier_a2.1));
        }
    }
    // Only the relative direction of a double bond's two carriers encodes
    // its geometry, so flipping both keeps it. Two double bonds that claim
    // the same bond with opposite values are reconciled by flipping one of
    // them (a two-colouring of the claim graph, as RDKit's direction
    // assignment does for conjugated polyenes). Only a claim graph with no
    // consistent colouring (an odd cycle) still rejects every double bond
    // in it.
    let mut parent: HashMap<BondIdx, (BondIdx, bool)> = HashMap::new();
    fn find(parent: &mut HashMap<BondIdx, (BondIdx, bool)>, x: BondIdx) -> (BondIdx, bool) {
        let (p, flip) = *parent.get(&x).unwrap_or(&(x, false));
        if p == x {
            return (x, false);
        }
        let (root, f) = find(parent, p);
        parent.insert(x, (root, flip ^ f));
        (root, flip ^ f)
    }
    // Constraints `flip(x) ^ flip(y) == want`:
    // * two double bonds claiming one bond must leave it one value;
    // * an alkene end whose two substituents both carry a marker (its own
    //   carrier and another double bond's) needs them on opposite sides, as
    //   every SMILES reader requires.
    let up_from = |bond: BondIdx, dir: BondOrder, end: AtomIdx| {
        (mol.bond(bond).atom1 == end) == (dir == BondOrder::Up)
    };
    let mut constraints: Vec<(BondIdx, BondIdx, bool)> = Vec::new();
    let mut claim_list: Vec<&Vec<(BondIdx, BondOrder)>> = claims.values().collect();
    claim_list.sort_by_key(|c| c.iter().map(|&(db, _)| db.0).min());
    for claimants in &claim_list {
        let (db0, v0) = claimants[0];
        for &(db, v) in &claimants[1..] {
            constraints.push((db0, db, v != v0));
        }
    }
    type Carrier = (BondIdx, BondOrder);
    let mut assigned_sorted: Vec<(BondIdx, Carrier, Carrier)> = outcomes
        .iter()
        .filter_map(|(&db, o)| match o {
            EzOutcome::Assigned {
                carrier_a1,
                carrier_a2,
            } => Some((db, *carrier_a1, *carrier_a2)),
            _ => None,
        })
        .collect();
    assigned_sorted.sort_by_key(|&(db, ..)| db.0);
    for &(db, c1, c2) in &assigned_sorted {
        let bond = mol.bond(db);
        for (end, (cx, vx)) in [(bond.atom1, c1), (bond.atom2, c2)] {
            for (_, m) in mol.neighbors(end) {
                if m == db || m == cx {
                    continue;
                }
                for &(other, vy) in claims.get(&m).map(Vec::as_slice).unwrap_or(&[]) {
                    if other != db {
                        let same = up_from(cx, vx, end) == up_from(m, vy, end);
                        constraints.push((db, other, same));
                    }
                }
            }
        }
    }
    let mut inconsistent: HashSet<BondIdx> = HashSet::new();
    for (db0, db, want) in constraints {
        {
            let (r0, f0) = find(&mut parent, db0);
            let (r1, f1) = find(&mut parent, db);
            if r0 == r1 {
                if f0 ^ f1 != want {
                    inconsistent.insert(r0);
                }
            } else {
                let (lo, hi, fl) = if r0.0 < r1.0 {
                    (r0, r1, f0 ^ f1 ^ want)
                } else {
                    (r1, r0, f0 ^ f1 ^ want)
                };
                parent.insert(hi, (lo, fl));
                if inconsistent.remove(&hi) {
                    inconsistent.insert(lo);
                }
            }
        }
    }
    let mut conflicted: HashSet<BondIdx> = HashSet::new();
    let mut flipped: HashSet<BondIdx> = HashSet::new();
    let assigned: Vec<BondIdx> = outcomes
        .iter()
        .filter(|(_, o)| matches!(o, EzOutcome::Assigned { .. }))
        .map(|(&db, _)| db)
        .collect();
    for db in assigned {
        let (root, flip) = find(&mut parent, db);
        if inconsistent.contains(&root) {
            conflicted.insert(db);
        } else if flip {
            flipped.insert(db);
        }
    }

    // Phase 3: commit. Only an Assigned outcome not caught by Phase 2 ever
    // writes to `mol`; everything else contributes at most one diagnostic.
    let mut diagnostics = Vec::with_capacity(double_bonds.len());
    for &bidx in &double_bonds {
        match outcomes.remove(&bidx).expect("classified in phase 1") {
            EzOutcome::NotRequested => {}
            EzOutcome::Rejected(reason) => {
                diagnostics.push(EzDirectionDiagnostic { bond: bidx, reason });
            }
            EzOutcome::Assigned {
                carrier_a1,
                carrier_a2,
            } => {
                if conflicted.contains(&bidx) {
                    diagnostics.push(EzDirectionDiagnostic {
                        bond: bidx,
                        reason: EzDirectionRejectionReason::CarrierConflict,
                    });
                } else {
                    let flip = |d: BondOrder| match (flipped.contains(&bidx), d) {
                        (true, BondOrder::Up) => BondOrder::Down,
                        (true, BondOrder::Down) => BondOrder::Up,
                        (_, d) => d,
                    };
                    mol.set_bond_direction(carrier_a1.0, flip(carrier_a1.1));
                    mol.set_bond_direction(carrier_a2.0, flip(carrier_a2.1));
                }
            }
        }
    }
    diagnostics
}

#[inline]
fn reference_classify_double_bonds(
    mol: &Molecule,
    coords: &[(f64, f64)],
    double_bonds: &[BondIdx],
    explicitly_unspecified: &HashSet<BondIdx>,
    aromaticity: &AromaticityModel,
) -> HashMap<BondIdx, EzOutcome> {
    let mut outcomes = HashMap::with_capacity(double_bonds.len());
    for &bidx in double_bonds {
        let outcome = classify_double_bond(mol, coords, bidx, explicitly_unspecified, aromaticity);
        outcomes.insert(bidx, outcome);
    }
    outcomes
}

fn compare(input: &Molecule, coords: &[(f64, f64)], unspecified: &HashSet<BondIdx>) {
    let mut expected = input.clone();
    let mut actual = input.clone();
    for _ in 0..2 {
        let old = reference_ez(&mut expected, coords, unspecified);
        let new = apply_ez_directions_from_2d_ex(&mut actual, coords, unspecified);
        assert_eq!(old, new);
        assert_eq!(
            expected.atoms().collect::<Vec<_>>(),
            actual.atoms().collect::<Vec<_>>()
        );
        assert_eq!(
            expected.bonds().collect::<Vec<_>>(),
            actual.bonds().collect::<Vec<_>>()
        );
        assert_eq!(expected.stereo_groups(), actual.stereo_groups());
        assert_eq!(expected.rdkit_bond_order(), actual.rdkit_bond_order());
        for (i, _) in expected.atoms() {
            assert_eq!(expected.neighbor_slice(i), actual.neighbor_slice(i));
            assert_eq!(
                expected.stereo_neighbor_order(i),
                actual.stereo_neighbor_order(i)
            );
            assert_eq!(expected.atom_tag(i), actual.atom_tag(i));
            assert_eq!(expected.r_group_label(i), actual.r_group_label(i));
        }
        for (i, _) in expected.bonds() {
            assert_eq!(expected.bond_direction(i), actual.bond_direction(i));
            assert_eq!(
                expected.bond_direction_anchor(i),
                actual.bond_direction_anchor(i)
            );
            assert_eq!(
                expected.smiles_ring_closure_metadata(i),
                actual.smiles_ring_closure_metadata(i)
            );
        }
        for i in 0..expected.stereo_groups().len() {
            assert_eq!(
                expected.stereo_group_bond_indices(i),
                actual.stereo_group_bond_indices(i)
            );
        }
    }
}

#[test]
fn inactive_double_gate_preserves_full_api_on_pinned_inputs() {
    let corpora = [
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"
        )),
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../validation/cosmolkit_comparison/corpora/stress-v1.smi"
        )),
    ];
    let mut rows = 0;
    let mut comparisons = 0;
    for corpus in corpora {
        for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
            let parsed =
                chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
            let kekule = chematic_core::kekulize(&parsed)
                .map(|k| chematic_core::apply_kekule(&parsed, &k))
                .unwrap_or_else(|_| parsed.clone());
            for mut mol in [parsed, kekule] {
                // Preserve annotations too; gate eligibility does not read them.
                if rows % 2 == 1 {
                    let ids = mol
                        .bonds()
                        .filter_map(|(bi, b)| {
                            (bi.0 % 3 == 0 && b.order != BondOrder::Double).then_some(bi)
                        })
                        .collect::<Vec<_>>();
                    for bi in ids {
                        mol.set_bond_direction(bi, BondOrder::Down);
                    }
                    if mol.atom_count() > 0 {
                        mol.set_tag(AtomIdx(0), Some(91));
                    }
                }
                let all = mol
                    .bonds()
                    .filter_map(|(bi, b)| (b.order == BondOrder::Double).then_some(bi))
                    .collect::<HashSet<_>>();
                let even = all
                    .iter()
                    .copied()
                    .filter(|bi| bi.0 % 2 == 0)
                    .collect::<HashSet<_>>();
                let normal = (0..mol.atom_count())
                    .map(|i| (i as f64, (i * i % 7) as f64))
                    .collect::<Vec<_>>();
                for coords in [
                    normal,
                    Vec::new(),
                    vec![(f64::NAN, f64::INFINITY); mol.atom_count()],
                ] {
                    for unspecified in [&HashSet::new(), &all, &even] {
                        compare(&mol, &coords, unspecified);
                        comparisons += 2;
                    }
                }
            }
            rows += 1;
        }
    }
    assert_eq!(rows, 14072);
    assert_eq!(comparisons, 506592);
}
