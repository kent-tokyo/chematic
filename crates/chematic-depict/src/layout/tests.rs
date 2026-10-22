use super::*;

/// Positive control for the group_ring_systems/place_chains determinism
/// fix: repeated calls in the SAME process on the SAME molecule must
/// produce bit-identical coordinates. This is a stricter, different bug
/// class than input-order dependence -- `group_ring_systems` (a plain
/// `HashMap`) and `place_chains` (iterating a `HashSet` without sorting)
/// used Rust's randomly-seeded default hasher, so a fresh `RandomState`
/// derived on every `HashMap`/`HashSet` construction could reorder
/// output even across calls with IDENTICAL input, within one run of one
/// binary -- undetectable by any test that only varies input spelling.
#[test]
fn compute_layout_is_deterministic_across_repeated_calls() {
    use chematic_smiles::parse;
    // Two separate ring systems joined by a chain, so both
    // group_ring_systems (2 union-find roots) and place_chains (a
    // multi-atom unplaced chain) have real tie material to shuffle.
    let mol = parse("c1ccccc1CCCCCc1ccccc1").unwrap();

    let first = compute_layout(&mol);
    for _ in 0..100 {
        let repeat = compute_layout(&mol);
        assert_eq!(
            repeat.coords, first.coords,
            "compute_layout produced different coordinates for the same molecule \
                 across repeated calls in one process"
        );
    }
}

// --- Ring-system-coincidence regression tests -----------------------
//
// Root cause: `place_regular_ring` always centered a new ring at the
// literal origin with no awareness of already-placed geometry, and was
// invoked unconditionally for every ring system -- so any two ring
// systems not fused/bridged into the same connected group (a chain- or
// bond-mediated substituent, or a spiro junction) landed on
// bit-for-bit identical coordinates. Fixed via `place_first_ring_anchored`
// plus a connectivity-driven growth pass (`grow_layout`) that discovers
// and anchors every ring system to its real attachment point.

/// (min non-bonded pairwise distance, min bonded distance, max bonded
/// distance) across every atom pair in `mol`'s `layout`.
fn layout_distance_summary(mol: &Molecule, layout: &Layout) -> (f64, f64, f64) {
    let n = mol.atom_count();
    let mut min_non_bonded = f64::MAX;
    let mut min_bonded = f64::MAX;
    let mut max_bonded = f64::MIN;
    for i in 0..n {
        for j in (i + 1)..n {
            let a = AtomIdx(i as u32);
            let b = AtomIdx(j as u32);
            let d = layout.get(a).dist(&layout.get(b));
            if mol.bond_between(a, b).is_some() {
                min_bonded = min_bonded.min(d);
                max_bonded = max_bonded.max(d);
            } else {
                min_non_bonded = min_non_bonded.min(d);
            }
        }
    }
    (min_non_bonded, min_bonded, max_bonded)
}

/// Asserts no exact/near coincidence (Tier A/B) and exact bond-length
/// fidelity (Tier C) -- the full set, for fixtures with no pre-existing
/// unrelated geometry bug to work around.
fn assert_layout_clean(smiles: &str) {
    use chematic_smiles::parse;
    let mol = parse(smiles).unwrap();
    let layout = compute_layout(&mol);
    let (min_non_bonded, min_bonded, max_bonded) = layout_distance_summary(&mol, &layout);
    assert!(
        min_non_bonded > BOND_LEN / 2.0,
        "{smiles}: non-bonded atoms too close (near-collision), min_non_bonded={min_non_bonded}"
    );
    assert!(
        (min_bonded - BOND_LEN).abs() < 1e-6,
        "{smiles}: bonded distance should equal BOND_LEN, min_bonded={min_bonded}"
    );
    assert!(
        (max_bonded - BOND_LEN).abs() < 1e-6,
        "{smiles}: bonded distance should equal BOND_LEN, max_bonded={max_bonded}"
    );
}

#[test]
fn ring_systems_joined_by_chain_do_not_collide() {
    // Simplest isolation of the bug: two benzene rings joined by a
    // plain chain, no shared atom, no direct bond. This SMILES is also
    // the determinism-test fixture above, which only ever checked
    // repeat-call stability -- it passed on top of the broken layout
    // for a long time.
    assert_layout_clean("c1ccccc1CCCCCc1ccccc1");
}

#[test]
fn ring_systems_joined_directly_do_not_collide() {
    // Direct-bond-mediated attachment (no intervening chain atoms) --
    // exercises the outer-worklist ring-discovery path in
    // `grow_layout`, distinct from the mid-chain discovery inside
    // `dfs_zigzag` the previous test exercises.
    assert_layout_clean("c1ccc(cc1)-c1ccccc1CC");
}

#[test]
fn three_substituent_rings_on_bridged_core_do_not_collide() {
    // The originally reported molecule: three separate phenyl-ring
    // substituents on a bridged bicyclic core, all landing on exactly
    // the same coordinates pre-fix (18 of 27 near-neighbor pairs at
    // distance 0.0). The bridged core shares a three-atom path between
    // rings, so placement must use all already-known ring atoms to choose
    // the correct side of the anchored regular polygon.
    use chematic_smiles::parse;
    let mol = parse("C1CC2CN(CC1N2c1ccccc1)c1cccc(c1)-c1ccccc1").unwrap();
    let layout = compute_layout(&mol);
    let (min_non_bonded, min_bonded, max_bonded) = layout_distance_summary(&mol, &layout);
    assert!(
        min_non_bonded > 1e-6,
        "no two atoms of unrelated ring systems should land on identical coordinates: \
             min_non_bonded={min_non_bonded}"
    );
    assert!(
        min_bonded > 30.0 && max_bonded < 50.0,
        "bridged-core bonds should stay near BOND_LEN: min={min_bonded}, max={max_bonded}; bonds={:?}",
        mol.bonds()
            .map(|(_, b)| (
                b.atom1.0,
                b.atom2.0,
                layout.get(b.atom1).dist(&layout.get(b.atom2))
            ))
            .collect::<Vec<_>>()
    );
}

#[test]
fn pure_chain_layout_unaffected() {
    assert_layout_clean("CCCCCCCC");
}

// --- Issue #347 regressions ------------------------------------------

#[test]
fn long_chain_does_not_wrap_onto_itself() {
    // Pre-fix, dfs_zigzag applied a constant -30°/bond drift instead of
    // alternating, so 12 consecutive bonds traced a full circle: a plain
    // 13-carbon chain's first and last atoms landed on identical
    // coordinates. 20 carbons is well past that wrap point -- a true
    // zigzag (ping-ponging between two directions) never revisits a
    // point no matter how long the chain runs, so this can't pass by
    // accident the way a shorter fixture might.
    assert_layout_clean(&"C".repeat(20));
}

#[test]
fn long_chain_bond_directions_genuinely_alternate() {
    // Direct angle check, not just non-collision: asserts the actual
    // zigzag pattern (ping-ponging between exactly two directions 30°
    // apart), not merely that nothing happened to collide.
    use chematic_smiles::parse;
    let mol = parse(&"C".repeat(10)).unwrap();
    let layout = compute_layout(&mol);
    let n = mol.atom_count();
    let dirs: Vec<f64> = (0..n - 1)
        .map(|i| {
            let a = layout.get(AtomIdx(i as u32));
            let b = layout.get(AtomIdx((i + 1) as u32));
            (b.y - a.y).atan2(b.x - a.x)
        })
        .collect();
    // Consecutive bond directions must differ by exactly 30° in
    // magnitude (a real zigzag turn), never 0° (collinear) or drifting
    // to some other value.
    for w in dirs.windows(2) {
        let mut turn = (w[1] - w[0]).abs();
        if turn > std::f64::consts::PI {
            turn = 2.0 * std::f64::consts::PI - turn;
        }
        assert!(
            (turn - std::f64::consts::PI / 6.0).abs() < 1e-6,
            "expected a 30° zigzag turn between consecutive bonds, got {turn} rad \
                 (dirs={dirs:?})"
        );
    }
    // And the pattern must actually alternate (ping-pong), not drift:
    // only 2 distinct directions should appear across the whole chain.
    let mut distinct: Vec<f64> = Vec::new();
    for &d in &dirs {
        if !distinct.iter().any(|&e: &f64| (e - d).abs() < 1e-6) {
            distinct.push(d);
        }
    }
    assert_eq!(
        distinct.len(),
        2,
        "expected exactly 2 distinct bond directions (ping-pong zigzag), got {distinct:?}"
    );
}

#[test]
fn ring_substituent_direction_bisects_the_open_angle() {
    // Issue #347 repro 2's exact scenario: a hexagon-ring atom with two
    // ring bonds at -30°/90° (120° apart, as any regular-hexagon vertex
    // is). The correct outward direction for an exocyclic substituent is
    // the bisector of the open 240° angle: 210°. Pre-fix,
    // best_outgoing_direction's coarse 60°-grid + broken tie-break
    // returned 180° here -- a real ~30° miss.
    let used_angles = [-std::f64::consts::PI / 6.0, std::f64::consts::PI / 2.0];
    let dir = best_direction_avoiding(&used_angles);
    let expected = 7.0 * std::f64::consts::PI / 6.0; // 210°
    let mut diff = (dir - expected).abs();
    if diff > std::f64::consts::PI {
        diff = 2.0 * std::f64::consts::PI - diff;
    }
    assert!(
        diff < 1e-6,
        "expected the 210° bisector, got {} rad ({} deg)",
        dir,
        dir.to_degrees()
    );
}

#[test]
fn ring_plus_exocyclic_branch_layout_unaffected() {
    // Full end-to-end version of issue #347 repro 2.
    assert_layout_clean("C1CCCC(C(=O)CC)C1");
}

#[test]
fn single_ring_layout_unaffected() {
    assert_layout_clean("c1ccccc1");
}

#[test]
fn fused_and_spiro_ring_systems_unaffected() {
    // Spiro exercises `place_ring_system`'s `already_placed.len() == 1`
    // fix directly; naphthalene/decalin (fused/bridged, 2-atom shared
    // edge) confirm `place_ring_anchored`'s existing, untouched path
    // still works after `place_ring_system`'s anchor-selection refactor.
    for smiles in ["C1CCC2(CC1)CCCCC2", "c1ccc2ccccc2c1", "C1CCC2CCCCC2C1"] {
        assert_layout_clean(smiles);
    }
}

#[test]
fn test_ring_radius_hexagon() {
    // Regular hexagon: all sides == BOND_LEN, circumradius == BOND_LEN.
    let r = ring_radius(6);
    let expected = BOND_LEN; // sin(PI/6) = 0.5, so BOND_LEN / (2*0.5) = BOND_LEN
    assert!((r - expected).abs() < 1e-9, "hexagon radius = {}", r);
}

// --- suggest_bond_direction tests ---

fn make_layout(coords: &[(f64, f64)]) -> Layout {
    Layout {
        coords: coords.iter().map(|&(x, y)| Point { x, y }).collect(),
    }
}

#[test]
fn test_suggest_direction_no_neighbors() {
    use chematic_smiles::parse;
    // Single atom: no neighbors → default 0.0 (pointing right).
    let mol = parse("C").unwrap();
    let layout = make_layout(&[(0.0, 0.0)]);
    let dir = suggest_bond_direction(&mol, AtomIdx(0), &layout);
    assert!((dir).abs() < 1e-9 || (dir - 2.0 * std::f64::consts::PI).abs() < 1e-9);
}

#[test]
fn test_suggest_direction_single_bond_avoids_existing() {
    use chematic_smiles::parse;
    // Ethane C-C: atom 0 at origin, atom 1 to the right (angle = 0°).
    // Suggested direction for a third atom from atom 0 should be far from 0°.
    let mol = parse("CC").unwrap();
    let layout = make_layout(&[(0.0, 0.0), (BOND_LEN, 0.0)]);
    let dir = suggest_bond_direction(&mol, AtomIdx(0), &layout);
    // Must be at least 90° away from 0° (the existing bond).
    let sep = min_angle_separation(dir, &[0.0_f64]);
    assert!(
        sep >= std::f64::consts::PI / 2.0,
        "suggested direction {dir:.3} should be ≥90° from existing bond, sep={sep:.3}"
    );
}

#[test]
fn test_suggest_direction_two_bonds_finds_gap() {
    use chematic_smiles::parse;
    // Three atoms: center C bonded to left (180°) and right (0°).
    // Suggested direction should be ≈ 90° or ≈ -90° (the open gap above/below).
    let mol = parse("CCC").unwrap();
    // atom 1 (center) has neighbors at atom 0 (left) and atom 2 (right).
    let layout = make_layout(&[(-BOND_LEN, 0.0), (0.0, 0.0), (BOND_LEN, 0.0)]);
    let dir = suggest_bond_direction(&mol, AtomIdx(1), &layout);
    // The gap is ≈ 90° (top) or ≈ 270° (bottom). Both have min-sep ≈ 90° from 0° and 180°.
    let sep = {
        let used = [0.0_f64, std::f64::consts::PI];
        min_angle_separation(dir, &used)
    };
    assert!(
        sep >= std::f64::consts::PI / 2.0 - 1e-6,
        "center atom: suggested direction {dir:.3} should be ~90° from both bonds, sep={sep:.3}"
    );
}

#[test]
fn test_suggest_direction_prefers_sp2_for_aromatic_ring() {
    use chematic_smiles::parse;
    // Benzene: use compute_layout to get real coordinates, then ask for
    // the exit direction from atom 0. It should be ≈120° from both ring bonds.
    let mol = parse("c1ccccc1").unwrap();
    let layout = compute_layout(&mol);
    let dir = suggest_bond_direction(&mol, AtomIdx(0), &layout);
    // Atom 0 has 2 ring neighbors. The best exit is ~120° from both.
    let p0 = layout.get(AtomIdx(0));
    let used: Vec<f64> = mol
        .neighbors(AtomIdx(0))
        .map(|(nb, _)| {
            let p = layout.get(nb);
            (p.y - p0.y).atan2(p.x - p0.x)
        })
        .collect();
    let sep = min_angle_separation(dir, &used);
    // Minimum separation from both ring bonds should be ≥ 60°.
    assert!(
        sep >= std::f64::consts::PI / 3.0 - 1e-6,
        "benzene exit direction should be ≥60° from ring bonds, sep={sep:.3}"
    );
}

/// Batch 14: bridged cages, quaternary ring atoms and crowded branches
/// are drawn with no atom within 0.4 bond lengths of a non-bonded atom,
/// and the uncaged ones with no crossing bond (each clashed or crossed
/// before).
#[test]
fn chain_forks_are_drawn_120_degrees_apart() {
    // A branch point's three bonds were drawn 150°/150°/60° (the ±30°
    // zigzag split): an ester's C=O against its C-C.
    use chematic_smiles::parse;
    for (smi, atom) in [("CCOC(=O)C", 3u32), ("CC(C)CC", 1), ("C=C(C)c1ccccc1", 1)] {
        let mol = parse(smi).unwrap();
        let layout = compute_layout(&mol);
        let c = layout.get(AtomIdx(atom));
        let mut angles: Vec<f64> = mol
            .neighbors(AtomIdx(atom))
            .map(|(nb, _)| {
                let p = layout.get(nb);
                (p.y - c.y)
                    .atan2(p.x - c.x)
                    .rem_euclid(std::f64::consts::TAU)
            })
            .collect();
        angles.sort_by(f64::total_cmp);
        for k in 0..3 {
            let gap = (angles[(k + 1) % 3] - angles[k]).rem_euclid(std::f64::consts::TAU);
            assert!((gap.to_degrees() - 120.0).abs() < 1e-6, "{smi}: {angles:?}");
        }
    }
}

#[test]
fn angle_separation_wraps_past_two_pi() {
    // 330° against -90° is 60° apart, not negative.
    let sep = min_angle_separation(330f64.to_radians(), &[(-90f64).to_radians()]);
    assert!(
        (sep.to_degrees() - 60.0).abs() < 1e-9,
        "{}",
        sep.to_degrees()
    );
}

#[test]
fn large_components_open_their_forks_too() {
    // A component of more than 60 atoms that needs narrow forks: its
    // ±45° drawing is relieved as well, and keeps no fork under 90°.
    use chematic_smiles::parse;
    let smi = "CC(C)Cc1cccc(-c2ccccc2C(C)C)c1OC1CC(CNC(=O)c2ccc(C=C3SC(=O)NC3=O)cc2)N(C(=O)c2ccccc2C(=O)c2ccc(F)cc2F)C1";
    let mol = parse(smi).unwrap();
    assert!(mol.atom_count() > 60);
    let layout = compute_layout(&mol);
    let ring_atoms: std::collections::HashSet<AtomIdx> = chematic_perception::find_sssr(&mol)
        .rings()
        .iter()
        .flatten()
        .copied()
        .collect();
    let mut narrow = 0;
    for i in 0..mol.atom_count() {
        let atom = AtomIdx(i as u32);
        if ring_atoms.contains(&atom) || mol.degree(atom) != 3 {
            continue;
        }
        let c = layout.get(atom);
        let mut angles: Vec<f64> = mol
            .neighbors(atom)
            .map(|(nb, _)| {
                let p = layout.get(nb);
                (p.y - c.y)
                    .atan2(p.x - c.x)
                    .rem_euclid(std::f64::consts::TAU)
            })
            .collect();
        angles.sort_by(f64::total_cmp);
        if (0..3).any(|k| {
            (angles[(k + 1) % 3] - angles[k]).rem_euclid(std::f64::consts::TAU)
                < 90f64.to_radians() - 1e-6
        }) {
            narrow += 1;
        }
    }
    assert_eq!(narrow, 0);
    assert!(detect_crossings(&layout, &mol).is_empty());
}

#[test]
fn crowded_forks_open_to_90_degrees_when_that_does_as_well() {
    // The carbonyl of a trityl phenyl ketone needs a narrower fork than
    // 120°; at ±30° its O ran 60° from a C-C, at ±45° 90°.
    use chematic_smiles::parse;
    let smi = "O=C(C1=CC=CC=C1)C(C2=CC=CC=C2)(C3=CC=CC=C3)C4=CC=CC=C4";
    let mol = parse(smi).unwrap();
    let layout = compute_layout(&mol);
    let c = layout.get(AtomIdx(1));
    let mut angles: Vec<f64> = mol
        .neighbors(AtomIdx(1))
        .map(|(nb, _)| {
            let p = layout.get(nb);
            (p.y - c.y)
                .atan2(p.x - c.x)
                .rem_euclid(std::f64::consts::TAU)
        })
        .collect();
    angles.sort_by(f64::total_cmp);
    for k in 0..3 {
        let gap = (angles[(k + 1) % 3] - angles[k]).rem_euclid(std::f64::consts::TAU);
        assert!(gap.to_degrees() > 90.0 - 1e-6, "{angles:?}");
    }
    assert!(detect_crossings(&layout, &mol).is_empty());
}

#[test]
fn bridged_and_crowded_layouts_have_no_clash_or_crossing() {
    use chematic_smiles::parse;
    for (smi, cage) in [
        ("CNC1CC2CCC1C2", true),        // norbornane
        ("C1C[S+]2CC[S+]1CC2", true),   // bicyclo[2.2.2]
        ("C1OCC2(CN3CCC2CC3)O1", true), // spiro quinuclidine
        // Tropane, a pinanol: the N-methyl and the bridge's methyls
        // swing into free faces instead of across a ring bond.
        ("CN1C2CCC1CC(O)C2", false),
        ("CN1C2CCC1CC(OC(=O)c1c[nH]c3ccccc13)C2", false),
        ("CC1(C)C2CCC(C)(O)C1C2", false),
        ("c1ccc(C(n2ccnc2)n2ccnc2)cc1", false), // two rings on one carbon
        ("CCC1(CC)C(=O)NC(=O)NC1=O", false),    // quaternary ring carbon
        ("OC(=O)C(F)(F)C(F)(F)C(O)=O", false),  // crowded chain
        ("c1ccc(C2(N3CCCCC3)CCCCC2)cc1", false),
        // A carboxyl stacked on a Kekulé-written ring: turned about the
        // ring atom, whose ring double bond sets no E/Z.
        (
            "OC(=O)C1=C(C=CC=C1)C(C2=CC=C(O)C(=C2)S(O)(=O)=O)=C3C=CC(=O)C(=C3)S(O)(=O)=O",
            false,
        ),
        // Aziridines on a cyclophosphazene's P (ring P=N).
        (
            "C1CCN(P2(N3CCCC3)=NP(N3CC3)(N3CC3)=NP(N3CC3)(N3CC3)=N2)C1",
            false,
        ),
        // Fused and spiro systems of three or more rings that overlap
        // as polygons (cephalotaxine core, spiro naphthodioxin).
        ("COC1=CC23CCCN2CCc2cc4c(cc2C3C1O)OCO4", false),
        ("O=C1C(O)=CC2(Oc3cccc4cccc(c34)O2)C23OC12C(O)CCC3O", false),
        // Gem-disubstituted ring entry atoms: phenytoin, a
        // 3,3-disubstituted glutarimide.
        ("O=C1NC(=O)C(c2ccccc2)(c2ccccc2)N1", false),
        ("O=C1CCC(c2ccccc2)(C2CCN(C)CC2)C(=O)N1", false),
        // Adamantane-type cages, drawn as RDKit's projection (one
        // crossing, as RDKit's), substituted on a bridgehead or a link.
        ("C1C2CC3CC1CC(C2)C3", true),
        ("C1N2CN3CN1CN(C2)C3", true),
        ("C=CC[N+]12CN3CN(CN(C3)C1)C2", true),
        ("NC12CC3CC(CC(C3)C1)C2", true),
        ("OC1C2CC3CC(C2)CC1C3", true),
        // Porphyrins (protoporphyrin IX, a zinc porphyrin): the pyrroles
        // were folded inside the macrocycle.
        (
            "CC1=C2NC(=C1CCC(O)=O)C=C3N=C(C=C4NC(=CC5=NC(=C2)C(=C5C)C=C)C(=C4C)C=C)C(=C3CCC(O)=O)C",
            false,
        ),
        (
            "C1=CC2=CC3=CC=C4N3[Zn]35N2C1=CC1=CC=C(N13)C=C1C=CC(=N15)C=4",
            false,
        ),
        // Perfluoroalkyl chains run straight through their CF2 atoms
        // instead of curling back (perfluorotributylamine; a
        // perfluorooctyl catechol).
        (
            "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)N(C(F)(F)C(F)(F)C(F)(F)C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
            false,
        ),
        (
            "Oc1cccc(CCCCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F)c1O",
            false,
        ),
        // An octahedral metal with four pyridines: their rings, 60°
        // apart, are hung one and a half bonds out.
        (
            "N#CS[Ni](SC#N)([N+]1=CC=CC=C1)([N+]2=CC=CC=C2)([N+]3=CC=CC=C3)[N+]4=CC=CC=C4",
            false,
        ),
        // A ring-fusion atom's substituent goes into the exterior gap
        // (it went into a ring when the gap's direction passed 2π).
        ("CC1(C)OC2C3C(COC2(COS(N)(=O)=O)O1)C3(Cl)Cl", false),
        // An aspartate carboxyl turned past 60° round its CH2, clear of
        // the neighbouring residue (three-atom groups take the wide turns).
        (
            "CCCN(NC(=O)C1CCCN1C(=O)C(NC(=O)C(NC(=O)C(CC(=O)O)NC(=O)C(CCC(=O)O)NC(=O)C(NC(=O)C(CC(=O)O)NC(C)=O)C(C)O)C(C)C)C(C)C)C(=O)c1cc(C(F)(F)F)cc(C(F)(F)F)c1",
            false,
        ),
        (
            "CC1(C)CCC2(C(=O)N3CCOCC3)CCC3(C)C(C(=O)C=C4C5(C)C=C(C#N)C(=O)C(C)(C)C5CCC43C)C2C1",
            false,
        ),
    ] {
        let mol = parse(smi).unwrap();
        let layout = compute_layout(&mol);
        let n = mol.atom_count();
        for i in 0..n {
            for j in i + 1..n {
                let (a, b) = (AtomIdx(i as u32), AtomIdx(j as u32));
                if mol.bond_between(a, b).is_none() {
                    let d = layout.get(a).dist(&layout.get(b));
                    assert!(d >= 0.4 * BOND_LEN, "{smi}: atoms {i} and {j} {d:.2} apart");
                }
            }
        }
        // A drawn cage may cross itself (the tropane bridge does).
        assert!(
            cage || detect_crossings(&layout, &mol).is_empty(),
            "{smi}: crossing bonds"
        );
    }
}
