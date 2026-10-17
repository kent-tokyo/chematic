use super::*;
use crate::canonical::canonical_smiles;
use chematic_core::{AtomIdx, BondIdx};

#[test]
fn bracket_atom_reads_two_digit_charge_and_h_count() {
    let a = parse("[C+10]").unwrap();
    assert_eq!(a.atom(AtomIdx(0)).charge, 10);
    let b = parse("[CH15]").unwrap();
    assert_eq!(b.atom(AtomIdx(0)).hydrogen_count, Some(15));
    let c = parse("[O-12]").unwrap();
    assert_eq!(c.atom(AtomIdx(0)).charge, -12);
    let d = parse("[NH4+]").unwrap();
    assert_eq!(d.atom(AtomIdx(0)).hydrogen_count, Some(4));
    assert_eq!(d.atom(AtomIdx(0)).charge, 1);
}
#[test]
fn test_parse_methane() {
    let mol = parse("C").unwrap();
    assert_eq!(mol.atom_count(), 1);
    assert_eq!(mol.bond_count(), 0);
}

#[test]
fn every_two_digit_ring_label_and_reuse_preserve_graph_and_stereo() {
    let expected = crate::canonical_smiles(&parse("N[C@]1(F)CCCC1").unwrap());
    for label in 0..100 {
        let token = format!("%{label:02}");
        let smiles = format!("N[C@]{token}(F)CCCC{token}");
        let mol = parse(&smiles).unwrap();
        assert_eq!(crate::canonical_smiles(&mol), expected, "{smiles}");
        let repeated = parse(&format!("C{token}CC{token}.C{token}CC{token}")).unwrap();
        assert_eq!(repeated.atom_count(), 6);
        assert_eq!(repeated.bond_count(), 6);
        assert!(matches!(parse(&format!("C{token}CC")),
            Err(SmilesError::UnmatchedRingClosure { ring_num, .. }) if ring_num == label));
    }
    for malformed in ["C%", "C%1", "C%1a", "C%a1", "C1C1"] {
        assert!(parse(malformed).is_err(), "{malformed}");
    }
    // Spill the inline close-partner buffer and reuse labels before a
    // chiral ring opening: slot ids must still resolve the right partner.
    let prefix = "C1CC1.".repeat(20);
    let mol = parse(&format!("{prefix}N[C@]1(F)CCCC1")).unwrap();
    assert_eq!(mol.atom_count(), 67);
    assert_eq!(mol.bond_count(), 67);
    let alternative = parse(&format!("{prefix}N[C@]%99(F)CCCC%99")).unwrap();
    assert_eq!(
        crate::canonical_smiles(&mol),
        crate::canonical_smiles(&alternative)
    );
}

#[test]
fn test_parse_ethane() {
    let mol = parse("CC").unwrap();
    assert_eq!(mol.atom_count(), 2);
    assert_eq!(mol.bond_count(), 1);
}

#[test]
fn test_parse_propane() {
    let mol = parse("CCC").unwrap();
    assert_eq!(mol.atom_count(), 3);
    assert_eq!(mol.bond_count(), 2);
}

#[test]
fn test_parse_isobutane() {
    // CC(C)C — branched structure
    let mol = parse("CC(C)C").unwrap();
    assert_eq!(mol.atom_count(), 4);
    assert_eq!(mol.bond_count(), 3);
}

#[test]
fn test_parse_double_bond() {
    let mol = parse("C=C").unwrap();
    assert_eq!(mol.bond_count(), 1);
    let (_, bond) = mol.bonds().next().unwrap();
    assert_eq!(bond.order, BondOrder::Double);
}

#[test]
fn test_parse_triple_bond() {
    let mol = parse("C#N").unwrap();
    let (_, bond) = mol.bonds().next().unwrap();
    assert_eq!(bond.order, BondOrder::Triple);
}

#[test]
fn test_parse_benzene_kekulized() {
    let mol = parse("C1=CC=CC=C1").unwrap();
    assert_eq!(mol.atom_count(), 6);
    assert_eq!(mol.bond_count(), 6);
}

#[test]
fn test_parse_benzene_aromatic() {
    let mol = parse("c1ccccc1").unwrap();
    assert_eq!(mol.atom_count(), 6);
    assert_eq!(mol.bond_count(), 6);
    for (_, atom) in mol.atoms() {
        assert!(atom.aromatic);
    }
    for (_, bond) in mol.bonds() {
        assert_eq!(bond.order, BondOrder::Aromatic);
    }
}

#[test]
fn test_parse_pyridine() {
    let mol = parse("c1ccncc1").unwrap();
    assert_eq!(mol.atom_count(), 6);
    assert_eq!(mol.bond_count(), 6);
}

#[test]
fn test_parse_naphthalene() {
    let mol = parse("c1ccc2ccccc2c1").unwrap();
    assert_eq!(mol.atom_count(), 10);
    assert_eq!(mol.bond_count(), 11);
}

#[test]
fn test_parse_bracket_water() {
    let mol = parse("[OH2]").unwrap();
    let atom = mol.atom(AtomIdx(0));
    assert_eq!(atom.element, Element::O);
    assert_eq!(atom.hydrogen_count, Some(2));
}

#[test]
fn test_parse_wildcard_atom_map() {
    let mol = parse("[*:17]C").unwrap();
    assert!(mol.atom(AtomIdx(0)).wildcard);
    assert_eq!(mol.atom(AtomIdx(0)).atom_map, Some(17));
}

#[test]
fn wildcard_in_aromatic_ring_has_aromatic_ring_bonds() {
    for smiles in ["c1cc*cc1", "c1cc[*]cc1", "c1cc-*-cc1", "c1cc*(F)cc1"] {
        let mol = parse(smiles).unwrap();
        let (idx, atom) = mol.atoms().find(|(_, atom)| atom.wildcard).unwrap();
        assert!(atom.aromatic, "{smiles}");
        assert_eq!(
            mol.neighbors(idx)
                .filter(|&(_, bond)| mol.bond(bond).order == BondOrder::Aromatic)
                .count(),
            2,
            "{smiles}"
        );
    }
}

#[test]
fn pendant_inter_ring_and_non_huckel_wildcards_remain_nonaromatic() {
    for smiles in [
        "c1cc(*)ccc1",
        "c1cc(*c2ccccc2)ccc1",
        "c1cc*ccc1",
        "c1cc*cccc1",
    ] {
        let mol = parse(smiles).unwrap();
        assert!(
            mol.atoms()
                .filter(|(_, atom)| atom.wildcard)
                .all(|(_, atom)| !atom.aromatic),
            "{smiles}"
        );
    }
}

#[test]
fn bare_wildcard_is_the_bracket_wildcard() {
    for (bare, bracket) in [("C*", "C[*]"), ("*C(*)=O", "[*]C([*])=O"), ("*", "[*]")] {
        let a = parse(bare).unwrap();
        let b = parse(bracket).unwrap();
        assert_eq!(
            crate::canonical_smiles(&a),
            crate::canonical_smiles(&b),
            "{bare}"
        );
        assert!(a.atoms().any(|(_, atom)| atom.wildcard));
    }
}

#[test]
fn wildcard_keeps_charge_isotope_and_h_count() {
    // RDKit keeps these on a dummy atom; they used to be dropped.
    let mol = parse("[13*H2-:3]").unwrap();
    let atom = mol.atom(AtomIdx(0));
    assert!(atom.wildcard);
    assert_eq!(
        (
            atom.isotope,
            atom.charge,
            atom.hydrogen_count,
            atom.atom_map
        ),
        (Some(13), -1, Some(2), Some(3))
    );
    assert_eq!(crate::canonical_smiles(&mol), "[13*H2-:3]");
    assert!(crate::canonical_smiles(&parse("C[*+]").unwrap()).contains("[*+]"));
}

#[test]
fn test_parse_ammonium() {
    let mol = parse("[NH4+]").unwrap();
    let atom = mol.atom(AtomIdx(0));
    assert_eq!(atom.charge, 1);
    assert_eq!(atom.hydrogen_count, Some(4));
}

#[test]
fn test_parse_13c() {
    let mol = parse("[13C]").unwrap();
    let atom = mol.atom(AtomIdx(0));
    assert_eq!(atom.isotope, Some(13));
}

#[test]
fn test_parse_ethanol() {
    let mol = parse("CCO").unwrap();
    assert_eq!(mol.atom_count(), 3);
    assert_eq!(mol.bond_count(), 2);
}

#[test]
fn test_parse_disconnected() {
    // Salt: sodium chloride
    let mol = parse("[Na+].[Cl-]").unwrap();
    assert_eq!(mol.atom_count(), 2);
    assert_eq!(mol.bond_count(), 0);
}

#[test]
fn test_parse_acetic_acid() {
    // CC(=O)O
    let mol = parse("CC(=O)O").unwrap();
    assert_eq!(mol.atom_count(), 4);
    assert_eq!(mol.bond_count(), 3);
}

#[test]
fn test_empty_smiles_error() {
    assert!(matches!(parse(""), Err(SmilesError::EmptyInput)));
}

#[test]
fn test_parse_cyclohexane() {
    let mol = parse("C1CCCCC1").unwrap();
    assert_eq!(mol.atom_count(), 6);
    assert_eq!(mol.bond_count(), 6);
}

#[test]
fn test_parse_percent_ring() {
    // %10 ring closure
    let mol = parse("C%10CCCCC%10").unwrap();
    assert_eq!(mol.atom_count(), 6);
    assert_eq!(mol.bond_count(), 6);
}

#[test]
fn test_parse_extended_smiles_plus_ring_label() {
    let mol = parse("C%(100)CC%(100)").unwrap();
    assert_eq!(mol.atom_count(), 3);
    assert_eq!(mol.bond_count(), 3);

    for malformed in [
        "C%()CC%()",
        "C%(99)CC%(99)",
        "C%(100",
        "C%(abc)CC%(abc)",
        "C%(4294967296)CC%(4294967296)",
        "C%100CC%100",
    ] {
        assert!(parse(malformed).is_err(), "{malformed}");
    }

    assert!(matches!(
        parse("C%(256)"),
        Err(SmilesError::UnmatchedExtendedRingClosure { ring_num: 256, .. })
    ));
    assert!(matches!(
        parse("C/%(256)CC/%(256)"),
        Err(SmilesError::ConflictingExtendedRingBond { ring_num: 256, .. })
    ));
    assert!(matches!(
        parse("C2.C1"),
        Err(SmilesError::UnmatchedRingClosure { ring_num: 1, .. })
    ));
}

#[test]
fn test_parse_chlorobenzene() {
    let mol = parse("c1ccccc1Cl").unwrap();
    assert_eq!(mol.atom_count(), 7);
    assert_eq!(mol.bond_count(), 7);
}

fn nested_branch_smiles(depth: usize) -> String {
    let mut smiles = String::from("C");
    for _ in 0..depth {
        smiles.push_str("(C");
    }
    for _ in 0..depth {
        smiles.push(')');
    }
    smiles
}

#[test]
fn test_branch_depth_limit_accepts_boundary() {
    let smiles = nested_branch_smiles(MAX_BRANCH_DEPTH);
    let mol = parse(&smiles).expect("maximum configured branch depth should parse");
    assert_eq!(mol.atom_count(), MAX_BRANCH_DEPTH + 1);
}

#[test]
fn test_branch_depth_limit_rejects_too_deep_input() {
    let smiles = nested_branch_smiles(MAX_BRANCH_DEPTH + 1);
    assert!(matches!(
        parse(&smiles),
        Err(SmilesError::NestingTooDeep { .. })
    ));
}

#[test]
fn test_malformed_ring_closures_return_errors() {
    assert!(matches!(
        parse("C1"),
        Err(SmilesError::UnmatchedRingClosure { ring_num: 1, .. })
    ));
    assert!(matches!(
        parse("C%"),
        Err(SmilesError::UnexpectedEnd { .. })
    ));
    assert!(matches!(
        parse("C%1"),
        Err(SmilesError::UnexpectedEnd { .. })
    ));
    assert!(matches!(
        parse("C%AA"),
        Err(SmilesError::UnexpectedEnd { .. })
    ));
    assert!(matches!(
        parse("C=1CC-1"),
        Err(SmilesError::ConflictingRingBond { ring_num: 1, .. })
    ));
}

#[test]
fn test_stereo_marker_between_aromatic_atoms_stays_aromatic() {
    // SMILES like `c1\c(O)c(O)\c1=N` use `/`/`\` to specify geometry of an
    // exocyclic double bond; the ring bonds themselves must remain Aromatic.
    // Regression: previously these were stored as Down/Up, breaking SMARTS (:a).
    use chematic_core::BondOrder;
    let mol = parse(r"N=c1\c(O)c(O)\c1=N").unwrap();
    let aromatic_ring_bonds: Vec<_> = mol
        .bonds()
        .filter(|(_, bond)| {
            mol.atom(bond.atom1).aromatic
                && mol.atom(bond.atom2).aromatic
                && bond.order == BondOrder::Aromatic
        })
        .collect();
    assert_eq!(
        aromatic_ring_bonds.len(),
        4,
        "all 4 ring bonds should be Aromatic"
    );
}

// -----------------------------------------------------------------
// Canonical-Stereo-D0: `resolve_aromatic_direction_stash` must be
// applied consistently at all three sites that can create a bond from a
// pending `/`/`\` marker -- the normal chain-edge path (already had the
// guard), the ring-closure path, and the branch-attachment path (the
// first bond of a `(...)` group) -- both of which silently stored the
// marker as the bond's own literal `Up`/`Down` order instead of
// stashing it. That inconsistency meant the SAME physical bond could be
// correctly stashed on one parse (when it happened to fall on a chain
// edge) and incorrectly stored literally on a later parse of chematic's
// own canonical output (when the SAME bond became a ring-closure or
// branch-attachment edge instead) -- silently changing which CIP
// elements were visible to `assign_ez` between the two parses, purely
// as an artifact of which of the three paths happened to run.
//
// These tests check the structural representation (bond order + stash)
// is correct and stable across a round trip for each of the three
// paths. They deliberately do NOT assert `assign_cip`'s E/Z output --
// `assign_ez`/`substituent_is_up` do not read `bond_direction` at all
// (a separate, pre-existing, real gap: this direction genuinely encodes
// legitimate exocyclic-imine E/Z per an independent RDKit check, not
// fabricated stereo -- see EZ-A0/EZ-S1). Pinning `[]` here as the
// expected CIP output would calcify that gap as intended behavior and
// get in the way of the follow-up that fixes it.
// -----------------------------------------------------------------

/// Assert `bidx` is `Aromatic` with a stashed direction, on a molecule
/// where both endpoints are aromatic -- the correct representation
/// regardless of which parser path created the bond.
#[track_caller]
fn assert_aromatic_with_stash(mol: &Molecule, bidx: BondIdx, path_label: &str) {
    let bond = mol.bond(bidx);
    assert!(
        mol.atom(bond.atom1).aromatic && mol.atom(bond.atom2).aromatic,
        "{path_label}: test setup sanity -- both endpoints must be aromatic"
    );
    assert_eq!(
        bond.order,
        BondOrder::Aromatic,
        "{path_label}: bond order must stay Aromatic, not the literal Up/Down marker"
    );
    assert!(
        mol.bond_direction(bidx).is_some(),
        "{path_label}: the direction must be stashed on the side channel"
    );
    assert_eq!(
        mol.bond_direction_anchor(bidx),
        Some(bond.atom1),
        "{path_label}: the stash anchor must follow the parser-side source endpoint"
    );
}

/// Re-parse `canonical_smiles(mol)` and assert the same structural
/// invariant holds again: some aromatic-aromatic bond is `Aromatic`
/// order with a stashed direction (not necessarily the same `BondIdx`,
/// since canonicalization renumbers atoms/bonds).
#[track_caller]
fn assert_round_trip_preserves_stash_representation(mol: &Molecule, path_label: &str) {
    let c1 = canonical_smiles(mol);
    let mol2 = parse(&c1).unwrap_or_else(|e| panic!("{path_label}: re-parse '{c1}': {e}"));
    let has_stashed_aromatic_bond = (0..mol2.bond_count()).any(|i| {
        let bidx = BondIdx(i as u32);
        let bond = mol2.bond(bidx);
        mol2.atom(bond.atom1).aromatic
            && mol2.atom(bond.atom2).aromatic
            && bond.order == BondOrder::Aromatic
            && mol2.bond_direction(bidx).is_some()
    });
    assert!(
        has_stashed_aromatic_bond,
        "{path_label}: round-tripped molecule '{c1}' must still represent the direction \
         as a stash on an Aromatic bond, not a literal Up/Down order"
    );
    let c2 = canonical_smiles(&mol2);
    assert_eq!(
        c1, c2,
        "{path_label}: canonical_smiles must be stable across a round trip"
    );
}

#[test]
fn direction_stash_normal_chain_edge() {
    // The `/`/`\` sits between the ring-opening atom and the very next
    // chain atom -- the plain tree-edge path (already had the guard
    // before this fix; kept as a path-1 regression pin, not a new fix).
    // Both imines carry markers at their N ends too: a marker on one
    // side only specifies nothing and canonical output drops it, as
    // RDKit does.
    let mol = parse(r"C/N=c1\c(O)c(O)\c1=N/C").unwrap();
    assert_aromatic_with_stash(&mol, BondIdx(2), "path1(chain-edge)");
    assert_round_trip_preserves_stash_representation(&mol, "path1(chain-edge)");
}

#[test]
fn direction_stash_ring_closure_edge() {
    // The `/` sits directly on the ring-closing bond itself (`...c/1`),
    // routed through `close_or_open_ring` -- previously stored as a
    // literal Up/Down order between two aromatic atoms.
    let mol = parse(r"C/N=c1ccccc/1").unwrap();
    assert_aromatic_with_stash(&mol, BondIdx(7), "path2(ring-closure)");
    assert_round_trip_preserves_stash_representation(&mol, "path2(ring-closure)");
}

#[test]
fn direction_stash_branch_attachment_edge() {
    // The `/` is the first bond inside a `(...)` branch, attaching the
    // branch's first atom to its parent -- routed through
    // `parse_chain`'s `attach_to` handling, not the plain chain-edge
    // loop. Previously stored as a literal Up/Down order between two
    // aromatic atoms (this is the real-world corpus mechanism: a
    // canonical DFS can route what was a chain edge on the original
    // parse through a branch attachment on the next).
    // (`O=c1cccc(/c1=N/C)`: the branch's first atom carries an exocyclic
    // imine whose N end is marked too, so the marker means something
    // and survives canonical output.)
    let mol = parse(r"O=c1cccc(/c1=N/C)").unwrap();
    assert_aromatic_with_stash(&mol, BondIdx(5), "path3(branch-attachment)");
    assert_round_trip_preserves_stash_representation(&mol, "path3(branch-attachment)");
}

#[test]
fn test_malformed_bracket_atoms_return_errors() {
    assert!(matches!(
        parse("[C"),
        Err(SmilesError::InvalidBracketAtom { .. })
    ));
    assert!(matches!(
        parse("[CH4"),
        Err(SmilesError::InvalidBracketAtom { .. })
    ));
    assert!(matches!(
        parse("[Q]"),
        Err(SmilesError::InvalidBracketAtom { .. })
    ));
    assert!(matches!(
        parse("[C@"),
        Err(SmilesError::InvalidBracketAtom { .. })
    ));
}

#[test]
fn parse_with_limits_reports_input_and_graph_limits() {
    let limits = SmilesParseLimits {
        max_input_bytes: 2,
        ..SmilesParseLimits::default()
    };
    assert!(matches!(
        parse_with_limits("CCO", &limits),
        Err(SmilesError::ResourceLimit {
            resource: "input bytes",
            ..
        })
    ));

    let limits = SmilesParseLimits {
        max_atoms: 2,
        ..SmilesParseLimits::default()
    };
    assert!(matches!(
        parse_with_limits("CCO", &limits),
        Err(SmilesError::ResourceLimit {
            resource: "atoms",
            ..
        })
    ));

    let limits = SmilesParseLimits {
        max_bonds: 1,
        ..SmilesParseLimits::default()
    };
    assert!(matches!(
        parse_with_limits("CCC", &limits),
        Err(SmilesError::ResourceLimit {
            resource: "bonds",
            ..
        })
    ));
}

#[test]
fn acyclic_parse_seeds_empty_ring_flags() {
    let molecule = parse("CCO.CN").unwrap();
    assert!(
        molecule
            .derived_if_computed::<bool>(chematic_core::DerivedSlot::RdkitModelGate)
            .is_some(),
        "the shared RDKit-model gate is computed once during parsing"
    );
    let flags = molecule
        .derived_if_computed::<Vec<bool>>(chematic_core::DerivedSlot::RingBondFlags)
        .expect("acyclic SMILES proves that every bond is non-cyclic");
    assert_eq!(&*flags, &vec![false; molecule.bond_count()]);

    let cyclic = parse("C1CCCCC1").unwrap();
    assert!(
        cyclic
            .derived_if_computed::<bool>(chematic_core::DerivedSlot::RdkitModelGate)
            .is_some(),
        "cyclic inputs also cache the cheap model-disagreement predicate"
    );
    assert!(
        cyclic
            .derived_if_computed::<Vec<bool>>(chematic_core::DerivedSlot::RingBondFlags)
            .is_none(),
        "cyclic inputs still require perception"
    );
}
