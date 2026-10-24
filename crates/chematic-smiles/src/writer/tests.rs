use super::*;

#[test]
fn vacant_nontetrahedral_centers_match_pinned_rdkit_after_reordering() {
    // RDKit 2025.09.3 canonical isomeric SMILES of the source structures.
    // Include root and non-root centers, three geometries, and tags for
    // which a reordered ligand list needs a different permutation number.
    for (input, expected) in [
        ("F[Pt@SP2](Cl)Br", "[F][Pt@SP2]([Cl])[Br]"),
        ("F[P@TB1](Cl)(Br)I", "F[P@TB1](Cl)(Br)I"),
        ("[P@TB19](F)(Cl)(Br)I", "F[P@TB16](Cl)(Br)I"),
        ("I[P@TB20](Br)(Cl)F", "F[P@TB3](Cl)(Br)I"),
        ("F[S@OH27](Cl)(Br)(I)N", "N[S@OH8](F)(Cl)(Br)I"),
        ("[S@OH30](F)(Cl)(Br)(I)N", "N[S@OH9](F)(Cl)(Br)I"),
        ("N[S@OH15](I)(Br)(Cl)F", "N[S@OH12](F)(Cl)(Br)I"),
    ] {
        let mol = crate::parse(input).unwrap();
        assert_eq!(
            crate::rdkit_canonical_smiles(&mol).unwrap(),
            expected,
            "{input}"
        );
        for written in [write(&mol), crate::canonical_smiles(&mol)] {
            let read_back = crate::parse(&written).unwrap();
            assert_eq!(
                crate::rdkit_canonical_smiles(&read_back).unwrap(),
                expected,
                "{input} -> {written}"
            );
        }
        for seed in 1..=16 {
            let written = crate::random_smiles(&mol, seed);
            let read_back = crate::parse(&written).unwrap();
            assert_eq!(
                crate::rdkit_canonical_smiles(&read_back).unwrap(),
                expected,
                "{input} -> {written} (seed {seed})"
            );
        }
    }
}

#[test]
fn three_ligand_square_planar_tags_survive_canonical_reordering() {
    // RDKit 2025.09.3 oracle: the canonical traversal starts at Cl and
    // swaps F/Cl relative to the input, while the vacant site stays put.
    for (input, expected) in [
        ("F[Pt@SP1](Cl)Br", "Cl[Pt@SP1](F)Br"),
        ("F[Pt@SP2](Cl)Br", "Cl[Pt@SP3](F)Br"),
        ("F[Pt@SP3](Cl)Br", "Cl[Pt@SP2](F)Br"),
    ] {
        let output = crate::canonical_smiles(&crate::parse(input).unwrap());
        assert_eq!(output, expected);
        assert_eq!(
            crate::canonical_smiles(&crate::parse(&output).unwrap()),
            output
        );
    }
}

#[test]
fn three_ligand_square_planar_remap_requires_same_distinct_ligands() {
    use chematic_core::SquarePlanarPermutation::SP1;
    assert_eq!(
        remap_three_ligand_square_planar(SP1, [1, 2, 3], [2, 1, 3]),
        Some(SP1)
    );
    assert_eq!(
        remap_three_ligand_square_planar(SP1, [1, 2, 3], [1, 2, 4]),
        None
    );
    assert_eq!(
        remap_three_ligand_square_planar(SP1, [1, 2, 3], [1, 1, 3]),
        None
    );
}

/// `write` re-expresses a stored `@`/`@@` for the order it emits: a
/// molecule whose recorded neighbour order differs from the DFS order
/// (here: atoms renumbered, as reaction products are) keeps its
/// configuration, and an organic-subset stereocentre is bracketed.
#[test]
fn write_keeps_configuration_when_dfs_order_differs() {
    for input in [
        "N[C@@H](C)C(=O)O",
        "O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O",
        "F[C@]1(Cl)CCC1",
        "C[C@]12CC[C@H]3[C@@H](CC=C4C[C@@H](O)CC[C@]34C)[C@@H]1CC[C@@H]2O",
        "F[Pt@SP2](Cl)Br",
        "F[P@TB1](Cl)(Br)I",
        "[P@TB19](F)(Cl)(Br)I",
        "F[S@OH27](Cl)(Br)(I)N",
        "[S@OH30](F)(Cl)(Br)(I)N",
    ] {
        let mol = crate::parse(input).unwrap();
        let want = crate::canonical_smiles(&mol);
        let n = mol.atom_count();
        // Reverse the atom order, carrying the stereo neighbour order.
        let mut b = chematic_core::MoleculeBuilder::new();
        for i in (0..n).rev() {
            let mut atom = mol.atom(AtomIdx(i as u32)).clone();
            if atom.chirality.nontetrahedral().is_none() {
                atom.hydrogen_count = None;
            }
            b.add_atom(atom);
        }
        let new = |i: u32| (n as u32 - 1) - i;
        for (_, bond) in mol.bonds() {
            b.add_bond(
                AtomIdx(new(bond.atom1.0)),
                AtomIdx(new(bond.atom2.0)),
                bond.order,
            )
            .unwrap();
        }
        for (i, _) in mol.atoms() {
            if let Some(order) = mol.stereo_neighbor_order(i) {
                let order = order
                    .iter()
                    .map(|&v| {
                        if matches!(v, STEREO_H_SENTINEL | STEREO_VACANT_SENTINEL) {
                            v
                        } else {
                            new(v)
                        }
                    })
                    .collect();
                b.set_stereo_neighbor_order(AtomIdx(new(i.0)), order);
            }
        }
        let reversed = b.build();
        let written = write(&reversed);
        assert_eq!(
            crate::canonical_smiles(&crate::parse(&written).unwrap()),
            want,
            "{input}: {written}"
        );
    }
}
use crate::parser::parse;
use chematic_core::{Atom, Element, MoleculeBuilder};

/// Parse → write → re-parse → verify atom/bond counts are preserved.
fn roundtrip(smiles: &str) {
    let mol1 = parse(smiles).expect(smiles);
    let out = write(&mol1);
    let mol2 = parse(&out).unwrap_or_else(|e| {
        panic!(
            "roundtrip failed for '{}': wrote '{}', error: {e}",
            smiles, out
        )
    });
    assert_eq!(
        mol1.atom_count(),
        mol2.atom_count(),
        "atom count mismatch: input='{}' output='{}'",
        smiles,
        out
    );
    assert_eq!(
        mol1.bond_count(),
        mol2.bond_count(),
        "bond count mismatch: input='{}' output='{}'",
        smiles,
        out
    );
}

#[test]
fn test_write_methane() {
    assert_eq!(write(&parse("C").unwrap()), "C");
}
#[test]
fn test_write_ethane() {
    assert_eq!(write(&parse("CC").unwrap()), "CC");
}

#[test]
fn writer_uses_extended_ring_labels_without_losing_closures() {
    // A 12×12 grid has cycle rank 121, forcing the DFS writer to allocate
    // labels beyond 99. It is intentionally a graph-built fixture because
    // legacy OpenSMILES cannot express the input without extended labels.
    let width = 12;
    let mut builder = MoleculeBuilder::with_capacity(width * width, 264);
    let atoms: Vec<_> = (0..width * width)
        .map(|_| builder.add_atom(Atom::new(Element::C)))
        .collect();
    for row in 0..width {
        for column in 0..width {
            let idx = row * width + column;
            if column + 1 < width {
                builder
                    .add_bond(atoms[idx], atoms[idx + 1], BondOrder::Single)
                    .unwrap();
            }
            if row + 1 < width {
                builder
                    .add_bond(atoms[idx], atoms[idx + width], BondOrder::Single)
                    .unwrap();
            }
        }
    }
    let mol = builder.build();
    let smiles = write(&mol);
    assert!(smiles.contains("%(100)"), "{smiles}");
    let reparsed = parse(&smiles).expect("extended ring labels must parse");
    assert_eq!(reparsed.atom_count(), mol.atom_count());
    assert_eq!(reparsed.bond_count(), mol.bond_count());
}

#[test]
fn test_roundtrip_propane() {
    roundtrip("CCC");
}
#[test]
fn test_roundtrip_isobutane() {
    roundtrip("CC(C)C");
}
#[test]
fn test_roundtrip_ethanol() {
    roundtrip("CCO");
}
#[test]
fn test_roundtrip_acetic_acid() {
    roundtrip("CC(=O)O");
}
#[test]
fn test_roundtrip_cyclohexane() {
    roundtrip("C1CCCCC1");
}
#[test]
fn test_roundtrip_benzene_kekule() {
    roundtrip("C1=CC=CC=C1");
}
#[test]
fn test_roundtrip_benzene_arom() {
    roundtrip("c1ccccc1");
}
#[test]
fn test_roundtrip_pyridine() {
    roundtrip("c1ccncc1");
}
#[test]
fn test_roundtrip_naphthalene() {
    roundtrip("c1ccc2ccccc2c1");
}
#[test]
fn test_roundtrip_chlorobenzene() {
    roundtrip("c1ccccc1Cl");
}
#[test]
fn test_roundtrip_13c() {
    roundtrip("[13C]");
}
#[test]
fn test_roundtrip_ammonium() {
    roundtrip("[NH4+]");
}
#[test]
fn test_roundtrip_disconnected() {
    roundtrip("[Na+].[Cl-]");
}
#[test]
fn test_roundtrip_aspirin() {
    roundtrip("CC(=O)Oc1ccccc1C(=O)O");
}
#[test]
fn test_roundtrip_caffeine() {
    roundtrip("Cn1cnc2c1c(=O)n(c(=O)n2C)C");
}

// Bracket atoms with `hydrogen_count: None` (implicit H left uninferred by the
// caller, e.g. after a programmatic charge/isotope/atom-map edit) must still get
// their implicit hydrogens written — regression tests for the bracket-H bug found
// via MRV oracle validation (isotope_0/2/3, charge_0/3, atom_map_0/1/2, disconnected_3
// in validation/mrv_io_parity_summary.json).
#[test]
fn test_bracket_implicit_h_ammonium_charge_only() {
    let mut b = MoleculeBuilder::new();
    let mut n = Atom::new(Element::N);
    n.charge = 1;
    b.add_atom(n);
    assert_eq!(write(&b.build()), "[NH4+]");
}

#[test]
fn test_bracket_implicit_h_isotope_only() {
    let mut b = MoleculeBuilder::new();
    let mut c = Atom::new(Element::C);
    c.isotope = Some(13);
    b.add_atom(c);
    assert_eq!(write(&b.build()), "[13CH4]");
}

#[test]
fn test_bracket_implicit_h_atom_map_only() {
    let mut b = MoleculeBuilder::new();
    let mut c0 = Atom::new(Element::C);
    c0.atom_map = Some(7);
    let c0 = b.add_atom(c0);
    let c1 = b.add_atom(Atom::new(Element::C));
    b.add_bond(c0, c1, BondOrder::Single).unwrap();
    assert_eq!(write(&b.build()), "[CH3:7]C");
}

#[test]
fn test_wildcard_roundtrip_preserves_atom_and_map() {
    let plain = parse("[*]C").unwrap();
    assert_eq!(write(&plain), "[*]C");

    let mapped = parse("[*:17]C").unwrap();
    assert_eq!(write(&mapped), "[*:17]C");
}

#[test]
fn test_bracket_implicit_h_isotope_and_atom_map() {
    let mut b = MoleculeBuilder::new();
    let mut c0 = Atom::new(Element::C);
    c0.isotope = Some(13);
    c0.atom_map = Some(7);
    let c0 = b.add_atom(c0);
    let c1 = b.add_atom(Atom::new(Element::C));
    b.add_bond(c0, c1, BondOrder::Single).unwrap();
    assert_eq!(write(&b.build()), "[13CH3:7]C");
}

#[test]
fn test_bracket_implicit_h_hydroxide_charge_only() {
    let mut b = MoleculeBuilder::new();
    let mut o = Atom::new(Element::O);
    o.charge = -1;
    b.add_atom(o);
    assert_eq!(write(&b.build()), "[OH-]");
}

// ── Standalone wedge/hash bond must not emit a meaningless SMILES
// directional token (docs/rfcs/stereo2d_reader_integration_rfc.md §3/§7,
// docs/rfcs/stereo2d_local_parity_calibration.md "Scope note") ────────────────
//
// `BondOrder::Up`/`Down` is overloaded for two unrelated concepts: 2D
// wedge/hash depiction (set by MDL/CDXML/MRV/KET-style readers on a
// stereocenter's substituent bond) and the OpenSMILES `/`/`\` E/Z
// directional-bond marker (only meaningful adjacent to a double bond).
// A wedge bond with no adjacent double bond has nothing to mark; writing
// `/`/`\` for it is self-consistent only within chematic's own
// round-trip and is silently dropped by any spec-compliant SMILES
// parser (e.g. RDKit) on re-parse.

/// Minimal repro from the RFC's own end-to-end trace (§3): a MOL V2000
/// wedge bond (`BondOrder::Up`) on a tetrahedral center's substituent,
/// no double bond anywhere in the molecule -- the RFC's own trace shows
/// naive `write()` producing `"C(F)(Cl)/Br"`, a token any RDKit re-parse
/// silently discards.
#[test]
fn test_standalone_solid_wedge_not_written_as_slash() {
    let mut b = MoleculeBuilder::new();
    let c = b.add_atom(Atom::new(Element::C));
    let f = b.add_atom(Atom::new(Element::F));
    let cl = b.add_atom(Atom::new(Element::CL));
    let br = b.add_atom(Atom::new(Element::BR));
    b.add_bond(c, f, BondOrder::Single).unwrap();
    b.add_bond(c, cl, BondOrder::Single).unwrap();
    let wedge_bond = b.add_bond(c, br, BondOrder::Up).unwrap();
    let mol = b.build();

    let out = write(&mol);
    assert!(
        !out.contains('/') && !out.contains('\\'),
        "standalone solid wedge (no adjacent double bond) must not be written \
             as a directional token: got '{out}'"
    );

    // The fix only changes what gets printed -- the bond's real order in
    // memory must be untouched (this is 2D depiction metadata, not
    // something the writer should mutate).
    assert_eq!(mol.bond(wedge_bond).order, BondOrder::Up);

    // Round-trip: connectivity survives re-parsing the token-free output.
    let mol2 = parse(&out).unwrap_or_else(|e| panic!("re-parse '{out}': {e}"));
    assert_eq!(mol2.atom_count(), mol.atom_count());
    assert_eq!(mol2.bond_count(), mol.bond_count());
}

/// Same shape as above but with a hash wedge (`BondOrder::Down`).
#[test]
fn test_standalone_hash_wedge_not_written_as_backslash() {
    let mut b = MoleculeBuilder::new();
    let c = b.add_atom(Atom::new(Element::C));
    let f = b.add_atom(Atom::new(Element::F));
    let cl = b.add_atom(Atom::new(Element::CL));
    let br = b.add_atom(Atom::new(Element::BR));
    b.add_bond(c, f, BondOrder::Single).unwrap();
    b.add_bond(c, cl, BondOrder::Single).unwrap();
    let wedge_bond = b.add_bond(c, br, BondOrder::Down).unwrap();
    let mol = b.build();

    let out = write(&mol);
    assert!(
        !out.contains('/') && !out.contains('\\'),
        "standalone hash wedge (no adjacent double bond) must not be written \
             as a directional token: got '{out}'"
    );
    assert_eq!(mol.bond(wedge_bond).order, BondOrder::Down);

    let mol2 = parse(&out).unwrap_or_else(|e| panic!("re-parse '{out}': {e}"));
    assert_eq!(mol2.atom_count(), mol.atom_count());
    assert_eq!(mol2.bond_count(), mol.bond_count());
}

/// A genuinely 4-heavy-atom stereocenter (no implicit/explicit H at all)
/// with a standalone wedge -- PR #130's `tetrahedral_4heavy_no_h`
/// fixture shape, added there specifically so it isn't conflated with
/// the 3-heavy+H case above.
#[test]
fn test_standalone_wedge_four_heavy_neighbors_no_h() {
    let mut b = MoleculeBuilder::new();
    let c = b.add_atom(Atom::new(Element::C));
    let f = b.add_atom(Atom::new(Element::F));
    let cl = b.add_atom(Atom::new(Element::CL));
    let br = b.add_atom(Atom::new(Element::BR));
    let i = b.add_atom(Atom::new(Element::I));
    b.add_bond(c, f, BondOrder::Single).unwrap();
    b.add_bond(c, cl, BondOrder::Single).unwrap();
    b.add_bond(c, br, BondOrder::Single).unwrap();
    b.add_bond(c, i, BondOrder::Up).unwrap();
    let mol = b.build();

    let out = write(&mol);
    assert!(
        !out.contains('/') && !out.contains('\\'),
        "4-heavy-neighbor standalone wedge must not be written as a \
             directional token: got '{out}'"
    );
}

/// A wedge bond that lands on a ring-closure edge (not a plain tree
/// edge) must be suppressed identically -- the fix touches both
/// emission sites in `write_chain`/`dfs_mark`.
#[test]
fn test_standalone_wedge_on_ring_closure_bond() {
    let mut b = MoleculeBuilder::new();
    let c1 = b.add_atom(Atom::new(Element::C));
    let c2 = b.add_atom(Atom::new(Element::C));
    let c3 = b.add_atom(Atom::new(Element::C));
    b.add_bond(c1, c2, BondOrder::Single).unwrap();
    b.add_bond(c2, c3, BondOrder::Single).unwrap();
    // Whichever of these two becomes the DFS back-edge, it carries a
    // wedge with no adjacent double bond anywhere in this molecule.
    b.add_bond(c1, c3, BondOrder::Up).unwrap();
    let mol = b.build();

    let out = write(&mol);
    assert!(
        !out.contains('/') && !out.contains('\\'),
        "standalone wedge on a ring-closure bond must not be written as \
             a directional token: got '{out}'"
    );
    let mol2 = parse(&out).unwrap_or_else(|e| panic!("re-parse '{out}': {e}"));
    assert_eq!(mol2.atom_count(), 3);
    assert_eq!(mol2.bond_count(), 3);
}

/// The fix must not disturb `Atom.chirality`/`Molecule::stereo_neighbor_order`
/// (the P1-S1a-core local-parity metadata) while suppressing a wedge
/// bond's spurious directional token on the same stereocenter -- these
/// are two unrelated mechanisms (bond-token emission vs `emit_atom`'s
/// chirality symbol) and must keep working independently.
#[test]
fn test_standalone_wedge_does_not_disturb_stereocenter_chirality() {
    let mut center = Atom::new(Element::C);
    center.chirality = chematic_core::Chirality::Clockwise;
    // Force bracket notation so the chirality symbol is actually
    // printed (`needs_bracket` doesn't key off chirality alone).
    center.hydrogen_count = Some(0);
    let mut b = MoleculeBuilder::new();
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

    let out = write(&mol);
    assert!(
        !out.contains('/') && !out.contains('\\'),
        "standalone wedge on a stereocenter's substituent must still be \
             suppressed: got '{out}'"
    );
    assert!(
        out.contains('@'),
        "the stereocenter's own chirality symbol must still be printed \
             alongside the (now-suppressed) wedge bond: got '{out}'"
    );

    // The fix must not have touched the stereocenter metadata itself --
    // only what gets printed for the wedge bond's token.
    assert_eq!(mol.atom(c).chirality, chematic_core::Chirality::Clockwise);
    assert_eq!(mol.bond(wedge_bond).order, BondOrder::Up);
    assert_eq!(
        mol.stereo_neighbor_order(c),
        Some([f.0, cl.0, br.0, iodine.0].as_slice())
    );
}

/// Round-trip check with charge and isotope present, per the required
/// test list: connectivity/charge/isotope must survive re-parsing the
/// (now token-free) output of a molecule carrying a standalone wedge.
#[test]
fn test_standalone_wedge_roundtrip_preserves_charge_and_isotope() {
    let mut b = MoleculeBuilder::new();
    let mut c = Atom::new(Element::C);
    c.isotope = Some(13);
    let c = b.add_atom(c);
    let f = b.add_atom(Atom::new(Element::F));
    let cl = b.add_atom(Atom::new(Element::CL));
    let mut o = Atom::new(Element::O);
    o.charge = -1;
    let o = b.add_atom(o);
    b.add_bond(c, f, BondOrder::Single).unwrap();
    b.add_bond(c, cl, BondOrder::Single).unwrap();
    let wedge_bond = b.add_bond(c, o, BondOrder::Up).unwrap();
    let mol = b.build();

    let out = write(&mol);
    assert!(
        !out.contains('/') && !out.contains('\\'),
        "standalone wedge must not be written as a directional token: got '{out}'"
    );

    let mol2 = parse(&out).unwrap_or_else(|e| panic!("re-parse '{out}': {e}"));
    assert_eq!(mol2.atom_count(), mol.atom_count());
    assert_eq!(mol2.bond_count(), mol.bond_count());
    assert!(
        (0..mol2.atom_count()).any(|i| mol2.atom(AtomIdx(i as u32)).isotope == Some(13)),
        "isotope must survive the round trip: '{out}'"
    );
    assert!(
        (0..mol2.atom_count()).any(|i| mol2.atom(AtomIdx(i as u32)).charge == -1),
        "charge must survive the round trip: '{out}'"
    );
    // The original wedge bond's order in memory is still untouched.
    assert_eq!(mol.bond(wedge_bond).order, BondOrder::Up);
}

/// Regression guard: a genuine E/Z directional marker (real double
/// bond, real geometry) must still be emitted -- this fix must not
/// suppress legitimate cases.
#[test]
fn test_genuine_ez_directional_bond_still_written() {
    let mol = parse("F/C=C/F").unwrap();
    let out = write(&mol);
    assert!(
        out.contains('/') || out.contains('\\'),
        "genuine double-bond-adjacent directional marker must survive \
             plain write(): got '{out}'"
    );
    let mol2 = parse(&out).unwrap_or_else(|e| panic!("re-parse '{out}': {e}"));
    let has_directional = (0..mol2.bond_count()).any(|i| {
        matches!(
            mol2.bond(BondIdx(i as u32)).order,
            BondOrder::Up | BondOrder::Down
        )
    });
    assert!(
        has_directional,
        "re-parsed molecule lost its E/Z directional bonds: '{out}'"
    );
}

// ── plain write() must read the bond_direction side channel too
// (previously canonical.rs-only; see docs/rfcs/stereo2d_reader_integration_rfc.md
// "both the non-canonical and canonical writer read the same effective
// direction") ─────────────────────────────────────────────────────────

/// F-C=C-F with the direction stashed via `Molecule::bond_direction`
/// instead of a literal `BondOrder::Up`/`Down` -- exactly the shape a
/// 2D-coordinate-derived E/Z direction (Track B) or the pre-existing
/// aromatic-bond stash produces. Before this fix, plain `write()` never
/// consulted the stash at all, so it silently emitted no `/`/`\` here.
#[test]
fn test_plain_write_reads_bond_direction_stash() {
    let mut b = MoleculeBuilder::new();
    let f1 = b.add_atom(Atom::new(Element::F));
    let c1 = b.add_atom(Atom::new(Element::C));
    let c2 = b.add_atom(Atom::new(Element::C));
    let f2 = b.add_atom(Atom::new(Element::F));
    let b1 = b.add_bond(f1, c1, BondOrder::Single).unwrap();
    b.add_bond(c1, c2, BondOrder::Double).unwrap();
    let b2 = b.add_bond(c2, f2, BondOrder::Single).unwrap();
    let mut mol = b.build();
    // Stash directions instead of literal Up/Down -- mirrors how the new
    // 2D E/Z perception stage writes (never mutating `bond.order`).
    mol.set_bond_direction(b1, BondOrder::Up);
    mol.set_bond_direction(b2, BondOrder::Up);

    let out = write(&mol);
    assert!(
        out.contains('/') || out.contains('\\'),
        "plain write() must surface a stashed bond_direction: got '{out}'"
    );
    // Round-trip through chematic's own parser must preserve real E/Z.
    let mol2 = parse(&out).unwrap();
    let has_directional = (0..mol2.bond_count()).any(|i| {
        matches!(
            mol2.bond(BondIdx(i as u32)).order,
            BondOrder::Up | BondOrder::Down
        )
    });
    assert!(has_directional, "stashed E/Z lost on round-trip: '{out}'");
}

/// A ring-closure bond is emitted from BOTH endpoints (the "open" atom
/// and, later, the "close" atom) but is one physical bond with one
/// `atom1`/`atom2` pair -- each side must independently re-derive its own
/// correct character rather than sharing one (necessarily backwards on
/// one side) value. This directly asserts the shared helper both
/// `dfs_mark` emission sites depend on: a bond declared `atom1 = c2`
/// must read as `Up` from `c2`'s side and `Down` from `c1`'s side.
#[test]
fn test_ring_closure_directional_bond_reoriented_per_endpoint() {
    let mut b = MoleculeBuilder::new();
    let c1 = b.add_atom(Atom::new(Element::C));
    let c2 = b.add_atom(Atom::new(Element::C));
    let bidx = b.add_bond(c2, c1, BondOrder::Up).unwrap(); // atom1=c2, atom2=c1
    let mol = b.build();

    // Seen from c2 (== atom1): Up stays Up.
    assert_eq!(direction_from(&mol, bidx, BondOrder::Up, c2), BondOrder::Up);
    // Seen from c1 (!= atom1): Up flips to Down.
    assert_eq!(
        direction_from(&mol, bidx, BondOrder::Up, c1),
        BondOrder::Down
    );
}

#[test]
fn test_direction_from_uses_aromatic_stash_anchor() {
    let mut b = MoleculeBuilder::new();
    let c1 = b.add_atom(Atom::aromatic(Element::C));
    let c2 = b.add_atom(Atom::aromatic(Element::C));
    let bidx = b.add_bond(c2, c1, BondOrder::Aromatic).unwrap();
    let mut mol = b.build();
    mol.set_bond_direction(bidx, BondOrder::Up);
    mol.set_bond_direction_anchor(bidx, c1);

    // The stored bond orientation is c2 -> c1, but the parser-side
    // marker was read at c1. Reorientation must follow that anchor.
    assert_eq!(direction_from(&mol, bidx, BondOrder::Up, c1), BondOrder::Up);
    assert_eq!(
        direction_from(&mol, bidx, BondOrder::Up, c2),
        BondOrder::Down
    );
}

/// A wedge whose direction is ALREADY suppressed (no adjacent double
/// bond) must stay suppressed even after being resolved through
/// `raw_bond_direction`/`direction_from` -- the reorientation step must
/// not accidentally resurrect a token `suppress_standalone_wedge` would
/// otherwise drop.
#[test]
fn test_standalone_wedge_still_suppressed_after_reorientation_helpers() {
    let mut b = MoleculeBuilder::new();
    let c = b.add_atom(Atom::new(Element::C));
    let f = b.add_atom(Atom::new(Element::F));
    let cl = b.add_atom(Atom::new(Element::CL));
    let br = b.add_atom(Atom::new(Element::BR));
    b.add_bond(c, f, BondOrder::Single).unwrap();
    b.add_bond(c, cl, BondOrder::Single).unwrap();
    b.add_bond(c, br, BondOrder::Up).unwrap();
    let mol = b.build();
    let out = write(&mol);
    assert!(!out.contains('/') && !out.contains('\\'), "got '{out}'");
}
