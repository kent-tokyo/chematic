//! Square-planar stereo (`@SP1`/`@SP2`/`@SP3`) MOL/SDF I/O integration
//! tests: read/write across V2000 and V3000, atom-renumbering invariance,
//! and the FlatZero round-trip limitation this design accepts.
//!
//! Design and RDKit oracle provenance:
//! `docs/rfcs/square_planar_mol_io_rfc.md`. Unit-level coverage of the
//! classifier/validator themselves lives in `mol2000.rs`'s own
//! `square_planar_tests` module; these tests cover the cross-cutting,
//! full-pipeline invariants that only exist once a real reader/writer round
//! trip is in the loop -- mirroring `stereo_reader_integration.rs`'s own
//! split between unit-level parity math and integration-level round trips.
//!
//! Every fixture here places the whole molecule off the z=0 plane (nonzero,
//! constant z) -- see the RFC §10 for why: a literal z=0 conformer is this
//! crate's pre-existing, load-bearing signal for "no real 3D data", and
//! using it here would produce a fixture whose tag silently fails to
//! round-trip, not a fixture that demonstrates round-tripping.

use chematic_core::{
    Atom, AtomIdx, BondOrder, Chirality, Coords3D, Element, MoleculeBuilder, Point3,
    SquarePlanarPermutation, remap_square_planar_tag,
};
use chematic_mol::mol2000::{
    MolFormat, MolMetadata, UnsupportedStereoReason, read_mol_with_diagnostics,
    validate_square_planar_for_write, write_mol_with_conformer, write_mol_with_conformer_checked,
    write_sdf_record_with_conformer_checked,
};
use chematic_mol::mol3000::{
    read_mol_v3000_with_diagnostics, write_mol_v3000_with_conformer_checked,
};
use chematic_mol::sdf::SdfRecordReader;
use chematic_smiles::{canonical_smiles, parse};

/// Ideal (undistorted) neighbor positions for `tag`, in `trans_pairs()` slot
/// order, all sharing `z` (nonzero -- see module docs). Reuses
/// `trans_pairs()` (the same oracle-verified source of truth
/// `chematic_core::remap_square_planar_tag` uses) rather than a hand-picked
/// per-tag layout, and applies `distortion_deg` symmetrically to both
/// members of each trans pair so the pair's angle moves off 180 degrees by
/// `2 * distortion_deg` -- used by the near-degenerate robustness test.
fn ideal_neighbor_positions(
    tag: SquarePlanarPermutation,
    z: f64,
    distortion_deg: f64,
) -> [Point3; 4] {
    let [(a, b), (c, d)] = tag.trans_pairs();
    let mut pts = [Point3::new(0.0, 0.0, z); 4];
    let at = |deg: f64| -> Point3 {
        let r = deg.to_radians();
        Point3::new(1.5 * r.cos(), 1.5 * r.sin(), z)
    };
    pts[a as usize] = at(45.0 + distortion_deg);
    pts[b as usize] = at(225.0 - distortion_deg);
    pts[c as usize] = at(135.0 - distortion_deg);
    pts[d as usize] = at(315.0 + distortion_deg);
    pts
}

/// One SMILES-declared square-planar fixture: the parsed molecule, a
/// synthetic `Coords3D` conformer whose real coordinates geometrically
/// encode the declared tag (using the parser's own `stereo_neighbor_order`,
/// not an assumed atom order), the declared tag itself, and its declared
/// neighbor order -- the latter is what makes the differential check below
/// frame-explicit rather than assuming the MOL round trip preserves
/// neighbor-listing order. The center atom's index is re-derived on demand
/// via `pt_atom_idx` rather than stored (this fixture family only ever has
/// one Pt-like center).
struct SquarePlanarFixture {
    mol: chematic_core::Molecule,
    conformer: Coords3D,
    tag: SquarePlanarPermutation,
    order: [u32; 4],
}

/// Panics (test setup, not the code under test) if `smiles` isn't a
/// square-planar center with exactly 4 explicit neighbors.
fn smiles_to_square_planar_mol_fixture(smiles: &str, z: f64) -> SquarePlanarFixture {
    let mol = parse(smiles).unwrap_or_else(|e| panic!("{smiles} must parse: {e}"));
    let (center, tag) = mol
        .atoms()
        .find_map(|(idx, atom)| match atom.chirality {
            Chirality::SquarePlanar(tag) => Some((idx, tag)),
            _ => None,
        })
        .unwrap_or_else(|| panic!("{smiles} must declare a square-planar center"));
    let order_vec = mol
        .stereo_neighbor_order(center)
        .unwrap_or_else(|| panic!("{smiles}'s center must have a stereo_neighbor_order"))
        .to_vec();
    let order: [u32; 4] = order_vec
        .try_into()
        .unwrap_or_else(|v: Vec<u32>| panic!("fixture must have 4 explicit neighbors, got {v:?}"));

    let neighbor_pts = ideal_neighbor_positions(tag, z, 0.0);
    let mut points = vec![Point3::new(0.0, 0.0, z); mol.atom_count()];
    for (slot, &atom_id) in order.iter().enumerate() {
        points[atom_id as usize] = neighbor_pts[slot];
    }
    SquarePlanarFixture {
        mol,
        conformer: Coords3D { points },
        tag,
        order,
    }
}

fn pt_atom_idx(mol: &chematic_core::Molecule) -> AtomIdx {
    mol.atoms()
        .find(|(_, a)| a.element == chematic_core::Element::PT)
        .map(|(idx, _)| idx)
        .expect("fixture must contain Pt")
}

/// The physical-arrangement-aware differential check the task asks for: not
/// "does the recovered tag's spelling match", but "does the recovered
/// `(tag, neighbor_order)` pair describe the *same physical arrangement* as
/// the original SMILES-declared `(tag, neighbor_order)` pair" -- computed via
/// `chematic_core::remap_square_planar_tag`, the same oracle-validated
/// (144-case table) function `chematic-smiles`'s own canonicalizer uses, not
/// a raw tag equality that would silently assume the MOL round trip
/// preserves `stereo_neighbor_order`'s exact listing order (it need not:
/// `perceive_square_planar_from_3d` derives its order from
/// `Molecule::neighbors`/bond-block order, which is not guaranteed to match
/// the original SMILES parser's branch-encounter order).
fn assert_same_physical_arrangement(
    original: &SquarePlanarFixture,
    read_mol: &chematic_core::Molecule,
    label: &str,
) {
    let read_center = pt_atom_idx(read_mol);
    let recovered_tag = match read_mol.atom(read_center).chirality {
        Chirality::SquarePlanar(tag) => tag,
        other => panic!("{label}: expected SquarePlanar, got {other:?}"),
    };
    let recovered_order: [u32; 4] = read_mol
        .stereo_neighbor_order(read_center)
        .unwrap_or_else(|| panic!("{label}: recovered center has no stereo_neighbor_order"))
        .to_vec()
        .try_into()
        .unwrap_or_else(|v: Vec<u32>| panic!("{label}: expected 4 neighbors, got {v:?}"));

    assert_eq!(
        remap_square_planar_tag(recovered_tag, recovered_order, original.order),
        Some(original.tag),
        "{label}: recovered ({recovered_tag:?}, {recovered_order:?}) does not describe the same \
         physical arrangement as the original ({:?}, {:?})",
        original.tag,
        original.order
    );
}

// ---------------------------------------------------------------------------
// Cisplatin/transplatin differential test (V2000 + V3000)
// ---------------------------------------------------------------------------

/// The same structures and tags as `chematic-smiles/tests/square_planar_stereo.rs`'s
/// `cisplatin_and_transplatin_have_distinct_canonical_identity` -- source of
/// truth for which tag means cis vs trans, reused here rather than
/// re-derived, exactly as the task requires: a differential test against the
/// already-oracle-validated SMILES path.
const CISPLATIN_SMILES: &str = "N->[Pt@SP1](<-N)(Cl)Cl";
const TRANSPLATIN_SMILES: &str = "N->[Pt@SP2](<-N)(Cl)Cl";

#[test]
fn cisplatin_transplatin_v2000_round_trip_recovers_distinct_tags() {
    let cis = smiles_to_square_planar_mol_fixture(CISPLATIN_SMILES, 1.5);
    let trans = smiles_to_square_planar_mol_fixture(TRANSPLATIN_SMILES, 1.5);
    assert_ne!(cis.tag, trans.tag, "fixture setup sanity check");

    let cis_block =
        write_mol_with_conformer_checked(&cis.mol, &MolMetadata::default(), &cis.conformer)
            .expect("cisplatin-shaped geometry must be representable");
    let trans_block =
        write_mol_with_conformer_checked(&trans.mol, &MolMetadata::default(), &trans.conformer)
            .expect("transplatin-shaped geometry must be representable");

    let cis_read = read_mol_with_diagnostics(&cis_block).expect("cisplatin V2000 parses");
    let trans_read = read_mol_with_diagnostics(&trans_block).expect("transplatin V2000 parses");

    // Differential check against the SMILES oracle -- frame-explicit via
    // `remap_square_planar_tag`, not a raw tag comparison (see the helper's
    // doc comment for why that would be unsound).
    assert_same_physical_arrangement(&cis, &cis_read.mol, "cisplatin V2000");
    assert_same_physical_arrangement(&trans, &trans_read.mol, "transplatin V2000");

    assert_ne!(
        cis_read.mol.atom(pt_atom_idx(&cis_read.mol)).chirality,
        trans_read.mol.atom(pt_atom_idx(&trans_read.mol)).chirality,
        "cisplatin and transplatin must stay distinct after a MOL round trip"
    );
    assert!(cis_read.square_planar_diagnostics.is_empty());
    assert!(trans_read.square_planar_diagnostics.is_empty());
}

#[test]
fn cisplatin_transplatin_v3000_round_trip_recovers_distinct_tags() {
    let cis = smiles_to_square_planar_mol_fixture(CISPLATIN_SMILES, 1.5);
    let trans = smiles_to_square_planar_mol_fixture(TRANSPLATIN_SMILES, 1.5);

    let cis_block =
        write_mol_v3000_with_conformer_checked(&cis.mol, &MolMetadata::default(), &cis.conformer)
            .expect("cisplatin-shaped geometry must be representable in V3000");
    let trans_block = write_mol_v3000_with_conformer_checked(
        &trans.mol,
        &MolMetadata::default(),
        &trans.conformer,
    )
    .expect("transplatin-shaped geometry must be representable in V3000");

    let cis_read = read_mol_v3000_with_diagnostics(&cis_block).expect("cisplatin V3000 parses");
    let trans_read =
        read_mol_v3000_with_diagnostics(&trans_block).expect("transplatin V3000 parses");

    assert_same_physical_arrangement(&cis, &cis_read.mol, "cisplatin V3000");
    assert_same_physical_arrangement(&trans, &trans_read.mol, "transplatin V3000");
    assert_ne!(
        cis_read.mol.atom(pt_atom_idx(&cis_read.mol)).chirality,
        trans_read.mol.atom(pt_atom_idx(&trans_read.mol)).chirality
    );

    // Both formats behave identically for this mechanism (checked, not
    // assumed -- RFC §5): the V2000 test above and this V3000 test recover
    // the same physical tags from the same source SMILES.

    // Full SMILES -> MOL -> SMILES closure (V3000 only -- V2000 collapses
    // Dative bonds to plain Single, which would perturb the canonical
    // string for reasons unrelated to stereo; see RFC/module docs).
    assert_eq!(
        canonical_smiles(&cis_read.mol),
        canonical_smiles(&cis.mol),
        "cisplatin: canonical SMILES must survive a V3000 MOL round trip"
    );
    assert_eq!(
        canonical_smiles(&trans_read.mol),
        canonical_smiles(&trans.mol),
        "transplatin: canonical SMILES must survive a V3000 MOL round trip"
    );
}

// ---------------------------------------------------------------------------
// Atom-renumbering invariance
// ---------------------------------------------------------------------------

/// Renumbering must preserve the recovered physical configuration. Per the
/// RFC (§ design notes) and this project's `stereo_geometry.rs`, a pure
/// atom relabeling does not require an explicit `remap_square_planar_tag`
/// call here: `Molecule`'s `stereo_neighbor_order` is itself defined
/// relative to atom ids, so renumbering the atoms and the conformer
/// together (the same ids used consistently on both sides) reproduces an
/// equivalent physical arrangement by construction -- `remap_square_planar_tag`
/// is exercised end to end via `canonical_smiles`'s own canonicalization
/// pipeline in the assertion below, not called directly by this test.
#[test]
fn atom_renumbering_preserves_square_planar_configuration() {
    let SquarePlanarFixture {
        mol,
        conformer: conf,
        tag,
        ..
    } = smiles_to_square_planar_mol_fixture(CISPLATIN_SMILES, 1.5);
    let n = mol.atom_count();

    // Reverse atom order: id `i` in the original becomes id `n-1-i`.
    let mut builder = chematic_core::MoleculeBuilder::new();
    let mut old_to_new = vec![0u32; n];
    for old in (0..n).rev() {
        let new_idx = builder.add_atom(mol.atom(AtomIdx(old as u32)).clone());
        old_to_new[old] = new_idx.0;
    }
    for (_, bond) in mol.bonds() {
        let a = AtomIdx(old_to_new[bond.atom1.0 as usize]);
        let b = AtomIdx(old_to_new[bond.atom2.0 as usize]);
        builder.add_bond(a, b, bond.order).unwrap();
    }
    // Remap stereo_neighbor_order (side channel, not carried by `add_atom`/
    // `Atom` itself) using the same old->new table.
    for (old_idx, atom) in mol.atoms() {
        if let Chirality::SquarePlanar(_) = atom.chirality
            && let Some(order) = mol.stereo_neighbor_order(old_idx)
        {
            let new_order: Vec<u32> = order.iter().map(|&v| old_to_new[v as usize]).collect();
            builder.set_stereo_neighbor_order(AtomIdx(old_to_new[old_idx.0 as usize]), new_order);
        }
    }
    let renumbered = builder.build();

    // Permute the conformer the same way.
    let mut renumbered_points = vec![Point3::new(0.0, 0.0, 0.0); n];
    for old in 0..n {
        renumbered_points[old_to_new[old] as usize] = conf.points[old];
    }
    let renumbered_conf = Coords3D {
        points: renumbered_points,
    };

    // Both orderings must write successfully and read back to the SAME
    // physical arrangement.
    let block_original = write_mol_with_conformer_checked(&mol, &MolMetadata::default(), &conf)
        .expect("original order writes");
    let block_renumbered =
        write_mol_with_conformer_checked(&renumbered, &MolMetadata::default(), &renumbered_conf)
            .expect("reversed order writes");

    let read_original = read_mol_with_diagnostics(&block_original).unwrap();
    let read_renumbered = read_mol_with_diagnostics(&block_renumbered).unwrap();

    assert_eq!(
        read_original
            .mol
            .atom(pt_atom_idx(&read_original.mol))
            .chirality,
        Chirality::SquarePlanar(tag)
    );
    assert_eq!(
        read_renumbered
            .mol
            .atom(pt_atom_idx(&read_renumbered.mol))
            .chirality,
        Chirality::SquarePlanar(tag),
        "atom-reversed molecule must recover the same tag as the original"
    );

    // Exercise `remap_square_planar_tag` end to end: canonical SMILES of
    // both orderings' read-back molecules must agree (same physical
    // identity), not merely "both happen to say SP1/SP2 in their own local
    // stereo_neighbor_order frame".
    assert_eq!(
        canonical_smiles(&read_original.mol),
        canonical_smiles(&read_renumbered.mol),
        "renumbering must not change canonical identity"
    );
}

// ---------------------------------------------------------------------------
// Near-degenerate / distorted geometry robustness
// ---------------------------------------------------------------------------

#[test]
fn moderately_distorted_geometry_still_recovers_the_correct_tag() {
    // 15 degrees of symmetric distortion per trans pair (well inside the
    // ~35-45 degree ambiguity band derived in the RFC/mol2000.rs comments).
    let mol = parse(CISPLATIN_SMILES).unwrap();
    let center_idx = mol
        .atoms()
        .find(|(_, a)| a.element == chematic_core::Element::PT)
        .unwrap()
        .0;
    let tag = match mol.atom(center_idx).chirality {
        Chirality::SquarePlanar(tag) => tag,
        other => panic!("expected square-planar, got {other:?}"),
    };
    let order = mol.stereo_neighbor_order(center_idx).unwrap().to_vec();

    let neighbor_pts = ideal_neighbor_positions(tag, 1.5, 15.0);
    let mut points = vec![Point3::new(0.0, 0.0, 1.5); mol.atom_count()];
    for (slot, &atom_id) in order.iter().enumerate() {
        points[atom_id as usize] = neighbor_pts[slot];
    }
    let conf = Coords3D { points };

    let block = write_mol_with_conformer_checked(&mol, &MolMetadata::default(), &conf)
        .expect("distorted-but-valid geometry must still be representable");
    let read = read_mol_with_diagnostics(&block).unwrap();
    assert_eq!(
        read.mol.atom(pt_atom_idx(&read.mol)).chirality,
        Chirality::SquarePlanar(tag),
        "15-degree distortion must not prevent tag recovery"
    );
}

// ---------------------------------------------------------------------------
// FlatZero: documented round-trip limitation (RFC §10)
// ---------------------------------------------------------------------------

#[test]
fn flat_z_conformer_is_a_documented_round_trip_limitation() {
    let SquarePlanarFixture {
        mol,
        conformer: mut conf,
        ..
    } = smiles_to_square_planar_mol_fixture(CISPLATIN_SMILES, 1.5);
    for p in &mut conf.points {
        p.z = 0.0;
    }

    // The checked writer must refuse (this is exactly the "writer says
    // encoded, reader says nothing" failure this design exists to prevent).
    let err = write_mol_with_conformer_checked(&mol, &MolMetadata::default(), &conf).unwrap_err();
    assert_eq!(
        err.reason,
        UnsupportedStereoReason::WholeMoleculeConformerFlat
    );
    assert!(matches!(
        validate_square_planar_for_write(&mol, Some(&conf), MolFormat::V2000),
        Err(e) if e.reason == UnsupportedStereoReason::WholeMoleculeConformerFlat
    ));

    // And even a hand-built, flat-z MOL block bypassing this crate's own
    // writer entirely (simulating an externally-produced file) does not
    // recover a tag on read -- documenting, not silently mishandling, the
    // limitation.
    let center = pt_atom_idx(&mol);
    let order = mol.stereo_neighbor_order(center).unwrap().to_vec();
    let neighbor_pts = ideal_neighbor_positions(
        match mol.atom(center).chirality {
            Chirality::SquarePlanar(t) => t,
            other => panic!("expected square-planar, got {other:?}"),
        },
        0.0,
        0.0,
    );
    let natoms = mol.atom_count();
    let mut lines = vec![
        String::new(),
        "     hand    3D".to_string(),
        String::new(),
        format!(
            "{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000",
            natoms,
            mol.bond_count()
        ),
    ];
    let mut xyz = vec![(0.0, 0.0, 0.0); natoms];
    for (slot, &atom_id) in order.iter().enumerate() {
        let p = neighbor_pts[slot];
        xyz[atom_id as usize] = (p.x, p.y, p.z);
    }
    for (idx, atom) in mol.atoms() {
        let (x, y, z) = xyz[idx.0 as usize];
        lines.push(format!(
            "{:>10.4}{:>10.4}{:>10.4} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0",
            x,
            y,
            z,
            atom.element.symbol()
        ));
    }
    for (_, bond) in mol.bonds() {
        lines.push(format!(
            "{:>3}{:>3}{:>3}{:>3}",
            bond.atom1.0 + 1,
            bond.atom2.0 + 1,
            1,
            0
        ));
    }
    lines.push("M  END".to_string());
    let hand_block = lines.join("\n") + "\n";

    let read = read_mol_with_diagnostics(&hand_block).expect("hand-built flat block still parses");
    assert_eq!(
        read.mol.atom(pt_atom_idx(&read.mol)).chirality,
        Chirality::None,
        "flat-z (z=0) square-planar geometry is a known, documented non-recoverable case"
    );
    assert!(read.conformer.is_none());
}

// ---------------------------------------------------------------------------
// Malformed / adversarial input never panics
// ---------------------------------------------------------------------------

#[test]
fn degenerate_and_truncated_conformers_never_panic() {
    let SquarePlanarFixture {
        mol,
        conformer: mut conf,
        ..
    } = smiles_to_square_planar_mol_fixture(CISPLATIN_SMILES, 1.5);

    // Every atom coincident at the same point (degenerate bond vectors).
    let center_pt = conf.points[pt_atom_idx(&mol).0 as usize];
    conf.points.fill(center_pt);
    let _ = write_mol_with_conformer_checked(&mol, &MolMetadata::default(), &conf);
    let _ = validate_square_planar_for_write(&mol, Some(&conf), MolFormat::V3000);

    // Truncated conformer (fewer points than atoms) -- must fail closed via
    // `MissingConformerPosition`, never index-panic.
    let short_conf = Coords3D {
        points: vec![Point3::new(0.0, 0.0, 1.5)],
    };
    let err =
        write_mol_with_conformer_checked(&mol, &MolMetadata::default(), &short_conf).unwrap_err();
    assert_eq!(
        err.reason,
        UnsupportedStereoReason::MissingConformerPosition
    );

    // Empty conformer entirely.
    let empty_conf = Coords3D { points: vec![] };
    let _ = write_mol_with_conformer_checked(&mol, &MolMetadata::default(), &empty_conf);

    // A MOL block whose Pt-adjacent atom coordinates are all identical
    // (degenerate) must parse and never panic, regardless of what stereo
    // (if any) gets perceived.
    let degenerate_block = "\
degenerate
  chematic          3D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    1.5000 Pt  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    1.5000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    1.5000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    1.5000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    1.5000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  1  5  1  0
M  END
";
    let read = read_mol_with_diagnostics(degenerate_block).expect("degenerate block still parses");
    assert_eq!(read.mol.atom(AtomIdx(0)).chirality, Chirality::None);
}

// ---------------------------------------------------------------------------
// Additional fixture coverage (per human-reviewer follow-up on this PR):
// SP3 + all-distinct ligands, a duplicate-ligand composition, degree-3/5
// negative controls, a malformed V3000 coordinate block, write-read-write
// stability, and the real (`$$$$`-terminated) SDF multi-record path.
// ---------------------------------------------------------------------------

/// Build a Pt-centered molecule with `ligand_elements.len()` explicit
/// single-bonded neighbors, no SMILES involved -- lets these tests probe
/// coordination numbers (3, 5) SMILES parsing wouldn't naturally produce,
/// and ligand compositions (4 distinct elements, 3+1 duplicate) SMILES
/// fixtures wouldn't cleanly express either.
fn build_pt_with_ligands(
    ligand_elements: &[Element],
) -> (chematic_core::Molecule, AtomIdx, Vec<AtomIdx>) {
    let mut b = MoleculeBuilder::new();
    let center = b.add_atom(Atom::new(Element::PT));
    let mut ligands = Vec::new();
    for &el in ligand_elements {
        let l = b.add_atom(Atom::new(el));
        b.add_bond(center, l, BondOrder::Single).unwrap();
        ligands.push(l);
    }
    (b.build(), center, ligands)
}

/// Declare `tag` on `center` with `ligands` as its neighbor order, and build
/// a matching `Coords3D` conformer (nonzero z -- see module docs).
fn declare_and_place(
    mol: &mut chematic_core::Molecule,
    center: AtomIdx,
    ligands: &[AtomIdx],
    tag: SquarePlanarPermutation,
) -> ([u32; 4], Coords3D) {
    let order: [u32; 4] = ligands
        .iter()
        .map(|a| a.0)
        .collect::<Vec<_>>()
        .try_into()
        .expect("exactly 4 ligands");
    mol.set_chirality(center, Chirality::SquarePlanar(tag));
    mol.set_stereo_neighbor_order(center, order.to_vec());

    let neighbor_pts = ideal_neighbor_positions(tag, 1.5, 0.0);
    let mut points = vec![Point3::new(0.0, 0.0, 1.5); mol.atom_count()];
    for (slot, &atom_id) in order.iter().enumerate() {
        points[atom_id as usize] = neighbor_pts[slot];
    }
    (order, Coords3D { points })
}

#[test]
fn all_four_ligands_distinct_round_trips_all_three_tags() {
    for tag in [
        SquarePlanarPermutation::SP1,
        SquarePlanarPermutation::SP2,
        SquarePlanarPermutation::SP3,
    ] {
        let (mut mol, center, ligands) =
            build_pt_with_ligands(&[Element::F, Element::CL, Element::BR, Element::I]);
        let (order, conf) = declare_and_place(&mut mol, center, &ligands, tag);

        let block = write_mol_with_conformer_checked(&mol, &MolMetadata::default(), &conf)
            .unwrap_or_else(|e| panic!("{tag:?}: must write: {e}"));
        let read = read_mol_with_diagnostics(&block).unwrap_or_else(|e| panic!("{tag:?}: {e}"));
        let read_center = pt_atom_idx(&read.mol);
        let recovered_tag = match read.mol.atom(read_center).chirality {
            Chirality::SquarePlanar(t) => t,
            other => panic!("{tag:?}: expected SquarePlanar, got {other:?}"),
        };
        let recovered_order: [u32; 4] = read
            .mol
            .stereo_neighbor_order(read_center)
            .unwrap()
            .to_vec()
            .try_into()
            .unwrap();
        assert_eq!(
            remap_square_planar_tag(recovered_tag, recovered_order, order),
            Some(tag),
            "{tag:?}: round trip must recover the same physical arrangement"
        );
    }

    // The 3 tags, applied to the SAME distinct-ligand composition, must
    // remain pairwise distinct after a round trip (not just individually
    // "recoverable") -- this is what makes them 3 genuinely different
    // stereoisomers for a 4-distinct-ligand square-planar complex, unlike
    // cisplatin/transplatin's 2xCl+2xN composition where SP1 and SP3
    // collapse to the same physical species (see
    // `chematic-core/src/stereo_geometry.rs`'s own duplicate-ligand note).
    let mut recovered_tags = Vec::new();
    for tag in [
        SquarePlanarPermutation::SP1,
        SquarePlanarPermutation::SP2,
        SquarePlanarPermutation::SP3,
    ] {
        let (mut mol, center, ligands) =
            build_pt_with_ligands(&[Element::F, Element::CL, Element::BR, Element::I]);
        let (_order, conf) = declare_and_place(&mut mol, center, &ligands, tag);
        let block = write_mol_with_conformer_checked(&mol, &MolMetadata::default(), &conf).unwrap();
        let read = read_mol_with_diagnostics(&block).unwrap();
        recovered_tags.push(canonical_smiles(&read.mol));
    }
    assert_ne!(recovered_tags[0], recovered_tags[1]);
    assert_ne!(recovered_tags[0], recovered_tags[2]);
    assert_ne!(recovered_tags[1], recovered_tags[2]);
}

#[test]
fn three_plus_one_duplicate_ligand_composition_round_trips() {
    // 3xCl + 1xN -- not itself a stereogenic arrangement chemically (no
    // isomerism when 3 of 4 ligands are identical), but the reader/writer
    // machinery must still round-trip whatever tag is declared without
    // crashing or silently corrupting it, the same "duplicate chemistry,
    // distinct ids" property `chematic-core`'s own unit tests check at the
    // geometry-module level.
    let (mut mol, center, ligands) =
        build_pt_with_ligands(&[Element::CL, Element::CL, Element::CL, Element::N]);
    let (order, conf) = declare_and_place(&mut mol, center, &ligands, SquarePlanarPermutation::SP1);

    let block = write_mol_with_conformer_checked(&mol, &MolMetadata::default(), &conf)
        .expect("duplicate-ligand geometry must still be representable");
    let read = read_mol_with_diagnostics(&block).expect("parses");
    let read_center = pt_atom_idx(&read.mol);
    let recovered_tag = match read.mol.atom(read_center).chirality {
        Chirality::SquarePlanar(t) => t,
        other => panic!("expected SquarePlanar, got {other:?}"),
    };
    let recovered_order: [u32; 4] = read
        .mol
        .stereo_neighbor_order(read_center)
        .unwrap()
        .to_vec()
        .try_into()
        .unwrap();
    assert_eq!(
        remap_square_planar_tag(recovered_tag, recovered_order, order),
        Some(SquarePlanarPermutation::SP1)
    );
}

#[test]
fn three_coordinate_center_is_never_perceived_as_square_planar() {
    let (mol, center, ligands) = build_pt_with_ligands(&[Element::CL, Element::CL, Element::N]);
    // Coplanar (real, nonzero-z, 3-around-a-center) positions -- the
    // perception candidate filter requires exactly 4 neighbors, so this
    // must never even be classified, regardless of how planar it looks.
    let mut points = vec![Point3::new(0.0, 0.0, 1.5); mol.atom_count()];
    let at = |deg: f64| -> Point3 {
        let r = deg.to_radians();
        Point3::new(1.5 * r.cos(), 1.5 * r.sin(), 1.5)
    };
    for (i, &l) in ligands.iter().enumerate() {
        points[l.0 as usize] = at(120.0 * i as f64);
    }
    let conf = Coords3D { points };

    let block = write_mol_with_conformer(&mol, &MolMetadata::default(), &conf);
    let read = read_mol_with_diagnostics(&block).expect("parses");
    assert_eq!(read.mol.atom(center).chirality, Chirality::None);
    assert!(read.square_planar_diagnostics.is_empty());
}

#[test]
fn five_coordinate_center_is_never_perceived_as_square_planar() {
    let (mol, center, ligands) = build_pt_with_ligands(&[
        Element::CL,
        Element::CL,
        Element::N,
        Element::N,
        Element::BR,
    ]);
    let mut points = vec![Point3::new(0.0, 0.0, 1.5); mol.atom_count()];
    let at = |deg: f64| -> Point3 {
        let r = deg.to_radians();
        Point3::new(1.5 * r.cos(), 1.5 * r.sin(), 1.5)
    };
    for (i, &l) in ligands.iter().enumerate() {
        points[l.0 as usize] = at(72.0 * i as f64);
    }
    let conf = Coords3D { points };

    let block = write_mol_with_conformer(&mol, &MolMetadata::default(), &conf);
    let read = read_mol_with_diagnostics(&block).expect("parses");
    assert_eq!(read.mol.atom(center).chirality, Chirality::None);
    assert!(read.square_planar_diagnostics.is_empty());
}

#[test]
fn truncated_v3000_coordinate_block_is_a_typed_error_not_a_panic() {
    // Same shape as `conformer_3d_io.rs`'s own `v3000_garbled_z_is_a_typed_error`,
    // applied to a square-planar-eligible (undefined-valence) center: the
    // malformed z field must fail during ordinary V3000 atom-line parsing,
    // before this PR's perception code ever runs -- confirming that new
    // code path doesn't need its own defense here, and that the existing
    // one still holds for this element class.
    let bad = "\
bad
  chematic

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 Pt 0.0 0.0 0.0 0
M  V30 2 Cl -1.0607 -1.0607 garbage 0
M  V30 3 Cl -1.0607 1.0607 1.5 0
M  V30 4 N 1.0607 1.0607 1.5 0
M  V30 5 N 1.0607 -1.0607 1.5 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 3
M  V30 3 1 1 4
M  V30 4 1 1 5
M  V30 END BOND
M  V30 END CTAB
M  END
";
    let result = read_mol_v3000_with_diagnostics(bad);
    assert!(matches!(
        result,
        Err(chematic_mol::MolParseError::InvalidAtomLine { .. })
    ));
}

#[test]
fn write_read_write_is_stable() {
    let cis = smiles_to_square_planar_mol_fixture(CISPLATIN_SMILES, 1.5);
    let block1 =
        write_mol_with_conformer_checked(&cis.mol, &MolMetadata::default(), &cis.conformer)
            .expect("first write");
    let read1 = read_mol_with_diagnostics(&block1).expect("first read");
    let conformer1 = read1.conformer.clone().expect("first read has a conformer");
    let block2 = write_mol_with_conformer_checked(&read1.mol, &MolMetadata::default(), &conformer1)
        .expect("second write");

    assert_eq!(
        block1, block2,
        "writing, reading back, and writing again must reproduce byte-identical output"
    );
}

#[test]
fn sdf_multi_record_round_trip_no_cross_contamination() {
    let cis = smiles_to_square_planar_mol_fixture(CISPLATIN_SMILES, 1.5);
    let trans = smiles_to_square_planar_mol_fixture(TRANSPLATIN_SMILES, 1.5);

    let mut cis_props = std::collections::HashMap::new();
    cis_props.insert("Name".to_string(), "cisplatin".to_string());
    let mut trans_props = std::collections::HashMap::new();
    trans_props.insert("Name".to_string(), "transplatin".to_string());

    let cis_record = write_sdf_record_with_conformer_checked(
        &cis.mol,
        &MolMetadata::default(),
        &cis.conformer,
        &cis_props,
    )
    .expect("cisplatin SD record must write");
    let trans_record = write_sdf_record_with_conformer_checked(
        &trans.mol,
        &MolMetadata::default(),
        &trans.conformer,
        &trans_props,
    )
    .expect("transplatin SD record must write");

    let sdf = format!("{cis_record}{trans_record}");
    let records: Vec<_> = SdfRecordReader::new(&sdf)
        .collect::<Result<Vec<_>, _>>()
        .expect("both records must parse");
    assert_eq!(
        records.len(),
        2,
        "must be exactly 2 records, not merged/split"
    );

    let rec0_tag = match records[0].mol.atom(pt_atom_idx(&records[0].mol)).chirality {
        Chirality::SquarePlanar(t) => t,
        other => panic!("record 0: expected SquarePlanar, got {other:?}"),
    };
    let rec1_tag = match records[1].mol.atom(pt_atom_idx(&records[1].mol)).chirality {
        Chirality::SquarePlanar(t) => t,
        other => panic!("record 1: expected SquarePlanar, got {other:?}"),
    };

    // Physical-arrangement-aware check, same as the top-level differential
    // tests -- not a raw tag comparison.
    let rec0_order: [u32; 4] = records[0]
        .mol
        .stereo_neighbor_order(pt_atom_idx(&records[0].mol))
        .unwrap()
        .to_vec()
        .try_into()
        .unwrap();
    let rec1_order: [u32; 4] = records[1]
        .mol
        .stereo_neighbor_order(pt_atom_idx(&records[1].mol))
        .unwrap()
        .to_vec()
        .try_into()
        .unwrap();
    assert_eq!(
        remap_square_planar_tag(rec0_tag, rec0_order, cis.order),
        Some(cis.tag),
        "record 0 must be cisplatin's arrangement"
    );
    assert_eq!(
        remap_square_planar_tag(rec1_tag, rec1_order, trans.order),
        Some(trans.tag),
        "record 1 must be transplatin's arrangement, not cisplatin's"
    );

    // Properties must attach to the correct record, not cross-contaminate.
    assert_eq!(
        records[0].properties.get("Name").map(|s| s.as_str()),
        Some("cisplatin")
    );
    assert_eq!(
        records[1].properties.get("Name").map(|s| s.as_str()),
        Some("transplatin")
    );
}
