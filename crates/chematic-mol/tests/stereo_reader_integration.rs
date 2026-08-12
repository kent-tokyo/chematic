//! Integration tests for 2D wedge/hash stereo perception wired into the
//! MOL V2000/V3000 readers and writers (see
//! `docs/rfcs/stereo2d_reader_integration_rfc.md` for the design background).
//!
//! Unit-level coverage of the parity math itself lives in
//! `chematic-perception`'s `stereo2d_local` tests; these tests cover the
//! cross-cutting invariants that only exist once a real reader/writer is in
//! the loop: V2000 vs. V3000 agreement, renumbering/reflection/rotation
//! invariance, and round-trip losslessness.

use chematic_core::{Atom, AtomIdx, BondOrder, Chirality, Element, MoleculeBuilder};
use chematic_mol::mol2000::{MolMetadata, write_mol_with_coords};
use chematic_mol::{parse_mol_v3000_with_coords, read_mol_with_diagnostics, write_mol_v3000};

/// Asymmetric, non-degenerate 4-position layout (matches
/// `chematic-perception`'s own `quad_positions()`) so no accidental
/// coplanarity or symmetry sneaks into these fixtures.
fn quad_positions() -> [(f64, f64); 4] {
    [(-1.0, 0.4), (0.9, 0.7), (-0.5, -1.1), (0.8, -0.6)]
}

/// Build a CHFClBr molecule with a wedge on C-F, and its matching coords.
fn chfclbr_wedge() -> (chematic_core::Molecule, Vec<(f64, f64)>) {
    let mut b = MoleculeBuilder::new();
    let c = b.add_atom(Atom::new(Element::C));
    let f = b.add_atom(Atom::new(Element::F));
    let cl = b.add_atom(Atom::new(Element::CL));
    let br = b.add_atom(Atom::new(Element::BR));
    let i = b.add_atom(Atom::new(Element::I));
    b.add_bond(c, f, BondOrder::Up).unwrap();
    b.add_bond(c, cl, BondOrder::Single).unwrap();
    b.add_bond(c, br, BondOrder::Single).unwrap();
    b.add_bond(c, i, BondOrder::Single).unwrap();
    let quad = quad_positions();
    let coords = vec![(0.0, 0.0), quad[0], quad[1], quad[2], quad[3]];
    (b.build(), coords)
}

#[test]
fn v2000_v3000_same_drawing_yield_identical_parity() {
    let (mol, coords) = chfclbr_wedge();

    let v2000_block = write_mol_with_coords(&mol, &MolMetadata::default(), &coords);
    let v3000_block = write_mol_v3000(&mol, &MolMetadata::default(), &coords);

    let (mol_v2000, _, _) = chematic_mol::parse_mol_with_coords(&v2000_block).expect("v2000");
    let (mol_v3000, _, _) = parse_mol_v3000_with_coords(&v3000_block).expect("v3000");

    let center = AtomIdx(0);
    assert_ne!(mol_v2000.atom(center).chirality, Chirality::None);
    assert_eq!(
        mol_v2000.atom(center).chirality,
        mol_v3000.atom(center).chirality
    );
    assert_eq!(
        mol_v2000.stereo_neighbor_order(center),
        mol_v3000.stereo_neighbor_order(center)
    );
}

#[test]
fn v3000_either_cfg_does_not_produce_a_defined_wedge() {
    // CFG=2 ("either"/unspecified) must not be treated as a definite wedge --
    // same policy as V2000's own code-4 handling.
    let block = "\
either
  chematic

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.0 0.4 0 0
M  V30 2 F 0.9 0.7 0 0
M  V30 3 Cl -0.5 -1.1 0 0
M  V30 4 Br 0.8 -0.6 0 0
M  V30 5 I 0.0 0.0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2 CFG=2
M  V30 2 1 1 3
M  V30 3 1 1 4
M  V30 4 1 1 5
M  V30 END BOND
M  V30 END CTAB
M  END
";
    let (mol, _, _) = parse_mol_v3000_with_coords(block).expect("v3000 parse");
    let bond = mol.bond_between(AtomIdx(0), AtomIdx(1)).unwrap().1;
    assert_eq!(bond.order, BondOrder::Single);
}

#[test]
fn atom_renumbering_preserves_physical_stereo() {
    // Same physical molecule (same real substituent positions, same wedge on
    // the same physical C-F bond), but atoms are added to the builder in the
    // opposite order. Renumbering changes stereo_neighbor_order's raw
    // indices by construction, so compare via canonical SMILES instead of
    // raw AtomIdx-relative chirality.
    let (mol_a, coords_a) = chfclbr_wedge();
    let mut mol_a = mol_a;
    chematic_perception::apply_local_parity_from_wedges(&mut mol_a, &coords_a);

    let mut b = MoleculeBuilder::new();
    let c = b.add_atom(Atom::new(Element::C));
    let i = b.add_atom(Atom::new(Element::I));
    let br = b.add_atom(Atom::new(Element::BR));
    let cl = b.add_atom(Atom::new(Element::CL));
    let f = b.add_atom(Atom::new(Element::F));
    b.add_bond(c, i, BondOrder::Single).unwrap();
    b.add_bond(c, br, BondOrder::Single).unwrap();
    b.add_bond(c, cl, BondOrder::Single).unwrap();
    b.add_bond(c, f, BondOrder::Up).unwrap();
    let quad = quad_positions();
    // coords[i] must match atom index i's real position: C, I, Br, Cl, F.
    let coords_b = vec![(0.0, 0.0), quad[3], quad[2], quad[1], quad[0]];
    let mut mol_b = b.build();
    chematic_perception::apply_local_parity_from_wedges(&mut mol_b, &coords_b);

    assert_ne!(mol_a.atom(AtomIdx(0)).chirality, Chirality::None);
    assert_ne!(mol_b.atom(AtomIdx(0)).chirality, Chirality::None);
    assert_eq!(
        chematic_smiles::canonical_smiles(&mol_a),
        chematic_smiles::canonical_smiles(&mol_b)
    );
}

#[test]
fn bond_declaration_order_reversal_preserves_parity() {
    // MOL-level analogue of chematic-perception's own
    // `bond_atom_order_inversion_flips_chirality`... except that test name
    // is about a DIFFERENT case (odd permutation of a 3-bond block, which
    // DOES flip). Here we reverse the FULL bond block for a 4-neighbor
    // center (an even permutation), which must preserve the sign -- same
    // physical molecule, same wedge, bonds just declared in the file in the
    // opposite order.
    let (mol_a, coords) = chfclbr_wedge();
    let mut mol_a = mol_a;
    chematic_perception::apply_local_parity_from_wedges(&mut mol_a, &coords);

    let mut b = MoleculeBuilder::new();
    let c = b.add_atom(Atom::new(Element::C));
    let f = b.add_atom(Atom::new(Element::F));
    let cl = b.add_atom(Atom::new(Element::CL));
    let br = b.add_atom(Atom::new(Element::BR));
    let i = b.add_atom(Atom::new(Element::I));
    // Same atoms, same coords, bonds declared in reverse (I, Br, Cl, F).
    b.add_bond(c, i, BondOrder::Single).unwrap();
    b.add_bond(c, br, BondOrder::Single).unwrap();
    b.add_bond(c, cl, BondOrder::Single).unwrap();
    b.add_bond(c, f, BondOrder::Up).unwrap();
    let mut mol_b = b.build();
    chematic_perception::apply_local_parity_from_wedges(&mut mol_b, &coords);

    assert_eq!(
        mol_a.atom(AtomIdx(0)).chirality,
        mol_b.atom(AtomIdx(0)).chirality
    );
}

#[test]
fn reflection_flips_parity() {
    let (mol, coords) = chfclbr_wedge();
    let mut mol_orig = mol.clone();
    chematic_perception::apply_local_parity_from_wedges(&mut mol_orig, &coords);

    let reflected: Vec<(f64, f64)> = coords.iter().map(|&(x, y)| (x, -y)).collect();
    let mut mol_reflected = mol;
    chematic_perception::apply_local_parity_from_wedges(&mut mol_reflected, &reflected);

    let center = AtomIdx(0);
    assert_ne!(mol_orig.atom(center).chirality, Chirality::None);
    assert_ne!(mol_reflected.atom(center).chirality, Chirality::None);
    assert_ne!(
        mol_orig.atom(center).chirality,
        mol_reflected.atom(center).chirality
    );
}

#[test]
fn rotation_and_translation_preserve_parity() {
    let (mol, coords) = chfclbr_wedge();
    let mut mol_orig = mol.clone();
    chematic_perception::apply_local_parity_from_wedges(&mut mol_orig, &coords);

    // 90-degree rotation + translation: (x, y) -> (-y + 5, x + 5).
    let transformed: Vec<(f64, f64)> = coords.iter().map(|&(x, y)| (-y + 5.0, x + 5.0)).collect();
    let mut mol_transformed = mol;
    chematic_perception::apply_local_parity_from_wedges(&mut mol_transformed, &transformed);

    let center = AtomIdx(0);
    assert_ne!(mol_orig.atom(center).chirality, Chirality::None);
    assert_eq!(
        mol_orig.atom(center).chirality,
        mol_transformed.atom(center).chirality
    );
}

#[test]
fn v2000_writer_roundtrip_preserves_chirality() {
    let (mol, coords) = chfclbr_wedge();
    let report = read_mol_with_diagnostics(&write_mol_with_coords(
        &mol,
        &MolMetadata::default(),
        &coords,
    ))
    .expect("parse original");
    assert!(report.stereo_diagnostics.is_empty());
    assert_ne!(report.mol.atom(AtomIdx(0)).chirality, Chirality::None);

    let rewritten = write_mol_with_coords(&report.mol, &MolMetadata::default(), &report.coords);
    let reparsed = read_mol_with_diagnostics(&rewritten).expect("re-parse");

    assert_eq!(
        report.mol.atom(AtomIdx(0)).chirality,
        reparsed.mol.atom(AtomIdx(0)).chirality
    );
    assert_eq!(
        report.mol.stereo_neighbor_order(AtomIdx(0)),
        reparsed.mol.stereo_neighbor_order(AtomIdx(0))
    );
}

#[test]
fn v3000_writer_roundtrip_preserves_chirality() {
    let (mol, coords) = chfclbr_wedge();
    let block = write_mol_v3000(&mol, &MolMetadata::default(), &coords);
    let report = parse_mol_v3000_with_coords(&block).expect("parse original");
    let (mol1, _, coords1) = report;
    assert_ne!(mol1.atom(AtomIdx(0)).chirality, Chirality::None);

    let rewritten = write_mol_v3000(&mol1, &MolMetadata::default(), &coords1);
    let (mol2, _, _) = parse_mol_v3000_with_coords(&rewritten).expect("re-parse");

    assert_eq!(
        mol1.atom(AtomIdx(0)).chirality,
        mol2.atom(AtomIdx(0)).chirality
    );
}

#[test]
fn mol_sourced_wedge_survives_to_canonical_smiles_without_meaningless_slash() {
    let (mol, coords) = chfclbr_wedge();
    let block = write_mol_with_coords(&mol, &MolMetadata::default(), &coords);
    let report = read_mol_with_diagnostics(&block).expect("parse");
    assert_ne!(report.mol.atom(AtomIdx(0)).chirality, Chirality::None);

    let smiles = chematic_smiles::canonical_smiles(&report.mol);
    assert!(
        smiles.contains('@'),
        "recovered chirality must survive to canonical SMILES: {smiles}"
    );
    assert!(
        !smiles.contains('/') && !smiles.contains('\\'),
        "a wedge with no adjacent double bond must never emit a directional \
         token: {smiles}"
    );
}

#[test]
fn contradictory_wedge_from_mol_yields_no_chirality_and_a_diagnostic() {
    let mut b = MoleculeBuilder::new();
    let c = b.add_atom(Atom::new(Element::C));
    let f = b.add_atom(Atom::new(Element::F));
    let cl = b.add_atom(Atom::new(Element::CL));
    let br = b.add_atom(Atom::new(Element::BR));
    let i = b.add_atom(Atom::new(Element::I));
    b.add_bond(c, f, BondOrder::Up).unwrap();
    b.add_bond(c, cl, BondOrder::Up).unwrap();
    b.add_bond(c, br, BondOrder::Single).unwrap();
    b.add_bond(c, i, BondOrder::Single).unwrap();
    let mol = b.build();
    let quad = quad_positions();
    let coords = vec![(0.0, 0.0), quad[0], quad[1], quad[2], quad[3]];
    let block = write_mol_with_coords(&mol, &MolMetadata::default(), &coords);

    let report = read_mol_with_diagnostics(&block).expect("parse");
    assert_eq!(report.mol.atom(AtomIdx(0)).chirality, Chirality::None);
    assert_eq!(report.stereo_diagnostics.len(), 1);
    assert_eq!(
        report.stereo_diagnostics[0].reason,
        chematic_perception::StereoRejectionReason::ContradictoryWedges
    );

    // The molecule as a whole must still be readable -- a malformed
    // stereocenter must never fail the entire parse.
    assert_eq!(report.mol.atom_count(), 5);
}
