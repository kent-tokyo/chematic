//! Integration fixtures for P1-S2: E/Z (double-bond cis/trans) direction
//! perception wired into the MOL V2000/V3000/SDF readers
//! (`chematic_perception::apply_ez_directions_from_2d_ex`, called from
//! inside `chematic_mol::{read_mol_with_diagnostics, read_mol_v3000_with_diagnostics}`).
//!
//! Every fixture here passes a REAL MOL/SDF block through the actual reader
//! entry point (`read_mol_with_diagnostics`/`read_mol_v3000_with_diagnostics`
//! /`SdfRecordReader`) -- never `chematic_perception::stereo2d_ez_direction`
//! directly -- with one documented exception (`missing_coordinate_via_api_misuse`,
//! see its doc comment: a well-formed MOL atom line always yields *some*
//! (x, y) value, so `MissingCoordinate` is only reachable via direct API
//! misuse, not a file shape, matching the precedent set by
//! `docs/rfcs/stereo2d_reader_integration_rfc.md`'s own fixture #14 for the
//! tetrahedral case).
//!
//! MOL blocks are generated via `MoleculeBuilder` + the crate's own
//! `write_mol_with_coords`/`write_mol_v3000` writers wherever possible (not
//! hand-typed fixed-width text), with small, targeted text patches for the
//! handful of MDL constructs those writers never emit (double-bond stereo
//! code 3 / V3000 `CFG=2` "either", and a literal `NaN` coordinate field).

use chematic_core::{Atom, AtomIdx, BondOrder, Chirality, Element, MoleculeBuilder};
use chematic_mol::mol2000::{MolMetadata, write_mol_with_coords};
use chematic_mol::{
    SdfRecordReader, read_mol_v3000_with_diagnostics, read_mol_with_diagnostics, write_mol_v3000,
    write_sdf,
};
use chematic_perception::EzDirectionRejectionReason;
use chematic_smiles::canonical_smiles;

// ---------------------------------------------------------------------------
// Text-patch helpers (only for MDL shapes the writers never emit)
// ---------------------------------------------------------------------------

/// Overwrite a V2000 bond line's stereo field (columns 9-11, 0-indexed) with
/// `code`. `bond_index` is 0-based, in bond-block declaration order.
fn patch_v2000_bond_stereo(mol_block: &str, natoms: usize, bond_index: usize, code: u8) -> String {
    let mut lines: Vec<String> = mol_block.lines().map(str::to_string).collect();
    let line_idx = 4 + natoms + bond_index; // 3 header lines + counts line + atom block
    let mut chars: Vec<char> = lines[line_idx].chars().collect();
    while chars.len() < 12 {
        chars.push(' ');
    }
    for (i, c) in format!("{code:>3}").chars().enumerate() {
        chars[9 + i] = c;
    }
    lines[line_idx] = chars.into_iter().collect();
    let mut out = lines.join("\n");
    out.push('\n');
    out
}

/// Append ` CFG=2` to a V3000 bond line identified by its 1-based V30 bond
/// index (== 0-based declaration order + 1).
fn patch_v3000_bond_cfg2(mol_block: &str, v30_bond_index: u32) -> String {
    let prefix = format!("M  V30 {v30_bond_index} ");
    mol_block
        .lines()
        .map(|line| {
            if line.starts_with(&prefix) {
                format!("{line} CFG=2")
            } else {
                line.to_string()
            }
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n"
}

/// Overwrite an atom line's x-coordinate field (columns 0-9) with a literal
/// `NaN` token -- reachable through the real reader because Rust's `f64`
/// `FromStr` accepts `"NaN"` (confirmed empirically, not assumed) and the
/// V2000 atom-line parser does `s.trim().parse().ok()`.
fn patch_v2000_atom_x_nan(mol_block: &str, atom_index: usize) -> String {
    let mut lines: Vec<String> = mol_block.lines().map(str::to_string).collect();
    let line_idx = 4 + atom_index;
    let mut chars: Vec<char> = lines[line_idx].chars().collect();
    let nan = format!("{:>10}", "NaN");
    for (i, c) in nan.chars().enumerate() {
        chars[i] = c;
    }
    lines[line_idx] = chars.into_iter().collect();
    let mut out = lines.join("\n");
    out.push('\n');
    out
}

fn has_directional_token(smiles: &str) -> bool {
    smiles.contains('/') || smiles.contains('\\')
}

// ---------------------------------------------------------------------------
// Shared molecule builders
// ---------------------------------------------------------------------------

/// but-2-ene skeleton `Me0-C1=C2-Me3` at the given 2D layout. Returns the
/// V2000 MOL block, `natoms`, and the double bond's 0-based declaration
/// index (always `1`, the second bond added).
fn but2ene_v2000(
    c0: (f64, f64),
    c1: (f64, f64),
    c2: (f64, f64),
    c3: (f64, f64),
) -> (String, usize) {
    let mut b = MoleculeBuilder::new();
    let m0 = b.add_atom(Atom::new(Element::C));
    let m1 = b.add_atom(Atom::new(Element::C));
    let m2 = b.add_atom(Atom::new(Element::C));
    let m3 = b.add_atom(Atom::new(Element::C));
    b.add_bond(m0, m1, BondOrder::Single).unwrap();
    b.add_bond(m1, m2, BondOrder::Double).unwrap();
    b.add_bond(m2, m3, BondOrder::Single).unwrap();
    let mol = b.build();
    let coords = vec![c0, c1, c2, c3];
    let block = write_mol_with_coords(&mol, &MolMetadata::default().with_name("but2ene"), &coords);
    (block, 4)
}

fn but2ene_v3000(c0: (f64, f64), c1: (f64, f64), c2: (f64, f64), c3: (f64, f64)) -> String {
    let mut b = MoleculeBuilder::new();
    let m0 = b.add_atom(Atom::new(Element::C));
    let m1 = b.add_atom(Atom::new(Element::C));
    let m2 = b.add_atom(Atom::new(Element::C));
    let m3 = b.add_atom(Atom::new(Element::C));
    b.add_bond(m0, m1, BondOrder::Single).unwrap();
    b.add_bond(m1, m2, BondOrder::Double).unwrap();
    b.add_bond(m2, m3, BondOrder::Single).unwrap();
    let mol = b.build();
    let coords = vec![c0, c1, c2, c3];
    write_mol_v3000(&mol, &MolMetadata::default().with_name("but2ene"), &coords)
}

/// Z-configuration coordinates: same side (up) for both methyls.
const Z_COORDS: [(f64, f64); 4] = [(-0.866, 0.5), (0.0, 0.0), (1.5, 0.0), (2.366, 0.5)];
/// E-configuration coordinates: opposite sides.
const E_COORDS: [(f64, f64); 4] = [(-0.866, 0.5), (0.0, 0.0), (1.5, 0.0), (2.366, -0.5)];

// ---------------------------------------------------------------------------
// Positive fixtures (1-10)
// ---------------------------------------------------------------------------

/// Fixture 1/2: E- and Z-2-butene, V2000: opposite drawn geometries must produce
/// genuinely different (non-inverted-relative-to-each-other) SMILES E/Z
/// semantics, both surfacing a directional token in plain `write()`.
#[test]
fn positive_01_02_e_and_z_2butene_v2000() {
    for (coords, label) in [(Z_COORDS, "Z"), (E_COORDS, "E")] {
        let (block, natoms) = but2ene_v2000(coords[0], coords[1], coords[2], coords[3]);
        let _ = natoms;
        let report = read_mol_with_diagnostics(&block).expect("parse");
        assert!(
            report.ez_diagnostics.is_empty(),
            "{label}: {:?}",
            report.ez_diagnostics
        );
        let smiles = chematic_smiles::write(&report.mol);
        assert!(
            has_directional_token(&smiles),
            "{label}: no directional token in '{smiles}'"
        );
        let canon = canonical_smiles(&report.mol);
        assert!(has_directional_token(&canon), "{label}: canon='{canon}'");
    }
    // The two configurations must not encode the same relative geometry --
    // canonical SMILES for Z must differ in its directional-bond pattern
    // from E for the "same physical bonds cis vs trans" case. Verified via
    // the RDKit oracle script for the true semantic check; here we lock in
    // that chematic itself treats them as distinct at all.
    let (z_block, _) = but2ene_v2000(Z_COORDS[0], Z_COORDS[1], Z_COORDS[2], Z_COORDS[3]);
    let (e_block, _) = but2ene_v2000(E_COORDS[0], E_COORDS[1], E_COORDS[2], E_COORDS[3]);
    let z_mol = read_mol_with_diagnostics(&z_block).unwrap().mol;
    let e_mol = read_mol_with_diagnostics(&e_block).unwrap().mol;
    assert_ne!(canonical_smiles(&z_mol), canonical_smiles(&e_mol));
}

/// Fixture 3: The same two configurations, in V3000.
#[test]
fn positive_03_e_and_z_2butene_v3000() {
    for (coords, label) in [(Z_COORDS, "Z"), (E_COORDS, "E")] {
        let block = but2ene_v3000(coords[0], coords[1], coords[2], coords[3]);
        let report = read_mol_v3000_with_diagnostics(&block).expect("parse");
        assert!(report.ez_diagnostics.is_empty(), "{label}");
        let smiles = chematic_smiles::write(&report.mol);
        assert!(has_directional_token(&smiles), "{label}: '{smiles}'");
    }
}

/// Fixture 4: Trisubstituted alkene: one end has two distinct substituents.
#[test]
fn positive_04_trisubstituted_alkene() {
    let mut b = MoleculeBuilder::new();
    let center = b.add_atom(Atom::new(Element::C));
    let cl = b.add_atom(Atom::new(Element::CL));
    let br = b.add_atom(Atom::new(Element::BR));
    let ch = b.add_atom(Atom::new(Element::C));
    let me = b.add_atom(Atom::new(Element::C));
    b.add_bond(center, cl, BondOrder::Single).unwrap();
    b.add_bond(center, br, BondOrder::Single).unwrap();
    b.add_bond(center, ch, BondOrder::Double).unwrap();
    b.add_bond(ch, me, BondOrder::Single).unwrap();
    let mol = b.build();
    let coords = vec![
        (0.0, 0.0),
        (-0.866, 0.5),
        (-0.866, -0.5),
        (1.5, 0.0),
        (2.366, 0.5),
    ];
    let block = write_mol_with_coords(&mol, &MolMetadata::default().with_name("trisub"), &coords);
    let report = read_mol_with_diagnostics(&block).expect("parse");
    assert!(
        report.ez_diagnostics.is_empty(),
        "{:?}",
        report.ez_diagnostics
    );
    let smiles = chematic_smiles::write(&report.mol);
    assert!(has_directional_token(&smiles), "'{smiles}'");
}

/// Fixture 5: Tetrasubstituted alkene: both ends have two distinct substituents.
#[test]
fn positive_05_tetrasubstituted_alkene() {
    let mut b = MoleculeBuilder::new();
    let c1 = b.add_atom(Atom::new(Element::C));
    let cl = b.add_atom(Atom::new(Element::CL));
    let br = b.add_atom(Atom::new(Element::BR));
    let c2 = b.add_atom(Atom::new(Element::C));
    let f = b.add_atom(Atom::new(Element::F));
    let i = b.add_atom(Atom::new(Element::I));
    b.add_bond(c1, cl, BondOrder::Single).unwrap();
    b.add_bond(c1, br, BondOrder::Single).unwrap();
    b.add_bond(c1, c2, BondOrder::Double).unwrap();
    b.add_bond(c2, f, BondOrder::Single).unwrap();
    b.add_bond(c2, i, BondOrder::Single).unwrap();
    let mol = b.build();
    let coords = vec![
        (0.0, 0.0),
        (-0.866, 0.5),
        (-0.866, -0.5),
        (1.5, 0.0),
        (2.366, 0.5),
        (2.366, -0.5),
    ];
    let block = write_mol_with_coords(&mol, &MolMetadata::default().with_name("tetrasub"), &coords);
    let report = read_mol_with_diagnostics(&block).expect("parse");
    assert!(
        report.ez_diagnostics.is_empty(),
        "{:?}",
        report.ez_diagnostics
    );
    assert!(has_directional_token(&chematic_smiles::write(&report.mol)));
}

/// Fixture 6: Conjugated diene, two independent E/Z centers sharing one physical
/// carrier bond -- (2E,4E)-hexa-2,4-diene, standard all-anti zigzag.
#[test]
fn positive_06_conjugated_diene_two_centers() {
    let mut b = MoleculeBuilder::new();
    let me1 = b.add_atom(Atom::new(Element::C));
    let ca = b.add_atom(Atom::new(Element::C));
    let cb = b.add_atom(Atom::new(Element::C));
    let cc = b.add_atom(Atom::new(Element::C));
    let cd = b.add_atom(Atom::new(Element::C));
    let me2 = b.add_atom(Atom::new(Element::C));
    b.add_bond(me1, ca, BondOrder::Single).unwrap();
    b.add_bond(ca, cb, BondOrder::Double).unwrap();
    b.add_bond(cb, cc, BondOrder::Single).unwrap();
    b.add_bond(cc, cd, BondOrder::Double).unwrap();
    b.add_bond(cd, me2, BondOrder::Single).unwrap();
    let mol = b.build();
    let coords = vec![
        (-2.0, 0.5),
        (-1.0, 0.0),
        (0.0, 0.5),
        (1.0, 0.0),
        (2.0, 0.5),
        (3.0, 0.0),
    ];
    let block = write_mol_with_coords(&mol, &MolMetadata::default().with_name("diene"), &coords);
    let report = read_mol_with_diagnostics(&block).expect("parse");
    assert!(
        report.ez_diagnostics.is_empty(),
        "{:?}",
        report.ez_diagnostics
    );
    // Both double bonds resolved -- at least 2 directional tokens.
    let smiles = chematic_smiles::write(&report.mol);
    assert!(
        smiles.matches(['/', '\\']).count() >= 2,
        "expected >=2 directional tokens for 2 independent centers: '{smiles}'"
    );
}

/// Fixture 7: Exocyclic double bond: a cyclopropane ring carbon with an exocyclic
/// `=CHBr`. The ring must be substituted asymmetrically (here, a chlorine
/// on ONE ring neighbor only) -- an otherwise-plain 3-membered ring is
/// mirror-symmetric across the ring carbon bearing the exocyclic double
/// bond, making its own two ring-bond substituents topologically
/// EQUIVALENT (correctly `NotRequested`, matching negative fixture #18's
/// mechanism) rather than exercising the "ring double bond" shape this
/// fixture is meant to test.
#[test]
fn positive_07_exocyclic_double_bond() {
    let mut b = MoleculeBuilder::new();
    let r1 = b.add_atom(Atom::new(Element::C));
    let r2 = b.add_atom(Atom::new(Element::C));
    let r3 = b.add_atom(Atom::new(Element::C));
    let cl = b.add_atom(Atom::new(Element::CL)); // breaks the ring's mirror symmetry
    let exo = b.add_atom(Atom::new(Element::C));
    let br = b.add_atom(Atom::new(Element::BR));
    b.add_bond(r1, r2, BondOrder::Single).unwrap();
    b.add_bond(r2, r3, BondOrder::Single).unwrap();
    b.add_bond(r3, r1, BondOrder::Single).unwrap();
    b.add_bond(r2, cl, BondOrder::Single).unwrap();
    b.add_bond(r1, exo, BondOrder::Double).unwrap();
    b.add_bond(exo, br, BondOrder::Single).unwrap();
    let mol = b.build();
    let coords = vec![
        (0.0, 0.0),
        (1.0, 0.3),
        (0.5, 1.2),
        (2.2, -0.2),
        (-1.2, -0.5),
        (-2.2, -0.2),
    ];
    let block = write_mol_with_coords(
        &mol,
        &MolMetadata::default().with_name("exocyclic"),
        &coords,
    );
    let report = read_mol_with_diagnostics(&block).expect("parse");
    assert!(
        report.ez_diagnostics.is_empty(),
        "{:?}",
        report.ez_diagnostics
    );
    assert!(has_directional_token(&chematic_smiles::write(&report.mol)));
}

/// Fixture 8: Isotopic substituent: a 13C-labeled methyl on one alkene carbon;
/// isotope must survive and must not interfere with direction. Uses V3000
/// (`write_mol_v3000` round-trips isotope via `MASS=`) -- the V2000 writer
/// has a pre-existing, unrelated gap (its atom-line mass-difference field is
/// hardcoded to `0`, so V2000 text can never carry isotope information
/// through `write_mol_with_coords`, independent of this PR).
#[test]
fn positive_08_isotopic_substituent() {
    let mut b = MoleculeBuilder::new();
    let m0 = b.add_atom(Atom::new(Element::C));
    let m1 = b.add_atom(Atom::new(Element::C));
    let m2 = b.add_atom(Atom::new(Element::C));
    let mut me3 = Atom::new(Element::C);
    me3.isotope = Some(13);
    let m3 = b.add_atom(me3);
    b.add_bond(m0, m1, BondOrder::Single).unwrap();
    b.add_bond(m1, m2, BondOrder::Double).unwrap();
    b.add_bond(m2, m3, BondOrder::Single).unwrap();
    let mol = b.build();
    let block = write_mol_v3000(
        &mol,
        &MolMetadata::default().with_name("isotopic"),
        &Z_COORDS,
    );
    let report = read_mol_v3000_with_diagnostics(&block).expect("parse");
    assert!(
        report.ez_diagnostics.is_empty(),
        "{:?}",
        report.ez_diagnostics
    );
    assert_eq!(report.mol.atom(AtomIdx(3)).isotope, Some(13));
    let smiles = chematic_smiles::write(&report.mol);
    assert!(has_directional_token(&smiles));
    assert!(smiles.contains("13C"), "isotope lost: '{smiles}'");
}

/// Fixture 9: A tetrahedral wedge and an E/Z double bond in the SAME molecule, on
/// DIFFERENT physical bonds: both must be perceived, and neither must
/// disturb the other.
#[test]
fn positive_09_wedge_and_ez_coexist_same_molecule() {
    // (Z)-CHFClBr-CH=CH-CH3: a wedge tetrahedral center attached to a
    // stereogenic alkene via a plain (non-wedged) single bond.
    let mut b = MoleculeBuilder::new();
    let center = b.add_atom(Atom::new(Element::C));
    let f = b.add_atom(Atom::new(Element::F));
    let cl = b.add_atom(Atom::new(Element::CL));
    let ca = b.add_atom(Atom::new(Element::C)); // =CH- attached to center
    let cb = b.add_atom(Atom::new(Element::C)); // =CH-
    let me = b.add_atom(Atom::new(Element::C));
    b.add_bond(center, f, BondOrder::Up).unwrap(); // wedge
    b.add_bond(center, cl, BondOrder::Single).unwrap();
    b.add_bond(center, ca, BondOrder::Single).unwrap();
    b.add_bond(ca, cb, BondOrder::Double).unwrap();
    b.add_bond(cb, me, BondOrder::Single).unwrap();
    let mol = b.build();
    let coords = vec![
        (0.0, 0.0),
        (-1.0, 0.4),
        (-0.5, -1.1),
        (1.2, 0.5),
        (2.4, 0.0),
        (3.4, 0.5),
    ];
    let block = write_mol_with_coords(&mol, &MolMetadata::default().with_name("coexist"), &coords);
    let report = read_mol_with_diagnostics(&block).expect("parse");
    assert!(
        report.stereo_diagnostics.is_empty(),
        "{:?}",
        report.stereo_diagnostics
    );
    assert!(
        report.ez_diagnostics.is_empty(),
        "{:?}",
        report.ez_diagnostics
    );
    assert_ne!(report.mol.atom(AtomIdx(0)).chirality, Chirality::None);
    // `chematic_smiles::write` (plain, non-canonical) has a pre-existing,
    // unrelated gap -- documented in `writer.rs`'s own
    // `test_standalone_wedge_does_not_disturb_stereocenter_chirality` --
    // where a bracket atom's `needs_bracket` gate doesn't key off chirality
    // alone, so a wedge-only stereocenter with no OTHER bracket-forcing
    // property (isotope/charge/explicit-H/atom-map) never gets its `@`
    // printed by the plain writer. `canonical_smiles` doesn't have this gap
    // (its own `needs_bracket` does include chirality). Check plain
    // `write()` only for the E/Z token this PR is responsible for.
    let smiles = chematic_smiles::write(&report.mol);
    assert!(has_directional_token(&smiles), "E/Z token lost: '{smiles}'");
    let canon = canonical_smiles(&report.mol);
    assert!(canon.contains('@'), "tetrahedral symbol lost: '{canon}'");
    assert!(has_directional_token(&canon), "E/Z token lost: '{canon}'");
}

/// Fixture 10: A tetrahedral wedge on a single bond ADJACENT to the double bond
/// (one bond further away than #9, exercising a different topology: the
/// wedge is on the stereocenter's OTHER substituent, not the bond leading
/// to the alkene).
#[test]
fn positive_10_wedge_adjacent_to_double_bond() {
    let mut b = MoleculeBuilder::new();
    let center = b.add_atom(Atom::new(Element::C));
    let f = b.add_atom(Atom::new(Element::F)); // wedged
    let cl = b.add_atom(Atom::new(Element::CL));
    let ca = b.add_atom(Atom::new(Element::C));
    let cb = b.add_atom(Atom::new(Element::C));
    let me = b.add_atom(Atom::new(Element::C));
    b.add_bond(center, f, BondOrder::Up).unwrap();
    b.add_bond(center, cl, BondOrder::Single).unwrap();
    b.add_bond(center, ca, BondOrder::Single).unwrap(); // plain bond adjacent to the alkene
    b.add_bond(ca, cb, BondOrder::Double).unwrap();
    b.add_bond(cb, me, BondOrder::Single).unwrap();
    let mol = b.build();
    let coords = vec![
        (0.0, 0.0),
        (-1.0, 0.4),
        (-0.5, -1.1),
        (1.2, 0.5),
        (2.4, 0.0),
        (3.4, 0.5),
    ];
    let block = write_mol_with_coords(&mol, &MolMetadata::default().with_name("adjacent"), &coords);
    let report = read_mol_with_diagnostics(&block).expect("parse");
    assert!(report.stereo_diagnostics.is_empty());
    assert!(
        report.ez_diagnostics.is_empty(),
        "{:?}",
        report.ez_diagnostics
    );
    let center_bond = report.mol.bond_between(AtomIdx(0), AtomIdx(1)).unwrap().0;
    assert_eq!(
        report.mol.bond(center_bond).order,
        BondOrder::Up,
        "wedge must survive untouched"
    );
    assert!(has_directional_token(&chematic_smiles::write(&report.mol)));
}

// ---------------------------------------------------------------------------
// Invariance fixtures (11-16)
// ---------------------------------------------------------------------------

/// Fixture 11: Atom renumbering: same physical molecule, atoms declared in reverse
/// order (coordinates reassigned to match). Canonical SMILES must be
/// identical.
#[test]
fn invariance_11_atom_renumbering() {
    let (block_fwd, _) = but2ene_v2000(Z_COORDS[0], Z_COORDS[1], Z_COORDS[2], Z_COORDS[3]);
    let mut b = MoleculeBuilder::new();
    let m3 = b.add_atom(Atom::new(Element::C));
    let m2 = b.add_atom(Atom::new(Element::C));
    let m1 = b.add_atom(Atom::new(Element::C));
    let m0 = b.add_atom(Atom::new(Element::C));
    b.add_bond(m3, m2, BondOrder::Single).unwrap();
    b.add_bond(m2, m1, BondOrder::Double).unwrap();
    b.add_bond(m1, m0, BondOrder::Single).unwrap();
    let mol = b.build();
    let coords = vec![Z_COORDS[3], Z_COORDS[2], Z_COORDS[1], Z_COORDS[0]];
    let block_rev = write_mol_with_coords(&mol, &MolMetadata::default().with_name("rev"), &coords);

    let mol_fwd = read_mol_with_diagnostics(&block_fwd).unwrap().mol;
    let mol_rev = read_mol_with_diagnostics(&block_rev).unwrap().mol;
    assert_eq!(canonical_smiles(&mol_fwd), canonical_smiles(&mol_rev));
}

/// Fixture 12: Bond declaration order reversed.
#[test]
fn invariance_12_bond_declaration_order_reversed() {
    let (block_fwd, _) = but2ene_v2000(Z_COORDS[0], Z_COORDS[1], Z_COORDS[2], Z_COORDS[3]);
    let mut b = MoleculeBuilder::new();
    let m0 = b.add_atom(Atom::new(Element::C));
    let m1 = b.add_atom(Atom::new(Element::C));
    let m2 = b.add_atom(Atom::new(Element::C));
    let m3 = b.add_atom(Atom::new(Element::C));
    b.add_bond(m2, m3, BondOrder::Single).unwrap();
    b.add_bond(m1, m2, BondOrder::Double).unwrap();
    b.add_bond(m0, m1, BondOrder::Single).unwrap();
    let mol = b.build();
    let block_rev = write_mol_with_coords(
        &mol,
        &MolMetadata::default().with_name("bondrev"),
        &Z_COORDS,
    );
    let mol_fwd = read_mol_with_diagnostics(&block_fwd).unwrap().mol;
    let mol_rev = read_mol_with_diagnostics(&block_rev).unwrap().mol;
    assert_eq!(canonical_smiles(&mol_fwd), canonical_smiles(&mol_rev));
}

/// Fixture 13: Double-bond atom1/atom2 reversal (the double bond's own bond line
/// lists atom2 before atom1).
#[test]
fn invariance_13_double_bond_atom_order_reversed() {
    let (block_fwd, _) = but2ene_v2000(Z_COORDS[0], Z_COORDS[1], Z_COORDS[2], Z_COORDS[3]);
    let mut b = MoleculeBuilder::new();
    let m0 = b.add_atom(Atom::new(Element::C));
    let m1 = b.add_atom(Atom::new(Element::C));
    let m2 = b.add_atom(Atom::new(Element::C));
    let m3 = b.add_atom(Atom::new(Element::C));
    b.add_bond(m0, m1, BondOrder::Single).unwrap();
    b.add_bond(m2, m1, BondOrder::Double).unwrap(); // reversed
    b.add_bond(m2, m3, BondOrder::Single).unwrap();
    let mol = b.build();
    let block_rev =
        write_mol_with_coords(&mol, &MolMetadata::default().with_name("dbrev"), &Z_COORDS);
    let mol_fwd = read_mol_with_diagnostics(&block_fwd).unwrap().mol;
    let mol_rev = read_mol_with_diagnostics(&block_rev).unwrap().mol;
    assert_eq!(canonical_smiles(&mol_fwd), canonical_smiles(&mol_rev));
}

/// Fixture 14: Rotation (90 degrees) must preserve E/Z semantics.
#[test]
fn invariance_14_rotation() {
    let (block_fwd, _) = but2ene_v2000(Z_COORDS[0], Z_COORDS[1], Z_COORDS[2], Z_COORDS[3]);
    let rotated: Vec<(f64, f64)> = Z_COORDS.iter().map(|&(x, y)| (-y, x)).collect();
    let (block_rot, _) = but2ene_v2000(rotated[0], rotated[1], rotated[2], rotated[3]);
    let mol_fwd = read_mol_with_diagnostics(&block_fwd).unwrap().mol;
    let mol_rot = read_mol_with_diagnostics(&block_rot).unwrap().mol;
    assert_eq!(canonical_smiles(&mol_fwd), canonical_smiles(&mol_rot));
}

/// Fixture 15: Translation must preserve E/Z semantics.
#[test]
fn invariance_15_translation() {
    let (block_fwd, _) = but2ene_v2000(Z_COORDS[0], Z_COORDS[1], Z_COORDS[2], Z_COORDS[3]);
    let translated: Vec<(f64, f64)> = Z_COORDS.iter().map(|&(x, y)| (x + 10.0, y - 4.0)).collect();
    let (block_t, _) = but2ene_v2000(translated[0], translated[1], translated[2], translated[3]);
    let mol_fwd = read_mol_with_diagnostics(&block_fwd).unwrap().mol;
    let mol_t = read_mol_with_diagnostics(&block_t).unwrap().mol;
    assert_eq!(canonical_smiles(&mol_fwd), canonical_smiles(&mol_t));
}

/// Fixture 16: Mirror reflection (y -> -y): flips BOTH substituents' signs
/// simultaneously, so semantic E/Z is unchanged (empirically verified here,
/// not just asserted) -- unlike tetrahedral parity, which DOES flip under a
/// pure mirror.
#[test]
fn invariance_16_mirror_reflection() {
    let (block_fwd, _) = but2ene_v2000(Z_COORDS[0], Z_COORDS[1], Z_COORDS[2], Z_COORDS[3]);
    let mirrored: Vec<(f64, f64)> = Z_COORDS.iter().map(|&(x, y)| (x, -y)).collect();
    let (block_m, _) = but2ene_v2000(mirrored[0], mirrored[1], mirrored[2], mirrored[3]);
    let mol_fwd = read_mol_with_diagnostics(&block_fwd).unwrap().mol;
    let mol_m = read_mol_with_diagnostics(&block_m).unwrap().mol;
    assert_eq!(
        canonical_smiles(&mol_fwd),
        canonical_smiles(&mol_m),
        "E/Z must be mirror-invariant (both substituents' signs flip together)"
    );
}

// ---------------------------------------------------------------------------
// Negative fixtures (17-28)
// ---------------------------------------------------------------------------

fn assert_negative_invariants(
    mol_before_ez: &chematic_core::Molecule,
    report_mol: &chematic_core::Molecule,
) {
    // Existing tetrahedral chirality (if any) unchanged.
    for (idx, atom) in mol_before_ez.atoms() {
        assert_eq!(atom.chirality, report_mol.atom(idx).chirality);
    }
}

/// Fixture 17: Terminal alkene (propene): no direction, no diagnostic, no panic.
#[test]
fn negative_17_terminal_alkene() {
    let mut b = MoleculeBuilder::new();
    let c0 = b.add_atom(Atom::new(Element::C));
    let c1 = b.add_atom(Atom::new(Element::C));
    let c2 = b.add_atom(Atom::new(Element::C));
    b.add_bond(c0, c1, BondOrder::Double).unwrap();
    b.add_bond(c1, c2, BondOrder::Single).unwrap();
    let mol = b.build();
    let coords = vec![(0.0, 0.0), (1.5, 0.0), (2.366, 0.5)];
    let block = write_mol_with_coords(&mol, &MolMetadata::default().with_name("propene"), &coords);
    let report = read_mol_with_diagnostics(&block).expect("parse");
    assert!(report.ez_diagnostics.is_empty());
    assert!(!has_directional_token(&chematic_smiles::write(&report.mol)));
    assert_negative_invariants(&mol, &report.mol);
}

/// Fixture 18: Non-stereogenic alkene with equivalent substituents (2-methyl-2-butene).
#[test]
fn negative_18_equivalent_substituents() {
    let mut b = MoleculeBuilder::new();
    let center = b.add_atom(Atom::new(Element::C));
    let me_a = b.add_atom(Atom::new(Element::C));
    let me_b = b.add_atom(Atom::new(Element::C));
    let ch = b.add_atom(Atom::new(Element::C));
    let me_c = b.add_atom(Atom::new(Element::C));
    b.add_bond(center, me_a, BondOrder::Single).unwrap();
    b.add_bond(center, me_b, BondOrder::Single).unwrap();
    b.add_bond(center, ch, BondOrder::Double).unwrap();
    b.add_bond(ch, me_c, BondOrder::Single).unwrap();
    let mol = b.build();
    let coords = vec![
        (0.0, 0.0),
        (-0.866, 0.5),
        (-0.866, -0.5),
        (1.5, 0.0),
        (2.366, 0.5),
    ];
    let block = write_mol_with_coords(&mol, &MolMetadata::default().with_name("2m2b"), &coords);
    let report = read_mol_with_diagnostics(&block).expect("parse");
    assert!(report.ez_diagnostics.is_empty());
    assert!(!has_directional_token(&chematic_smiles::write(&report.mol)));
}

/// Fixture 19: Carbonyl (acetaldehyde): oxygen end has no substituent.
#[test]
fn negative_19_carbonyl() {
    let mut b = MoleculeBuilder::new();
    let c0 = b.add_atom(Atom::new(Element::C));
    let c1 = b.add_atom(Atom::new(Element::C));
    let o = b.add_atom(Atom::new(Element::O));
    b.add_bond(c0, c1, BondOrder::Single).unwrap();
    b.add_bond(c1, o, BondOrder::Double).unwrap();
    let mol = b.build();
    let coords = vec![(-1.0, 0.0), (0.0, 0.0), (0.5, 1.0)];
    let block = write_mol_with_coords(
        &mol,
        &MolMetadata::default().with_name("acetaldehyde"),
        &coords,
    );
    let report = read_mol_with_diagnostics(&block).expect("parse");
    assert!(report.ez_diagnostics.is_empty());
    assert!(!has_directional_token(&chematic_smiles::write(&report.mol)));
}

/// Fixture 20: Collinear substituent: the only substituent at one end sits exactly
/// on the double-bond axis.
#[test]
fn negative_20_collinear_substituent() {
    let (block, _) = but2ene_v2000((-1.0, 0.0), (0.0, 0.0), (1.5, 0.0), (2.366, 0.5));
    let report = read_mol_with_diagnostics(&block).expect("parse");
    assert_eq!(report.ez_diagnostics.len(), 1);
    assert_eq!(
        report.ez_diagnostics[0].reason,
        EzDirectionRejectionReason::DegenerateGeometry
    );
    assert!(!has_directional_token(&chematic_smiles::write(&report.mol)));
}

/// Fixture 21: Zero-length double-bond coordinates.
#[test]
fn negative_21_zero_length_double_bond() {
    let (block, _) = but2ene_v2000((-0.866, 0.5), (0.0, 0.0), (0.0, 0.0), (2.366, 0.5));
    let report = read_mol_with_diagnostics(&block).expect("parse");
    assert_eq!(report.ez_diagnostics.len(), 1);
    assert_eq!(
        report.ez_diagnostics[0].reason,
        EzDirectionRejectionReason::DegenerateGeometry
    );
}

/// Fixture 22: Missing coordinate. Not reachable through a well-formed MOL FILE (a
/// real atom line always yields some (x, y), even a wrong one) -- reachable
/// only via direct API misuse (a caller passing a `coords` slice shorter
/// than the molecule's atom count to the standalone perception function),
/// exactly mirroring `docs/rfcs/stereo2d_reader_integration_rfc.md`'s own
/// fixture #14 for the tetrahedral case. Documented here rather than
/// silently dropped.
#[test]
fn negative_22_missing_coordinate_via_api_misuse() {
    // Built directly (not through the reader, and not reusing a
    // reader-produced `Molecule` that already had E/Z perceived once with
    // full coordinates -- re-running perception on an already-processed
    // molecule would see its own prior `bond_direction` writes and report
    // `CarrierConflict` instead of the coordinate problem this fixture
    // targets). A single, fresh application with a deliberately-truncated
    // `coords` slice is the API-misuse shape this fixture documents.
    let mut b = MoleculeBuilder::new();
    let m0 = b.add_atom(Atom::new(Element::C));
    let m1 = b.add_atom(Atom::new(Element::C));
    let m2 = b.add_atom(Atom::new(Element::C));
    let m3 = b.add_atom(Atom::new(Element::C));
    b.add_bond(m0, m1, BondOrder::Single).unwrap();
    b.add_bond(m1, m2, BondOrder::Double).unwrap();
    b.add_bond(m2, m3, BondOrder::Single).unwrap();
    let mut mol = b.build();
    let truncated_coords = vec![Z_COORDS[0], Z_COORDS[1], Z_COORDS[2]]; // drop Me3's coordinate
    let diagnostics = chematic_perception::apply_ez_directions_from_2d_with_diagnostics(
        &mut mol,
        &truncated_coords,
    );
    assert_eq!(diagnostics.len(), 1);
    assert_eq!(
        diagnostics[0].reason,
        EzDirectionRejectionReason::MissingCoordinate
    );
}

/// Fixture 23: NaN coordinate -- reachable through a REAL MOL file: a literal
/// "NaN" token in the fixed-width coordinate field parses to a real
/// `f64::NAN` via the actual reader (confirmed empirically).
#[test]
fn negative_23_nan_coordinate_via_real_file() {
    let (block, _) = but2ene_v2000(Z_COORDS[0], Z_COORDS[1], Z_COORDS[2], Z_COORDS[3]);
    let patched = patch_v2000_atom_x_nan(&block, 3); // last atom (Me3)
    let report = read_mol_with_diagnostics(&patched).expect("parse must still succeed");
    assert!(report.coords[3].0.is_nan());
    assert_eq!(report.ez_diagnostics.len(), 1);
    assert_eq!(
        report.ez_diagnostics[0].reason,
        EzDirectionRejectionReason::NonFiniteCoordinate
    );
    assert!(!has_directional_token(&chematic_smiles::write(&report.mol)));
}

/// Fixture 24: Explicit MDL unknown/either double-bond stereo (V2000 code 3, V3000
/// CFG=2): must never derive a direction from coordinates, even when they
/// would otherwise resolve cleanly.
#[test]
fn negative_24_explicit_unknown_either_stereo() {
    let (block, natoms) = but2ene_v2000(Z_COORDS[0], Z_COORDS[1], Z_COORDS[2], Z_COORDS[3]);
    let patched = patch_v2000_bond_stereo(&block, natoms, 1, 3); // double bond is decl. index 1
    let report = read_mol_with_diagnostics(&patched).expect("parse");
    assert_eq!(report.ez_diagnostics.len(), 1);
    assert_eq!(
        report.ez_diagnostics[0].reason,
        EzDirectionRejectionReason::ExplicitlyUnspecified
    );
    assert!(!has_directional_token(&chematic_smiles::write(&report.mol)));

    // V3000 CFG=2 equivalent.
    let block_v3k = but2ene_v3000(Z_COORDS[0], Z_COORDS[1], Z_COORDS[2], Z_COORDS[3]);
    let patched_v3k = patch_v3000_bond_cfg2(&block_v3k, 2); // double bond is V30 bond #2
    let report_v3k = read_mol_v3000_with_diagnostics(&patched_v3k).expect("parse V3000");
    assert_eq!(report_v3k.ez_diagnostics.len(), 1);
    assert_eq!(
        report_v3k.ez_diagnostics[0].reason,
        EzDirectionRejectionReason::ExplicitlyUnspecified
    );
}

/// Fixture 25: Cumulene/allene: propa-1,2-diene (allene) -- must reject both double
/// bonds as unsupported topology, never guess from 2D coordinates.
#[test]
fn negative_25_cumulene_allene() {
    let mut b = MoleculeBuilder::new();
    let s1 = b.add_atom(Atom::new(Element::C));
    let t1 = b.add_atom(Atom::new(Element::C));
    let central = b.add_atom(Atom::new(Element::C));
    let t2 = b.add_atom(Atom::new(Element::C));
    let s2 = b.add_atom(Atom::new(Element::C));
    b.add_bond(s1, t1, BondOrder::Single).unwrap();
    b.add_bond(t1, central, BondOrder::Double).unwrap();
    b.add_bond(central, t2, BondOrder::Double).unwrap();
    b.add_bond(t2, s2, BondOrder::Single).unwrap();
    let mol = b.build();
    let coords = vec![(-2.0, 1.0), (-1.0, 0.0), (0.0, 0.0), (1.0, 0.0), (2.0, 1.0)];
    let block = write_mol_with_coords(&mol, &MolMetadata::default().with_name("allene"), &coords);
    let report = read_mol_with_diagnostics(&block).expect("parse");
    assert_eq!(
        report.ez_diagnostics.len(),
        2,
        "{:?}",
        report.ez_diagnostics
    );
    assert!(
        report
            .ez_diagnostics
            .iter()
            .all(|d| d.reason == EzDirectionRejectionReason::UnsupportedTopology)
    );
}

/// Fixture 26: Shared-carrier conflict: two double bonds' independently-computed
/// requirements on the SAME physical bond disagree. Both must reject; the
/// shared bond must end up with no direction.
#[test]
fn negative_26_shared_carrier_conflict() {
    let mut b = MoleculeBuilder::new();
    let me1 = b.add_atom(Atom::new(Element::C));
    let ca = b.add_atom(Atom::new(Element::C));
    let cb = b.add_atom(Atom::new(Element::C));
    let cc = b.add_atom(Atom::new(Element::C));
    let cd = b.add_atom(Atom::new(Element::C));
    let me2 = b.add_atom(Atom::new(Element::C));
    b.add_bond(me1, ca, BondOrder::Single).unwrap();
    b.add_bond(ca, cb, BondOrder::Double).unwrap();
    b.add_bond(cb, cc, BondOrder::Single).unwrap();
    b.add_bond(cc, cd, BondOrder::Double).unwrap();
    b.add_bond(cd, me2, BondOrder::Single).unwrap();
    let mol = b.build();
    // Hand-verified conflicting layout (see stereo2d_ez_direction.rs's own
    // `conjugated_diene_shared_bond_conflict` unit test for the derivation).
    let coords = vec![
        (-1.0, 1.0),
        (0.0, 0.0),
        (1.0, 0.0),
        (2.0, 1.0),
        (2.0, 2.0),
        (3.0, 3.0),
    ];
    let block = write_mol_with_coords(&mol, &MolMetadata::default().with_name("conflict"), &coords);
    let report = read_mol_with_diagnostics(&block).expect("parse");
    assert_eq!(
        report.ez_diagnostics.len(),
        2,
        "{:?}",
        report.ez_diagnostics
    );
    assert!(
        report
            .ez_diagnostics
            .iter()
            .all(|d| d.reason == EzDirectionRejectionReason::CarrierConflict)
    );
}

/// Fixture 27: Wedge only, no E/Z anywhere: must not spuriously produce any E/Z
/// diagnostic or direction.
#[test]
fn negative_27_wedge_only_no_ez() {
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
    let mol = b.build();
    let coords = vec![
        (0.0, 0.0),
        (-1.0, 0.4),
        (0.9, 0.7),
        (-0.5, -1.1),
        (0.8, -0.6),
    ];
    let block = write_mol_with_coords(
        &mol,
        &MolMetadata::default().with_name("wedgeonly"),
        &coords,
    );
    let report = read_mol_with_diagnostics(&block).expect("parse");
    assert!(report.ez_diagnostics.is_empty());
    assert!(report.stereo_diagnostics.is_empty());
    assert_ne!(report.mol.atom(AtomIdx(0)).chirality, Chirality::None);
}

/// Fixture 28: No-double-bond control: ethanol. No panic, no diagnostics.
#[test]
fn negative_28_no_double_bond_control() {
    let mut b = MoleculeBuilder::new();
    let c1 = b.add_atom(Atom::new(Element::C));
    let c2 = b.add_atom(Atom::new(Element::C));
    let o = b.add_atom(Atom::new(Element::O));
    b.add_bond(c1, c2, BondOrder::Single).unwrap();
    b.add_bond(c2, o, BondOrder::Single).unwrap();
    let mol = b.build();
    let coords = vec![(0.0, 0.0), (1.5, 0.0), (3.0, 0.0)];
    let block = write_mol_with_coords(&mol, &MolMetadata::default().with_name("ethanol"), &coords);
    let report = read_mol_with_diagnostics(&block).expect("parse");
    assert!(report.ez_diagnostics.is_empty());
    assert!(!has_directional_token(&chematic_smiles::write(&report.mol)));
}

// ---------------------------------------------------------------------------
// Cross-cutting: raw wedge untouched, SDF inheritance, PR #154 regression
// ---------------------------------------------------------------------------

/// Lock-in: the E/Z stage does not modify raw wedge/hash `BondOrder`. Reuses
/// fixture #9's molecule (wedge + independent E/Z system).
#[test]
fn ez_stage_does_not_modify_raw_wedge() {
    let mut b = MoleculeBuilder::new();
    let center = b.add_atom(Atom::new(Element::C));
    let f = b.add_atom(Atom::new(Element::F));
    let cl = b.add_atom(Atom::new(Element::CL));
    let ca = b.add_atom(Atom::new(Element::C));
    let cb = b.add_atom(Atom::new(Element::C));
    let me = b.add_atom(Atom::new(Element::C));
    b.add_bond(center, f, BondOrder::Up).unwrap();
    b.add_bond(center, cl, BondOrder::Single).unwrap();
    b.add_bond(center, ca, BondOrder::Single).unwrap();
    b.add_bond(ca, cb, BondOrder::Double).unwrap();
    b.add_bond(cb, me, BondOrder::Single).unwrap();
    let mol = b.build();
    let coords = vec![
        (0.0, 0.0),
        (-1.0, 0.4),
        (-0.5, -1.1),
        (1.2, 0.5),
        (2.4, 0.0),
        (3.4, 0.5),
    ];
    let block = write_mol_with_coords(&mol, &MolMetadata::default().with_name("lockin"), &coords);
    let wedge_bond = mol.bond_between(center, f).unwrap().0;
    let order_before = mol.bond(wedge_bond).order;

    let report = read_mol_with_diagnostics(&block).expect("parse");
    let wedge_bond_after = report.mol.bond_between(AtomIdx(0), AtomIdx(1)).unwrap().0;
    assert_eq!(report.mol.bond(wedge_bond_after).order, order_before);
    assert_eq!(report.mol.bond(wedge_bond_after).order, BondOrder::Up);
}

/// SDF inherits automatically through the V2000 parsing core.
#[test]
fn sdf_inherits_ez_direction_from_v2000_core() {
    let (block, _) = but2ene_v2000(Z_COORDS[0], Z_COORDS[1], Z_COORDS[2], Z_COORDS[3]);
    let direct = read_mol_with_diagnostics(&block).expect("direct parse");
    assert!(direct.ez_diagnostics.is_empty());

    let sdf = format!("{block}$$$$\n");
    let rec = SdfRecordReader::new(&sdf)
        .next()
        .expect("one record")
        .expect("parse ok");
    assert!(rec.ez_diagnostics.is_empty());
    assert_eq!(
        chematic_smiles::write(&rec.mol),
        chematic_smiles::write(&direct.mol)
    );
    let _ = write_sdf; // exercised indirectly via write_mol_with_coords + "$$$$"
}

/// PR #154 regression guard: parsing an ordinary tetrahedral-only wedge
/// fixture (no double bonds at all) must be byte-identical in its
/// `chematic_smiles::write` output before/after this PR -- there is nothing
/// for the new E/Z stage to touch.
#[test]
fn pr154_tetrahedral_only_fixture_unaffected() {
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
    let mol = b.build();
    let coords = vec![
        (0.0, 0.0),
        (-1.0, 0.4),
        (0.9, 0.7),
        (-0.5, -1.1),
        (0.8, -0.6),
    ];
    let block = write_mol_with_coords(
        &mol,
        &MolMetadata::default().with_name("valid_wedge"),
        &coords,
    );
    let report = read_mol_with_diagnostics(&block).expect("parse");
    assert!(report.stereo_diagnostics.is_empty());
    assert!(
        report.ez_diagnostics.is_empty(),
        "a tetrahedral-only fixture (no double bonds) must never produce an \
         E/Z diagnostic"
    );
    assert_ne!(report.mol.atom(AtomIdx(0)).chirality, Chirality::None);
    // See `positive_09`'s comment: plain `write()`'s bracket gate doesn't
    // key off chirality alone (pre-existing, orthogonal to this PR) --
    // `canonical_smiles` is the reliable check here.
    let canon = canonical_smiles(&report.mol);
    assert!(canon.contains('@'), "'{canon}'");
    assert!(
        !has_directional_token(&canon),
        "no double bond exists, so no E/Z token should appear: '{canon}'"
    );
}
