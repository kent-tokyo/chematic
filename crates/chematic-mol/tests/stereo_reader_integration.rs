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

// Open Babel 3.2.1 writes its 0D tetrahedral stereo in V2000 atom parity.
// Its own MOL re-reader returns F[C@](Cl)(Br)I for parity 2 and the opposite
// enantiomer for parity 1. No coordinates or wedge bonds are present.
const OPENBABEL_0D_PARITY_2: &str = r#"
 OpenBabel10152614202D

  5  4  0  0  1  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  2  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  2  4  1  0  0  0  0
  2  5  1  0  0  0  0
M  END
"#;

#[test]
fn openbabel_0d_atom_parity_preserves_both_enantiomers() {
    for (parity, expected_smiles) in [(2, "F[C@](Cl)(Br)I"), (1, "F[C@@](Cl)(Br)I")] {
        let block =
            OPENBABEL_0D_PARITY_2.replace(" C   0  0  2  0", &format!(" C   0  0  {parity}  0"));
        let read = read_mol_with_diagnostics(&block).expect("0D V2000 parses");
        let expected = chematic_smiles::parse(expected_smiles).expect("reference SMILES parses");
        assert_eq!(
            chematic_smiles::canonical_smiles(&read.mol),
            chematic_smiles::canonical_smiles(&expected),
            "atom parity {parity}"
        );
    }
}

// Open Babel keeps a leading explicit hydrogen as atom 1, but its 0D parity
// orders that hydrogen after the other substituents of the tetrahedral center.
const OPENBABEL_EXPLICIT_H_FIRST_0D: &str = r#"
 OpenBabel10152619472D

  5  4  0  0  1  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  2  4  1  0  0  0  0
  2  5  1  0  0  0  0
M  END
"#;

#[test]
fn openbabel_0d_atom_parity_places_explicit_hydrogen_last() {
    for (block, expected_smiles) in [
        (
            OPENBABEL_EXPLICIT_H_FIRST_0D.to_string(),
            "[H][C@](F)(Cl)Br",
        ),
        (
            OPENBABEL_EXPLICIT_H_FIRST_0D.replace("M  END", "M  ISO  1   1   2\nM  END"),
            "[2H][C@](F)(Cl)Br",
        ),
    ] {
        let read = read_mol_with_diagnostics(&block).expect("0D V2000 parses");
        let expected = chematic_smiles::parse(expected_smiles).expect("reference SMILES parses");
        assert_eq!(
            chematic_smiles::canonical_smiles(&read.mol),
            chematic_smiles::canonical_smiles(&expected),
            "explicit hydrogen in {expected_smiles}"
        );
    }
}

const OPENBABEL_0D_SULFOXIDE: &str = r#"
 OpenBabel10152619492D

  9  9  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 S   0  0  2  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  2  4  1  0  0  0  0
  4  9  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  7  1  0  0  0  0
  7  8  2  0  0  0  0
  8  9  1  0  0  0  0
M  END
"#;

#[test]
fn openbabel_0d_atom_parity_preserves_sulfoxide_stereo() {
    for (parity, expected_smiles) in [(2, "C1=CC=CC=C1[S@@](=O)C"), (1, "C1=CC=CC=C1[S@](=O)C")] {
        let block =
            OPENBABEL_0D_SULFOXIDE.replace(" S   0  0  2  0", &format!(" S   0  0  {parity}  0"));
        let read = read_mol_with_diagnostics(&block).expect("0D sulfoxide parses");
        assert_eq!(
            chematic_smiles::canonical_smiles(&read.mol),
            expected_smiles,
            "sulfoxide parity {parity}"
        );
    }
}

const OPENBABEL_ALANINE_0D: &str = r#"
 OpenBabel10152614232D

  6  5  0  0  1  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  2  4  1  0  0  0  0
  4  5  2  0  0  0  0
  4  6  1  0  0  0  0
M  END
"#;

#[test]
fn openbabel_0d_atom_parity_with_implicit_hydrogen() {
    for (parity, expected_smiles) in [(1, "N[C@@H](C)C(=O)O"), (2, "N[C@H](C)C(=O)O")] {
        let block =
            OPENBABEL_ALANINE_0D.replace(" C   0  0  1  0", &format!(" C   0  0  {parity}  0"));
        let read = read_mol_with_diagnostics(&block).expect("0D V2000 parses");
        let expected = chematic_smiles::parse(expected_smiles).expect("reference SMILES parses");
        assert_eq!(
            chematic_smiles::canonical_smiles(&read.mol),
            chematic_smiles::canonical_smiles(&expected),
            "atom parity {parity}"
        );
    }
}

#[test]
fn atom_parity_is_not_inferred_from_nonzero_coordinates() {
    let block = OPENBABEL_0D_PARITY_2.replace(
        "    0.0000    0.0000    0.0000 F  ",
        "    1.0000    0.0000    0.0000 F  ",
    );
    let read = read_mol_with_diagnostics(&block).expect("nonzero-coordinate MOL parses");
    assert_eq!(read.mol.atom(AtomIdx(1)).chirality, Chirality::None);
}

#[test]
fn atom_parity_does_not_invent_chirality_with_duplicate_substituents() {
    let block = OPENBABEL_0D_PARITY_2.replace(" I   0  0  0  0", " F   0  0  0  0");
    let read = read_mol_with_diagnostics(&block).expect("parity on achiral atom parses");
    assert_eq!(read.mol.atom(AtomIdx(1)).chirality, Chirality::None);
}

#[test]
fn atom_parity_does_not_override_2d_wedge_geometry() {
    let original = read_mol_with_diagnostics(RDKIT_L_ALANINE_2D_V2000).expect("2D alanine parses");
    for parity in [1, 2] {
        let block =
            RDKIT_L_ALANINE_2D_V2000.replace(" C   0  0  0  0", &format!(" C   0  0  {parity}  0"));
        let read = read_mol_with_diagnostics(&block).expect("2D alanine with atom parity parses");
        assert_eq!(
            chematic_smiles::canonical_smiles(&read.mol),
            chematic_smiles::canonical_smiles(&original.mol),
            "2D wedges take precedence over parity {parity}"
        );
    }
}

const OPENBABEL_0D_ALKENE: &str = r#"
 OpenBabel10152614252D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  1  0  0  0
  2  3  2  0  0  0  0
  3  4  1  6  0  0  0
M  END
"#;

const OPENBABEL_0D_CONJUGATED_DIENE: &str = r#"
 OpenBabel10152619522D

  6  5  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0
  1  2  1  1  0  0  0
  2  3  2  0  0  0  0
  3  4  1  6  0  0  0
  4  5  2  0  0  0  0
  5  6  1  1  0  0  0
M  END
"#;

#[test]
fn openbabel_0d_shared_direction_between_double_bonds_is_reversed() {
    for (middle_marker, expected_smiles) in [(6, "C/C=C/C=C/C"), (1, "C/C=C\\C=C/C")] {
        let block = OPENBABEL_0D_CONJUGATED_DIENE
            .replace("  3  4  1  6  0", &format!("  3  4  1  {middle_marker}  0"));
        let read = read_mol_with_diagnostics(&block).expect("0D diene parses");
        let expected = chematic_smiles::parse(expected_smiles).expect("reference SMILES parses");
        assert_eq!(
            chematic_smiles::canonical_smiles(&read.mol),
            chematic_smiles::canonical_smiles(&expected),
            "middle marker {middle_marker}"
        );
    }
}

// Two substituents at C3 and C8 each carry an Open Babel direction marker.
// Skipping those double-bond ends inverted the final imine's E/Z assignment.
const OPENBABEL_0D_POLYIMINE: &str = r#"
 OpenBabel10152620252D

 11 10  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0
  1  2  1  1  0  0  0
  2  3  2  0  0  0  0
  3  4  1  1  0  0  0
  3  5  1  6  0  0  0
  5  6  1  1  0  0  0
  6  7  2  0  0  0  0
  7  8  1  6  0  0  0
  8  9  1  1  0  0  0
  8 10  2  0  0  0  0
 10 11  1  1  0  0  0
M  END
"#;

#[test]
fn openbabel_0d_multiple_markers_at_one_double_bond_end() {
    let read = read_mol_with_diagnostics(OPENBABEL_0D_POLYIMINE).expect("0D polyimine parses");
    let expected = chematic_smiles::parse("C/N=C(\\S)N/N=C/C(C)=N/N").expect("source parses");
    assert_eq!(
        chematic_smiles::canonical_smiles(&read.mol),
        chematic_smiles::canonical_smiles(&expected)
    );
}

#[test]
fn openbabel_0d_alkene_directions_preserve_e_and_z() {
    for (last_marker, expected_smiles) in [(6, "F/C=C/F"), (1, "F/C=C\\F")] {
        let block =
            OPENBABEL_0D_ALKENE.replace("  3  4  1  6", &format!("  3  4  1  {last_marker}"));
        let read = read_mol_with_diagnostics(&block).expect("0D alkene parses");
        let expected = chematic_smiles::parse(expected_smiles).expect("reference SMILES parses");
        assert_eq!(
            chematic_smiles::canonical_smiles(&read.mol),
            chematic_smiles::canonical_smiles(&expected),
            "last marker {last_marker}"
        );
    }
}

const OPENBABEL_0D_BRANCHED_ALKENE: &str = r#"
 OpenBabel10152614292D

  7  6  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  1  0  0  0
  2  3  2  0  0  0  0
  3  4  1  1  0  0  0
  3  5  1  6  0  0  0
  5  6  2  0  0  0  0
  5  7  1  0  0  0  0
M  END
"#;

const OPENBABEL_0D_AROMATIC_ALKENE: &str = r#"
 OpenBabel10152614292D

  9  9  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  9  2  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  1  0  0  0
  3  7  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  1  0  0  0
  7  8  2  0  0  0  0
  8  9  1  0  0  0  0
M  END
"#;

#[test]
fn openbabel_0d_alkene_directions_ignore_carbonyl_and_aromatic_doubles() {
    for (block, smiles) in [
        (OPENBABEL_0D_BRANCHED_ALKENE, "F/C=C(\\C)C(=O)N"),
        (OPENBABEL_0D_AROMATIC_ALKENE, "c1cc(/C=C\\F)ccc1"),
    ] {
        let read = read_mol_with_diagnostics(block).expect("0D alkene parses");
        let expected = chematic_smiles::parse(smiles).expect("reference SMILES parses");
        let aromatic_normalized = chematic_perception::apply_aromaticity(&read.mol);
        assert_eq!(
            chematic_smiles::canonical_smiles(&aromatic_normalized),
            chematic_smiles::canonical_smiles(&expected),
            "{smiles}"
        );
    }
}

const RDKIT_L_ALANINE_2D_V2000: &str = r#"
     RDKit          2D

  6  5  0  0  0  0  0  0  0  0999 V2000
    1.5000   -1.2990    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000   -1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  6
  2  4  1  0
  4  5  2  0
  4  6  1  0
M  END
"#;

// RDKit 2025.09.3 2D V2000 output for a four-center bridged-ring example.
// Keep this format-level fixture separate from its source SMILES: the point
// is to prove that the wedge on the three-heavy-neighbor center survives the
// reader, before the SMILES writer's traversal-dependent parity remapping is
// involved.
const RDKIT_BRIDGED_RING_2D_V2000: &str = r#"
     RDKit          2D

 28 32  0  0  0  0  0  0  0  0999 V2000
    1.0462   -1.7455    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0125   -0.6828    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4621   -1.0683    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5207   -0.0056    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9703   -0.3910    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0290    0.6716    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.4786    0.2862    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.5372    1.3489    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.9869    0.9634    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.3779   -0.4847    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -8.3192   -1.5474    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.8696   -1.1620    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1297    1.4426    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6801    1.8280    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.3785    0.7653    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8282    1.1508    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8868    0.0881    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4958   -1.3601    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.3364    0.4735    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.3951   -0.5891    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.8447   -0.2037    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.2357    1.2444    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.1771    2.3071    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.6186   -0.5877    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0041   -2.0373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.4537   -1.6518    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9454   -0.9746    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.7274    1.9217    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  1  0
  6  7  1  0
  7  8  1  0
  8  9  1  0
  9 10  1  0
 10 11  1  0
 11 12  1  0
  4 13  1  0
 13 14  1  0
 14 15  2  0
 15 16  1  0
 16 17  1  0
 17 18  2  0
 19 17  1  6
 19 20  1  0
 21 20  1  6
 21 22  1  0
 23 22  1  6
 23 24  1  0
 25 24  1  1
 25 26  1  0
 25 27  1  0
 23 28  1  0
 15  2  1  0
 27 19  1  0
 12  7  1  0
 28 19  1  0
 26 21  1  0
M  END
"#;

#[test]
fn rdkit_v2000_implicit_h_wedge_keeps_semantic_smiles_parity() {
    // RDKit 2025.09.3 2D output for N[C@@H](C)C(=O)O. Its canonical
    // isomeric SMILES is C[C@H](N)C(=O)O. Keep the expected form through
    // chematic's own SMILES parser so this fixture asserts the physical
    // configuration rather than an arbitrary writer spelling.
    let read = read_mol_with_diagnostics(RDKIT_L_ALANINE_2D_V2000).expect("RDKit V2000 parses");
    let expected = chematic_smiles::canonical_smiles(
        &chematic_smiles::parse("C[C@H](N)C(=O)O").expect("reference SMILES parses"),
    );
    assert_eq!(chematic_smiles::canonical_smiles(&read.mol), expected);
}

#[test]
fn rdkit_bridged_ring_implicit_h_wedge_is_perceived_before_smiles_writing() {
    let read = read_mol_with_diagnostics(RDKIT_BRIDGED_RING_2D_V2000)
        .expect("RDKit bridged-ring V2000 parses");
    // Atom 21 in the V2000 record is a three-heavy-neighbor carbon with an
    // implicit H. Its `6` bond stereo flag must produce a local parity, even
    // though the later canonical SMILES traversal currently has a separate
    // complex-ring residual to resolve.
    assert_ne!(
        read.mol.atom(AtomIdx(20)).chirality,
        Chirality::None,
        "implicit H={} diagnostics={:?}",
        chematic_core::implicit_hcount(&read.mol, AtomIdx(20)),
        read.stereo_diagnostics
    );
    assert!(
        !read
            .stereo_diagnostics
            .iter()
            .any(|diagnostic| diagnostic.atom == AtomIdx(20)),
        "the T-shaped implicit-H center must not be rejected: {:?}",
        read.stereo_diagnostics
    );
}

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
