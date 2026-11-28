//! The plain (non-canonical) writer keeps stereo tokens and charges that
//! the canonical paths already cover: allene, three-ligand tetrahedral and
//! square-planar centres, and multiply charged atoms.

use chematic_smiles::{canonical_smiles, parse, write};

fn assert_round_trip(input: &str) -> String {
    let mol = parse(input).unwrap();
    let written = write(&mol);
    let reparsed = parse(&written).unwrap_or_else(|e| panic!("{input} -> {written}: {e}"));
    assert_eq!(
        canonical_smiles(&reparsed),
        canonical_smiles(&mol),
        "{input} -> {written}"
    );
    written
}

#[test]
fn allene_tags_are_written_as_axial_tokens() {
    assert!(assert_round_trip("CC=[C@AL1]=CC").contains("@AL1"));
    assert!(assert_round_trip("CC=[C@AL2]=CC").contains("@AL2"));
}

#[test]
fn three_ligand_centres_keep_their_parity() {
    for input in [
        "C[S@](=O)CC",
        "O=[S@@](C)CC",
        "[S@](=O)(C)CC",
        "F[Pt@SP1](Cl)Br",
        "[Pt@SP2](F)(Cl)Br",
    ] {
        assert_round_trip(input);
    }
}

#[test]
fn multiple_charges_are_written_with_their_magnitude() {
    assert!(assert_round_trip("[O-2]").contains("-2"));
    assert!(assert_round_trip("[Fe+3]").contains("+3"));
}
