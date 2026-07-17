//! invert_stereocenter-correctness-P0: `chematic_chem::invert_stereocenter`'s
//! chemical output, verified through `chematic_inchi::inchi` -- a code path
//! that never reads `Molecule::stereo_neighbor_order` at all (confirmed in
//! the Kekule-S0 InChI suite), so agreement here is an independent
//! cross-check on top of the CIP- and canonical-SMILES-based checks in
//! `chematic-chem`'s own unit tests, not a re-statement of them.
//!
//! Ground truth is always an independently-parsed SMILES string (the known
//! mirror form), never a comparison of `@`/`@@` characters in the input --
//! per review, "chemical identity, not string comparison."

use chematic_chem::invert_stereocenter;
use chematic_core::{AtomIdx, Chirality, Molecule};
use chematic_inchi::inchi;
use chematic_smiles::parse;

fn find_by_map(m: &Molecule, map_num: u16) -> AtomIdx {
    m.atoms()
        .find(|(_, a)| a.atom_map == Some(map_num))
        .map(|(idx, _)| idx)
        .expect("atom map tag not found")
}

#[test]
fn invert_stereocenter_matches_independent_enantiomer_inchi() {
    let l_ala = parse("N[C@@H](C)C(=O)O").expect("valid SMILES");
    let target = l_ala
        .atoms()
        .find(|(_, a)| a.chirality != Chirality::None)
        .expect("test setup sanity: has a stereocenter")
        .0;

    let inverted = invert_stereocenter(&l_ala, target);
    let d_ala = parse("N[C@H](C)C(=O)O").expect("valid SMILES");

    assert_eq!(
        inchi(&inverted),
        inchi(&d_ala),
        "inverting L-alanine's stereocenter must produce the same InChI as \
         independently-parsed D-alanine, including the /t stereo layer"
    );
}

#[test]
fn invert_stereocenter_flips_only_target_center_inchi() {
    // Two independent, both-resolvable stereocenters (each carries an
    // implicit H). Inverting only the map-1 center must: (a) match an
    // independently-parsed molecule with only that center's @/@@ flipped,
    // and (b) leave the non-stereo InChI layers (connectivity/H counts)
    // byte-identical to the uninverted original -- proving nothing but the
    // one center's configuration changed.
    let base = "C[C@H:1](O)[C@@H:2](N)C(=O)O";
    let m = parse(base).expect("valid SMILES");
    let target = find_by_map(&m, 1);

    let inverted = invert_stereocenter(&m, target);
    let expected = parse("C[C@@H:1](O)[C@@H:2](N)C(=O)O").expect("valid SMILES");

    assert_eq!(
        inchi(&inverted),
        inchi(&expected),
        "inverting only the mapped-1 center must match an independently-parsed \
         molecule with only that center's @/@@ flipped"
    );

    fn non_stereo_layers(inchi_str: &str) -> &str {
        inchi_str.find("/t").map_or(inchi_str, |i| &inchi_str[..i])
    }
    assert_eq!(
        non_stereo_layers(&inchi(&m)),
        non_stereo_layers(&inchi(&inverted)),
        "connectivity/H-count layers must be unaffected by inverting one stereocenter"
    );
}

#[test]
fn invert_stereocenter_double_inversion_matches_original_inchi() {
    let m = parse("N[C@@H](C)C(=O)O").expect("valid SMILES");
    let target = m
        .atoms()
        .find(|(_, a)| a.chirality != Chirality::None)
        .expect("test setup sanity: has a stereocenter")
        .0;

    let once = invert_stereocenter(&m, target);
    let twice = invert_stereocenter(&once, target);

    assert_eq!(
        inchi(&twice),
        inchi(&m),
        "inverting the same center twice must reproduce the original molecule's \
         full InChI string (not just its stereo layer)"
    );
}
