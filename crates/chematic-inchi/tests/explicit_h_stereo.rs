//! Regression tests for two bugs in explicit-graph-hydrogen handling
//! (`crates/chematic-inchi/src/native/convert.rs`):
//!
//! 1. Isotope-tallying bug: an explicit graph H atom (`[2H]`, `[3H]`, as opposed
//!    to bracket-H notation like `[C@H]`) had its `isotope` field ignored when
//!    tallying `num_iso_h` -- every explicit H was counted as ordinary H
//!    regardless of isotope, so no `/i` layer was ever produced for D/T.
//! 2. Stereo0D-drop bug: when a tetrahedral stereocentre's 4th substituent is a
//!    real explicit graph H atom (any isotope, including plain `[H]`), the
//!    whole stereo descriptor silently disappeared -- no `/t` layer for either
//!    enantiomer, making them InChI-indistinguishable.
//!
//! Expected values below were independently verified against RDKit
//! (`rdkit==2026.03.3`, isolated venv) via `Chem.MolToInchi`, confirmed
//! byte-exact -- see the PR body for the generation script and the venv path.

#![cfg(feature = "native-inchi")]

use chematic_inchi::standard_inchi;
use chematic_smiles::parse;

// (label, SMILES, expected standard InChI -- RDKit 2026.03.3 byte-exact)
const FIXTURES: &[(&str, &str, &str)] = &[
    (
        "plain_H_R",
        "[C@](Br)(Cl)(F)[H]",
        "InChI=1S/CHBrClF/c2-1(3)4/h1H/t1-/m1/s1",
    ),
    (
        "plain_H_S",
        "[C@@](Br)(Cl)(F)[H]",
        "InChI=1S/CHBrClF/c2-1(3)4/h1H/t1-/m0/s1",
    ),
    (
        "D_R",
        "[C@](Br)(Cl)(F)[2H]",
        "InChI=1S/CHBrClF/c2-1(3)4/h1H/t1-/m1/s1/i1D",
    ),
    (
        "D_S",
        "[C@@](Br)(Cl)(F)[2H]",
        "InChI=1S/CHBrClF/c2-1(3)4/h1H/t1-/m0/s1/i1D",
    ),
    (
        "T_R",
        "[C@](Br)(Cl)(F)[3H]",
        "InChI=1S/CHBrClF/c2-1(3)4/h1H/t1-/m1/s1/i1T",
    ),
    (
        "T_S",
        "[C@@](Br)(Cl)(F)[3H]",
        "InChI=1S/CHBrClF/c2-1(3)4/h1H/t1-/m0/s1/i1T",
    ),
    // bracket-H control: already-working common case, must stay unchanged.
    (
        "bracketH_control",
        "[C@H](Br)(Cl)F",
        "InChI=1S/CHBrClF/c2-1(3)4/h1H/t1-/m0/s1",
    ),
    (
        "bracketH_control_mirror",
        "[C@@H](Br)(Cl)F",
        "InChI=1S/CHBrClF/c2-1(3)4/h1H/t1-/m1/s1",
    ),
    // heavy-atom isotope control: different code path (isotopic_mass on a
    // heavy atom), must stay unchanged.
    (
        "heavy_isotope_control",
        "[13CH3]C",
        "InChI=1S/C2H6/c1-2/h1-2H3/i1+1",
    ),
    // combined charge + isotope + stereo.
    (
        "combined_charge_isotope_stereo_R",
        "[NH3+][C@](Br)([2H])C(=O)[O-]",
        "InChI=1S/C2H4BrNO2/c3-1(4)2(5)6/h1H,4H2,(H,5,6)/t1-/m1/s1/i1D",
    ),
    (
        "combined_charge_isotope_stereo_S",
        "[NH3+][C@@](Br)([2H])C(=O)[O-]",
        "InChI=1S/C2H4BrNO2/c3-1(4)2(5)6/h1H,4H2,(H,5,6)/t1-/m0/s1/i1D",
    ),
    // Atom-renumbering invariance: same absolute configuration as plain_H_R
    // (RDKit-confirmed: full reversal of a 4-element neighbor list is an even
    // permutation, so the descriptor is preserved), atoms declared in a
    // different order (H first, Cl/Br/F re-slotted).
    (
        "plain_H_R_renumbered",
        "[H][C@](Cl)(Br)F",
        "InChI=1S/CHBrClF/c2-1(3)4/h1H/t1-/m1/s1",
    ),
    // Bond-declaration-order-reversal invariance: full reversal of the
    // substituent list order (Br,Cl,F,H -> H,F,Cl,Br), also RDKit-confirmed
    // to be the same configuration as plain_H_R.
    (
        "plain_H_R_bond_reversed",
        "[C@]([H])(F)(Cl)Br",
        "InChI=1S/CHBrClF/c2-1(3)4/h1H/t1-/m1/s1",
    ),
    // Explicit protium ([1H]): a distinct isotope tag from plain (untagged) H,
    // even though 1 is hydrogen's natural/most-abundant mass number -- InChI
    // still emits a dedicated /i layer for an explicitly-declared isotope.
    (
        "protium_R",
        "[C@](Br)(Cl)(F)[1H]",
        "InChI=1S/CHBrClF/c2-1(3)4/h1H/t1-/m1/s1/i1H",
    ),
    (
        "protium_S",
        "[C@@](Br)(Cl)(F)[1H]",
        "InChI=1S/CHBrClF/c2-1(3)4/h1H/t1-/m0/s1/i1H",
    ),
    // Two independent stereocentres in one molecule, each with its own
    // explicit-graph-H substituent of a different isotope (D and T) -- checks
    // that the two manufactured Stereo0D stand-in atoms don't collide/cross-talk.
    (
        "two_stereocenters_D_T",
        "[C@](Br)(Cl)([2H])[C@@]([3H])(F)I",
        "InChI=1S/C2H2BrClFI/c3-1(4)2(5)6/h1-2H/t1-,2-/m0/s1/i1D,2T",
    ),
    (
        "two_stereocenters_D_T_mirror",
        "[C@@](Br)(Cl)([2H])[C@]([3H])(F)I",
        "InChI=1S/C2H2BrClFI/c3-1(4)2(5)6/h1-2H/t1-,2-/m1/s1/i1D,2T",
    ),
];

#[test]
fn explicit_h_fixtures_match_rdkit_byte_exact() {
    let mut failures = Vec::new();
    for &(label, smiles, expected) in FIXTURES {
        let mol = parse(smiles).unwrap_or_else(|e| panic!("parse {smiles:?}: {e}"));
        let got =
            standard_inchi(&mol).unwrap_or_else(|e| panic!("standard_inchi({smiles:?}): {e}"));
        if got != expected {
            failures.push(format!(
                "[{label}] SMILES {smiles:?}\n  got:      {got}\n  expected: {expected}"
            ));
        }
    }
    if !failures.is_empty() {
        panic!(
            "{}/{} fixtures mismatched RDKit reference:\n{}",
            failures.len(),
            FIXTURES.len(),
            failures.join("\n\n")
        );
    }
}

/// Core bug-2 acceptance criterion: the enantiomer pairs must differ. Before
/// the fix, ALL of these collapsed to the identical InChI string (no /t layer
/// at all), making enantiomers indistinguishable.
#[test]
fn explicit_h_enantiomer_pairs_differ() {
    let pairs: &[(&str, &str)] = &[
        ("[C@](Br)(Cl)(F)[H]", "[C@@](Br)(Cl)(F)[H]"),
        ("[C@](Br)(Cl)(F)[2H]", "[C@@](Br)(Cl)(F)[2H]"),
        ("[C@](Br)(Cl)(F)[3H]", "[C@@](Br)(Cl)(F)[3H]"),
        (
            "[NH3+][C@](Br)([2H])C(=O)[O-]",
            "[NH3+][C@@](Br)([2H])C(=O)[O-]",
        ),
    ];
    for &(r_smi, s_smi) in pairs {
        let r_mol = parse(r_smi).unwrap();
        let s_mol = parse(s_smi).unwrap();
        let r_inchi = standard_inchi(&r_mol).unwrap();
        let s_inchi = standard_inchi(&s_mol).unwrap();
        assert_ne!(
            r_inchi, s_inchi,
            "enantiomers {r_smi:?} / {s_smi:?} must produce different InChI strings"
        );
        // And both must actually carry a /t (tetrahedral stereo) layer -- the
        // failure mode was silently DROPPING the layer, not producing a wrong one.
        assert!(
            r_inchi.contains("/t"),
            "expected a /t layer in {r_inchi:?} (from {r_smi:?})"
        );
        assert!(
            s_inchi.contains("/t"),
            "expected a /t layer in {s_inchi:?} (from {s_smi:?})"
        );
    }
}

/// Core bug-1 acceptance criterion: D/T must produce an /i layer distinct
/// from plain H and from each other.
#[test]
fn explicit_h_isotope_layers_distinguished() {
    let plain = standard_inchi(&parse("[C@](Br)(Cl)(F)[H]").unwrap()).unwrap();
    let d = standard_inchi(&parse("[C@](Br)(Cl)(F)[2H]").unwrap()).unwrap();
    let t = standard_inchi(&parse("[C@](Br)(Cl)(F)[3H]").unwrap()).unwrap();
    assert!(
        !plain.contains("/i"),
        "plain H must have no /i layer: {plain}"
    );
    assert!(d.contains("/i1D"), "D must have an /i1D layer: {d}");
    assert!(t.contains("/i1T"), "T must have an /i1T layer: {t}");
    assert_ne!(d, t);
    assert_ne!(d, plain);
    assert_ne!(t, plain);
}
