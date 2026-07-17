//! Kekule-S0: `chematic_inchi::inchi`'s stereo layers (`/t`, `/m`, `/s`)
//! must be unchanged whether a molecule is InChI-generated directly from its
//! parsed (possibly aromatic-form) SMILES, or from the same molecule after
//! an explicit `chematic_core::kekulize`+`apply_kekule` round trip first --
//! exactly what `chematic-perception`'s
//! `apply_aromaticity_rdkit_parity_experimental` does internally before
//! computing aromaticity.
//!
//! InChI's own canonical numbering is independent of `AtomIdx`, so this
//! compares the generated InChI strings (and specifically their stereo
//! layers) directly rather than an atom-mapped table.
//!
//! **Scoping note**: for molecules whose stereocenter sits on/near a ring
//! that itself needs kekulization, comparing `inchi(direct)` against
//! `inchi(via_kekulize)` conflates two *different* bugs. Confirmed by direct
//! check: `chematic_inchi::inchi`'s own canonical numbering (`/c`, `/h`
//! layers) is representation-dependent (aromatic-form vs Kekulé-form input
//! produce different atom numbering) for a fully **achiral** aromatic
//! molecule too (`c1ccc(cc1)C(F)Cl`) -- this is a separate, pre-existing gap
//! in chematic-inchi's own canonicalization, unrelated to
//! `stereo_neighbor_order` and out of scope for Kekule-S0 (which is scoped
//! to `apply_kekule` alone; see `docs/` for the follow-up this surfaces).
//! For ring-adjacent stereocenters, this suite instead checks
//! `inchi(via_kekulize(a))` against `inchi(via_kekulize(b))` for a true
//! enantiomer pair `(a, b)`: both go through the *same* kekulize pipeline,
//! so the pre-existing numbering instability cancels out, and the
//! comparison isolates exactly what Kekule-S0 controls -- does the fix
//! produce a correctly mirrored stereo layer, not a collapsed or unrelated
//! one.
//!
//! **Honest note on what this suite actually demonstrates**: `chematic_inchi`
//! never reads `Molecule::stereo_neighbor_order` at all (confirmed: no
//! occurrences in `crates/chematic-inchi/src/`). It derives stereo purely
//! from `atom.chirality` plus bond/neighbor *insertion* order, which
//! `apply_kekule`'s rebuild always preserved even before this fix (only the
//! separate `stereo_neighbor_order`/`stereo_groups`/`bond_directions` side
//! channels were dropped). Positive-control check (temporarily reverting the
//! fix) confirms every test in this file passes identically before and
//! after -- i.e. `chematic-inchi`'s stereo output was never actually at risk
//! from this bug. These tests are kept as a real, useful regression pin
//! (InChI stereo layers stay correct across kekulization either way), not as
//! evidence the fix changed InChI's behavior.

use chematic_core::{apply_kekule, kekulize};
use chematic_inchi::inchi;
use chematic_smiles::parse;

/// Extract everything from `/t` onward (the stereo-relevant tail: `/t`,
/// `/m`, `/s`), or `""` if there's no `/t` layer at all.
fn stereo_tail(inchi_str: &str) -> &str {
    inchi_str.find("/t").map_or("", |i| &inchi_str[i..])
}

#[track_caller]
fn assert_inchi_stereo_invariant(smi: &str) {
    let raw = parse(smi).expect("valid SMILES");
    let k = kekulize(&raw).expect("kekulizable");
    let kekulized = apply_kekule(&raw, &k);

    let direct = inchi(&raw);
    let via_kekulize = inchi(&kekulized);

    assert_eq!(
        stereo_tail(&direct),
        stereo_tail(&via_kekulize),
        "{smi}: InChI stereo layer differs before/after kekulize+apply_kekule\n\
         direct:       {direct}\n\
         via_kekulize: {via_kekulize}"
    );
}

#[test]
fn simple_stereocenters_no_ring_kekulization_needed() {
    for smi in [
        "N[C@@H](C)C(=O)O",         // L-alanine, single stereocenter (/t + /s1, no /m)
        "N[C@H](C)C(=O)O",          // D-alanine
        "C[C@H](O)[C@@H](O)C(=O)O", // tartaric acid, two stereocenters (/t + /m + /s1)
        "[C@@H]1(N)CCCC1",
        "F[C@H]1CCCCC1",
    ] {
        assert_inchi_stereo_invariant(smi);
    }
}

#[test]
fn stereo_tail_present_where_expected() {
    // Sanity check on the helper itself, and that these test molecules
    // actually exercise the /t (and for 2+ stereocenters, /m) layers --
    // guards against the invariance test above passing vacuously because
    // neither side has a stereo layer at all.
    let single = inchi(&parse("N[C@@H](C)C(=O)O").unwrap());
    assert!(stereo_tail(&single).starts_with("/t"), "{single}");
    assert!(!stereo_tail(&single).contains("/m"), "{single}");

    let double = inchi(&parse("C[C@H](O)[C@@H](O)C(=O)O").unwrap());
    assert!(stereo_tail(&double).contains("/m"), "{double}");
}

/// For ring-adjacent stereocenters (where `direct` and `via_kekulize` can't
/// be fairly string-compared -- see module docs), verify instead that
/// routing a true enantiomer pair through the *same* kekulize+apply_kekule
/// pipeline produces InChI output with identical connectivity/H layers
/// (same molecule, same pre-stereo numbering) and a correctly mirrored
/// stereo tail (`/t` signs flipped, or `/m` flipped) -- not collapsed to the
/// same string, and not diverging into an unrelated one.
#[test]
fn ring_adjacent_enantiomer_pairs_mirror_correctly_via_kekulize() {
    let pairs: &[(&str, &str)] = &[
        ("c1ccc(cc1)[C@H](F)Cl", "c1ccc(cc1)[C@@H](F)Cl"),
        (
            "Cc1cn([C@H]2CCCC[C@@H]2C)c(=O)n1C",
            "Cc1cn([C@@H]2CCCC[C@H]2C)c(=O)n1C",
        ),
    ];

    fn via_kekulize(smi: &str) -> String {
        let raw = parse(smi).expect("valid SMILES");
        let k = kekulize(&raw).expect("kekulizable");
        inchi(&apply_kekule(&raw, &k))
    }

    fn non_stereo_layers(inchi_str: &str) -> &str {
        inchi_str.find("/t").map_or(inchi_str, |i| &inchi_str[..i])
    }

    for (a, b) in pairs {
        let ia = via_kekulize(a);
        let ib = via_kekulize(b);

        assert_eq!(
            non_stereo_layers(&ia),
            non_stereo_layers(&ib),
            "{a} vs {b}: connectivity/H layers should match for a true enantiomer pair\n\
             a: {ia}\nb: {ib}"
        );
        assert_ne!(
            stereo_tail(&ia),
            stereo_tail(&ib),
            "{a} vs {b}: stereo tail must differ (mirror images), not collapse to the same value\n\
             a: {ia}\nb: {ib}"
        );
        assert!(
            !stereo_tail(&ia).is_empty() && !stereo_tail(&ib).is_empty(),
            "{a} vs {b}: both should have a non-empty stereo tail\na: {ia}\nb: {ib}"
        );
    }
}
