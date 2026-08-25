//! Regression tests for issue #399: several `standardize.rs` functions
//! rebuilt a molecule via a bare `MoleculeBuilder` without carrying
//! `stereo_neighbor_order`/`bond_directions`/`stereo_groups` forward, so any
//! surviving stereocenter/declared-E-Z-bond silently lost its recorded
//! order -- `remove_hydrogens`'s adjacency-based fallback then reconstructed
//! a *transposed* order for ring-opening stereocenters, flipping `@`/`@@` in
//! the canonical writer, on a subsequent canonicalize pass. Each test below
//! exercises one of the 8 fixed functions (or, for `disconnect_metals`, via
//! `standardize()` itself since it is a private pipeline-only stage).
//!
//! Note on methodology: a naive "canonical form must be textually unchanged
//! before/after `f`" check is wrong whenever `f` performs a real chemical
//! transformation (nitro normalization, deprotonation, isotope removal,
//! fragment/metal disconnection) -- those legitimately shift canonical rank
//! and can legitimately flip a *written* `@`/`@@` character while preserving
//! the same physical configuration, and a multi-fragment molecule's surviving
//! organic component can legitimately get different ring-closure numbering
//! than the same fragment parsed standalone. The property that actually
//! matters for #399 is idempotency: `canonicalize(f(x))` must equal
//! `canonicalize(f(reparse(canonicalize(f(x)))))` -- i.e. running `f` again
//! on the reparsed *output* of `f` must reproduce the same result. That is
//! what every test here checks.

use chematic_chem::{
    StandardizeOptions, neutralize_charges, normalize_groups, normalize_zwitterion, prefer_organic,
    reionize, remove_isotopes, standardize, uncharge,
};
use chematic_core::Molecule;
use chematic_smiles::{canonical_smiles, parse};

#[track_caller]
fn assert_f_idempotent(label: &str, smi: &str, f: impl Fn(&Molecule) -> Molecule) {
    let mol = parse(smi).unwrap_or_else(|e| panic!("{label}: parse '{smi}' failed: {e}"));
    let once = canonical_smiles(&f(&mol));
    let reparsed =
        parse(&once).unwrap_or_else(|e| panic!("{label}: re-parse '{once}' failed: {e}"));
    let twice = canonical_smiles(&f(&reparsed));
    assert_eq!(
        once, twice,
        "{label}: '{smi}' not idempotent through f+canonicalize: once='{once}' twice='{twice}'"
    );
}

#[test]
fn normalize_groups_idempotent_with_unrelated_stereocenter() {
    assert_f_idempotent(
        "normalize_groups",
        "O=[N+]([O-])c1ccc([C@@H](N)C)cc1",
        normalize_groups,
    );
}

#[test]
fn normalize_zwitterion_active_path_idempotent_with_unrelated_stereocenter() {
    // A genuine zwitterion (ammonium/carboxylate) with an unrelated ring
    // stereocenter elsewhere -- exercises the active proton-transfer path,
    // not the has_zwitterion() early-return.
    assert_f_idempotent(
        "normalize_zwitterion",
        "[NH3+][C@@H](CCCC(=O)[O-])c1ccccc1",
        normalize_zwitterion,
    );
}

#[test]
fn neutralize_charges_noop_preserves_stereocenter() {
    // The confirmed minimal repro for issue #399: zero charged atoms at all,
    // so neutralize_charges makes no modifications -- a true no-op, so
    // canonical form must be exactly unchanged (not just idempotent).
    let smi = "CN1CCC[C@H]1c1cccnc1";
    let mol = parse(smi).expect("parse");
    let before = canonical_smiles(&mol);
    let after = canonical_smiles(&neutralize_charges(&mol));
    assert_eq!(
        before, after,
        "neutralize_charges: true no-op case changed canonical form ({before} -> {after})"
    );
}

#[test]
fn neutralize_charges_idempotent_with_unrelated_charge() {
    assert_f_idempotent(
        "neutralize_charges (with unrelated charge)",
        "[O-]C(=O)c1ccc([C@@H](N)C)cc1",
        neutralize_charges,
    );
}

#[test]
fn remove_isotopes_idempotent_with_unrelated_stereocenter() {
    assert_f_idempotent(
        "remove_isotopes",
        "[13CH3][C@@H](N)c1ccccc1",
        remove_isotopes,
    );
}

#[test]
fn reionize_idempotent_with_unrelated_stereocenter() {
    assert_f_idempotent("reionize", "OC(=O)c1ccc([C@@H](N)C)cc1", reionize);
}

#[test]
fn uncharge_idempotent_with_unrelated_stereocenter() {
    assert_f_idempotent("uncharge", "[NH3+]CCc1ccc([C@@H](N)C)cc1", uncharge);
}

#[test]
fn prefer_organic_idempotent_with_kept_fragment_stereocenter() {
    assert_f_idempotent(
        "prefer_organic",
        "[Cl-].[NH3+][C@@H](C)c1ccccc1",
        prefer_organic,
    );
}

#[test]
fn disconnect_metals_path_idempotent_with_organic_stereocenter() {
    // disconnect_metals fires unconditionally inside standardize() whenever
    // any atom is a metal; the organic stereocenter it doesn't touch must
    // still survive a second standardize+canonicalize pass unperturbed.
    let opts = StandardizeOptions::default();
    let f = |m: &Molecule| standardize(m, &opts);
    assert_f_idempotent(
        "disconnect_metals (via standardize)",
        "[Na+].[O-]C(=O)[C@@H](N)c1ccccc1",
        f,
    );
}

/// Standardize+canonicalize round trip must stay idempotent regardless of
/// whether the stereocenter plays the ring-*closing* role (its ring bond is
/// added to physical adjacency synchronously at parse time) or the
/// ring-*opening* role (the bond is deferred until the matching digit is
/// parsed) in the original SMILES text -- issue #399's defect 2 mechanism
/// specifically depended on this role differing between the original text
/// and the canonical rewrite.
#[test]
fn standardize_idempotent_for_both_ring_closing_and_ring_opening_stereocenter_roles() {
    let opts = StandardizeOptions::default();
    let cases = [
        // Ring-closing role: N opens ring 1, [C@H] consumes (closes) it.
        ("ring-closing", "CN1CCC[C@H]1c1cccnc1"),
        // Ring-opening role: the stereocenter itself opens ring 1; N closes it.
        ("ring-opening", "[C@@H]1(c2cccnc2)CCCN1C"),
    ];
    for (label, smi) in cases {
        let mol = parse(smi).unwrap_or_else(|e| panic!("{label}: parse '{smi}' failed: {e}"));
        let once = canonical_smiles(&standardize(&mol, &opts));
        let reparsed =
            parse(&once).unwrap_or_else(|e| panic!("{label}: re-parse '{once}' failed: {e}"));
        let twice = canonical_smiles(&standardize(&reparsed, &opts));
        assert_eq!(
            once, twice,
            "{label} ('{smi}'): standardize+canonicalize not idempotent: once='{once}' twice='{twice}'"
        );
    }
}

/// Mirror-image distinctness must survive `standardize()`: a genuine (R)/(S)
/// pair (four distinct substituents -- unlike e.g. monosubstituted
/// cyclopentane, whose two ring arms are homotopic and so is not actually a
/// stereocenter, a pre-existing, unrelated characteristic of that fixture,
/// not a regression) must never collapse to the same canonical form.
#[test]
fn standardize_preserves_mirror_image_distinctness() {
    let opts = StandardizeOptions::default();
    let pairs = [
        ("N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O"), // L-/D-alanine
        ("CN1CCC[C@H]1c1cccnc1", "CN1CCC[C@@H]1c1cccnc1"),
        (
            "[NH3+][C@@H](CCCC(=O)[O-])c1ccccc1",
            "[NH3+][C@H](CCCC(=O)[O-])c1ccccc1",
        ),
    ];
    for (r, s) in pairs {
        let mr = parse(r).expect("parse R");
        let ms = parse(s).expect("parse S");
        let cr = canonical_smiles(&standardize(&mr, &opts));
        let cs = canonical_smiles(&standardize(&ms, &opts));
        assert_ne!(
            cr, cs,
            "mirror images '{r}' and '{s}' collapsed to the same canonical form after standardize: {cr}"
        );
    }
}

/// Atom-order permutation invariance through `standardize()`: the same
/// molecule written starting from different atoms (different input
/// traversal order) must converge to the same canonical form post-standardize.
///
/// Both alternative spellings below are independently verified (RDKit
/// `MolToInchi` on both spellings agrees) rather than hand-derived --
/// deriving stereo-tag parity for a rewritten traversal order by hand is
/// genuinely error-prone (an earlier draft of this test got the second
/// pair's tag backwards and would have asserted equality of two *opposite*
/// enantiomers had it not been checked against RDKit first).
#[test]
fn standardize_atom_order_permutation_invariance() {
    let opts = StandardizeOptions::default();
    let pairs = [
        ("CN1CCC[C@H]1c1cccnc1", "c1ccncc1[C@@H]1CCCN1C"),
        (
            "[NH3+][C@@H](CCCC(=O)[O-])c1ccccc1",
            "c1ccccc1[C@@H]([NH3+])CCCC(=O)[O-]",
        ),
    ];
    for (a, b) in pairs {
        let ma = parse(a).expect("parse a");
        let mb = parse(b).expect("parse b");
        let ca = canonical_smiles(&standardize(&ma, &opts));
        let cb = canonical_smiles(&standardize(&mb, &opts));
        assert_eq!(
            ca, cb,
            "permutation invariance: '{a}' and '{b}' (same molecule, different traversal) diverged after standardize: {ca} vs {cb}"
        );
    }
}
