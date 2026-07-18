//! EZ-S1 sibling gap: `standard_inchi()` (native-inchi feature) must read
//! the same stashed aromatic-bond direction `chematic_chem::assign_cip`'s
//! `substituent_is_up` was fixed to read (see `crates/chematic-chem/src/cip.rs`).
//! `crates/chematic-inchi/src/native/convert.rs`'s `find_stereo_sub` closure
//! mirrored the pre-fix bug exactly: it read only a bond's own `order`,
//! never `Molecule::bond_direction`, for InChI `Stereo0D` double-bond
//! descriptors -- so a mancude ring's exocyclic double bond geometry was
//! silently dropped from the InChI `/b` layer.
//!
//! `apply_kekule` (run on any molecule with aromatic bonds before InChI
//! conversion) preserves `bond_direction` verbatim on the same bond index
//! while only ever updating `order` -- so the stash is still the only place
//! the real direction lives even after kekulization turns the bond's own
//! `order` from `Aromatic` into `Single`/`Double`.
//!
//! Unlike `chematic-chem`'s own CIP E/Z label, InChI's Stereo0D format
//! doesn't need a CIP-priority-highest substituent -- any one determinate
//! substituent per alkene end is sufficient, since the parity is defined
//! relative to whichever specific neighbor is fed in (`inchi_api.h`'s 0D
//! stereo notes). So only the stash-read fix applies here; the
//! CIP-priority-fallback and tie-guard bugs found and fixed alongside
//! `substituent_is_up` in `chematic-chem` don't have an analogue in this
//! code (confirmed by reasoning about the format, not assumed).
//!
//! Full 5,000-molecule corpus diff (`ez_stash_inchi_snapshot` example):
//! before vs after this fix, exactly 12/5,000 molecules changed, all 12 a
//! pure `/b`-layer gain (0 lost, 0 changed, 0 non-`/b`-layer side effects),
//! and all 12 post-fix outputs match RDKit's `Chem.MolToInchi` byte-for-byte
//! -- confirmed by direct comparison, not assumed from the chematic-chem
//! corpus numbers.
//!
//! Reference InChI strings below were independently verified against
//! `Chem.MolToInchi` (RDKit, which wraps the same underlying IUPAC InChI
//! reference C library `standard_inchi` links against).

#![cfg(feature = "native-inchi")]

use chematic_inchi::standard_inchi;
use chematic_smiles::parse;

#[track_caller]
fn assert_inchi(smi: &str, expected: &str) {
    let mol = parse(smi).unwrap_or_else(|e| panic!("parse {smi:?}: {e}"));
    let got = standard_inchi(&mol).unwrap_or_else(|e| panic!("standard_inchi {smi:?}: {e:?}"));
    assert_eq!(got, expected, "InChI mismatch for {smi:?}");
}

#[test]
fn ez_stash_ring_gains_b_layer() {
    // A squarate-diimine mancude ring (both exocyclic imines substituted,
    // unlike the bare `=NH` D0 fixture, which has no real substituent to
    // be stereogenic about at all). RDKit-confirmed `/b8-4-,9-5-`.
    assert_inchi(
        r"CC/N=c1\c(O)c(O)\c1=N/C",
        "InChI=1S/C7H10N2O2/c1-3-9-5-4(8-2)6(10)7(5)11/h10-11H,3H2,1-2H3/b8-4-,9-5-",
    );
}

#[test]
fn ez_stash_ring_inverted_pair_differs() {
    // Same structure as above with just one ring-imine's slash inverted --
    // RDKit-confirmed the b-layer digit for that bond flips sign (`9-5+`
    // instead of `9-5-`), while the other bond's descriptor is unchanged.
    // A "the two just differ" check alone would also pass under a global
    // sign-inversion bug, so both this and the previous test pin specific
    // values.
    assert_inchi(
        r"CC/N=c1/c(O)c(O)\c1=N/C",
        "InChI=1S/C7H10N2O2/c1-3-9-5-4(8-2)6(10)7(5)11/h10-11H,3H2,1-2H3/b8-4-,9-5+",
    );
    assert_inchi(
        r"CC/N=c1\c(O)c(O)\c1=N\C",
        "InChI=1S/C7H10N2O2/c1-3-9-5-4(8-2)6(10)7(5)11/h10-11H,3H2,1-2H3/b8-4+,9-5-",
    );
}

#[test]
fn ez_stash_mobile_h_tautomer_ring_stays_undefined() {
    // A pyridine-imine mancude ring with a ring N-H: standard InChI's
    // mobile-H tautomer normalization legitimately merges this bond's
    // geometry into the `(H,7,8)` tautomer layer -- RDKit produces the
    // identical InChI (no `/b` layer) for both slash directions, confirmed
    // independently via `Chem.MolToInchi`. This is not something this fix
    // changes or should change; pinned so a future, unrelated change to
    // tautomer handling doesn't silently start fabricating a `/b` layer
    // InChI itself doesn't define for this structural class.
    let no_b = "InChI=1S/C6H8N2/c1-7-6-4-2-3-5-8-6/h2-5H,1H3,(H,7,8)";
    assert_inchi(r"C/N=c1\cccc[nH]1", no_b);
    assert_inchi(r"C/N=c1/cccc[nH]1", no_b);
}
