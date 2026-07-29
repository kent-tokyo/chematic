//! Dative (`->` / `<-`) bond direction, end to end (issue #194).
//!
//! `BondOrder::Dative` stores donor→acceptor as `atom1`→`atom2`. Two
//! independent things have to respect that:
//!
//! * the **parser**, which must map `A<-B` to `atom1 = B` (B is the donor),
//!   not to the same bond as `A->B`; and
//! * the **writers**, which must emit the whole two-character token *and*
//!   pick `->` vs `<-` according to which endpoint the DFS happens to write
//!   first.
//!
//! Expectations below are pinned to hand-written SMILES literals and to
//! molecules built through `MoleculeBuilder` (never to this crate's own
//! output), so a writer and a parser that were wrong in the *same* direction
//! could not make them pass.

use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};
use chematic_smiles::{canonical_smiles, parse, write};

/// Index of the single dative bond, plus its stored (donor, acceptor)
/// elements.
fn dative_donor_acceptor(mol: &chematic_core::Molecule) -> (Element, Element) {
    let bidx = (0..mol.bond_count() as u32)
        .map(chematic_core::BondIdx)
        .find(|&b| mol.bond(b).order == BondOrder::Dative)
        .expect("molecule has a dative bond");
    let bond = mol.bond(bidx);
    (mol.atom(bond.atom1).element, mol.atom(bond.atom2).element)
}

fn atom(element: Element) -> Atom {
    Atom::new(element)
}

/// N donates to Fe: `atom1` = N, `atom2` = Fe.
fn n_to_fe() -> chematic_core::Molecule {
    let mut b = MoleculeBuilder::new();
    let n = b.add_atom(atom(Element::N));
    let fe = b.add_atom(atom(Element::FE));
    b.add_bond(n, fe, BondOrder::Dative).unwrap();
    b.build()
}

/// Same graph, opposite stored direction: Fe donates to N.
fn fe_to_n() -> chematic_core::Molecule {
    let mut b = MoleculeBuilder::new();
    let n = b.add_atom(atom(Element::N));
    let fe = b.add_atom(atom(Element::FE));
    b.add_bond(fe, n, BondOrder::Dative).unwrap();
    b.build()
}

// ── parser ──────────────────────────────────────────────────────────────

/// `N->[Fe]` and `N<-[Fe]` describe *opposite* donor/acceptor relationships
/// and must not parse to the same bond. Before the fix both arrows produced
/// `atom1 = N`, silently reversing the second one.
#[test]
fn parser_distinguishes_arrow_direction() {
    let forward = parse("N->[Fe]").expect("N->[Fe]");
    let backward = parse("N<-[Fe]").expect("N<-[Fe]");

    assert_eq!(
        dative_donor_acceptor(&forward),
        (Element::N, Element::FE),
        "`N->[Fe]`: nitrogen donates to iron"
    );
    assert_eq!(
        dative_donor_acceptor(&backward),
        (Element::FE, Element::N),
        "`N<-[Fe]`: iron donates to nitrogen"
    );
    assert_ne!(
        dative_donor_acceptor(&forward),
        dative_donor_acceptor(&backward),
        "the two arrows must not collapse to the same bond"
    );
}

/// Same check one level in: the arrow inside a branch.
#[test]
fn parser_arrow_direction_in_branch() {
    let forward = parse("[Fe](Cl)(Cl)<-N").expect("<- in chain after branches");
    assert_eq!(
        dative_donor_acceptor(&forward),
        (Element::N, Element::FE),
        "`<-N` written from Fe means N is the donor"
    );
    let branch = parse("[Fe](<-N)Cl").expect("<- as the first token of a branch");
    assert_eq!(
        dative_donor_acceptor(&branch),
        (Element::N, Element::FE),
        "`(<-N)` written from Fe means N is the donor"
    );
}

/// A dative bond spelled as a ring closure, in both directions, marked at
/// the opening occurrence.
#[test]
fn parser_arrow_direction_at_ring_closure() {
    // Ring: N-C-C-[Fe], closed by a dative bond between N (opener) and Fe.
    let donor_first = parse("N->1CC[Fe]1").expect("N->1CC[Fe]1");
    assert_eq!(
        dative_donor_acceptor(&donor_first),
        (Element::N, Element::FE),
        "`N->1 … [Fe]1`: the arrow at the opener points N→Fe"
    );
    let acceptor_first = parse("N<-1CC[Fe]1").expect("N<-1CC[Fe]1");
    assert_eq!(
        dative_donor_acceptor(&acceptor_first),
        (Element::FE, Element::N),
        "`N<-1 … [Fe]1`: the arrow at the opener points Fe→N"
    );
}

/// The marker written at the *closing* occurrence reads close→open, so it
/// must be flipped before it is compared with the opener's marker. Both
/// spellings of the same ring bond must agree rather than raise
/// `ConflictingRingBond`.
#[test]
fn parser_arrow_direction_at_ring_close_occurrence() {
    let at_close = parse("N1CC[Fe]<-1").expect("N1CC[Fe]<-1");
    assert_eq!(
        dative_donor_acceptor(&at_close),
        (Element::N, Element::FE),
        "`[Fe]<-1` closing back to N means N is the donor"
    );
    let both_ends = parse("N->1CC[Fe]<-1").expect("consistent markers at both ring ends");
    assert_eq!(
        dative_donor_acceptor(&both_ends),
        (Element::N, Element::FE),
        "`->` at the opener and `<-` at the closer describe the same N→Fe bond"
    );

    // …and two markers that genuinely contradict each other are rejected,
    // exactly as a `/` vs `/` ring-bond conflict already is, instead of one
    // silently winning.
    assert!(
        parse("N->1CC[Fe]->1").is_err(),
        "`->` at both ring ends claims both atoms are the donor"
    );
}

// ── writers ─────────────────────────────────────────────────────────────

/// The whole `->`/`<-` token survives both writers. Before the fix,
/// `canonical_smiles` truncated it to a bare `-` (issue #194) and `write`
/// emitted `->` regardless of direction.
#[test]
fn writers_emit_the_full_two_character_token() {
    for mol in [n_to_fe(), fe_to_n()] {
        for out in [write(&mol), canonical_smiles(&mol)] {
            assert!(
                out.contains("->") || out.contains("<-"),
                "dative token truncated or dropped: {out}"
            );
            // A bare `-` adjacent to nothing else would be the truncation bug.
            let arrows = out.matches("->").count() + out.matches("<-").count();
            assert_eq!(arrows, 1, "expected exactly one dative arrow in {out}");
        }
    }
}

/// The critical case: two molecules that differ *only* in which endpoint of
/// the dative bond is stored as `atom1`. Their canonical rank vectors are
/// identical (same elements, same connectivity — nothing in the ranking is
/// dative-direction-aware), so the canonical writer reaches the same atom
/// first in both. It picks the iron; for `n_to_fe` that is the *acceptor*,
/// so the arrow has to be written backwards relative to storage order.
///
/// Both expected strings are derived from `BondOrder::Dative`'s own
/// definition rather than from what the writer happens to produce: reading
/// `[Fe]<-N` left to right, the arrow points from N into Fe, i.e. N is the
/// donor — which is `n_to_fe`, whose `atom1` is N even though it is written
/// second.
#[test]
fn canonical_writer_flips_the_arrow_when_the_acceptor_is_written_first() {
    let a = canonical_smiles(&n_to_fe());
    let b = canonical_smiles(&fe_to_n());

    assert_eq!(a, "[Fe]<-N", "acceptor written first ⇒ reversed arrow");
    assert_eq!(b, "[Fe]->N", "donor written first ⇒ forward arrow");

    // And each one still means what it meant before it was written.
    assert_eq!(
        dative_donor_acceptor(&parse(&a).unwrap()),
        (Element::N, Element::FE)
    );
    assert_eq!(
        dative_donor_acceptor(&parse(&b).unwrap()),
        (Element::FE, Element::N)
    );
}

/// Full round trip through both writers, for both stored directions.
#[test]
fn dative_direction_round_trips_through_both_writers() {
    for (mol, expected) in [
        (n_to_fe(), (Element::N, Element::FE)),
        (fe_to_n(), (Element::FE, Element::N)),
    ] {
        assert_eq!(
            dative_donor_acceptor(&mol),
            expected,
            "builder precondition"
        );
        for out in [write(&mol), canonical_smiles(&mol)] {
            let reparsed = parse(&out).unwrap_or_else(|e| panic!("re-parse of {out}: {e}"));
            assert_eq!(
                dative_donor_acceptor(&reparsed),
                expected,
                "donor/acceptor lost writing {out}"
            );
        }
    }
}

/// Same round trip for a dative bond that is a ring-closure back-edge
/// rather than a tree edge, in both stored directions.
#[test]
fn dative_ring_closure_round_trips_through_both_writers() {
    for donor_is_n in [true, false] {
        let mut b = MoleculeBuilder::new();
        let n = b.add_atom(atom(Element::N));
        let c1 = b.add_atom(atom(Element::C));
        let c2 = b.add_atom(atom(Element::C));
        let fe = b.add_atom(atom(Element::FE));
        b.add_bond(n, c1, BondOrder::Single).unwrap();
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        b.add_bond(c2, fe, BondOrder::Single).unwrap();
        let (donor, acceptor) = if donor_is_n { (n, fe) } else { (fe, n) };
        b.add_bond(donor, acceptor, BondOrder::Dative).unwrap();
        let mol = b.build();

        let expected = if donor_is_n {
            (Element::N, Element::FE)
        } else {
            (Element::FE, Element::N)
        };
        assert_eq!(
            dative_donor_acceptor(&mol),
            expected,
            "builder precondition"
        );

        for out in [write(&mol), canonical_smiles(&mol)] {
            // Both writers really do route this bond through a ring-closure
            // digit rather than a tree edge, and both ends of it print: the
            // donor end `->`, the acceptor end `<-`, each immediately before
            // its ring number.
            assert!(
                out.contains("->1") && out.contains("<-1"),
                "ring-closure dative token truncated or not on the ring digit: {out}"
            );
            let reparsed = parse(&out).unwrap_or_else(|e| panic!("re-parse of {out}: {e}"));
            assert_eq!(
                dative_donor_acceptor(&reparsed),
                expected,
                "donor/acceptor lost writing ring closure {out}"
            );
        }
    }
}

/// The crate's central invariant (see the `lib.rs` doctest): re-parsing a
/// canonical SMILES and writing it again reproduces the same string.
/// Nothing else covers it for dative bonds — the ChEMBL corpus fixtures
/// contain none — and the ring spelling (`N->1CC[Fe]<-1`, arrows on both
/// ring digits) is one nothing had round-tripped before this change.
#[test]
fn canonical_smiles_is_stable_for_dative_molecules() {
    for smiles in [
        "N->[Fe]",
        "[Fe]<-N",
        "N->1CC[Fe]<-1",
        "N<-1CC[Fe]->1",
        "N1CC[Fe]<-1",
        "[Fe](Cl)(Cl)<-N",
    ] {
        let canon = canonical_smiles(&parse(smiles).unwrap());
        assert_eq!(
            canonical_smiles(&parse(&canon).unwrap()),
            canon,
            "canonical SMILES of {smiles} is not stable"
        );
        assert!(
            canon.contains("->") || canon.contains("<-"),
            "dative bond lost canonicalizing {smiles}: {canon}"
        );
    }
}

/// Guard against misreading the test above: a molecule built through
/// `MoleculeBuilder` (rather than parsed) can shift once on its first
/// canonicalization, because the builder and the parser do not populate
/// every atom field identically. That is pre-existing and has nothing to do
/// with dative bonds — the same ring with an ordinary single bond in place
/// of the dative one shifts in exactly the same way — so it must not be
/// "fixed" by reshaping the dative writer.
#[test]
fn builder_first_shift_is_not_dative_specific() {
    let ring = |closing: BondOrder| {
        let mut b = MoleculeBuilder::new();
        let n = b.add_atom(atom(Element::N));
        let c1 = b.add_atom(atom(Element::C));
        let c2 = b.add_atom(atom(Element::C));
        let fe = b.add_atom(atom(Element::FE));
        b.add_bond(n, c1, BondOrder::Single).unwrap();
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        b.add_bond(c2, fe, BondOrder::Single).unwrap();
        b.add_bond(n, fe, closing).unwrap();
        let mol = b.build();
        let first = canonical_smiles(&mol);
        let second = canonical_smiles(&parse(&first).unwrap());
        (first, second)
    };

    let (single_first, single_second) = ring(BondOrder::Single);
    let (dative_first, dative_second) = ring(BondOrder::Dative);
    assert_ne!(
        single_first, single_second,
        "control: the shift is expected without any dative bond"
    );
    assert_ne!(dative_first, dative_second, "same shift, same cause");
    // And in both cases it settles after one parse.
    assert_eq!(
        canonical_smiles(&parse(&single_second).unwrap()),
        single_second
    );
    assert_eq!(
        canonical_smiles(&parse(&dative_second).unwrap()),
        dative_second
    );
}
