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

/// `n_to_fe`/`fe_to_n` differ *only* in which endpoint of the dative bond is
/// stored as `atom1`, but their canonical rank vectors are no longer
/// identical: the donor's implicit-H count is part of its ranking
/// invariant (`initial_invariant` in `canonical.rs`), and a donor's
/// implicit H count now correctly excludes the dative bond's own
/// contribution (see `chematic_core::valence::valence_inferred_hcount` --
/// found and fixed via the platinum coordination-chemistry benchmark,
/// `validation/platinum/FEASIBILITY.md`: an un-bracketed dative donor like
/// bare `N` must still mean NH3, not NH2). That changed N's invariant in
/// `n_to_fe` (N is the donor there) but not in `fe_to_n` (N is the
/// *acceptor* there, untouched by the donor-side-only fix) -- so `a`'s
/// canonical form flipped which atom is written first, while `b`'s did not.
///
/// Both expected strings are derived from `BondOrder::Dative`'s own
/// definition rather than from what the writer happens to produce.
#[test]
fn canonical_writer_orders_dative_endpoints_by_current_rank() {
    let a = canonical_smiles(&n_to_fe());
    let b = canonical_smiles(&fe_to_n());

    assert_eq!(a, "N->[Fe]", "donor (N) now ranks first ⇒ forward arrow");
    assert_eq!(b, "[Fe]->N", "donor (Fe) written first ⇒ forward arrow");

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

/// The arrow-flip-on-acceptor-first path (the actual code under test by the
/// name of the test above, before O/N/Fe's specific ranks shifted it away
/// from N/Fe) is still exercised here with an O donor instead: an oxygen
/// donor's implicit-H count is *also* correctly donor-exempted by the same
/// fix, but O still ranks below Fe, so Fe (the acceptor) is written first
/// and the writer must still emit a reversed `<-` arrow to keep the
/// donor/acceptor pair intact.
#[test]
fn canonical_writer_flips_the_arrow_when_the_acceptor_is_written_first() {
    let mut b = MoleculeBuilder::new();
    let o = b.add_atom(atom(Element::O));
    let fe = b.add_atom(atom(Element::FE));
    b.add_bond(o, fe, BondOrder::Dative).unwrap();
    let mol = b.build();

    let out = canonical_smiles(&mol);
    assert_eq!(out, "[Fe]<-O", "acceptor written first ⇒ reversed arrow");
    assert_eq!(
        dative_donor_acceptor(&parse(&out).unwrap()),
        (Element::O, Element::FE)
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
///
/// This molecule has two paths between N and Fe (the 3-single-bond chain,
/// and the direct dative bond), so which one the DFS picks as the ring-
/// closure "back edge" vs. a tree edge is incidental to canonical atom
/// ordering, not a fixed property of the graph. `write` (plain, non-
/// canonical, insertion-order-driven) always routes the dative bond through
/// the ring closure for this builder-constructed molecule, so it keeps the
/// stricter ring-closure-digit assertion. `canonical_smiles` legitimately
/// changed which edge becomes the ring closure once `initial_invariant`'s
/// explicit/implicit-H unification fix (issue #205) corrected this
/// molecule's Morgan ranks -- the dative bond is now a tree edge, written
/// inline right after Fe, rather than at the ring digit. Only the invariant
/// that actually matters -- donor/acceptor direction survives the round
/// trip, and the dative token itself is never dropped -- is asserted for
/// `canonical_smiles`, not the incidental tree-edge-vs-ring-closure choice.
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

        let plain = write(&mol);
        assert!(
            plain.contains("->1") && plain.contains("<-1"),
            "plain writer must still route this dative bond through the ring-closure \
             digit (insertion-order-driven, unaffected by canonical ranking): {plain}"
        );

        let canon = canonical_smiles(&mol);
        assert!(
            canon.contains("->") || canon.contains("<-"),
            "dative bond token lost entirely canonicalizing {canon}"
        );

        for out in [plain, canon] {
            let reparsed = parse(&out).unwrap_or_else(|e| panic!("re-parse of {out}: {e}"));
            assert_eq!(
                dative_donor_acceptor(&reparsed),
                expected,
                "donor/acceptor lost writing {out}"
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
/// `MoleculeBuilder` (rather than parsed) used to be able to shift once on
/// its first canonicalization, because the builder and the parser did not
/// populate every atom field identically -- specifically, `hydrogen_count`
/// (`None` from a builder-constructed atom vs. an explicit `Some(_)` after a
/// bracket atom like `[Fe]` is written and re-parsed). That was pre-existing
/// and had nothing to do with dative bonds -- the same ring with an ordinary
/// single bond in place of the dative one shifted in exactly the same way --
/// so it must never be "fixed" by reshaping the dative writer specifically.
///
/// Fixed at the root (issue #205): `initial_invariant`'s Morgan-rank seed
/// and `emit_atom`'s bracket-necessity check now both go through
/// `implicit_hcount`/`valence_inferred_hcount` instead of raw
/// `atom.hydrogen_count`, so an explicit H count that merely repeats what
/// valence inference would already give (e.g. bracket `[Fe]` with H=0,
/// implicit-organic-subset atoms aside since Fe isn't organic-subset) no
/// longer changes ranking or bracket emission. Builder- and parser-
/// constructed molecules of the same graph now canonicalize identically on
/// the very first call -- this test now asserts convergence, not the shift.
#[test]
fn builder_first_canonicalization_matches_parser_no_shift() {
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
    assert_eq!(
        single_first, single_second,
        "builder- and parser-constructed molecules must canonicalize identically \
         (no more first-call shift) -- not dative-specific, same fix for both"
    );
    assert_eq!(dative_first, dative_second, "same invariant, same cause");
    // Stable under repeated round-trips too.
    assert_eq!(
        canonical_smiles(&parse(&single_second).unwrap()),
        single_second
    );
    assert_eq!(
        canonical_smiles(&parse(&dative_second).unwrap()),
        dative_second
    );
}
