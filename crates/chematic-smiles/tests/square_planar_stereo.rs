//! Square-planar coordination stereochemistry (`@SP1`/`@SP2`/`@SP3`), end to
//! end (P0 gap measured, not fixed, by `validation/platinum/FEASIBILITY.md`;
//! design in `docs/rfcs/square_planar_stereo_rfc.md`).
//!
//! Scope: SMILES parsing, canonicalization, and writing only -- MOL/SDF I/O,
//! `chematic-3d` embedding, and CIP/R-S-analog labeling for square-planar
//! centers are all explicitly out of scope for this PR (see the RFC).
//!
//! The trans-pairing rule per tag (SP1: (0,2)+(1,3), SP2: (0,1)+(2,3), SP3:
//! (0,3)+(1,2)) and the full permutation-remap table were derived
//! empirically against RDKit 2026.03.3, not hand-derived from the
//! OpenSMILES spec's unreproduced diagram -- see
//! `scripts/square_planar_permutation_oracle.py` (144/144 cases, 0
//! mismatches, re-runnable independently of this test file).

use chematic_core::{
    Atom, BondOrder, Chirality, Element, MoleculeBuilder, SquarePlanarPermutation,
};
use chematic_smiles::{canonical_smiles, parse};

/// The killer-benchmark condition this whole platinum-complex project has
/// been chasing: cisplatin and transplatin must NOT collapse to the same
/// canonical identity once their declared stereochemistry is given.
#[test]
fn cisplatin_and_transplatin_have_distinct_canonical_identity() {
    // Same donor/acceptor grammar as validation/platinum/pt_corpus.jsonl's
    // own (untagged) cisplatin/transplatin entries, with an @SP tag added.
    let cisplatin = parse("N->[Pt@SP1](<-N)(Cl)Cl").expect("cisplatin parses");
    let transplatin = parse("N->[Pt@SP2](<-N)(Cl)Cl").expect("transplatin parses");

    let cis_smi = canonical_smiles(&cisplatin);
    let trans_smi = canonical_smiles(&transplatin);
    assert_ne!(
        cis_smi, trans_smi,
        "cisplatin and transplatin must have distinct canonical identity once @SP-tagged"
    );

    // Same physical molecules re-parsed from their own canonical form must
    // still disagree with each other (canonicalization didn't accidentally
    // erase the distinction it just proved existed).
    assert_ne!(
        canonical_smiles(&parse(&cis_smi).unwrap()),
        canonical_smiles(&parse(&trans_smi).unwrap()),
    );
}

/// `(tag name, enum value, trans-pair-of-positions template)` for all 3
/// square-planar tags -- shared fixture data for both
/// `oracle_verified_permutation_table_matches_end_to_end` (end-to-end
/// through the parser/writer) and `remap_square_planar_tag_matches_oracle_table_directly`
/// (direct unit-level call), so the 144-case sweep isn't hand-transcribed
/// twice.
type TagEntry<'a> = (&'a str, SquarePlanarPermutation, [(u8, u8); 2]);
const TAGS: [TagEntry<'static>; 3] = [
    ("SP1", SquarePlanarPermutation::SP1, [(0u8, 2u8), (1, 3)]),
    ("SP2", SquarePlanarPermutation::SP2, [(0, 1), (2, 3)]),
    ("SP3", SquarePlanarPermutation::SP3, [(0, 3), (1, 2)]),
];

fn build(order: [usize; 4], ligands: &[&str; 4], tag: &str) -> String {
    format!(
        "{}[Pt@{tag}]({})({}){}",
        ligands[order[0]], ligands[order[1]], ligands[order[2]], ligands[order[3]]
    )
}

/// Predict which tag describes the same physical arrangement as
/// `(tag_pairs, order)` when re-expressed against the reference order
/// `[0,1,2,3]`. `order[i]` names which ligand id sits at original SMILES
/// slot `i`. The reference target order IS the identity `[0,1,2,3]`, so
/// "the position within the reference order of ligand id v" is simply `v`
/// itself -- no lookup needed (unlike the general case `remap_square_planar_tag`
/// handles, which looks a real `canonical` sequence up).
fn predict(tag_pairs: [(u8, u8); 2], order: [usize; 4]) -> &'static str {
    let mut new_pairs: Vec<(u8, u8)> = tag_pairs
        .iter()
        .map(|&(i, j)| {
            let (a, b) = (order[i as usize] as u8, order[j as usize] as u8);
            (a.min(b), a.max(b))
        })
        .collect();
    new_pairs.sort_unstable();
    for &(name, _, pairs) in &TAGS {
        let mut t = pairs.to_vec();
        t.sort_unstable();
        if t == new_pairs {
            return name;
        }
    }
    unreachable!("one of the 3 tags must match")
}

/// Every `@SP1`/`@SP2`/`@SP3` RDKit-oracle-verified permutation case
/// (`scripts/square_planar_permutation_oracle.py`, 144/144, 0 mismatches),
/// re-verified here end to end through chematic's own parser + canonical
/// writer: a molecule built with neighbors in a permuted order and tag
/// `tag` must canonicalize identically to the same molecule built with
/// neighbors in reference order 0,1,2,3 and the *predicted* new tag.
#[test]
fn oracle_verified_permutation_table_matches_end_to_end() {
    // C[Pt@tag](F)(Cl)[H] -- same shape as the oracle script's
    // "simple_4_distinct_ligands", ligand order [C, F, Cl, H].
    let ligands = ["C", "F", "Cl", "[H]"];

    let mut checked = 0;
    for order in permutations_of_4() {
        for &(tag_name, _, tag_pairs) in &TAGS {
            let smi = build(order, &ligands, tag_name);
            let mol = parse(&smi).unwrap_or_else(|e| panic!("{smi} failed to parse: {e}"));
            let canon = canonical_smiles(&mol);

            let predicted = predict(tag_pairs, order);
            let reference_smi = build([0, 1, 2, 3], &ligands, predicted);
            let reference_canon = canonical_smiles(&parse(&reference_smi).unwrap());

            assert_eq!(
                canon, reference_canon,
                "order={order:?} tag={tag_name} ({smi}) should canonicalize the same as \
                 reference-order+{predicted} ({reference_smi})"
            );
            checked += 1;
        }
    }
    assert_eq!(checked, 24 * 3, "must check all 24 permutations x 3 tags");
}

/// Direct unit-level counterpart to the end-to-end test above: calls
/// `chematic_core::remap_square_planar_tag` (the generalized-stereo-geometry
/// replacement for the removed `remap_square_planar`) directly against the
/// same 24-permutations x 3-tags sweep and the same `predict` fixture data
/// (not hand-transcribed a second time), skipping the SMILES parser/writer
/// entirely -- proves the new module reproduces the oracle-verified
/// `trans_pairs()` semantics at the function level, not just "some writer
/// output looks unchanged."
#[test]
fn remap_square_planar_tag_matches_oracle_table_directly() {
    let reference: [u32; 4] = [0, 1, 2, 3];
    let mut checked = 0;
    for order in permutations_of_4() {
        let original: [u32; 4] = [
            order[0] as u32,
            order[1] as u32,
            order[2] as u32,
            order[3] as u32,
        ];
        for &(tag_name, tag, tag_pairs) in &TAGS {
            let predicted_name = predict(tag_pairs, order);
            let predicted_tag = TAGS
                .iter()
                .find(|&&(name, _, _)| name == predicted_name)
                .map(|&(_, t, _)| t)
                .unwrap();

            assert_eq!(
                chematic_core::remap_square_planar_tag(tag, original, reference),
                Some(predicted_tag),
                "order={order:?} tag={tag_name}"
            );
            checked += 1;
        }
    }
    assert_eq!(checked, 24 * 3, "must check all 24 permutations x 3 tags");
}

fn permutations_of_4() -> Vec<[usize; 4]> {
    let mut out = Vec::with_capacity(24);
    for a in 0..4 {
        for b in 0..4 {
            if b == a {
                continue;
            }
            for c in 0..4 {
                if c == a || c == b {
                    continue;
                }
                for d in 0..4 {
                    if d == a || d == b || d == c {
                        continue;
                    }
                    out.push([a, b, c, d]);
                }
            }
        }
    }
    out
}

/// Several differently-ordered, differently-tagged SMILES for the *same*
/// physical arrangement must reach the same canonical identity.
#[test]
fn reordered_and_retagged_smiles_for_same_arrangement_reach_same_identity() {
    // All 3 describe the same physical cis-[PtCl2(NH3)2] arrangement.
    let variants = [
        "N->[Pt@SP1](<-N)(Cl)Cl",
        "N->[Pt@SP3](Cl)(<-N)Cl", // swap N and Cl branch order; predicted SP3
        "[NH3]->[Pt@SP1]([Cl])([Cl])<-[NH3]",
    ];
    let canon: Vec<String> = variants
        .iter()
        .map(|s| canonical_smiles(&parse(s).unwrap()))
        .collect();
    for (i, c) in canon.iter().enumerate().skip(1) {
        assert_eq!(
            &canon[0], c,
            "variant {i} ({}) should match variant 0 ({})",
            variants[i], variants[0]
        );
    }
}

/// Different tags for the *same* SMILES-text neighbor order must NOT
/// collapse to the same canonical identity (no false merge).
#[test]
fn different_tags_same_order_stay_distinct() {
    let sp1 = canonical_smiles(&parse("C[Pt@SP1](F)(Cl)[H]").unwrap());
    let sp2 = canonical_smiles(&parse("C[Pt@SP2](F)(Cl)[H]").unwrap());
    let sp3 = canonical_smiles(&parse("C[Pt@SP3](F)(Cl)[H]").unwrap());
    assert_ne!(sp1, sp2);
    assert_ne!(sp2, sp3);
    assert_ne!(sp1, sp3);
}

/// `canonical(parse(canonical(parse(s)))) == canonical(parse(s))` -- needed
/// in addition to the round-trip checks above, since a writer/parser
/// convention mismatch on a 3-state tag produces a *different valid tag*,
/// not a parse error, so it would NOT be caught by round-trip identity
/// alone.
#[test]
fn canonicalization_is_idempotent_for_square_planar_fixtures() {
    for smi in [
        "N->[Pt@SP1](<-N)(Cl)Cl",
        "N->[Pt@SP2](<-N)(Cl)Cl",
        "C[Pt@SP1](F)(Cl)[H]",
        "C[Pt@SP3](F)(Cl)[H]",
    ] {
        let once = canonical_smiles(&parse(smi).unwrap());
        let twice = canonical_smiles(&parse(&once).unwrap());
        assert_eq!(once, twice, "canonicalization must be idempotent for {smi}");
    }
}

/// A malformed/duplicate-neighbor square-planar record, built directly via
/// `MoleculeBuilder` (bypassing the parser, which itself never produces
/// this) must drop to unspecified on write, never emit a
/// plausible-but-wrong tag.
#[test]
fn malformed_duplicate_neighbor_order_drops_to_unspecified() {
    let mut b = MoleculeBuilder::new();
    let pt = b.add_atom(Atom::bracket(
        Element::PT,
        None,
        Chirality::SquarePlanar(SquarePlanarPermutation::SP1),
        0,
        0,
        None,
    ));
    let c = b.add_atom(Atom::new(Element::C));
    let f = b.add_atom(Atom::new(Element::F));
    let cl = b.add_atom(Atom::new(Element::CL));
    let h = b.add_atom(Atom::new(Element::H));
    b.add_bond(pt, c, BondOrder::Single).unwrap();
    b.add_bond(pt, f, BondOrder::Single).unwrap();
    b.add_bond(pt, cl, BondOrder::Single).unwrap();
    b.add_bond(pt, h, BondOrder::Single).unwrap();
    // Duplicate id (c.0 twice) -- a data-integrity problem, not a real
    // 4-distinct-neighbor order.
    b.set_stereo_neighbor_order(pt, vec![c.0, c.0, cl.0, h.0]);
    let mol = b.build();

    let smi = canonical_smiles(&mol);
    let reparsed = parse(&smi).expect("still parses (as unspecified stereo)");
    assert_eq!(
        reparsed.atom(pt).chirality,
        Chirality::None,
        "malformed stereo_neighbor_order must drop to unspecified, not emit a guessed tag: {smi}"
    );
}

/// A plain, untagged 4-coordinate Pt SMILES must stay `Chirality::None` --
/// never auto-promoted to a declared square-planar identity just because
/// the atom happens to have 4 neighbors.
#[test]
fn untagged_four_coordinate_pt_is_never_auto_promoted() {
    let mol = parse("N[Pt](N)(Cl)Cl").expect("plain untagged Pt center parses");
    let pt = (0..mol.atom_count())
        .map(|i| chematic_core::AtomIdx(i as u32))
        .find(|&i| mol.atom(i).element == Element::PT)
        .expect("has a Pt atom");
    assert_eq!(mol.atom(pt).chirality, Chirality::None);

    let smi = canonical_smiles(&mol);
    assert!(
        !smi.contains("@SP"),
        "untagged input must never gain a square-planar tag on write: {smi}"
    );
}

/// A chelate/ring-closure-shaped fixture, matching
/// `validation/platinum/pt_corpus.jsonl`'s own carboplatin grammar
/// (`N->[Pt]1(<-N)OC(=O)...O1`) with an `@SP` tag added -- exercises the
/// ring-closure `stereo_neighbor_order` path end to end. Verified directly
/// against chematic's own real parser/writer, not a separately re-modeled
/// Python approximation of ring-closure encounter order (see
/// `scripts/square_planar_permutation_oracle.py`'s scope note for why that
/// case is checked here instead of there).
#[test]
fn chelate_ring_closure_shaped_fixture_round_trips() {
    let carboplatin_like = "N->[Pt@SP1]1(<-N)OC(=O)C2(CCC2)C(=O)O1";
    let mol = parse(carboplatin_like).expect("carboplatin-shaped SP1 fixture parses");
    let canon = canonical_smiles(&mol);
    let reparsed = parse(&canon).expect("its own canonical form re-parses");
    assert_eq!(
        canonical_smiles(&reparsed),
        canon,
        "chelate-shaped square-planar fixture must round-trip stably"
    );
    assert!(canon.contains("@SP"), "stereo must survive: {canon}");
}
