//! Accurate CIP on a ring centre with three single bonds and a lone pair
//! (bridgehead amine, cyclic phosphine, arsine, sulfonium, selenonium).
//!
//! RDKit's SMILES parser (`chiralAtomNeedsTagInversion`) inverts `@`/`@@` on
//! such a centre when exactly one ring-closure digit is written on it, so the
//! same stereoisomer has spellings that RDKit and an OpenSMILES reader read
//! as mirror images. chematic parses as OpenSMILES; the accurate engine
//! reports RDKit's label for the spelling parsed and marks it
//! spelling-dependent. Expected labels are RDKit 2026.03.6 `rdCIPLabeler` on
//! `MolFromSmiles` of each spelling: the spellings of one molecule were
//! written by RDKit (`MolToSmiles(doRandom=True)`), so they all read as one
//! stereoisomer there and carry one label. 520 such spellings (13
//! molecules x 40) agree; the cases below are a sample with zero, one and
//! two ring-closure digits on the centre.

use chematic_chem::{CipMode, assign_cip_with_mode, cip_label_depends_on_smiles_spelling};
use chematic_core::{AtomIdx, CipCode};
use chematic_smiles::parse;

fn accurate(smiles: &str, atom: u32) -> Option<CipCode> {
    let mol = parse(smiles).expect("parses");
    assign_cip_with_mode(&mol, CipMode::Accurate)
        .expect("accurate CIP")
        .assignments
        .iter()
        .find(|(i, _)| i.0 == atom)
        .map(|&(_, c)| c)
}

const CASES: &[(&str, u32, CipCode)] = &[
    // aziridine N
    ("C[N@]1C[C@H]1C", 1, CipCode::R),
    ("C1[C@@H](C)[N@]1C", 3, CipCode::R),
    ("[C@H]1(C[N@@]1C)C", 2, CipCode::R),
    ("C[N@@]1C[C@H]1C", 1, CipCode::S),
    ("[N@]1([C@@H](C1)C)C", 0, CipCode::S),
    ("[N@@]1(C[C@H]1C)C", 0, CipCode::S),
    ("c1(ccccc1)[C@@H]1[N@@](C1)CC", 7, CipCode::R),
    ("c1ccc(cc1)[C@H]1C[N@@]1CC", 8, CipCode::R),
    // phospholane P
    ("C[C@@H]1CC[P@](c2ccccc2)C1", 4, CipCode::S),
    ("[P@@]1(c2ccccc2)C[C@H](C)CC1", 0, CipCode::S),
    ("C1[P@@](c2ccccc2)CC[C@H]1C", 1, CipCode::S),
    ("C1C[C@H](C[P@]1c1ccccc1)C", 4, CipCode::R),
    // thiolanium S+
    ("C[C@H]1CC[S@+](C)C1", 4, CipCode::S),
    ("[C@@H]1(C)C[S@@+](C)CC1", 3, CipCode::S),
    ("[S@@+]1(CC)CC[C@@H](C1)C", 0, CipCode::R),
    // quinuclidine bridgehead N (two ring-closure digits on the centre)
    ("FC1C[N@@]2C[C@H](Cl)C1CC2", 3, CipCode::S),
    ("C12CC[N@](C[C@@H]2Cl)CC1F", 3, CipCode::S),
    ("[N@]12CC(C(CC1)[C@H](C2)Cl)F", 0, CipCode::S),
    // arsolane As, selenolanium Se+
    ("C1[As@@](CC[C@@H]1C)C", 1, CipCode::R),
    ("[As@]1(CC[C@@H](C1)C)C", 0, CipCode::R),
    ("[Se@@+]1(C[C@H](CC1)C)C", 0, CipCode::R),
    ("C1C[Se@@+](C)C[C@H]1C", 2, CipCode::R),
];

#[test]
fn ring_lone_pair_centre_gets_rdkits_label_for_the_spelling() {
    for &(smiles, atom, want) in CASES {
        assert_eq!(accurate(smiles, atom), Some(want), "{smiles} atom {atom}");
        let mol = parse(smiles).unwrap();
        assert!(
            cip_label_depends_on_smiles_spelling(&mol, AtomIdx(atom)),
            "{smiles} atom {atom} is spelling-dependent"
        );
    }
}

#[test]
fn other_centres_are_not_spelling_dependent() {
    for (smiles, atom) in [
        ("C[S@@](=O)CC", 1),         // acyclic sulfoxide
        ("C[C@H]1CC[S@@](=O)C1", 4), // ring sulfoxide: RDKit reads it as written
        ("C[C@H](N)O", 1),           // carbon
        ("C[P@](CC)c1ccccc1", 1),    // acyclic phosphine
    ] {
        let mol = parse(smiles).unwrap();
        assert!(
            !cip_label_depends_on_smiles_spelling(&mol, AtomIdx(atom)),
            "{smiles} atom {atom}"
        );
    }
}

#[test]
fn pseudoasymmetric_ring_lone_pair_centres_get_rdkits_r_s() {
    // RDKit 2026.03.6 rdCIPLabeler; rule 5 between the two enantiomorphic
    // ring arms, with RDKit's reading of the spelling (one ring-closure
    // digit on the aziridine N). 1,080 of 1,120 RDKit spellings of ten such
    // molecules agree; the other 40 are the quinuclidinone carbon below.
    for (smiles, atom, want) in [
        ("C[C@@H]1C[P@](C)C[C@H](C)C1", 3u32, CipCode::LowerS),
        ("C[C@@H]1C[P@@](C)C[C@H](C)C1", 3, CipCode::LowerR),
        ("C[C@@H]1[C@H](C)[N@]1C", 4, CipCode::LowerS),
        ("C[C@@H]1[C@H](C)[N@@]1CC", 4, CipCode::LowerR),
        ("C[C@@H]1C[S@+](C)C[C@H](C)C1", 3, CipCode::LowerS),
        ("C[C@@H]1CC[P@](C)CC[C@H](C)C1", 4, CipCode::LowerS),
    ] {
        assert_eq!(accurate(smiles, atom), Some(want), "{smiles} atom {atom}");
    }
}

#[test]
fn carbon_whose_rule5_needs_a_lone_pair_centre_gets_rdkits_label() {
    // 3-quinuclidinone-like (RDKit: N s, C s): the N's arms are told apart
    // through the carbon and the carbon's through the N's lone pair, which
    // rule 5 now takes as an auxiliary descriptor (lone pair lowest, RDKit's
    // reading of the spelling).
    let mol = parse("O=C1C[N@]2CC[C@H]1CC2").unwrap();
    let r = assign_cip_with_mode(&mol, CipMode::Accurate).unwrap();
    assert_eq!(r.get(AtomIdx(3)), Some(CipCode::LowerS));
    assert_eq!(r.get(AtomIdx(6)), Some(CipCode::LowerS));
}
