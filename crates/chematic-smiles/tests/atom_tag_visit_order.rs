//! Regression coverage for molecule atom tags against SMILES visit-order helpers.
//!
//! Tags must not change written / canonical SMILES. Visit order from
//! `write_with_atom_order` / `canonical_smiles_with_atom_order` remaps tags
//! across write/parse.

use std::collections::{BTreeMap, HashSet};
use std::num::NonZeroU16;

use chematic_core::{AtomIdx, Molecule};
use chematic_smiles::{
    canonical_smiles, canonical_smiles_with_atom_order, parse, write, write_with_atom_order,
};

const CORPUS: &[&str] = &[
    "C",
    "O",
    "CCO",
    "CC(=O)O",
    "c1ccccc1",
    "c1ccncc1",
    "C1CCCCC1",
    "C/C=C/C",
    "C/C=C\\C",
    "[NH4+].[Cl-]",
    "CC(=O)Oc1ccccc1C(=O)O",      // aspirin
    "Cn1cnc2c1c(=O)n(c(=O)n2C)C", // caffeine
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O", // ibuprofen
    "c1ccc2ccccc2c1",             // naphthalene
    "NCC(=O)O",
];

fn stamp_tags(mol: &mut Molecule) {
    for i in 0..mol.atom_count() {
        // 1-based NonZero ids; stable for a given molecule size.
        mol.set_tag(AtomIdx(i as u32), Some((i as u16) + 1));
    }
}

fn tag_element_map(mol: &Molecule) -> BTreeMap<u16, u8> {
    (0..mol.atom_count())
        .filter_map(|i| {
            let atom = mol.atom(AtomIdx(i as u32));
            mol.atom_tag(AtomIdx(i as u32))
                .map(|t| (t.get(), atom.element.atomic_number()))
        })
        .collect()
}

fn assert_permutation(order: &[AtomIdx], n: usize) {
    assert_eq!(order.len(), n);
    let mut seen: Vec<_> = order.iter().map(|a| a.0).collect();
    seen.sort_unstable();
    assert_eq!(
        seen,
        (0..n as u32).collect::<Vec<_>>(),
        "visit order must be a permutation of 0..{n}"
    );
}

fn remap_via_order(mol: &Molecule, order: &[AtomIdx], fresh: &mut Molecule) {
    assert_eq!(fresh.atom_count(), order.len());
    for (new_i, &old_idx) in order.iter().enumerate() {
        let tag = mol.atom_tag(old_idx).map(NonZeroU16::get);
        fresh.set_tag(AtomIdx(new_i as u32), tag);
    }
}

#[test]
fn tags_do_not_change_write_or_canonical_across_corpus() {
    for smi in CORPUS {
        let plain = parse(smi).unwrap_or_else(|e| panic!("parse {smi}: {e}"));
        let mut tagged = plain.clone();
        stamp_tags(&mut tagged);

        assert_eq!(
            write(&tagged),
            write(&plain),
            "write changed by tags for {smi}"
        );
        assert_eq!(
            canonical_smiles(&tagged),
            canonical_smiles(&plain),
            "canonical_smiles changed by tags for {smi}"
        );

        let out = write(&tagged);
        assert!(
            !out.chars().any(|c| c == ':'),
            "tag leaked into SMILES for {smi}: {out}"
        );
    }
}

#[test]
fn tagged_molecule_atom_order_apis_match_untagged() {
    for smi in CORPUS {
        let mut mol = parse(smi).unwrap_or_else(|e| panic!("parse {smi}: {e}"));
        let (plain_w, plain_wo) = write_with_atom_order(&mol);
        let (plain_c, plain_co) = canonical_smiles_with_atom_order(&mol);
        stamp_tags(&mut mol);
        let (tag_w, tag_wo) = write_with_atom_order(&mol);
        let (tag_c, tag_co) = canonical_smiles_with_atom_order(&mol);
        assert_eq!(tag_w, plain_w, "write string changed by tags for {smi}");
        assert_eq!(tag_wo, plain_wo, "write order changed by tags for {smi}");
        assert_eq!(tag_c, plain_c, "canonical string changed by tags for {smi}");
        assert_eq!(
            tag_co, plain_co,
            "canonical order changed by tags for {smi}"
        );
        assert_permutation(&tag_wo, mol.atom_count());
        assert_permutation(&tag_co, mol.atom_count());
    }
}

#[test]
fn write_parse_remaps_tags_by_atom_order() {
    for smi in CORPUS {
        let mut mol = parse(smi).unwrap_or_else(|e| panic!("parse {smi}: {e}"));
        if mol.atom_count() == 0 {
            continue;
        }
        stamp_tags(&mut mol);
        let before = tag_element_map(&mol);

        let (written, order) = write_with_atom_order(&mol);
        let mut fresh = parse(&written).unwrap_or_else(|e| panic!("reparse {written}: {e}"));
        remap_via_order(&mol, &order, &mut fresh);

        assert_eq!(
            tag_element_map(&fresh),
            before,
            "tag→element map lost on write/parse for {smi}"
        );
        let tags: HashSet<_> = (0..fresh.atom_count())
            .filter_map(|i| fresh.atom_tag(AtomIdx(i as u32)).map(NonZeroU16::get))
            .collect();
        assert_eq!(tags.len(), fresh.atom_count());
    }
}

#[test]
fn canonical_write_parse_remaps_tags_by_atom_order() {
    for smi in CORPUS {
        let mut mol = parse(smi).unwrap_or_else(|e| panic!("parse {smi}: {e}"));
        if mol.atom_count() == 0 {
            continue;
        }
        stamp_tags(&mut mol);
        let before = tag_element_map(&mol);

        let (csmi, order) = canonical_smiles_with_atom_order(&mol);
        let mut fresh = parse(&csmi).unwrap_or_else(|e| panic!("reparse {csmi}: {e}"));
        remap_via_order(&mol, &order, &mut fresh);

        assert_eq!(
            tag_element_map(&fresh),
            before,
            "tag→element map lost on canonical write/parse for {smi}"
        );
    }
}

#[test]
fn reparse_does_not_invent_tags() {
    for smi in CORPUS {
        let mut tagged = parse(smi).unwrap();
        stamp_tags(&mut tagged);
        let out = write(&tagged);
        let re = parse(&out).unwrap();
        assert!(re.atoms().all(|(i, _)| re.atom_tag(i).is_none()));
    }
}

#[test]
fn tags_survive_fragments() {
    let mut mol = parse("CCO.O").unwrap();
    stamp_tags(&mut mol);
    let before = tag_element_map(&mol);
    let frags = mol.fragments_with_source_atoms();
    assert!(frags.len() >= 2);
    let mut recovered = BTreeMap::new();
    for (frag, source) in &frags {
        for (i, &src) in source.iter().enumerate() {
            let tag = frag.atom_tag(AtomIdx(i as u32)).map(NonZeroU16::get);
            assert_eq!(
                tag,
                mol.atom_tag(src).map(NonZeroU16::get),
                "fragment tag mismatch at source {src:?}"
            );
            if let Some(t) = tag {
                recovered.insert(t, frag.atom(AtomIdx(i as u32)).element.atomic_number());
            }
        }
    }
    assert_eq!(recovered, before);
}

#[test]
fn tags_do_not_change_atom_equality() {
    let mut a = parse("C").unwrap();
    let b = parse("C").unwrap();
    assert_eq!(a.atom(AtomIdx(0)), b.atom(AtomIdx(0)));
    a.set_tag(AtomIdx(0), Some(7));
    assert_eq!(a.atom(AtomIdx(0)), b.atom(AtomIdx(0)));
    a.set_tag(AtomIdx(0), None);
    assert_eq!(a.atom(AtomIdx(0)), b.atom(AtomIdx(0)));
}

#[test]
fn write_does_not_emit_atom_tag() {
    let mut mol = parse("C").unwrap();
    mol.set_tag(AtomIdx(0), Some(42));
    assert_eq!(write(&mol), "C");
    assert_eq!(mol.atom_tag(AtomIdx(0)).map(NonZeroU16::get), Some(42));
}
