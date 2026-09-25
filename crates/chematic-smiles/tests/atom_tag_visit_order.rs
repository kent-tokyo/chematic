//! Regression coverage for `Atom.tag` + SMILES visit-order helpers.
//!
//! Guarantees:
//! - tags never change written / canonical SMILES
//! - `*_with_order` strings match the non-order APIs
//! - visit order is a permutation usable to remap tags across write/parse

use std::collections::{BTreeMap, HashSet};
use std::num::NonZeroU16;

use chematic_core::{AtomIdx, Molecule};
use chematic_smiles::{
    canonical_smiles, canonical_smiles_with_order, parse, write, write_with_order,
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
    "CC(=O)Oc1ccccc1C(=O)O", // aspirin
    "Cn1cnc2c1c(=O)n(c(=O)n2C)C", // caffeine
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O", // ibuprofen
    "c1ccc2ccccc2c1", // naphthalene
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
            atom.tag
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
        let tag = mol.atom(old_idx).tag.map(NonZeroU16::get);
        fresh.set_tag(AtomIdx(new_i as u32), tag);
    }
}

#[test]
fn empty_molecule_order_apis() {
    let mol = chematic_core::MoleculeBuilder::new().build();
    assert_eq!(write(&mol), "");
    assert_eq!(write_with_order(&mol), (String::new(), Vec::new()));
    assert_eq!(canonical_smiles(&mol), "");
    assert_eq!(
        canonical_smiles_with_order(&mol),
        (String::new(), Vec::new())
    );
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

        // Tags must not appear as atom-map `:n` in the string.
        let out = write(&tagged);
        assert!(
            !out.chars().any(|c| c == ':'),
            "tag leaked into SMILES for {smi}: {out}"
        );
    }
}

#[test]
fn with_order_strings_match_plain_apis_across_corpus() {
    for smi in CORPUS {
        let mol = parse(smi).unwrap_or_else(|e| panic!("parse {smi}: {e}"));
        let (w, w_order) = write_with_order(&mol);
        assert_eq!(w, write(&mol), "write_with_order string mismatch for {smi}");
        assert_permutation(&w_order, mol.atom_count());

        let (c, c_order) = canonical_smiles_with_order(&mol);
        assert_eq!(
            c,
            canonical_smiles(&mol),
            "canonical_smiles_with_order string mismatch for {smi}"
        );
        assert_permutation(&c_order, mol.atom_count());
    }
}

#[test]
fn write_parse_remaps_tags_by_visit_order() {
    for smi in CORPUS {
        let mut mol = parse(smi).unwrap_or_else(|e| panic!("parse {smi}: {e}"));
        if mol.atom_count() == 0 {
            continue;
        }
        stamp_tags(&mut mol);
        let before = tag_element_map(&mol);

        let (written, order) = write_with_order(&mol);
        let mut fresh = parse(&written).unwrap_or_else(|e| panic!("reparse {written}: {e}"));
        remap_via_order(&mol, &order, &mut fresh);

        assert_eq!(
            tag_element_map(&fresh),
            before,
            "tag→element map lost on write/parse for {smi}"
        );
        // Every stamped tag is unique and present.
        let tags: HashSet<_> = (0..fresh.atom_count())
            .filter_map(|i| fresh.atom(AtomIdx(i as u32)).tag.map(NonZeroU16::get))
            .collect();
        assert_eq!(tags.len(), fresh.atom_count());
    }
}

#[test]
fn canonical_write_parse_remaps_tags_by_visit_order() {
    for smi in CORPUS {
        let mut mol = parse(smi).unwrap_or_else(|e| panic!("parse {smi}: {e}"));
        if mol.atom_count() == 0 {
            continue;
        }
        stamp_tags(&mut mol);
        let before = tag_element_map(&mol);

        let (csmi, order) = canonical_smiles_with_order(&mol);
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
fn tagged_and_untagged_molecules_roundtrip_identically() {
    for smi in CORPUS {
        let plain = parse(smi).unwrap();
        let mut tagged = plain.clone();
        stamp_tags(&mut tagged);

        let plain_out = write(&plain);
        let tagged_out = write(&tagged);
        assert_eq!(plain_out, tagged_out);

        let plain_re = parse(&plain_out).unwrap();
        let tagged_re = parse(&tagged_out).unwrap();
        assert_eq!(plain_re.atom_count(), tagged_re.atom_count());
        assert_eq!(plain_re.bond_count(), tagged_re.bond_count());
        // Reparse never invents tags.
        assert!(tagged_re.atoms().all(|(_, a)| a.tag.is_none()));
    }
}

#[test]
fn disconnected_fragments_visit_order_covers_all_atoms() {
    let mol = parse("[Na+].[Cl-].CCO").unwrap();
    let (smi, order) = write_with_order(&mol);
    assert_eq!(smi, write(&mol));
    assert_permutation(&order, mol.atom_count());
    assert!(smi.contains('.'));

    let (csmi, c_order) = canonical_smiles_with_order(&mol);
    assert_eq!(csmi, canonical_smiles(&mol));
    assert_permutation(&c_order, mol.atom_count());
}

#[test]
fn clearing_tags_restores_default_atom_equality() {
    let mut a = parse("C").unwrap();
    let b = parse("C").unwrap();
    assert_eq!(a.atom(AtomIdx(0)), b.atom(AtomIdx(0)));
    a.set_tag(AtomIdx(0), Some(7));
    assert_ne!(a.atom(AtomIdx(0)), b.atom(AtomIdx(0)));
    a.set_tag(AtomIdx(0), None);
    assert_eq!(a.atom(AtomIdx(0)), b.atom(AtomIdx(0)));
}
