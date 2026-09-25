//! Issue #650: SMILES atom output order.
//!
//! For `write_with_atom_order` and `canonical_smiles_with_atom_order`:
//!
//! 1. the string equals `write` / `canonical_smiles`;
//! 2. `order` is a permutation of the atom indices;
//! 3. parsing the string back yields a molecule whose atom `k` corresponds to
//!    `mol.atom(order[k])` (element, charge, isotope, aromaticity, atom map)
//!    and whose bonds are exactly the images of `mol`'s bonds under that
//!    correspondence (plain bond order; E/Z direction markers are spelling).
//!
//! Set `CHEMATIC_ATOM_ORDER_CORPUS` to a colon-separated list of `.smi`
//! files to run the same check over whole corpora.

use chematic_core::{AtomIdx, BondOrder, Molecule};
use chematic_smiles::{
    canonical_smiles, canonical_smiles_with_atom_order, parse, write, write_with_atom_order,
};

const SAMPLES: &[&str] = &[
    "C",
    "CCO",
    "OCC",
    "c1ccccc1",
    "c1ccc2ccccc2c1",
    "c1ccc2[nH]ccc2c1",
    "CC(C)(C)c1ccc(O)cc1",
    "N[C@@H](C)C(=O)O",
    "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
    "F/C=C/F",
    "F/C=C\\F",
    "C(/F)=C/Cl",
    "C/C=C/C=C/C",
    "[Na+].[Cl-]",
    "[K+].[O-]S(=O)(=O)[O-]",
    "[13CH3]C(=O)[O-]",
    "[2H]C([2H])([2H])O",
    "[CH3:1][OH:2]",
    "C1CC2CCC1CC2",
    "C12C3C4C1C5C2C3C45",
    "[H]C([H])([H])[H]",
    "*C(=O)O",
    "C[N+](C)(C)C.[Br-]",
    "CC(=O)Nc1ccc(O)cc1.CCO.O",
    "c1ccc2c(c1)-c1ccccc1-2",
    "N->[Co]",
    "O=C1NC(=O)c2ccccc21",
    "CC(C)C[C@H](NC(=O)[C@@H](Cc1ccccc1)NC(=O)c1cnccn1)B(O)O",
    "C1=CC=CC=C1",
    "c1cc/c(=C/C)cc1",
];

fn plain(order: BondOrder) -> BondOrder {
    match order {
        BondOrder::Up | BondOrder::Down => BondOrder::Single,
        other => other,
    }
}

fn check(mol: &Molecule, smiles: &str, order: &[AtomIdx], label: &str, input: &str) {
    let n = mol.atom_count();
    assert_eq!(order.len(), n, "{label} {input}: order length");
    let mut inv = vec![usize::MAX; n];
    for (k, a) in order.iter().enumerate() {
        let a = a.0 as usize;
        assert!(
            a < n && inv[a] == usize::MAX,
            "{label} {input}: not a permutation"
        );
        inv[a] = k;
    }
    let reparsed =
        parse(smiles).unwrap_or_else(|e| panic!("{label} {input}: reparse {smiles}: {e}"));
    assert_eq!(
        reparsed.atom_count(),
        n,
        "{label} {input} -> {smiles}: atom count"
    );
    assert_eq!(
        reparsed.bond_count(),
        mol.bond_count(),
        "{label} {input} -> {smiles}: bond count"
    );
    for (k, &a) in order.iter().enumerate() {
        let (x, y) = (mol.atom(a), reparsed.atom(AtomIdx(k as u32)));
        assert_eq!(
            (
                x.element, x.charge, x.isotope, x.aromatic, x.atom_map, x.wildcard
            ),
            (
                y.element, y.charge, y.isotope, y.aromatic, y.atom_map, y.wildcard
            ),
            "{label} {input} -> {smiles}: atom {k} (source {a:?})"
        );
    }
    for (_, bond) in mol.bonds() {
        let (a, b) = (inv[bond.atom1.0 as usize], inv[bond.atom2.0 as usize]);
        let Some((_, image)) = reparsed.bond_between(AtomIdx(a as u32), AtomIdx(b as u32)) else {
            panic!("{label} {input} -> {smiles}: missing bond {a}-{b}");
        };
        if !matches!(bond.order, BondOrder::Dative) {
            assert_eq!(
                plain(image.order),
                plain(bond.order),
                "{label} {input} -> {smiles}: bond {a}-{b} order"
            );
        }
    }
}

fn check_input(input: &str) {
    let Ok(mol) = parse(input) else { return };
    let (s, order) = write_with_atom_order(&mol);
    assert_eq!(s, write(&mol), "write {input}");
    check(&mol, &s, &order, "write", input);
    let (c, order) = canonical_smiles_with_atom_order(&mol);
    assert_eq!(c, canonical_smiles(&mol), "canonical {input}");
    check(&mol, &c, &order, "canonical", input);
}

#[test]
fn atom_order_matches_reparsed_atoms() {
    for input in SAMPLES {
        check_input(input);
    }
}

#[test]
fn write_order_is_visit_order_not_canonical_rank() {
    // "OCC" is written in input order; the canonical string starts elsewhere.
    let mol = parse("OCC").unwrap();
    let (s, order) = write_with_atom_order(&mol);
    assert_eq!(s, "OCC");
    assert_eq!(order, vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)]);
    let (c, corder) = canonical_smiles_with_atom_order(&mol);
    assert_eq!(c, "C(C)O");
    assert_eq!(corder, vec![AtomIdx(1), AtomIdx(2), AtomIdx(0)]);
}

#[test]
fn empty_molecule_has_empty_order() {
    let mol = chematic_core::MoleculeBuilder::new().build();
    assert_eq!(write_with_atom_order(&mol), (String::new(), Vec::new()));
    assert_eq!(
        canonical_smiles_with_atom_order(&mol),
        (String::new(), Vec::new())
    );
}

#[test]
fn atom_order_over_corpus_from_env() {
    let Ok(paths) = std::env::var("CHEMATIC_ATOM_ORDER_CORPUS") else {
        return;
    };
    let mut rows = 0usize;
    for path in paths.split(':').filter(|p| !p.is_empty()) {
        let text = std::fs::read_to_string(path).expect("corpus file");
        for line in text.lines() {
            if let Some(smi) = line.split_whitespace().next() {
                check_input(smi);
                rows += 1;
            }
        }
    }
    eprintln!("atom order checked on {rows} corpus rows");
}
