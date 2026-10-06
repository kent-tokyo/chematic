//! Canonical SMILES of an explicit-H reaction product in the atom order the
//! reaction left it in (a BioTransformer secondary-carbon hydroxylation of
//! tri(dodecylthio)methane, 117 atoms, read from a MOL block that keeps that
//! order). The automorphism search behind canonical ranking used to take
//! vertices in index order and spent tens of milliseconds per product here
//! (175 s over the rule's 4,320 product sets); it now extends from the
//! mapped region. The answer must not depend on the atom order.

#[test]
fn product_order_and_smiles_order_give_one_canonical_smiles() {
    let block = include_str!("data/tri-dodecylthio-methanol-product-explicit-h.mol");
    let (mol, _) = chematic_mol::parse_mol(block).expect("MOL block");
    assert_eq!(mol.atom_count(), 117);
    let t = std::time::Instant::now();
    let smiles = chematic_smiles::canonical_smiles(&mol);
    let elapsed = t.elapsed();
    let reparsed = chematic_smiles::parse(&smiles).expect("round trip");
    assert_eq!(chematic_smiles::canonical_smiles(&reparsed), smiles);
    // About 16 ms in a debug build on Linux x86-64 (24 us optimised; the old
    // index-order search took 34 ms optimised).
    assert!(elapsed.as_millis() < 500, "{elapsed:?}");
}
