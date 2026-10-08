//! Compare `rdkit_canonical_smiles(add_hydrogens(parse(s)))` with RDKit's
//! `MolToSmiles(AddHs(MolFromSmiles(s)))` from a TSV `input<TAB>expected`.
//! (chematic-smiles cannot call chematic-chem's `add_hydrogens`, so the
//! example adds the hydrogens the same way: every implicit H as a new
//! atom, heavy atom by heavy atom.)

use std::io::BufRead;

use chematic_core::{Atom, AtomIdx, BondIdx, BondOrder, Element, MoleculeBuilder, implicit_hcount};

fn add_hs(mol: &chematic_core::Molecule) -> chematic_core::Molecule {
    let mut b = MoleculeBuilder::new();
    for i in 0..mol.atom_count() {
        let mut a = mol.atom(AtomIdx(i as u32)).clone();
        a.hydrogen_count = Some(0);
        b.add_atom(a);
    }
    for i in 0..mol.bond_count() {
        let bond = mol.bond(BondIdx(i as u32));
        let _ = b.add_bond(bond.atom1, bond.atom2, bond.order);
    }
    b.copy_stereo_from(mol);
    b.copy_bond_directions_from(mol);
    let mut sentinel = Vec::new();
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        for _ in 0..implicit_hcount(mol, idx) {
            let h = b.add_atom(Atom::new(Element::H));
            let _ = b.add_bond(idx, h, BondOrder::Single);
            sentinel.push((idx, h));
        }
    }
    for (idx, h) in sentinel {
        if let Some(order) = mol.stereo_neighbor_order(idx) {
            let o: Vec<u32> = order
                .iter()
                .map(|&v| {
                    if v == chematic_core::STEREO_H_SENTINEL {
                        h.0
                    } else {
                        v
                    }
                })
                .collect();
            b.set_stereo_neighbor_order(idx, o);
        }
    }
    b.build()
}

fn main() {
    let path = std::env::args().nth(1).expect("FILE.tsv");
    let max_shown: usize = std::env::args().nth(2).map_or(20, |s| s.parse().unwrap());
    let (mut ok, mut bad, mut shown) = (0, 0, 0);
    for line in std::io::BufReader::new(std::fs::File::open(path).unwrap()).lines() {
        let line = line.unwrap();
        let Some((input, expected)) = line.split_once('\t') else {
            continue;
        };
        let got = chematic_smiles::parse(input)
            .map_err(|e| e.to_string())
            .and_then(|m| {
                chematic_smiles::rdkit_canonical_smiles(&add_hs(&m)).map_err(|e| e.to_string())
            });
        if got.as_deref() == Ok(expected) {
            ok += 1;
        } else {
            bad += 1;
            if shown < max_shown {
                shown += 1;
                println!("MISMATCH {input}\n  rdkit {expected}\n  ours  {got:?}");
            }
        }
    }
    println!("match {ok} mismatch {bad}");
}
