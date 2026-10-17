use chematic_core::AtomIdx;
use chematic_mol::{Mol2Error, parse_mol2_record, write_mol2_record};

fn input(ids: &[u32], bonds: &[(u32, u32)]) -> String {
    use std::fmt::Write;
    let mut text = format!(
        "@<TRIPOS>MOLECULE\nID order\n{} {} 0 0 0\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n",
        ids.len(),
        bonds.len()
    );
    for (i, id) in ids.iter().enumerate() {
        writeln!(text, "{id} C{i} 0 0 0 C.3 7 RES 0 status").unwrap();
    }
    text.push_str("@<TRIPOS>BOND\n");
    for (i, (a, b)) in bonds.iter().enumerate() {
        writeln!(text, "{} {a} {b} 1 bond-status", i + 1).unwrap();
    }
    text
}

#[test]
fn sparse_reverse_ids_keep_atom_order_metadata_and_bond_endpoints() {
    let ids: Vec<_> = (0..96).map(|i| u32::MAX - i * 17).collect();
    let text = input(&ids, &[(ids[80], ids[3]), (ids[3], ids[0])]);
    let record = parse_mol2_record(&text).unwrap();
    let bonds: Vec<_> = record
        .molecule
        .bonds()
        .map(|(_, b)| (b.atom1, b.atom2))
        .collect();
    assert_eq!(bonds, [(AtomIdx(80), AtomIdx(3)), (AtomIdx(3), AtomIdx(0))]);
    assert_eq!(
        record.atoms.iter().map(|a| a.atom_id).collect::<Vec<_>>(),
        ids
    );
    assert_eq!(record.atoms[80].status_bits, ["status"]);
    assert_eq!(record.bonds[0].status_bits, ["bond-status"]);
    let rewritten = write_mol2_record(&record);
    let roundtrip = parse_mol2_record(&rewritten).unwrap();
    assert_eq!(record.atoms, roundtrip.atoms);
    assert_eq!(record.bonds, roundtrip.bonds);
    assert_eq!(record.coords, roundtrip.coords);
}

#[test]
fn duplicate_and_zero_ids_keep_the_first_source_atom() {
    for n in [32, 33, 80] {
        let mut ids: Vec<_> = (0..n).collect();
        ids[n as usize - 1] = 0;
        ids[n as usize - 2] = 1;
        let record = parse_mol2_record(&input(&ids, &[(0, 1), (1, 0)])).unwrap();
        let (_, bond) = record.molecule.bonds().next().unwrap();
        assert_eq!((bond.atom1, bond.atom2), (AtomIdx(0), AtomIdx(1)));
        assert_eq!(record.duplicate_bonds.len(), 1);
        assert_eq!(record.duplicate_bonds[0].atom1_id, 1);
        assert_eq!(record.duplicate_bonds[0].atom2_id, 0);
    }
}

#[test]
fn missing_endpoint_keeps_original_error_line_and_precedence() {
    let ids: Vec<_> = (1..=64).rev().collect();
    let text = input(&ids, &[(900, 901)]);
    assert_eq!(
        parse_mol2_record(&text).err(),
        Some(Mol2Error::InvalidBondLine {
            line: text.lines().count(),
            detail: "atom_id 900 not found".into(),
        })
    );
}
