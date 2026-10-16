use chematic_core::AtomIdx;
use chematic_mol::{
    Mol2SectionOrder, parse_mol2_record, write_mol2, write_mol2_record,
    write_mol2_record_with_section_order, write_mol2_with_section_order,
};

fn markers(text: &str) -> Vec<&str> {
    text.lines()
        .filter(|line| line.starts_with("@<TRIPOS>"))
        .collect()
}

#[test]
fn charged_graph_writers_emit_atom_bond_then_unity_even_with_no_bonds() {
    for smiles in ["[O-2]", "[Na+]", "C[N+](C)(C)C", "[NH3+]CC[O-]"] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let output = write_mol2_with_section_order(&mol, &[], Mol2SectionOrder::AtomBondAdjacent);
        let default = write_mol2(&mol, &[]);
        assert!(default.find("@<TRIPOS>UNITY_ATOM_ATTR") < default.find("@<TRIPOS>BOND"));
        assert_eq!(
            default,
            write_mol2_with_section_order(&mol, &[], Mol2SectionOrder::default())
        );
        assert_eq!(
            &markers(&output)[..4],
            [
                "@<TRIPOS>MOLECULE",
                "@<TRIPOS>ATOM",
                "@<TRIPOS>BOND",
                "@<TRIPOS>UNITY_ATOM_ATTR"
            ]
        );
        let parsed = parse_mol2_record(&output).unwrap();
        assert_eq!(parsed.molecule.atom_count(), mol.atom_count());
        assert_eq!(parsed.molecule.bond_count(), mol.bond_count());
        for (idx, atom) in mol.atoms() {
            assert_eq!(parsed.molecule.atom(idx).charge, atom.charge);
        }
        assert!(
            parsed
                .atoms
                .iter()
                .all(|atom| atom.partial_charge == Some(0.0))
        );
    }
}

#[test]
fn editable_record_retains_attributes_ids_and_partial_charges() {
    let input = "@<TRIPOS>MOLECULE\noriginal\n1 0 0 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n17 nitrogen 1.25 -2.5 3.75 N.4 2 LIG -0.7 ACTIVE\n@<TRIPOS>UNITY_ATOM_ATTR\n17 2\ncharge 1\ncolor red\n@<TRIPOS>BOND\n@<TRIPOS>VENDOR\nopaque retained\n";
    let record = parse_mol2_record(input).unwrap();
    let output = write_mol2_record_with_section_order(&record, Mol2SectionOrder::AtomBondAdjacent);
    let default = write_mol2_record(&record);
    assert!(default.find("@<TRIPOS>UNITY_ATOM_ATTR") < default.find("@<TRIPOS>BOND"));
    assert_eq!(
        markers(&output),
        [
            "@<TRIPOS>MOLECULE",
            "@<TRIPOS>ATOM",
            "@<TRIPOS>BOND",
            "@<TRIPOS>UNITY_ATOM_ATTR",
            "@<TRIPOS>VENDOR"
        ]
    );
    let after = parse_mol2_record(&output).unwrap();
    assert_eq!(after.atoms, record.atoms);
    assert_eq!(after.coords, record.coords);
    assert_eq!(after.unity_atom_attributes, record.unity_atom_attributes);
    assert_eq!(after.opaque_sections, record.opaque_sections);
    assert_eq!(after.molecule.atom(AtomIdx(0)).charge, 1);
}

#[test]
fn neutral_graph_writer_bytes_are_unchanged() {
    let mol = chematic_smiles::parse("CO").unwrap();
    assert_eq!(
        write_mol2(&mol, &[]),
        include_str!(
            "../../../validation/openbabel_issues/mol2-section-order/neutral-expected.mol2"
        )
    );
}

#[test]
fn immutable_source_writer_does_not_reorder_original_text() {
    let input =
        include_str!("../../../validation/openbabel_issues/mol2-section-order/charged-source.mol2");
    for text in [input.to_string(), input.replace('\n', "\r\n")] {
        let record = chematic_mol::parse_mol2_source_record(&text).unwrap();
        assert_eq!(chematic_mol::write_mol2_source_record(&record), text);
    }
}
