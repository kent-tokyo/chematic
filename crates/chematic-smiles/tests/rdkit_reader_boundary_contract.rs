use chematic_smiles::{
    RdkitSmilesError, rdkit_mol_from_mol2_block, rdkit_mol_from_pdb_block, rdkit_mol_from_xyz_block,
};

fn assert_reader_error(result: Result<impl Sized, RdkitSmilesError>, needle: &str) {
    let error = result
        .err()
        .unwrap_or_else(|| panic!("fixture must be rejected ({needle})"));
    assert!(
        error.to_string().contains(needle),
        "expected {needle:?}, got {error}"
    );
}

#[test]
fn xyz_reader_covers_numeric_eof_and_count_boundaries() {
    let (molecule, coordinates) =
        rdkit_mol_from_xyz_block("2\nno final newline\nC +.5 0 -1\nCL 1,5 2 3")
            .expect("valid XYZ without a final newline");
    assert_eq!(molecule.atom_count(), 2);
    assert_eq!(coordinates, [[0.5, 0.0, -1.0], [1.0, 2.0, 3.0]]);
    let (empty, coordinates) =
        rdkit_mol_from_xyz_block("   \ncomment\n").expect("blank count follows RDKit's zero");
    assert_eq!(empty.atom_count(), 0);
    assert!(coordinates.is_empty());

    for (text, needle) in [
        ("bad\ncomment\n", "recognize"),
        ("100001\ncomment\n", "maximum"),
        ("2\ncomment", "EOF"),
        ("2\ncomment\nC 0 0 0\n", "Missing coordinates"),
        ("1\ncomment\nC x 0 0\n", "Cannot convert"),
        ("1\ncomment\nXx 0 0 0\n", "not found"),
        ("1\ncomment\nC 0 0 0\nextra\n", "More lines"),
    ] {
        assert_reader_error(rdkit_mol_from_xyz_block(text), needle);
    }
}

fn pdb_atom(serial: i32, name: &str, residue: &str, alt: char, x: &str, element: &str) -> String {
    format!(
        "ATOM  {serial:>5} {name:<4}{alt}{residue:>3} A   1    {x:>8}{:>8}{:>8}  1.00 20.00          {element:>2}  ",
        "0.000", "0.000"
    )
}

#[test]
fn pdb_reader_skips_rdkit_ignored_records_and_rejects_bad_fields() {
    let text = [
        pdb_atom(1, "C", "GLY", ' ', "0.000", "C"),
        pdb_atom(2, "C", "GLY", 'B', "1.000", "C"),
        pdb_atom(3, " Q", "GLY", ' ', "2.000", "C"),
        pdb_atom(4, "DU", "DUM", ' ', "3.000", "C"),
        "END".into(),
    ]
    .join("\n");
    let parsed = rdkit_mol_from_pdb_block(&text, false, false, 0, false)
        .expect("ignored records are not parse errors")
        .expect("one ordinary atom remains");
    assert_eq!(parsed.molecule.atom_count(), 2);
    assert!(parsed.smiles.is_empty());

    let wildcard = format!("{}\nEND\n", pdb_atom(1, "*", "LIG", ' ', "0.000", "*"));
    assert_reader_error(
        rdkit_mol_from_pdb_block(&wildcard, false, false, 1, false),
        "determine element",
    );

    let bad_element = format!("{}\nEND\n", pdb_atom(1, "XX", "LIG", ' ', "0.000", "Xx"));
    assert_reader_error(
        rdkit_mol_from_pdb_block(&bad_element, true, true, 1, false),
        "not found",
    );
    let bad_coordinate = format!("{}\nEND\n", pdb_atom(1, "C", "LIG", ' ', "oops", "C"));
    assert_reader_error(
        rdkit_mol_from_pdb_block(&bad_coordinate, true, true, 1, false),
        "coordinate",
    );
}

#[test]
fn mol2_reader_reports_truncated_query_and_index_failures() {
    let header = "@<TRIPOS>MOLECULE\nfixture\n1 0 0 0 0\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n";
    for (atom, needle) in [
        ("", "no info"),
        ("1", "premature"),
        ("1 C x 0 0 C.3 1 RES 0", "coordinates"),
        ("1 Q 0 0 0 ANY 1 RES 0", "query atom"),
        ("1 X 0 0 0 Xx 1 RES 0", "not found"),
    ] {
        assert_reader_error(
            rdkit_mol_from_mol2_block(&format!("{header}{atom}\n"), true, true, true),
            needle,
        );
    }

    let with_bond = "@<TRIPOS>MOLECULE\nfixture\n2 1 0 0 0\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3 1 RES 0\n2 C2 1 0 0 C.3 1 RES 0\n@<TRIPOS>BOND\n1 1 3 1\n";
    assert_reader_error(
        rdkit_mol_from_mol2_block(with_bond, true, true, true),
        "index mismatch",
    );
    assert_reader_error(
        rdkit_mol_from_mol2_block(&with_bond.replace("1 1 3 1", "1 1 2 du"), true, true, true),
        "unspecified",
    );
}
