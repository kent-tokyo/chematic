//! Edge cases of the RDKit-ported PDB and Mol2 readers that follow directly
//! from RDKit's reader code (`PDBParser.cpp`, `Mol2FileParser.cpp`).
use chematic_smiles::{rdkit_mol_from_mol2_block, rdkit_mol_from_pdb_block};

const C1: &str = "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n";

fn pdb_atom_count(text: &str, flavor: u32) -> usize {
    rdkit_mol_from_pdb_block(text, false, false, flavor, false)
        .unwrap()
        .map_or(0, |m| m.molecule.atom_count())
}

#[test]
fn pdb_placeholder_coordinates_are_skipped_unless_flavor_bit_one() {
    let placeholder =
        "HETATM    2  C2  LIG A   1    9999.0009999.0009999.000  1.00  0.00           C\n";
    let text = format!("{C1}{placeholder}END\n");
    assert_eq!(pdb_atom_count(&text, 0), 1);
    assert_eq!(pdb_atom_count(&text, 1), 2);
}

#[test]
fn pdb_second_alternate_location_is_skipped_unless_flavor_bit_one() {
    let alt_b = "HETATM    2  C2 BLIG A   1       1.500   0.000   0.000  0.50  0.00           C\n";
    let text = format!("{C1}{alt_b}END\n");
    assert_eq!(pdb_atom_count(&text, 0), 1);
    assert_eq!(pdb_atom_count(&text, 1), 2);
}

fn mol2(extra: &str) -> String {
    format!(
        "@<TRIPOS>MOLECULE\nfixture\n2 1 0 0 0\nSMALL\nUSER_CHARGES\n\n\
         @<TRIPOS>ATOM\n\
         1 N1 0.000 0.000 0.000 N.4 1 RES 0\n\
         2 C1 1.500 0.000 0.000 C.3 1 RES 0\n\
         @<TRIPOS>BOND\n1 1 2 1\n{extra}"
    )
}

#[test]
fn mol2_unity_atom_attr_skips_non_charge_attributes() {
    let text = mol2("@<TRIPOS>UNITY_ATOM_ATTR\n1 1\ncharge 1\n2 1\nAtomExpr 0\n");
    let parsed = rdkit_mol_from_mol2_block(&text, false, false, false).unwrap();
    assert_eq!(parsed.molecule.atom_count(), 2);
}

#[test]
fn mol2_unity_atom_attr_charge_is_applied() {
    let text = mol2("@<TRIPOS>UNITY_ATOM_ATTR\n1 1\nAtomExpr 1\n@<TRIPOS>SUBSTRUCTURE\n");
    let parsed = rdkit_mol_from_mol2_block(&text, false, false, false).unwrap();
    let charges: Vec<i8> = parsed.molecule.atoms().map(|(_, a)| a.charge).collect();
    assert_eq!(charges, vec![1, 0]);
}

#[test]
fn mol2_unity_atom_attr_rejects_out_of_range_atom() {
    let text = mol2("@<TRIPOS>UNITY_ATOM_ATTR\n3 1\nAtomExpr 1\n");
    let error = rdkit_mol_from_mol2_block(&text, false, false, false)
        .err()
        .unwrap();
    assert!(
        error.to_string().contains("atom index out of range"),
        "{error}"
    );
}

#[test]
fn mol2_reads_only_the_first_molecule() {
    let first = mol2("");
    let text = format!("{first}{}", first.replace("fixture", "second"));
    let parsed = rdkit_mol_from_mol2_block(&text, false, false, false).unwrap();
    assert_eq!(parsed.molecule.atom_count(), 2);
}
