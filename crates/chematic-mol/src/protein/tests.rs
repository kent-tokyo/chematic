use super::*;

const PDB: &str = concat!(
    "MODEL        1\n",
    "ATOM      1  N   ALA A   1      11.104  13.207   8.678  1.00 20.00           N  \n",
    "ATOM      2  CA  ALA A   1      12.560  13.300   8.500  1.00 20.00           C  \n",
    "HETATM    3 SE   MSE A   2      13.000  14.000   9.000  1.00 21.00          SE  \n",
    "HETATM    4  O   HOH A   3      10.000  10.000  10.000  1.00 22.00           O  \n",
    "ENDMDL\n",
    "MODEL        2\n",
    "ATOM      5  N  AGLY B   7A     15.000  16.000  17.000  0.50 23.00           N  \n",
    "ENDMDL\n",
);

const MMCIF: &str = r#"data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . ALA A 1 1 ? 0 0 0 1.0 10.0 1
HETATM 2 Se SE . MSE A 2 2 ? 1 0 0 1.0 11.0 1
HETATM 3 O O . HOH A 3 . ? 2 0 0 1.0 12.0 1
ATOM 4 C CA . DA B 1 1 ? 3 0 0 1.0 13.0 1
"#;

#[test]
fn pdb_projection_keeps_models_chains_modifications_and_altloc() {
    let projection = project_protein_from_pdb(PDB).unwrap();
    assert_eq!(projection.model_count(), 2);
    assert_eq!(projection.chain_count(), 2);
    assert_eq!(projection.residue_count(), 3);
    assert_eq!(projection.atom_count(), 4);
    assert_eq!(projection.chains[0].sequence(), "AM");
    assert_eq!(projection.chains[1].sequence(), "G");
    assert!(projection.chains[0].residues[1].modified);
    assert_eq!(
        projection.chains[0].residues[1].atoms[0].element,
        Element::SE
    );
    assert_eq!(projection.chains[1].residues[0].insertion_code, Some('A'));
    assert_eq!(projection.chains[1].residues[0].atoms[0].alt_loc, Some('A'));
}

#[test]
fn mmcif_projection_excludes_water_and_nucleic_acids() {
    let projection = project_protein_from_mmcif(MMCIF).unwrap();
    assert_eq!(projection.model_count(), 1);
    assert_eq!(projection.chain_count(), 1);
    assert_eq!(projection.residue_count(), 2);
    assert_eq!(projection.atom_count(), 2);
    assert_eq!(projection.chains[0].sequence(), "AM");
}

#[test]
fn ligand_only_input_returns_an_empty_projection() {
    let pdb = "HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00 20.00           C  \n";
    let projection = project_protein_from_pdb(pdb).unwrap();
    assert_eq!(projection, ProteinProjection::default());
}

#[test]
fn invalid_explicit_pdb_element_is_rejected_instead_of_becoming_carbon() {
    let pdb = "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00          ZZ  \n";
    assert!(matches!(
        project_protein_from_pdb(pdb),
        Err(ProteinProjectionError::UnknownPdbElement { serial: 1, .. })
    ));
}
