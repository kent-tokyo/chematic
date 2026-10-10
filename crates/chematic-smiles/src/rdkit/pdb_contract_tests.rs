//! PDB and expected canonical graphs captured with RDKit 2026.03.1:
//! MolFromSequence -> MolToPDBBlock -> MolFromPDBBlock(proximityBonding=False).
use super::rdkit_mol_from_pdb_block;

#[test]
fn protein_and_nucleic_acid_pdb_graphs_match_rdkit_reference() {
    let rows: serde_json::Value = serde_json::from_str(include_str!(
        "../../../../validation/rdkit-2026.03.1-pdb-residue-contract.json"
    ))
    .unwrap();
    for row in rows.as_array().unwrap() {
        let pdb = rdkit_mol_from_pdb_block(row["pdb"].as_str().unwrap(), true, true, 0, false)
            .unwrap()
            .unwrap();
        assert_eq!(
            pdb.smiles,
            row["smiles"].as_str().unwrap(),
            "{}",
            row["sequence"]
        );
        assert_eq!(
            pdb.molecule.atom_count(),
            row["atoms"].as_u64().unwrap() as usize
        );
        assert_eq!(
            pdb.molecule.bond_count(),
            row["bonds"].as_u64().unwrap() as usize
        );
        assert_eq!(pdb.coords.len(), pdb.molecule.atom_count());
    }
}
