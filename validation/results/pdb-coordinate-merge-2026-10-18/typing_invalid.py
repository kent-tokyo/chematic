from chematic import merge_pdb_coordinate_blocks
merge_pdb_coordinate_blocks([b"PDB"])
merge_pdb_coordinate_blocks(["PDB"], chain_maps=[{"A": 1}])
merge_pdb_coordinate_blocks(["PDB"], allow_metadata_loss="yes")
wrong: int = merge_pdb_coordinate_blocks(["PDB"]).pdb
