from chematic import split_pdb_models

split_pdb_models(b"END\n")
split_pdb_models("END\n", prefix=1)
split_pdb_models("END\n", max_models=1.5)
split_pdb_models("END\n").files[0].filename = "changed.pdb"
