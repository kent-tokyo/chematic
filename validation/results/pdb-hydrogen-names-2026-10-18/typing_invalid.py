from chematic import name_pdb_hydrogens_checked
name_pdb_hydrogens_checked(b"PDB")
name_pdb_hydrogens_checked("", [1.2])
name_pdb_hydrogens_checked("", {1})
r = name_pdb_hydrogens_checked("")
r.pdb = "changed"
r.assignments[0].new_name = " H1 "
