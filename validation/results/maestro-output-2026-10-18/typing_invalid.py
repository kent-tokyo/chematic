from chematic import from_smiles, write_mae_checked, v2000_to_mae_checked
m = from_smiles("[H][H]")
write_mae_checked("wrong", [])
write_mae_checked(m, [["x", 0., 0.]])
write_mae_checked(m, [], title=1)
v2000_to_mae_checked(b"wrong")
r = write_mae_checked(m, [])
r.text = "changed"
