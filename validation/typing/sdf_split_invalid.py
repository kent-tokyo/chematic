from chematic import plan_sdf_split

plan_sdf_split("input.sdf", 1.5)
plan_sdf_split("input.sdf", 2, prefix=b"out")
plan_sdf_split("input.sdf", 2).iter_part("0")
plan_sdf_split("input.sdf", 2).parts[0].filename = "changed.sdf"
