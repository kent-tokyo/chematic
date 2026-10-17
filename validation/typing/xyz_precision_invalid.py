import chematic

chematic.write_xyz_checked([1], [[0.0, 0.0, 0.0]])
chematic.write_xyz_checked(["O"], [[0.0, 0.0, 0.0]], notation="general")
chematic.write_xyz_checked(["O"], [[0.0, 0.0, 0.0]], decimal_places=1.5)
chematic.write_xyz_checked(["O"], [[0.0, 0.0, 0.0]], min_width="20")
