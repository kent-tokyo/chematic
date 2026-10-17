import chematic


def geometry() -> str:
    positions = ((x, 0.0, 0.0) for x in [0.1, 1.2, -0.3])
    return chematic.write_xyz_checked(
        iter(["O", "H", "H"]),
        positions,
        notation="scientific",
        decimal_places=16,
        min_width=0,
    )
