from chematic import parse_pymatgen_structures

parse_pymatgen_structures({})
parse_pymatgen_structures(["structure"])
parse_pymatgen_structures(None)
parse_pymatgen_structures(0)
parse_pymatgen_structures("{}")[0].to_periodic_structure(allow_metadata_loss="yes")

bad_geometry: str = parse_pymatgen_structures("{}")[0].to_periodic_structure()
