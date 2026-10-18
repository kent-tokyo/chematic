from chematic import read_cube_datasets
read_cube_datasets(3)
read_cube_datasets("text", max_atoms=1.5)
read_cube_datasets("text").select("0")
wrong: int = read_cube_datasets("text").select(0).text
