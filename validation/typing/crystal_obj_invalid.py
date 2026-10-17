import chematic as ch

source = ch.PeriodicStructure.from_cif("example")
ch.write_crystal_obj_checked("example", [.5])
ch.write_crystal_obj_checked(source, [".5"])
ch.write_crystal_obj_checked(source, [.5], subdivisions="2")
result = ch.write_crystal_obj_checked(source, [.5])
wrong_text: int = result.text
result.face_count = 0
