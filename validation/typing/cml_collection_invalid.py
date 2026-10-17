from chematic import read_cml_collection, write_cml_collection

read_cml_collection({})
read_cml_collection(None)
write_cml_collection([1])
write_cml_collection([None])
wrong_return: int = write_cml_collection([])
