from chematic import nucleic_acid_linear_svg

nucleic_acid_linear_svg({})
nucleic_acid_linear_svg(None)
nucleic_acid_linear_svg("{}", antiparallel="yes")
nucleic_acid_linear_svg("{}", base_pairs=[(1, 2)])
wrong_return: int = nucleic_acid_linear_svg("{}")
