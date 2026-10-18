from chematic.rdkit_compat import rdInfoTheory as I
generator = I.BitCorrMatGenerator()
generator.SetBitList([1.5])
generator.SetBitList(3)
generator.CollectVotes([1, 3])
wrong: list[float] = generator.GetCorrMatrix()
