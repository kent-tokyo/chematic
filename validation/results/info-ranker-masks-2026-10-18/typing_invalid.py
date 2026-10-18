from chematic.rdkit_compat import rdInfoTheory as I
r = I.InfoBitRanker(4, 2)
r.SetMaskBits([1.5])
r.SetMaskBits(["1"])
r.SetMaskBits(None)
r.SetMaskBits(1)
