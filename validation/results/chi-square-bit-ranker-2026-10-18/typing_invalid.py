from chematic.rdkit_compat import rdInfoTheory
r = rdInfoTheory.InfoBitRanker(8, 2, rdInfoTheory.InfoType.CHISQUARE)
rdInfoTheory.InfoBitRanker(8, 2, "CHISQUARE")
r.AccumulateVotes([1, 2], 0)
r.AccumulateVotes(None, 0)
r.GetTopN(1.0)
