from chematic.rdkit_compat import rdInfoTheory

ranker = rdInfoTheory.InfoBitRanker(4, 2)
ranker.AccumulateVotes([0, 1], 0)
ranker.AccumulateVotes("0101", 1)
ranker.GetTopN("2")
rdInfoTheory.InfoBitRanker("4", 2)
