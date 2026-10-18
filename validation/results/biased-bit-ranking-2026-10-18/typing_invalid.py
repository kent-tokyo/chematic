"""Four deliberate errors, expected to be rejected by strict mypy."""

from chematic.rdkit_compat import rdInfoTheory as I

ranker = I.InfoBitRanker(3, 2, I.InfoType.BIASCHISQUARE)
ranker.SetBiasList(1)
ranker.SetBiasList([1.0])
ranker.SetBiasList(["0"])
ranker.AccumulateVotes([0, 1], 0)
