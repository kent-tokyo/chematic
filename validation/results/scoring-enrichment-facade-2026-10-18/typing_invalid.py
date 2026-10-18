"""Four deliberately invalid public argument types."""

from chematic.rdkit_compat import Scoring

Scoring.CalcEnrichment(1, 1, [0.5])
Scoring.CalcEnrichment([1], 0, [0.5])
Scoring.CalcEnrichment([(1.0, 1)], 1.5, [0.5])
Scoring.CalcEnrichment([(1.0, 1)], 1, ["0.5"])
