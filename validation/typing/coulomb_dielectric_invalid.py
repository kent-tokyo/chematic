"""Four intentional public type errors; do not execute."""

from chematic.rdkit_compat import rdMIF

rdMIF.CoulombDielectric([1.0], [(0, 0, 0)], epsilon="water")
rdMIF.CoulombDielectric([1.0], (row for row in [(0, 0, 0)]))
rdMIF.CoulombDielectric(["1"], [(0, 0, 0)])
rdMIF.CoulombDielectric([1.0], [(0, 0, 0)])(0, 0, 0, "all")
