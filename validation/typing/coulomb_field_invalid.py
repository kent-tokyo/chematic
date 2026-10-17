"""Expected four public consumer errors; never execute this file."""

from chematic.rdkit_compat import rdMIF

rdMIF.Coulomb((x for x in [1.0]), [(0.0, 0.0, 0.0)])
rdMIF.Coulomb(["charge"], [(0.0, 0.0, 0.0)])
rdMIF.Coulomb([1.0], [0.0])
rdMIF.Coulomb([1.0], [(0.0, 0.0, 0.0)], probeCharge="1")
