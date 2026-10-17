"""Public typed dielectric field usage."""

from chematic.rdkit_compat import rdMIF

field = rdMIF.CoulombDielectric([1.0], [(0.0, 0.0, 0.0)], epsilon=80.0, xi=4.0)
energy: float = field(2.0, 0.0, 0.0, 25.0)
count: int = field.site_count
soft = rdMIF.CoulombDielectric((1, -1), ((0, 0, 0), (1, 0, 0)), softcoreParam=0.2)
