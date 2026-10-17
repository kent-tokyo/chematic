"""Public consumer: copied explicit-charge Coulomb field signatures."""

from chematic.rdkit_compat import rdMIF

field = rdMIF.Coulomb([1.0, -1.0], [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)])
energy: float = field(2.0, 0.0, 0.0, 25.0)
count: int = field.site_count
soft = rdMIF.Coulomb((1, -1), ((0, 0, 0), (1, 0, 0)), softcoreParam=0.5)
soft_energy: float = soft(2, 0, 0, 25)
