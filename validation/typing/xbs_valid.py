"""Public typed scene output consumer."""

from chematic import XbsSpecies, XbsBondRule, write_xbs

style = XbsSpecies("C", 1.0, 0.5)
rule = XbsBondRule("C", "C", 0.0, 2.0, 0.1, 0.7)
text: str = write_xbs(["C"], [(0.0, 0.0, 0.0)], species=[style], bonds=[rule])
