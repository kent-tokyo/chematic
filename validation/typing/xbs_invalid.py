"""Four intentional type errors; do not execute."""

from chematic import XbsSpecies, XbsBondRule, write_xbs

XbsSpecies("C", "large", 0.5)
XbsBondRule("C", "C", 0, 2, 0.1, "gray")
write_xbs(["C"], [(0, 0, 0)], species=["C"])
write_xbs(["C"], ["xyz"], species=[XbsSpecies("C", 1, 0.5)])
