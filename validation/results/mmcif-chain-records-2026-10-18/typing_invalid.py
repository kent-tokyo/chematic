"""Four deliberate invalid public input types; mypy must reject all four."""

from pathlib import Path
from chematic import read_mmcif_chains

read_mmcif_chains(b"data_x")
read_mmcif_chains(Path("input.cif"))
read_mmcif_chains("data_x", max_atoms=2.5)
read_mmcif_chains("data_x", max_chains="2")
