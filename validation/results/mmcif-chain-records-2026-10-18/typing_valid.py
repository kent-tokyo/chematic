"""Type-check a public consumer of immutable model/chain records."""

from pathlib import Path
import numpy as np
from chematic import MmcifChain, MmcifChainRead, read_mmcif_chains

text = Path("input.cif").read_text()
result: MmcifChainRead = read_mmcif_chains(text, max_chains=np.int64(32))
chain: MmcifChain = result.chains[0]
label: str = chain.label_asym_id
position: tuple[float, float, float] = chain.atoms[0].position_angstrom
component: str | None = chain.atoms[0].residue.label_comp_id
