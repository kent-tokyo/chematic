from typing import Sequence
from chematic import nucleic_acid_linear_svg

pairs: Sequence[tuple[str, str]] = [("a", "b")]
svg: str = nucleic_acid_linear_svg("{}", antiparallel=True, base_pairs=pairs)
svg_bytes: str = nucleic_acid_linear_svg(b"{}")
