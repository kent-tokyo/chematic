# RDKit ring-count version contract

`rdkit_num_rings` is a bounded compatibility implementation, not a claim that
every RDKit release returns the same ring set. RDKit 2026.03.6 and 2026.09.1
disagree on at least four pinned stress inputs because their default ring
backends differ.

The machine-readable contract is
[`validation/rdkit_ring_count_profiles_v1.json`](../validation/rdkit_ring_count_profiles_v1.json).
It records each input, both oracle values, the current chematic default, and the
classification `oracle_version_difference`. CI recomputes the chematic value
from the stored SMILES and rejects accidental drift.

The current default deliberately retains a bounded fallback for highly
repetitive cyclic oligosaccharides. Therefore it agrees with the 2026.03.6
values on these four boundary cases and must not be described as universal
2026.09.1 `CalcNumRings` parity. Reports must name the RDKit version and keep
version differences separate from parser rejection, typed abstention, notation
difference, and ordinary numeric mismatch.
