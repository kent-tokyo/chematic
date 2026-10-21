# RDKit Reionize source follow-up (2026-10-21)

The source at `d4667e16` uses the ordered acid/base catalog already used by
Cleanup for the public Reionize operation. Isolated-atom charge correction
remains a separate first step. The baseline is the 2026-10-20 source record;
its source-file hashes match the pre-change `main` checkout.

| Operation | Before | After | Change |
| --- | ---: | ---: | ---: |
| Normalize | 4,069/4,069 | 4,069/4,069 | 0 |
| Cleanup | 4,053/4,069 | 4,053/4,069 | 0 |
| Reionize | 4,039/4,069 | 4,064/4,069 | +25, no regressed rows |

The corpus is `validation/cosmolkit_comparison/corpora/stress-v1.smi`
(SHA-256 `879166ab5481cac805ecd2d77f37d628383ec124085b90e2ee1659d55ce76277`).
The comparator is RDKit 2026.09.1. Three of the 4,072 input rows are excluded
from the shared-input exact-output denominator. The raw Reionize count is
4,067/4,072. The five remaining shared-input mismatches are rows 1432, 1878,
2914, 3153 and 4023. In each case RDKit's in-memory Reionize product writes a
tetrahedral tag on a three-coordinate carbon anion, while reparsing that
product SMILES removes the tag. The chematic result retains an in-memory
tetrahedral tag, but the RDKit-compatible writer sanitizes its graph before
writing. These five rows remain exact-output mismatches.

The [machine record](../validation/results/rdkit-standardization-source-2026-10-21-reionize-catalog.json)
contains source and corpus hashes, row indices, classification and the full
comparison. This is a source-only development corpus. It does not establish
published-package accuracy, independent chemical correctness, or overall
superiority to RDKit. Reaction and R-group algorithms were not changed in this
follow-up.
