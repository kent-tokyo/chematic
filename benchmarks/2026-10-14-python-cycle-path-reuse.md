# Intermediate Python cycle fixed-path reuse (2026-10-14)

This retained intermediate source candidate caches fixed-tree branches and edge vectors. It predates the final root-pruning candidate and is not the working-tree source.

Baseline is main `cac480fb23cd8adeec6e5e9bed1f18677f87b392`. On unfiltered offset-0 rows (1,000 molecules), 20 rotating four-arm fresh-process blocks demonstrated a first-use ring-count improvement over baseline of **1.473x**, paired 95% CI **[1.430, 1.495]**, **20/20** faster blocks. No operation in this one-slice matrix had its entire baseline-ratio interval below 1. Only offset 0 was measured for this intermediate candidate; there is no offset-5000 result or toolkit-wide claim.

- [All raw blocks, output gates and comparisons](2026-10-14-python-cycle-path-reuse-slice0-paired20.json)
- [Frozen source manifest](2026-10-14-python-cycle-path-reuse-source.json)
- [Frozen candidate patch, including tests](2026-10-14-python-cycle-path-reuse-candidate.patch)
- [14,072-input byte-exact descriptor/Morgan parity hashes](2026-10-14-python-cycle-path-reuse-parity.json)

First use excludes parsing; the raw matrix separately records actual parse-inclusive and repeated-use timing. First-use RDKit ring info is eagerly prepared by its parser. The strict criterion is CI lower endpoint above 1 plus all 20 faster blocks. Unequal comparator outputs remain unranked. The final candidate additionally avoids roots that cannot close a cycle and allocates path rows only for reached 2-core vertices.
