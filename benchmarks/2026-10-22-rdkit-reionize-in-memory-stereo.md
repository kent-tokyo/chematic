# RDKit Reionize in-memory stereo follow-up (2026-10-22)

The source candidate based on `296c6dc7` retains a tetrahedral tag when
Reionize removes hydrogen from a tagged carbon anion. The RDKit-compatible
writer now distinguishes that edited, in-memory atom from a freshly parsed
three-coordinate `[C@-]` atom using the recorded hydrogen stereo ligand.
Reparsing the product SMILES still removes the tag, as RDKit does.

The corpus is `validation/cosmolkit_comparison/corpora/stress-v1.smi`
(SHA-256 `879166ab5481cac805ecd2d77f37d628383ec124085b90e2ee1659d55ce76277`).
The comparator is the saved RDKit 2026.09.1 output in the
[previous source record](../validation/results/rdkit-standardization-source-2026-10-21-reionize-catalog.json).
Both source versions were run on all 4,072 input rows with
`cargo run -p chematic-chem --example cleanup_rdkit_dump`.

| Operation | Before | Candidate | Changed rows |
| --- | ---: | ---: | ---: |
| Normalize | 4,069/4,069 | 4,069/4,069 | 0 |
| Cleanup | 4,053/4,069 | 4,053/4,069 | 0 |
| Reionize | 4,064/4,069 | 4,069/4,069 | 5 |

The only changed rows were 1432, 1878, 2914, 3153 and 4023. All five now
match the saved RDKit Reionize strings. Their molecular graphs and charges
already agreed; this change addresses serialization of the in-memory stereo
tag. The remaining 4,067 outputs were byte-identical to the source baseline.
The [machine record](2026-10-22-rdkit-reionize-in-memory-stereo.json) includes
the dump hashes and changed indices.

The Rust `chematic-chem` and `chematic-smiles` unit suites passed. A Python
binding regression test was added for `mol.reionize().rdkit_smiles`; it was not
run locally against a newly built extension. This is source-level agreement
on one development corpus, using saved comparator output. It does not establish
published-package parity, independent chemical correctness, or superiority
over RDKit. Reaction and R-group accuracy were not assessed in this change.
