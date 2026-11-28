# Selenium-containing development crystals: 3D stereo diagnostics

This development audit selects every selenium-containing molecule from the
existing COD development partition. That partition contains 1,345 records;
all parsed in this selection and 18 contained selenium. The original split,
source and development-input checksums are verified by the selector. No
holdout rows were scored. These structures have already been exposed during
prior development and are not an independent holdout. None contains the newly
supported neutral three-coordinate selenium stereo motif, so this is a
surrounding-case regression check, not a direct validation of that motif.

The script reads the 18 selected records' SDF coordinates, removes explicit
hydrogens, and maps each SMILES heavy-atom graph to the reference heavy-atom
graph using a chirality-aware RDKit substructure match. All 18 mapped; all
had at least two possible matches under a search capped at two. The script
records the first mapping, so this comparison is sensitive to mapping
ambiguity and the chirality-aware choice. Coordinates are the selected SDF
coordinates, not a new conformer-generation result. For each, the script
runs chematic's per-element diagnostics and RDKit 2026.09.2 3D CIP
reassignment on the same mapped coordinates and on their reflection.

Across 18 original and 18 reflected inputs, the verifier reported ten atom
stereo checks as satisfied in the originals and violated in the reflections;
one square-planar center was unevaluable in both. Ten molecules had no
reported declared-stereo elements. RDKit produced a comparable CIP assignment
for eight of the ten tetrahedral centers: all 16 original/reflected
comparisons matched chematic's satisfied/violated status. The remaining two
centers had no comparable RDKit CIP assignment and were not counted as
agreement. This checks only reported centers, not all possible stereochemistry,
physical geometry, or coordinate-generation accuracy.

The former and changed source wheels were also run on the identical 36
coordinate sets using the pre-existing boolean stereo gate. All 36 results
matched exactly: 27 true and nine false. A true value includes cases with no
reported stereo and must not be interpreted as complete stereo correctness.
The probe specifically selected no instances of the newly supported neutral
selenium motif, and thus does not measure its recall. The narrow
[neutral-selenium development comparison](2026-11-27-neutral-selenium-stereo-development.md)
contains the separately authored parity fixtures.

The host used macOS 27.0.1 arm64, Python 3.13.6 and RDKit 2026.09.2, with
`RAYON_NUM_THREADS=1` for the selection run. The source-built, nonpublished
wheel SHA-256 was
`6fdfac5a97b2d61531bfbffad95ed8c3caae8ac3d8d965f3154d861ef1f38587`.
The prior wheel used for the boolean comparison had SHA-256
`65740f3f6d3d596c26001977d55f0078f98135f88c70b70979702ae75cd0ac1d`.
The SDF SHA-256 was `e39322d1eba60bc041fa380f043146e95247d4209d000374591582e6e2acd117`;
the split SHA-256 was `e13f3e18f11eeb75ccd9dd3e96a28426e47ef1d8bc56fa28c53b02ae9bfc8da8`;
the development input SHA-256 was
`e3c2cf4aa5c7c5740fd3d23c670e26535515b97f16a02ee2bf387661513726a1`.
The [selection and diagnostic script](../scripts/validation/selenium_crystal_stereo_development.py)
SHA-256 was `69ca530eb75a3ded1866317a4082b43ea714bf0718e23f55c80671f215da4185`;
the [gate comparison script](../scripts/validation/selenium_crystal_gate_compare.py)
SHA-256 was `63039b61fc72789b20b0efb5216e97154e3b9d66b1042205b761824a35054a5c`.
The resulting mapped-coordinate packet SHA-256 was
`44a64e20b7b9ebc7a26067af0378f9c03a350a660556c6b6e18ba370468884dc`.
Selected SMILES, mappings, coordinates, diagnostic results and gate records
are in `validation/results/2026-11-27-selenium-crystal-stereo/`, alongside the
complete development input used for selection. Regenerating the coordinate
packet also requires the hash-matched source SDF and split manifest.
