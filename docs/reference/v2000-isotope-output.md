# Explicit isotopes in regenerated V2000 output

Regenerated MOL/SDF now carries every positive isotope label in `M  ISO`
records. An explicitly labelled most-common isotope, such as `[12CH4]`, stays
distinct from an unspecified natural-abundance carbon atom. Previously its
atom-line mass difference was zero and the label disappeared.

The legacy mass-difference field now contains only -3 through +4. Larger
differences go through `M  ISO`, avoiding overflow into neighboring charge or
atom-map columns. All labelled atoms appear in the property block, including
those with a nonzero legacy difference, because `M  ISO` supersedes atom-block
isotope information. Records split after eight pairs.
These mechanisms follow the
[CTfile atom and property records](https://mdli.com/downloads/literature/ctfile.pdf).

Both the 2D and conformer-based native V2000 writers use this property writer,
and their SDF wrappers inherit it. The reader also recognizes adjacent
four-column isotope pairs: `M  ISO  1   19999` assigns mass 9999 to atom 1.
Historical whitespace-separated property input remains accepted.

Python's `Mol.to_sdf_with_charges` now delegates its molecular block to the
native conformer writer. Its former independent formatter dropped isotope,
formal charge and atom mapping. It retains the supplied Cartesian values,
the existing four-decimal `PARTIAL_CHARGES` data field and record delimiter,
while the header now correctly identifies the coordinates as 3D.
Public signatures are unchanged.

## Same-input observation, 2026-10-24

The fixed baseline is `342d1d5f9ad09aa6ade4e3079a902cba856209e2`; the candidate
branch starts at `90014be569aa3f2e2fc8f490b3f723afc4c7d32b`. The four existing
production files involved in this change are identical between those two
commits. Previously built, separately retained baseline Python and Node WASM
artifacts were reused; the candidate artifacts were rebuilt. Reports hash
the production files, packet, harness, runtime artifacts and external readers.
Measurement paths describe a local run and need not remain available.

The original synthetic packet contains **15 ordinary isotope cases** and
**4 separate hypothetical isotope-storage controls**. Ordinary cases include
unspecified atoms, explicit most-common masses, C-13, heavy water, U-235,
mixed labelled/unlabelled atoms, a charged salt, mapping and nine labels in
one record. The storage controls use carbon mass numbers 1, 112, 999 and
9999; their preservation does not imply those are physically valid isotopes.

| Route, ordinary cases | Baseline | Candidate |
|---|---|---|
| Python direct/coordinate/common MOL/common SDF, each | 5/15 | 15/15 |
| Python partial-charge SDF | 1/15 | 15/15 |
| WASM direct/common MOL/common SDF/array SDF, each | 5/15 | 15/15 |
| All listed Python routes jointly | 1/15 | 15/15 |
| All listed WASM routes jointly | 5/15 | 15/15 |

All four storage controls change from 0/4 to 4/4 on every listed Python and
WASM route. Native integration tests additionally exercise both 2D and 3D
writers, reusable buffers, report/checked wrappers and SDF record/multi-record
writers. The CLI test checks its `mol` output and existing `sdf` alias: that
alias currently emits a MOL block, so it is not credited as SDF framing.

RDKit **2025.09.3** reads every output with `sanitize=False`, `removeHs=False`
and `strictParsing=True`. Its indexed atomic numbers, formal charges, isotope
mass numbers, atom maps and sorted bond endpoints/orders must match the
authored oracle. A separate narrow layout gate checks the legacy mass-difference
range and fixed-width `M  ISO` assignments. Neither reader acceptance nor that
layout gate establishes full CTfile or chemical validity.

Open Babel **3.2.1** receives the identical SMILES in two profiles:
`obabel -ismi -omol` and the same command with `-xa`, its documented
[MOL atom-class output option](https://openbabel.org/docs/FileFormats/MDL_MOL_format.html).

| Open Babel profile | Ordinary | Storage controls |
|---|---|---|
| Default | 14/15 | 2/4 |
| Atom classes (`-xa`) | 15/15 | 4/4 |

Default mismatches lose mapping. With the documented option, Open Babel already
preserves all observed fields. This fix closes a chematic preservation gap;
it does **not** establish superiority over properly configured Open Babel.
Each external process has a 10-second and combined 1-MiB output bound.

## Limits and evidence

The checked packet uses ordinary fixed-width V2000 counts and maps and positive
isotope numbers fitting the four-column property field (1..9999). The legacy
infallible writers do not validate every possible core-model value: this change
does not establish lossless V2000 output for larger isotope numbers, oversized
counts/maps, arbitrary charges, query attributes or stereochemistry. Use V3000
for fields that exceed V2000's representation. Hydrogen counts, aromatic or
stereochemical equivalence, arbitrary geometry precision and speed are outside
this comparison. The Python data/geometry regression separately checks two
specified coordinate rows and their partial-charge values.

This workflow regenerates output; it is not a source-byte-retention result.
Six observer controls detect changed fields, unreadable output, shifted fields,
bad property counts/targets and duplicate isotope assignments. They run without
RDKit in CI. The native, Python, Node and CLI tests run through the existing
test discovery jobs.

The [fixture packet and reproduction commands](https://github.com/kent-tokyo/chematic/tree/main/validation/v2000_isotope_output)
link both compressed observations. These are source-build results on a small
authored packet, not a published-release, holdout-corpus or broad-format claim.
