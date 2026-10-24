# V2000 property precedence

MOL/SDF reading now treats `M  CHG` and `M  ISO` as replacements for their
corresponding **whole atom-block channels**, rather than updates to only the
listed atoms. This follows the
[MDL CTfile property rules, pages 2-9 and 2-10](https://mdli.com/downloads/literature/ctfile.pdf).

For example, an old atom-line `+1` charge on N must become neutral when a
`M  CHG` record lists only O. Similarly, an atom-line N-15 mass difference must
become unspecified when `M  ISO` lists only a different atom. The previous
reader incorrectly retained those old values.

The implementation tracks property presence separately for charge and isotope,
resets each present channel once, then applies all records in that channel.
Multiple property lines accumulate; a later line does not erase an earlier
line's assignments. Absent channels retain their atom-block values. Parsing
validation still runs before any property update is accepted. The change is
shared by the diagnostic reader and the graph-only SDF fast reader; public
signatures are unchanged.

## Fixed packet, 2026-10-24

The original synthetic packet contains five controls and nine precedence
reproductions, with a spec-authored oracle for atomic number, formal charge,
isotope, atom mapping and single/double bonds. It exercises sparse assignments,
both charge signs, several property lines, independent channels and reversed
CHG/ISO ordering. It is not a physical-validity or holdout corpus.

The fixed baseline is `342d1d5f9ad09aa6ade4e3079a902cba856209e2`; the candidate
starts at `506a52197463fa621516696ad6a29a8fd66426a7`. The affected reader file is
identical between those commits. Separately retained baseline Python and Node
WASM artifacts were reused; candidate artifacts were rebuilt. Reports hash
the reader, packet, harness, runtime artifacts and external tools. Filesystem
paths describe this local run and need not remain available.

| Interpretation of source input | Controls | Precedence reproductions |
|---|---|---|
| Baseline Python, all six routes | 5/5 | 0/9 |
| Candidate Python, all six routes | 5/5 | 9/9 |
| Baseline WASM, all three routes | 5/5 | 0/9 |
| Candidate WASM, all three routes | 5/5 | 9/9 |
| Open Babel 3.2.1 | 5/5 | 9/9 |
| RDKit 2025.09.3 direct-input diagnostic | 5/5 | 4/9 |

Python routes are direct MOL, coordinate-returning MOL, common MOL-to-SMILES,
string SDF iteration, file `SDMolSupplier` fast reading and forward-stream SDF
reading. WASM routes are direct MOL, common MOL-to-SMILES and SDF-to-SMILES.
All listed routes must satisfy the oracle for a case to count.

Open Babel receives the identical input in
`obabel -imol -osmi -xn -xa`: omit names and retain atom classes, using its
documented [SMILES output options](https://openbabel.org/docs/FileFormats/SMILES_format.html).
The selected profile preserves every observed field on all 14 inputs. This
change closes a chematic gap; it does not establish superiority over Open Babel.
Each external process is limited to 10 seconds and 1 MiB of combined output.

## Source interpretation and output observation differ

RDKit 2025.09.3 reads **emitted SMILES** with sanitization disabled and
hydrogens not removed. The five observed fields must match the authored
oracle. To account for canonical SMILES reordering, the packet requires exactly
one atom per atomic number and compares atom labels and bond endpoints through
that unique-element correspondence. The observer explicitly rejects repeated
elements; it is not a general graph-isomorphism implementation. Atom-map values
must remain on their corresponding atoms.

Direct RDKit **MOL input parsing** is a separate diagnostic, with sanitization
disabled, hydrogens retained and strict parsing enabled. In five isotope
precedence cases it retains an atom-block isotope on an atom omitted from
`M  ISO`, unlike the authored MDL rule and measured Open Babel behavior.
That diagnostic is not used to generate the oracle or gate candidate success.
This also documents an intentional compatibility difference for MOL/SDF APIs
whose names resemble RDKit's. No general claim about either toolkit follows.

Native tests check indexed fields on LF/CRLF input through six reader paths.
Normal and fast file readers, Python file/forward suppliers and Node SDF
reading also process all cases followed by their reverse order. Property-like
lines inside SD field values after `M  END` remain opaque and never affect the
graph; property presence never leaks between records. Native tests separately
cover already-accepted empty property records; those legacy zero-entry inputs
are outside the scored standard packet. Six observer controls test reordered
atoms, changed labels/maps/bonds, malformed observations and stale isotope
values without requiring RDKit in CI.

## Boundaries

This change concerns `M  CHG` and `M  ISO` on the packet's ordinary heavy-atom
graphs. It does not add `M  RAD`/spin-state support, interpret aliases or query
properties, change isotope-output generation, or validate every CTfile field.
Repeated assignments to the same atom keep their prior behavior. Hydrogen,
aromaticity, stereo, geometry, source-byte retention and performance are outside
the comparison. The separate previous isotope-output work is not included here.

The [packet and reproduction commands](https://github.com/kent-tokyo/chematic/tree/main/validation/v2000_property_precedence)
link complete baseline/candidate observations. These are bounded source-build
results, not a published-release or broad-format superiority claim.
