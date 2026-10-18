# PDB coordinate blocks without reclassification

Open [Open Babel #2181](https://github.com/openbabel/openbabel/issues/2181)
reports protein/ligand PDB joining that changes ATOM records to HETATM,
changes residue identity and omits the separating TER. A local Open Babel 3.2.1
`-j` control reproduces the record-kind, residue-name and TER changes on an
original three-atom metadata fixture. The Windows GUI from the report and its
attached structure were not used.

Python now provides `merge_pdb_coordinate_blocks`. It joins validated source
coordinate records without molecular graph projection or residue/element
perception. It reuses the native [coordinate selector](pdb-atom-selection.md)
to validate every source and the complete merged output. This adds the missing
assembly boundary; native selection and molecular PDB writers are unchanged.

```python
from chematic import merge_pdb_coordinate_blocks

# Original synthetic fragments, not a complete protein/physical model.
protein = """ATOM     20  CA AALA A   7I      0.000   0.000   0.000  0.25 10.00           C1+
ATOM     35  N   ALA A   7I      0.500  -1.250   0.125  0.50 10.50           N
END
"""
ligand = """HETATM   20  C1  LIG A   7      20.125   0.000   1.500  0.25 10.00           C
END
"""
result = merge_pdb_coordinate_blocks(
    [protein, ligand], chain_maps=[{}, {"A": "L"}]
)
assert [a.record_type for a in result.atoms] == ["ATOM", "ATOM", "HETATM"]
assert [a.new_serial for a in result.atoms] == [1, 2, 4]
assert result.inserted_ter_after_sources == (0, 1)
assert result.source_texts == (protein, ligand)
```

## Identity and references

Output follows source and atom order. Atom and existing TER serials are assigned
consecutive positive decimal values starting at one, including inserted TER
records. Every ATOM/HETATM keeps its record kind, raw four-column atom name,
alternate location, residue name/number/insertion code, numeric spelling,
occupancy, B factor and element/charge text. Only serial and explicitly mapped
chain fields change. A ` CA ` atom with element C remains carbon; names are
never used to assign elements.

ANISOU receives its associated atom's new serial and chain; tensor fields and
all other text remain unchanged. Existing TER records remain at their source
positions with updated serial/chain fields. If a source has unterminated atoms
after its last TER, a source-boundary TER is added after its last atom/ANISOU.
This also applies to HET-only sources: the separator does not assert polymer
classification. Internal chain boundaries without TER are not reconstructed.

CONECT source/target serials are remapped within their original source. Target
order and duplicate-target multiplicity are retained, with all CONECT records
moved after all coordinate records. References between different sources are
not invented. One final END is emitted. Output uses LF, omits blank lines and
normalizes CONECT spacing; complete original text remains available separately.
The field handling follows [wwPDB coordinates](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html)
and [CONECT](https://www.wwpdb.org/documentation/file-format-content/format33/sect10.html).

`PdbMergeResult` is immutable and contains:

- `pdb`: merged coordinate-block text; the API writes no files.
- `source_texts`: exact originals, including CRLF and excluded metadata.
- `atoms`: zero-based source/atom indices, old/new serial, old/new chain and
  unchanged ATOM/HETATM kind for every input atom.
- `losses`: source index, original one-based physical line and excluded record
  name, when metadata-loss consent was supplied.
- `inserted_ter_after_sources`: zero-based sources receiving a new final TER.

## Chain namespaces, metadata and bounds

With no mapping, chain IDs are retained. Chains belonging to different inputs
must remain disjoint, even when their residue numbers do not overlap. Supply
one `chain_maps` mapping per source to resolve collisions. Unmapped chains
retain their ID; unknown keys and noninjective mappings refuse. IDs are one
ASCII letter, digit or space; the blank chain is `" "`. Source chains that are
already internally ambiguous are not repaired.

Each source must contain atoms in one implicit model. Explicit MODEL/ENDMDL
wrappers, malformed coordinates/ANISOU/TER/CONECT, unknown/self references,
hybrid-36 and five-column serial overflow refuse. Whole-entry metadata such
as HEADER, SEQRES, LINK, SSBOND and CRYST1 refuses by default. Explicit
`allow_metadata_loss=True` permits exclusion with a loss sidecar; it never
overrides malformed supported records or the model/chain restrictions.

Hard caps are 32 sources, 8 MiB per source, 16 MiB aggregate input/output,
50,000 atoms, 250,000 physical input/output records and 1,024 bytes per line.
Source/mapping iteration consumes at most one excess item to detect caps.
Caller iterator/mapping errors propagate, without exposing partial results.
There is no timeout for arbitrary caller methods. Wrong types raise TypeError;
invalid content, counts, mappings or bounds raise ValueError.

The caller supplies coordinates already in a common angstrom frame. There is
no alignment, bond/hydrogen inference, residue-template validation, crystallographic
assembly, alternate-location selection, minimization or cross-source chemistry.
This is a coordinate block, not a complete deposition-ready PDB entry. Richer
identifiers and complete metadata should remain in mmCIF/source sidecars.

## Validation

231 Python tests pass (43 new plus existing selection/compatibility regressions).
Tests exercise source exceptions, metadata consent, reference multiplicity,
immutable results and actual 32-source/50000-atom/250000-record boundaries.
The five-column serial limit is tested at 99999 and refusal of 100000.

Gemmi 0.7.3 independently cross-reads 96 assemblies from 382 source blocks:
1,338 atoms, 382 ANISOU records and 1,146 CONECT target references. Residue,
atom, alternate-location, element/charge, occupancy/B-factor, tensor and explicit
connectivity identities agree after the declared serial/chain mapping. Coordinate
difference is zero. Raw coordinate/ANISOU fields also match outside the two
allowed serial/chain ranges.

Open Babel's separate three-atom join control changes two ATOM kinds, three
residue names and three chains, and emits no TER; the new API retains the
source identity and emits two boundary TER records. This is one observed CLI
path, not an all-version/API or complete protein benchmark.

Evidence and reproduction are under
`validation/results/pdb-coordinate-merge-2026-10-18/`. Tests reuse a verified
extension whose relevant native source matches the fixture revision. No fresh
current-main native build, published wheel or performance claim is made.

## Python reference

::: chematic._pdb_merge
