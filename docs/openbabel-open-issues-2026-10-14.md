# Open Babel open issue implementation audit — 2026-10-14

Inspected the 100 most recently updated open issues returned by GitHub's
`repo:openbabel/openbabel is:issue is:open` search. This is a prioritized
sample, not an audit of every open issue. Source baseline: `ffe682c3` plus
the pre-existing working-tree changes, including the separate RDKit MOL2
charge work. Existing changes were retained.

| Open Babel issue | chematic source finding | Decision |
| --- | --- | --- |
| [#2481: deletion changes remaining PDB identity](https://github.com/openbabel/openbabel/issues/2481) | Reproduced on Open Babel 3.2.1 with an original six-atom fixture; existing PDB graph projection has no source-preserving atom selector. | [Added checked PDB coordinate selection](pdb-atom-selection.md), exact atom/ANISOU lines, filtered CONECT, validated TER updates and explicit unknown-record loss consent; this row is a 2026-10-15 follow-up. |
| [#2995: CIF Cartesian coordinates labeled as fractional](https://github.com/openbabel/openbabel/issues/2995) | `write_cif(..., None)` had the same labeling defect and produced output that `parse_cif` rejected for missing cell parameters. | Implemented here, with coordinate-family and finite-value reader checks. |
| [#2655: explicit CML hydrogens gain implicit hydrogens](https://github.com/openbabel/openbabel/issues/2655) | A follow-up corrected explicit-H inference and total `hydrogenCount` reconciliation in both reader and writer. | [Implemented contract, compatibility and validation](cml-hydrogen-accounting.md): CH2, mixed H representations, isotopes, invalid totals, and a live Open Babel cross-read. |
| [#2818: singular unit cell produces NaN](https://github.com/openbabel/openbabel/issues/2818) | The unchecked plain writer accepted invalid `UnitCell` values. | [Implemented a checked Rust/Python CIF writer](cif-checked-writer.md), shared plain-reader cell validation, finite/count checks and precise float serialization. No POSCAR support; external CIF cross-read remains inconclusive on this host. |
| [#2845: PQR charge/radius lost](https://github.com/openbabel/openbabel/issues/2845) | `PqrAtomRecord` already stores charge and radius separately. `to_molecule()` is a graph projection that does not carry those fields. | [Implemented a checked rigid PQR-to-PDBQT adapter](pqr-pdbqt-charge-preservation.md): caller-supplied types, source charges/identity, radius-loss consent and a numeric-rounding sidecar. Open Babel 3.2.1 charge loss reproduced on two synthetic ions. |
| [#2817: GRO two-letter element truncation](https://github.com/openbabel/openbabel/issues/2817) | No native GRO reader existed; atom names alone cannot establish element identity. | [Added a typed single-frame reader](gro-element-resolution.md): unresolved or explicit caller elements, preserved names, nm/nm-per-ps fields, variable precision and raw boxes; 720 field/element checks and an Open Babel 3.2.1 diagnostic. |
| [#2794: GZMAT precision control](https://github.com/openbabel/openbabel/issues/2794) | Existing Gaussian input support uses Cartesian coordinates; no native Cartesian-to-Z-matrix conversion existed. Open Babel 3.2.1 produces identical rounded GZMAT with and without `-x8` on an original synthetic probe. | [Added checked caller-referenced Z-matrix conversion and writing](gzmat-precision.md), explicit 1–17 significant digits, signed torsions, singular/rounding validation and independent RDKit/Open Babel geometry checks. No Gaussian or ORCA job-execution claim. |
| [#2681: RDFile support](https://github.com/openbabel/openbabel/issues/2681) | No RDFile reader existed; native MOL/RXN parsers can supply embedded graph projections. | [Added checked typed RDFile records](rdfile-records.md), registry-only slots, ordered/duplicate data, source preservation and explicit SMILES projection; 288 two-way Indigo cross-read records with an explicit unwrapped-data profile. Default folding interoperability is reported separately. |
| [#2821: QCSchema Molecule](https://github.com/openbabel/openbabel/issues/2821) | Native QCSchema types, readers, converters, and tests already exist. | No duplicate implementation; broaden interoperability fixtures separately. |
| [#3027: five-character residue names truncated in PDB](https://github.com/openbabel/openbabel/issues/3027) | mmCIF has full residue-name fields; the protein projection is deliberately limited to amino-acid residues. | [Implemented a checked coordinate-block adapter](mmcif-pdb-checked.md): explicit namespace and metadata-loss policy, no truncated identifiers, collision/width checks and original-input sidecar. Model 1 only; no topology or full PDB deposition export. |
| [#3028: branched glycan residues merged](https://github.com/openbabel/openbabel/issues/3028) | mmCIF retains per-atom residue fields, but the protein projection excludes glycans. | [Implemented a general atom_site residue reader](mmcif-residue-identity.md) with distinct auth/label identities and strict ambiguity errors; 720 row-order permutations and an Open Babel 3.2.1 control. Glycan connectivity and PDB export remain separate. |

## Implemented coordinate contract

`write_cif` takes Cartesian coordinates in Å:

- With no cell, output uses `_atom_site_Cartn_x/y/z` and keeps the values
  unchanged to the writer's five-decimal precision. It declares no invented cell.
- With a cell, output continues to use `_atom_site_fract_x/y/z`; tests cover
  orthogonal and triclinic Cartesian-to-fractional-to-Cartesian round trips.
- The reader selects a complete coordinate family. A complete fractional
  triple takes precedence; otherwise a complete Cartesian triple is used.
  It never combines fractional and Cartesian components.
- Non-finite input coordinates and non-finite results of fractional-to-Cartesian
  conversion return `CifError::InvalidCoordinate`.

The tag semantics follow the
[IUCr Cartesian coordinate definition](https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_Cartn_.html).
This is the small-molecule CIF atom/coordinate adapter; mmCIF uses a separate
adapter and tag family.

## Reproduction and validation

The independent two-atom fixture is
`validation/openbabel_issues/cif/cartesian-no-cell.xyz`.
The installed `/opt/homebrew/bin/obabel -V` reported **Open Babel 3.2.1**.

```sh
obabel -ixyz validation/openbabel_issues/cif/cartesian-no-cell.xyz -ocif
```

It emitted `_atom_site_fract_x/y/z` without cell parameters, with the unchanged
Cartesian values `(2.5, -3.25, 4.75)` and `(-1.125, 0.625, 8.5)`.
This reproduces the reported labeling defect without 3D-generation randomness.

Before the fix, four of the five new tests failed: no-cell output labels,
mixed incomplete coordinate families, a complete Cartesian set accompanied
by a stray fractional column, and non-finite coordinates. The existing
cell-based round trip passed before and after the fix.

```sh
cargo test -p chematic-mol --test cif_coordinate_contract --offline --quiet
cargo test -p chematic-mol --offline --quiet
cargo clippy -p chematic-mol --all-targets --offline --no-deps -- -D warnings
```

The regression suite is `crates/chematic-mol/tests/cif_coordinate_contract.rs`.
Local validation passed 980 tests, including five new regression tests, with
zero failures or ignored tests. Clippy with warnings denied and formatting
checks passed for the changed Rust files.
Validation concerns local source; it does not constitute package publication
or an upstream Open Babel change.

## October 15 follow-up: HIN aromatic bonds

[Open Babel #3000](https://github.com/openbabel/openbabel/issues/3000) remains
open. An independently constructed six-carbon aromatic HIN ring reproduces
MDL query type 5 output on Open Babel 3.2.1. The new checked HIN reader maps
`a` to the existing aromatic graph bond, so chematic emits MDL type 4.
See the [bounded HIN profile and validation](openbabel-hin-aromatic-bonds.md).

## October 16 follow-up: Chemical Component Dictionary

[Open Babel #1516](https://github.com/openbabel/openbabel/issues/1516) remains
open. The existing mmCIF adapter reads `_atom_site`, not CCD molecular templates.
A [single-component CCD reader](ccd-components.md) now reads `_chem_comp_atom`
and `_chem_comp_bond` in Rust/Python/WASM. It retains stereo annotations and
requires explicit consent for an unassigned graph projection. Bulk dictionary
streaming and CCD writing are outside this implementation.
