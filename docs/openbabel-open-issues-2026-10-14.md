# Open Babel open issue implementation audit — 2026-10-14

Inspected the 100 most recently updated open issues returned by GitHub's
`repo:openbabel/openbabel is:issue is:open` search. This is a prioritized
sample, not an audit of every open issue. Source baseline: `ffe682c3` plus
the pre-existing working-tree changes, including the separate RDKit MOL2
charge work. Existing changes were retained.

| Open Babel issue | chematic source finding | Decision |
| --- | --- | --- |
| [#2669: GRO coordinates and cell scale inconsistently](https://github.com/openbabel/openbabel/issues/2669) | An original Open Babel 3.2.1 probe scales coordinates but replaces a 2×3×4 nm box with a 10 Å cube. Native GRO and LAMMPS documents had no checked unit adapter. | [Added explicit GRO→LAMMPS atomic geometry](gro-lammps-units.md) for Rust/Python: caller types/masses, real/metal units, cell/tilt/velocity scaling and exact source/ID sidecars. No molecular topology inference or LAMMPS engine execution claim. |
| [#3030: auto-detected gzip state leaks into later inputs](https://github.com/openbabel/openbabel/issues/3030) | Reproduced plain string/stream replacement failures on Open Babel 3.2.1. Native checked SDF files had no gzip selection; no analogous shared conversion session exists. | [Added per-input SDF compression](sdf-gzip-input.md) in Rust/Python: explicit plain/gzip/auto, magic detection, concatenated members, terminal integrity errors and decoded-data limits. This adds a missing capability, not a demonstrated chematic sticky-state fix. |
| [#3032: Clear retains the previous title](https://github.com/openbabel/openbabel/issues/3032) | A four-cycle probe reproduces the title leak in Open Babel 3.2.1; chematic had no full reset entry point. | [Added complete molecule reset](molecule-clear.md): Rust graph/metadata/cache clearing and native/compatibility Python title/property reset, preserving snapshots. This row is a 2026-10-16 follow-up. |
| [#2811: dielectric changes leave prepared electrostatics stale](https://github.com/openbabel/openbabel/issues/2811) | Reproduced on Open Babel 3.2.1 with an original two-ion probe. Prepared native MMFF fields had no dielectric-update API. | [Added checked constant-dielectric updates](mmff-dielectric-updates.md) to Rust/Python MMFF94/MMFF94s fields, with immediate energy/gradient/minimizer use and atomic refusal of invalid settings. UFF/GAFF and solvent modelling are outside this scope. |
| [#2101: PDBQT coordinate columns overflow](https://github.com/openbabel/openbabel/issues/2101) | Reproduced on Open Babel 3.2.1 with authored XYZ; the legacy native graph writer also treats format widths as minimum widths. | [Added checked rigid PDBQT output](pdbqt-checked-output.md) in Rust/Python/WASM, with exact counts, post-rounding widths, finite values, serial limits and header validation. Existing writer behavior is retained. |
| [#2668: Kekulé bonds in MOL2 output](https://github.com/openbabel/openbabel/issues/2668) | Graph writers and exact source round-trips retain `ar` bond tokens; no explicit bond-encoding choice existed. | [Added a bounded source-token conversion](mol2-kekule-bond-output.md) in Rust/Python/WASM, retaining atom types and all other source bytes; neutral C/N/O/S cycles only, with matching, H/valence and aromaticity checks. |
| [#2671: RMSD for all MOL2 reference conformations](https://github.com/openbabel/openbabel/issues/2671) | Native single-record readers rejected concatenated MOL2; no explicit all-reference-versus-one RMSD API existed. Open Babel's default pairs nth records; reversed inputs with `-f` work, so the one-line diagnostic alone does not establish a parser bug. | [Added bounded source-preserving record iteration and checked all-reference RMSD](mol2-reference-rmsd.md) in Rust/Python, with stable source/indexed-graph correspondence, optional proper rigid alignment and errors rather than skipped records. |
| [#2350: order of MOL2 sections](https://github.com/openbabel/openbabel/issues/2350) | Default editable/generated output also places UNITY between ATOM and BOND. Moving it globally causes Open Babel 3.2.1 to lose formal charges in six original graph controls. | [Added an explicit adjacent-core output option](mol2-section-order.md); default bytes and source preservation remain unchanged. This row is a 2026-10-16 follow-up; Antechamber was not executed. |
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

## October 16 follow-up: explicit CIF molecular bonds

[#2726](https://github.com/openbabel/openbabel/issues/2726) remains open, but
the comments explain that `-xg` already enables Open Babel bond output; the
reported `-ag` selects a read option. Do not count it as missing upstream
functionality. Added the previously missing chematic
[checked molecular bond writer](cif-molecular-bonds.md) in Rust/Python, with
explicit graph edges, S/D/T/A types and finite positive Cartesian distances.
The positions-only reader and periodic adapter retain their existing scope.
## October 16 follow-up: XYZ element symbols

[#2333](https://github.com/openbabel/openbabel/issues/2333) remains open.
Open Babel 3.2.1 maps original CL/BR/CA/FE XYZ probes to dummy atoms; the
previous chematic readers refused them. Added
[format-local uppercase species support](xyz-element-symbols.md) to native
XYZ/extxyz and their stream readers. Core lookup, SMILES grammar and PDB atom
names remain case-sensitive; coordinates and opaque properties are unaffected.
## October 16 follow-up: total hydrogen accounting

[#3008](https://github.com/openbabel/openbabel/issues/3008) and
[#2624](https://github.com/openbabel/openbabel/issues/2624) remain open.
Open Babel 3.2.1 reproduces formula-H loss on original explicit-H/component
probes. Chematic formulas already pass, but its separate total-H descriptor
double-counted bracket H and missed H atom nodes. [Fix the counter and WASM
atom information](hydrogen-count-accounting.md), keeping formula and
serialization behavior as controls.

## October 17 follow-up: named CML atoms

[#2120](https://github.com/openbabel/openbabel/issues/2120) remains open. Added
a [checked Rust/Python writer](cml-named-atoms.md) with exact caller-supplied
atom IDs, consistent bond references and optional Cartesian 3D coordinates.
The documented MOL2 recipe preserves names in this graph export; source-only
MOL2 metadata and unsupported stereo remain outside the contract.
