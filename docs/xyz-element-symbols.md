# XYZ uppercase element symbols

Native XYZ and Extended XYZ readers now accept canonical element symbols and
uppercase two-letter aliases: `CL` → Cl, `BR` → Br, `CA` → Ca, `FE` → Fe, and
the corresponding aliases for all 118 supported elements.

[Open Babel #2333](https://github.com/openbabel/openbabel/issues/2333) reports
chlorine becoming a dummy atom when an XYZ file uses `CL`. The issue remains
open as inspected on 2026-10-16. Original one-atom probes reproduce the
identity loss on Open Babel 3.2.1 for CL, BR, CA and FE. The previous chematic
reader instead rejected those aliases with an unknown-element error.

## Format-local contract

The shared `chematic_3d::xyz::element_from_xyz_symbol` helper accepts existing
canonical tokens first, then normalizes exactly two uppercase ASCII letters.
It is used only by the native 3D XYZ reader and the molecular XYZ/extxyz
species readers. Trajectory and buffered readers inherit the same rule.

- `CA` is calcium **in the XYZ species column**. PDB atom names, such as an
  alpha-carbon name `CA`, are not processed by this helper.
- Core `Element::from_symbol` and SMILES remain case-sensitive. Aromatic `c`
  still has its existing SMILES meaning; `CL` is not made valid SMILES.
- Lowercase/noncanonical mixed-case aliases, suffixed atom labels (`CL1`),
  atomic-number tokens, D/T aliases, unknown symbols and Unicode lookalikes
  retain their previous rejection behavior. Unknown input never becomes `*`.
- Atom count, row order, coordinates, explicit H, comments, lattice metadata
  and non-species property values are retained under their existing contracts.
  A string property `name=CA` or frame value `label=CL` is not normalized.
- Writers emit canonical symbols. Raw source spelling is not preserved by
  the typed XYZ atom record. Existing error variants and limits are unchanged.

## Existing entry points

Rust: `chematic_3d::parse_xyz`, `chematic_mol::parse_xyz`, `parse_extxyz`, and
their trajectory/buffered readers. Python: `from_xyz`, `from_extxyz` and
`from_extxyz_all`. WASM: `mol_from_xyz`, `mol_from_extxyz`, extxyz frame JSON
and the existing XYZ/extxyz trajectory paths.

```python
import chematic

mol, coords = chematic.from_xyz("1\nchlorine\nCL 0 0 0\n")
print(mol.to_xyz(coords))  # writes Cl
```

XYZ carries positions and elements, not bonds. These entry points do not
infer connectivity or reproduce Open Babel's full XYZ-to-SMILES workflow.
All explicit H rows remain atoms; downstream chemical perception is separate.
The WASM topology-only XYZ helper still discards coordinates; use the extxyz
frame JSON accessor when coordinates are required. The separately versioned
RDKit-compatible XYZ path retains its existing behavior.

## Reproduction and verification

The original five-atom fixture is
`validation/openbabel_issues/xyz-elements/chloromethane-uppercase.xyz`.
`crates/chematic-mol/tests/xyz_symbol_case.rs` exercises both native parsers,
all 118 elements in two casing modes, streaming, canonical writing, metadata,
strict invalid-token handling and the core/SMILES boundary.

`scripts/validation/validate_xyz_symbol_case.py` uses RDKit 2026.09.1's independent
periodic table for 236 element/casing cases (222 distinct species tokens).
It checks Python and optional real Node/WASM, diagnoses Open Babel with
canonical and uppercase controls, and optionally records a retained pre-fix
runtime in a separate process. Its compressed packet contains every input
and external output; the summary records source/runtime hashes and separate
canonical/uppercase mismatch counts. This is element/coordinate identity
evidence for a source candidate, not a published-package or general chemical
conversion benchmark.

The 2026-10-16 source run matched all 236 cases through each Python reader
(472 checks); real Node/WASM matched another 472 reader checks and all 236
extxyz coordinate checks. The retained pre-fix source runtime rejected 104 uppercase
aliases. Open Babel 3.2.1 preserved all 118 canonical controls, but replaced
104 uppercase two-letter species with dummy atoms. One-letter symbols are
identical in the two casing modes. These counts describe this corpus and
installed runtime; they do not measure bond perception or other formats.

Local checks passed: 74 Rust tests (66 existing XYZ tests, seven new
regressions and one existing public-API compatibility test), 20 Python
tests, the Node/WASM regression script, Clippy with warnings denied,
formatting, the shared binding/format contracts and strict MkDocs build.
The summary and compressed packet are retained under
`validation/results/openbabel-xyz-elements-2026-10-16*`.
