# Component serialization and isolated-ion MMFF coverage

This change exercises existing coverage targets without changing workflow scope,
ignores, thresholds, or exclusions. Reference data use the pinned native RDKit
2026.03.1 wheel; no RDKit installation or network access is needed in Rust CI.

## Behavioral correction

MMFF assigned covalent halogen types to isolated halide anions, losing their
ionic formal charge. It also rejected supported isolated metal cations despite
already containing their numeric type registry, charge, and van der Waals tables.
The typer now selects these existing rows for isolated, correctly charged ions:
F-/Cl-/Br-, Li+/Na+/K+, Mg2+/Ca2+/Zn2+, Fe2+/Fe3+, and Cu+/Cu2+.
Other charge states and bonded metals retain the unsupported-type error.
The reference implementation is RDKit's
[`AtomTyper.cpp`, Release_2026_03_1](https://github.com/rdkit/rdkit/blob/Release_2026_03_1/Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp).

The nonbonded fixture checks 13 ion pairs at separations 0, 1e-8, 1 and 4 angstrom,
with interfragment interactions enabled and disabled: 104 cases. Every case
compares atom types, partial charges, energy, and all gradient components with
native RDKit. At nonzero separations of at least 1 angstrom, an independent
central finite difference also checks the gradient. Exactly coincident atoms
exercise RDKit's deterministic collision guard, not a physically meaningful
minimum or a differentiable geometry. Pair gradients must be antisymmetric;
ignoring interfragment interactions must give exactly zero energy and gradient.

## Additional reference contracts

- 196 rooted SMARTS and 17 CXSMARTS outputs cover disconnected salts, isotope
  and atom-map metadata, tetrahedral and double-bond stereo, dative bonds,
  radicals, and the empty graph. An empty `MoleculeBuilder` graph is used for
  writer tests because the SMILES parser intentionally rejects empty input.
- 4,698 branched-polyene spellings compare canonical parsing/writing with native
  RDKit. Reference generation reparses each actual serialized input before
  freezing its expected output. Native ambiguity cleanup can remove conflicting
  direction assignments; these cases do not assert physical feasibility or
  preservation of every externally assigned stereo tag.
- 25 InChI inputs come from unusual-valence cleanup families documented by
  RDKit's `External/INCHI-API/inchi.cpp`. The raw unsanitized structures are only
  generator inputs. Rust receives the generated identifier and must agree with
  native InChI reading or rejection. This does not assert that the raw source
  encodings are valid sanitized molecules.
- 78 default 12-iteration Gasteiger cases compare implicit and explicit hydrogen
  representations. Unsupported parameter families retain native NaN/infinity
  classification rather than being presented as usable finite charges. Finite
  fully explicit charge sums must equal the molecular formal charge.

- 106 PDB fixed-width and CONECT cases compare native reading/rejection,
  element-name inference, charge fields, short records, line endings, and
  repeated/reversed bonds. Native zero-order inter-residue bonds and radicals
  check the documented explicit unsupported-feature error instead of silently
  discarding these features.
- 171 SMARTS queries and 11 reaction SMARTS inputs check Boolean precedence,
  atom/bond negation, ranges, isotopes, atom maps, grouped components, directional
  bonds, and syntax rejection against native serialization.
- 13 raw cleanup inputs compare actual native charge/bond/dative-donor edits
  before perception. Six inputs change and seven are unchanged. The public
  inspection API must report edits without modifying its input graph.
- A deliberately overvalent public graph checks that standardization retains
  input and per-stage warnings, distinguishes unchanged from changed-with-warning
  reports, and keeps the source graph intact. It is a validation test, not a
  chemically valid-molecule reference.

## Reproduction

Set `PYTHONPATH` to a directory containing the pinned `rdkit==2026.03.1` wheel,
then run the eight `scripts/generate_*_boundary_fixtures.py` generators named for
`smarts_components`, `branched_directions`, `inchi_unusual`, `gasteiger`, and
`mmff_nonbonded`, `pdb_records`, `query_serialization`, and `cleanup_edits`. Each generator asserts the RDKit version. Regeneration must
produce byte-identical JSON files under `validation/`.

Native InChI reference generation and polyene parsing emit diagnostic warnings
for the deliberately unusual input encodings. Expected results are measured
from the actual native reader, including rejection and ambiguity cleanup.

Fixture SHA-256 values:

| Fixture family | SHA-256 |
| --- | --- |
| smarts-components | `0c56d91cadaf30b24b011f9f922423dcaca5e29f0915f87ea8e29953af9851ea` |
| branched-directions | `2037144ce6b1ffc07a9ecbcd70d555447e28eac2c5ddcc80e9106a1674186195` |
| inchi-unusual | `c3527ead8cd502b53f4cb8503f858a9a76e38fa8781e1dc08f61e5203febaec8` |
| gasteiger | `d7f8c0c8a67931ea995881ba3d7ee98a73cac2df852468e097313d3fc22ae067` |
| mmff-nonbonded | `5c8d6b65629437837bda880307d697005928d80a5486b9139d8b6f4208f6b078` |
| pdb-records | `3562cb2991f1cd2e1cc9801ed477d1020e808172edba2e5010de12b6a69152b5` |
| query-serialization | `23d16093c4404523f02fd42d4b57739c86fde382eb9286f7f8b0f1e9194d83ff` |
| cleanup-edits | `c65a25b1be158ffde477eb4635ba2d7ac79373cbe84dce163199f18ac8f8fd81` |

## Validation result

Clean local coverage with the unchanged workflow command passed 5,384 tests
(9 ignored) across 101 test executables. LCOV records 155,174 of 161,869 lines
in 340 files: 95.863939% locally. This is a macOS source measurement; the Linux
workflow upload and processed Codecov main report remain authoritative for the
published percentage. Workspace all-target Clippy (`-D warnings`), formatting,
and `git diff --check` passed. All eight fixture families regenerated with
byte-identical SHA-256 values.

The earlier attempt exhausted disk space and failed to save some LLVM profiles.
Its incomplete report was discarded. The numbers above come from a fresh
`cargo llvm-cov clean --workspace` followed by the complete successful run.
